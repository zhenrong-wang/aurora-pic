#include "pic/Simulation2D.hpp"
#include "pic/Convergence.hpp"
#include "pic/ParticleState.hpp"
#include "pic/Pusher.hpp"
#include "pic/Runtime.hpp"
#include "pic/Units.hpp"
#include "pic/VTKWriter.hpp"
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <stdexcept>
#include <string>
#include <utility>

namespace pic {
namespace {
constexpr const char* kCheckpointMagicV1 = "AuroraPIC-checkpoint-v1";
constexpr const char* kCheckpointMagicV2 = "AuroraPIC-checkpoint-v2";
constexpr const char* kCheckpointMagicV3 = "AuroraPIC-checkpoint-v3";

double wrap_periodic(double value, double length) {
    return std::fmod(std::fmod(value, length) + length, length);
}

ParticleBoundary resolve_particle_boundary(ParticleBoundary configured, Boundary field_boundary) {
    if (configured != ParticleBoundary::Auto) return configured;
    return field_boundary == Boundary::Periodic ? ParticleBoundary::Periodic : ParticleBoundary::Absorbing;
}

void reflect_coordinate(double& coordinate, double& velocity, double length) {
    while (coordinate < 0.0 || coordinate > length) {
        if (coordinate < 0.0) {
            coordinate = -coordinate;
            velocity = -velocity;
        } else {
            coordinate = 2.0 * length - coordinate;
            velocity = -velocity;
        }
    }
}

bool apply_lower_boundary(double& coordinate,
                          double& velocity,
                          double length,
                          ParticleBoundary boundary,
                          std::size_t& absorbed_count) {
    if (coordinate >= 0.0) return true;
    switch (boundary) {
        case ParticleBoundary::Absorbing:
            ++absorbed_count;
            return false;
        case ParticleBoundary::Reflecting:
            reflect_coordinate(coordinate, velocity, length);
            return true;
        case ParticleBoundary::Periodic:
            coordinate = wrap_periodic(coordinate, length);
            return true;
        case ParticleBoundary::Auto:
            throw std::logic_error("unresolved lower 2D particle boundary policy");
    }
    return true;
}

bool apply_upper_boundary(double& coordinate,
                          double& velocity,
                          double length,
                          ParticleBoundary boundary,
                          std::size_t& absorbed_count) {
    if (coordinate <= length) return true;
    switch (boundary) {
        case ParticleBoundary::Absorbing:
            ++absorbed_count;
            return false;
        case ParticleBoundary::Reflecting:
            reflect_coordinate(coordinate, velocity, length);
            return true;
        case ParticleBoundary::Periodic:
            coordinate = wrap_periodic(coordinate, length);
            return true;
        case ParticleBoundary::Auto:
            throw std::logic_error("unresolved upper 2D particle boundary policy");
    }
    return true;
}

std::filesystem::path checkpoint_path_for_step(const Simulation2DConfig& cfg, std::size_t step) {
    if (!cfg.checkpoint_path.empty()) return cfg.checkpoint_path;
    return cfg.output_dir / ("checkpoint_" + std::to_string(step) + ".apc");
}

void ensure_parent_directory(const std::filesystem::path& path) {
    const auto parent = path.parent_path();
    if (!parent.empty()) std::filesystem::create_directories(parent);
}

template <typename T>
void require_stream(T& stream, const std::string& message) {
    if (!stream) throw std::runtime_error(message);
}

bool has_magnetic_field(const Simulation2DConfig& cfg) {
    return cfg.magnetic_field_profile.has_value() ||
           cfg.magnetic_field_x != 0.0 ||
           cfg.magnetic_field_y != 0.0 ||
           cfg.magnetic_field_z != 0.0;
}

Vec3 magnetic_field(
    const Simulation2DConfig& cfg,
    Vec2 position) {
    if (cfg.magnetic_field_profile) {
        return cfg.magnetic_field_profile->evaluate(
            {position.x, position.y, 0.0});
    }
    return {
        cfg.magnetic_field_x,
        cfg.magnetic_field_y,
        cfg.magnetic_field_z};
}

void initialize_particle_pusher(Particle2D& particle, Vec2 electric, double charge_to_mass, const Simulation2DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        initialize_boris_half_step(
            particle, electric, magnetic_field(cfg, particle.position),
            charge_to_mass, cfg.dt);
    } else {
        initialize_leapfrog_half_step(particle, electric, charge_to_mass, cfg.dt);
    }
}

void kick_particle(Particle2D& particle, Vec2 electric, double charge_to_mass, const Simulation2DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        kick_boris(
            particle, electric, magnetic_field(cfg, particle.position),
            charge_to_mass, cfg.dt);
    } else {
        kick_leapfrog(particle, electric, charge_to_mass, cfg.dt);
    }
}

void synchronize_particle(Particle2D& particle, Vec2 electric, double charge_to_mass, const Simulation2DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        synchronize_boris(
            particle, electric, magnetic_field(cfg, particle.position),
            charge_to_mass, cfg.dt);
    } else {
        synchronize_leapfrog(particle, electric, charge_to_mass, cfg.dt);
    }
}

void write_vtk_outputs(const Mesh2D& mesh, const std::filesystem::path& output_dir, std::size_t step, VTKOutputFormat format) {
    const auto stem = output_dir / ("fields_" + std::to_string(step));
    if (format == VTKOutputFormat::Legacy || format == VTKOutputFormat::Both) {
        write_legacy_vtk(mesh, stem.string() + ".vtk");
    }
    if (format == VTKOutputFormat::Xml || format == VTKOutputFormat::Both) {
        write_vtk_xml(mesh, stem.string() + ".vts");
    }
}
} // namespace

Simulation2D::Simulation2D(Simulation2DConfig cfg)
    : cfg_(std::move(cfg)),
      mesh_(cfg_.nx, cfg_.ny, cfg_.length_x, cfg_.length_y,
            cfg_.boundary, cfg_.boundary_config),
      solver_(cfg_.units.permittivity()),
      rng_(cfg_.seed) {
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) cfg_.checkpoint_interval = cfg_.output_interval;
    if (!cfg_.restart_path.empty() &&
        !cfg_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "2D restart_path and initial_state_path are mutually exclusive");
    }
    if (cfg_.initial_state_signature &&
        cfg_.initial_state_path.empty()) {
        throw std::invalid_argument(
            "2D initial_state_signature requires initial_state_path");
    }
    validate_initialization_acceptance(
        cfg_.initialization_acceptance,
        "2D initialization acceptance config");
    if (!std::isfinite(cfg_.dt) || cfg_.dt <= 0.0) throw std::invalid_argument("2D simulation dt must be positive and finite");
    if (cfg_.output_interval == 0) throw std::invalid_argument("2D output_interval must be positive");
    if (cfg_.particle_output_stride == 0) throw std::invalid_argument("2D particle_output_stride must be positive");
    if (!std::isfinite(cfg_.magnetic_field_x) ||
        !std::isfinite(cfg_.magnetic_field_y) ||
        !std::isfinite(cfg_.magnetic_field_z)) {
        throw std::invalid_argument(
            "2D magnetic_field components must be finite");
    }
    if (cfg_.magnetic_field_profile) {
        if (cfg_.magnetic_field_x != 0.0 ||
            cfg_.magnetic_field_y != 0.0 ||
            cfg_.magnetic_field_z != 0.0) {
            throw std::invalid_argument(
                "2D uniform magnetic_field components and magnetic_field_profile are mutually exclusive");
        }
        cfg_.magnetic_field_profile->validate_domain(
            {0.0, 0.0, 0.0},
            {cfg_.length_x, cfg_.length_y, 0.0},
            "2D simulation");
    }
    if (cfg_.mode == RunMode::SteadyState) {
        if (cfg_.max_steps == 0) throw std::invalid_argument("2D max_steps must be positive for steady-state mode");
        if (cfg_.steady_window == 0) throw std::invalid_argument("2D steady_window must be positive");
        if (!std::isfinite(cfg_.steady_tolerance) || cfg_.steady_tolerance <= 0.0) {
            throw std::invalid_argument("2D steady_tolerance must be positive and finite");
        }
    }
    validate_runtime_policy(cfg_.runtime);
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) {
        throw std::invalid_argument("2D checkpoint_interval must be positive when checkpoint_output is enabled");
    }
    cfg_.particle_boundary_config.left = resolve_particle_boundary(cfg_.particle_boundary_config.left, cfg_.boundary);
    cfg_.particle_boundary_config.right = resolve_particle_boundary(cfg_.particle_boundary_config.right, cfg_.boundary);
    cfg_.particle_boundary_config.bottom = resolve_particle_boundary(cfg_.particle_boundary_config.bottom, cfg_.boundary);
    cfg_.particle_boundary_config.top = resolve_particle_boundary(cfg_.particle_boundary_config.top, cfg_.boundary);
    for (const auto& sc : cfg_.species) species_.emplace_back(sc);
    if (species_.empty()) species_.emplace_back(Species2DConfig{});
}

void Simulation2D::initialize() {
    time_ = 0.0;
    step_ = 0;
    boundary_losses_ = {};
    if (cfg_.initial_state_path.empty()) {
        for (auto& sp : species_) sp.initialize(mesh_, rng_);
    } else {
        std::vector<ExternalSpeciesExpectation> expected;
        expected.reserve(species_.size());
        for (const auto& species : species_) {
            expected.push_back({
                species.name(),
                species.config().particles});
        }
        for (auto& species : species_) {
            species.particles().resize(
                species.config().particles);
        }
        initial_state_metadata_ =
            load_validated_external_particle_state_bounded(
                cfg_.initial_state_path, 2,
                cfg_.units.system, expected,
                "2D simulation",
                [&](std::size_t species_index,
                    std::size_t record_index,
                    const ExternalParticleRecord& record) {
                    auto& species =
                        species_.at(species_index);
                    const auto& species_config =
                        species.config();
                    const double maximum_x =
                        species_config.init_x_max < 0.0
                            ? mesh_.length_x()
                            : species_config.init_x_max;
                    const double maximum_y =
                        species_config.init_y_max < 0.0
                            ? mesh_.length_y()
                            : species_config.init_y_max;
                    const bool outside =
                        record.position.x <
                            species_config.init_x_min ||
                        record.position.x > maximum_x ||
                        record.position.y <
                            species_config.init_y_min ||
                        record.position.y > maximum_y ||
                        (mesh_.boundary() == Boundary::Periodic &&
                         (record.position.x == mesh_.length_x() ||
                          record.position.y == mesh_.length_y()));
                    if (outside) {
                        throw std::runtime_error(
                            "external particle for species '" +
                            species.name() +
                            "' lies outside the 2D domain");
                    }
                    auto& particle =
                        species.particles().at(record_index);
                    particle.position = {
                        record.position.x,
                        record.position.y};
                    particle.velocity = {
                        record.velocity.x,
                        record.velocity.y};
                    particle.velocity_z =
                        record.velocity.z;
                    particle.velocity_half =
                        particle.velocity;
                    particle.velocity_half_z =
                        particle.velocity_z;
                    particle.alive = true;
                },
                cfg_.initial_state_signature);
    }
    deposit_and_solve();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& particle = particles[particle_id];
            if (particle.alive) initialize_particle_pusher(particle, interpolate_electric(mesh_, particle.position), qm, cfg_);
        });
    }
    initialized_ = true;
}

void Simulation2D::deposit_and_solve() {
    mesh_.clear_charge();
    for (const auto& sp : species_) sp.deposit_charge(mesh_);
    solver_.solve(mesh_);
}

void Simulation2D::apply_particle_boundaries(Particle2D& particle) {
    if (!apply_lower_boundary(particle.position.x, particle.velocity_half.x, mesh_.length_x(),
                              cfg_.particle_boundary_config.left, boundary_losses_.absorbed_left)) {
        particle.alive = false;
        return;
    }
    if (!apply_upper_boundary(particle.position.x, particle.velocity_half.x, mesh_.length_x(),
                              cfg_.particle_boundary_config.right, boundary_losses_.absorbed_right)) {
        particle.alive = false;
        return;
    }
    if (!apply_lower_boundary(particle.position.y, particle.velocity_half.y, mesh_.length_y(),
                              cfg_.particle_boundary_config.bottom, boundary_losses_.absorbed_bottom)) {
        particle.alive = false;
        return;
    }
    if (!apply_upper_boundary(particle.position.y, particle.velocity_half.y, mesh_.length_y(),
                              cfg_.particle_boundary_config.top, boundary_losses_.absorbed_top)) {
        particle.alive = false;
        return;
    }
}

void Simulation2D::step() {
    if (!initialized_) initialize();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        for (auto& particle : sp.particles()) {
            if (!particle.alive) continue;
            kick_particle(particle, interpolate_electric(mesh_, particle.position), qm, cfg_);
            drift_leapfrog(particle, cfg_.dt);
            apply_particle_boundaries(particle);
        }
    }

    deposit_and_solve();
    for (auto& sp : species_) {
        const double qm = sp.charge() / sp.mass();
        auto& particles = sp.particles();
        runtime_parallel_for(std::size_t{0}, particles.size(), cfg_.runtime, [&](std::size_t particle_id) {
            auto& particle = particles[particle_id];
            if (particle.alive) synchronize_particle(particle, interpolate_electric(mesh_, particle.position), qm, cfg_);
        });
    }
    ++step_;
    time_ += cfg_.dt;
}

void Simulation2D::save_checkpoint(const std::filesystem::path& path) const {
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open 2D checkpoint for writing: " + path.string());
    out << std::setprecision(17);
    out << kCheckpointMagicV3 << '\n';
    out << "dimension 2\n";
    out << "units " << to_string(cfg_.units.system) << ' '
        << cfg_.units.relative_permittivity << ' '
        << cfg_.units.permittivity() << "\n";
    out << "step " << step_ << "\n";
    out << "time " << time_ << "\n";
    out << "boundary_losses " << boundary_losses_.absorbed_left << ' ' << boundary_losses_.absorbed_right << ' '
        << boundary_losses_.absorbed_bottom << ' ' << boundary_losses_.absorbed_top << "\n";
    out << "species_count " << species_.size() << "\n";
    out << "rng " << rng_ << "\n";
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& sp = species_[species_id];
        out << "species " << species_id << ' ' << sp.name() << ' ' << sp.particles().size() << "\n";
        for (const auto& p : sp.particles()) {
            out << p.position.x << ' ' << p.position.y << ' '
                << p.velocity.x << ' ' << p.velocity.y << ' '
                << p.velocity_z << ' '
                << p.velocity_half.x << ' ' << p.velocity_half.y << ' '
                << p.velocity_half_z << ' '
                << (p.alive ? 1 : 0) << "\n";
        }
    }
    require_stream(out, "failed while writing 2D checkpoint: " + path.string());
}

void Simulation2D::load_checkpoint(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open 2D checkpoint for reading: " + path.string());
    std::string magic;
    std::getline(in, magic);
    const bool checkpoint_v1 = magic == kCheckpointMagicV1;
    const bool checkpoint_v2 = magic == kCheckpointMagicV2;
    const bool checkpoint_v3 = magic == kCheckpointMagicV3;
    if (!checkpoint_v1 && !checkpoint_v2 && !checkpoint_v3) {
        throw std::runtime_error("invalid checkpoint magic in: " + path.string());
    }

    std::string key;
    unsigned dimension = 0;
    in >> key >> dimension;
    if (key != "dimension" || dimension != 2) throw std::runtime_error("checkpoint dimension does not match 2D simulation");
    in >> key;
    if (key == "units") {
        std::string unit_system;
        double relative_permittivity = 0.0;
        double permittivity = 0.0;
        in >> unit_system >> relative_permittivity >> permittivity;
        if ((!checkpoint_v2 && !checkpoint_v3) ||
            unit_system != to_string(cfg_.units.system) ||
            relative_permittivity != cfg_.units.relative_permittivity ||
            permittivity != cfg_.units.permittivity()) {
            throw std::runtime_error(
                "checkpoint unit system does not match 2D config");
        }
        in >> key;
    } else if (checkpoint_v2 || checkpoint_v3 ||
               cfg_.units.system != UnitSystem::Normalized ||
               cfg_.units.relative_permittivity != 1.0) {
        throw std::runtime_error(
            "legacy checkpoint without unit metadata requires normalized units");
    }
    in >> step_;
    if (key != "step") throw std::runtime_error("checkpoint missing step");
    in >> key >> time_;
    if (key != "time") throw std::runtime_error("checkpoint missing time");
    in >> key >> boundary_losses_.absorbed_left >> boundary_losses_.absorbed_right
       >> boundary_losses_.absorbed_bottom >> boundary_losses_.absorbed_top;
    if (key != "boundary_losses") throw std::runtime_error("checkpoint missing 2D boundary loss counters");
    std::size_t species_count = 0;
    in >> key >> species_count;
    if (key != "species_count" || species_count != species_.size()) throw std::runtime_error("checkpoint species count does not match 2D config");
    in >> key;
    if (key != "rng") throw std::runtime_error("checkpoint missing rng state");
    in >> rng_;

    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        std::size_t stored_species_id = 0;
        std::string stored_name;
        std::size_t particle_count = 0;
        in >> key >> stored_species_id >> stored_name >> particle_count;
        if (key != "species" || stored_species_id != species_id || stored_name != species_[species_id].name()) {
            throw std::runtime_error("checkpoint species metadata does not match 2D config");
        }
        auto& particles = species_[species_id].particles();
        particles.resize(particle_count);
        for (auto& p : particles) {
            int alive = 0;
            in >> p.position.x >> p.position.y
               >> p.velocity.x >> p.velocity.y;
            if (checkpoint_v3) {
                in >> p.velocity_z;
            } else {
                p.velocity_z = 0.0;
            }
            in >> p.velocity_half.x >> p.velocity_half.y;
            if (checkpoint_v3) {
                in >> p.velocity_half_z;
            } else {
                p.velocity_half_z = 0.0;
            }
            in >> alive;
            p.alive = alive != 0;
        }
    }
    require_stream(in, "failed while reading 2D checkpoint: " + path.string());
    deposit_and_solve();
    initialized_ = true;
}

RunSummary2D Simulation2D::run() {
    if (!cfg_.restart_path.empty()) load_checkpoint(cfg_.restart_path);
    else initialize();

    write_unit_metadata(cfg_.output_dir, cfg_.units, 2);
    if (!cfg_.initial_state_path.empty()) {
        write_external_particle_state_metadata(
            cfg_.output_dir /
                "initial_state_metadata.txt",
            cfg_.initial_state_path,
            initial_state_metadata_,
            cfg_.initial_state_signature);
    }
    std::vector<InitializationSpeciesMoments> initialization_moments;
    initialization_moments.reserve(species_.size());
    for (const auto& species : species_) {
        initialization_moments.push_back(
            summarize_initialization(species));
    }
    write_initialization_report(
        cfg_.output_dir / "initialization.csv", 2,
        !cfg_.restart_path.empty()
            ? "restart"
            : (!cfg_.initial_state_path.empty()
                   ? "external"
                   : "generated"),
        initialization_moments);
    const auto initialization_acceptance =
        assess_initialization_acceptance(
            cfg_.initialization_acceptance,
            initialization_moments, 3);
    write_initialization_acceptance_report(
        cfg_.output_dir / "initialization_acceptance.csv",
        initialization_acceptance);
    enforce_initialization_acceptance(
        initialization_acceptance);
    Diagnostics2D diag(
        cfg_.output_dir, species_, cfg_.units.permittivity());
    diag.write_header();
    auto s0 = diag.sample(step_, time_, mesh_, species_, boundary_losses_);
    diag.write_sample(s0);
    if (cfg_.vtk_output) write_vtk_outputs(mesh_, cfg_.output_dir, step_, cfg_.vtk_format);
    if (cfg_.particle_output) {
        diag.write_particle_sample(step_, species_, cfg_.particle_output_stride, cfg_.particle_sample_count);
    }
    if (cfg_.checkpoint_output) save_checkpoint(checkpoint_path_for_step(cfg_, step_));

    const std::size_t particle_interval = cfg_.particle_output_interval == 0 ? cfg_.output_interval : cfg_.particle_output_interval;
    const std::size_t limit = cfg_.mode == RunMode::SteadyState ? cfg_.max_steps : cfg_.steps;
    RunSummary2D summary;
    summary.final_sample = s0;
    while (step_ < limit) {
        step();
        bool reached_steady = false;
        if (step_ % cfg_.output_interval == 0 || step_ == limit) {
            auto s = diag.sample(step_, time_, mesh_, species_, boundary_losses_);
            diag.write_sample(s);
            summary.final_sample = s;
            if (cfg_.vtk_output) {
                write_vtk_outputs(mesh_, cfg_.output_dir, step_, cfg_.vtk_format);
            }
            reached_steady = cfg_.mode == RunMode::SteadyState &&
                             adjacent_energy_windows_converged(diag.history(), cfg_.steady_window, cfg_.steady_tolerance);
            if (reached_steady) summary.steady_state_reached = true;
        }
        if (cfg_.particle_output && (step_ % particle_interval == 0 || step_ == limit || reached_steady)) {
            diag.write_particle_sample(step_, species_, cfg_.particle_output_stride, cfg_.particle_sample_count);
        }
        if (cfg_.checkpoint_output &&
            (step_ % cfg_.checkpoint_interval == 0 || step_ == limit || reached_steady)) {
            save_checkpoint(checkpoint_path_for_step(cfg_, step_));
        }
        if (reached_steady) break;
    }
    summary.steps_completed = step_;
    summary.final_time = time_;
    if (summary.final_sample.step != step_) summary.final_sample = sample();
    return summary;
}

DiagnosticSample2D Simulation2D::sample() const {
    DiagnosticSample2D s;
    s.step = step_;
    s.time = time_;
    s.boundary_losses = boundary_losses_;
    s.live_particles_by_species.reserve(species_.size());
    for (const auto& sp : species_) {
        s.kinetic_energy += sp.kinetic_energy();
        const auto live = sp.live_count();
        s.live_particles += live;
        s.live_particles_by_species.push_back(live);
    }
    for (std::size_t j = 0; j < mesh_.ny(); ++j) {
        for (std::size_t i = 0; i < mesh_.nx(); ++i) {
            const auto idx = mesh_.index(i, j);
            const double e2 = mesh_.electric_x()[idx] * mesh_.electric_x()[idx] + mesh_.electric_y()[idx] * mesh_.electric_y()[idx];
            const double area = mesh_.node_area(i, j);
            s.field_energy +=
                0.5 * cfg_.units.permittivity() * e2 * area;
            s.charge_l1 += std::abs(mesh_.rho()[idx]) * area;
        }
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    return s;
}
}
