#include "pic/Simulation3D.hpp"
#include "pic/Convergence.hpp"
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
            throw std::logic_error("unresolved lower 3D particle boundary policy");
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
            throw std::logic_error("unresolved upper 3D particle boundary policy");
    }
    return true;
}

std::filesystem::path checkpoint_path_for_step(const Simulation3DConfig& cfg, std::size_t step) {
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

bool has_magnetic_field(const Simulation3DConfig& cfg) {
    return cfg.magnetic_field.x != 0.0 || cfg.magnetic_field.y != 0.0 || cfg.magnetic_field.z != 0.0;
}

void initialize_particle_pusher(Particle3D& particle, Vec3 electric, double charge_to_mass, const Simulation3DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        initialize_boris_half_step(particle, electric, cfg.magnetic_field, charge_to_mass, cfg.dt);
    } else {
        initialize_leapfrog_half_step(particle, electric, charge_to_mass, cfg.dt);
    }
}

void kick_particle(Particle3D& particle, Vec3 electric, double charge_to_mass, const Simulation3DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        kick_boris(particle, electric, cfg.magnetic_field, charge_to_mass, cfg.dt);
    } else {
        kick_leapfrog(particle, electric, charge_to_mass, cfg.dt);
    }
}

void synchronize_particle(Particle3D& particle, Vec3 electric, double charge_to_mass, const Simulation3DConfig& cfg) {
    if (has_magnetic_field(cfg)) {
        synchronize_boris(particle, electric, cfg.magnetic_field, charge_to_mass, cfg.dt);
    } else {
        synchronize_leapfrog(particle, electric, charge_to_mass, cfg.dt);
    }
}

void write_vtk_outputs(const Mesh3D& mesh, const std::filesystem::path& output_dir, std::size_t step, VTKOutputFormat format) {
    const auto stem = output_dir / ("fields_" + std::to_string(step));
    if (format == VTKOutputFormat::Legacy || format == VTKOutputFormat::Both) {
        write_legacy_vtk(mesh, stem.string() + ".vtk");
    }
    if (format == VTKOutputFormat::Xml || format == VTKOutputFormat::Both) {
        write_vtk_xml(mesh, stem.string() + ".vts");
    }
}
} // namespace

Simulation3D::Simulation3D(Simulation3DConfig cfg)
    : cfg_(std::move(cfg)),
      mesh_(cfg_.nx, cfg_.ny, cfg_.nz, cfg_.length_x,
            cfg_.length_y, cfg_.length_z, cfg_.boundary),
      solver_(cfg_.units.permittivity()),
      rng_(cfg_.seed) {
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) cfg_.checkpoint_interval = cfg_.output_interval;
    if (!std::isfinite(cfg_.dt) || cfg_.dt <= 0.0) throw std::invalid_argument("3D simulation dt must be positive and finite");
    if (cfg_.output_interval == 0) throw std::invalid_argument("3D output_interval must be positive");
    if (cfg_.particle_output_stride == 0) throw std::invalid_argument("3D particle_output_stride must be positive");
    if (!std::isfinite(cfg_.magnetic_field.x) || !std::isfinite(cfg_.magnetic_field.y) || !std::isfinite(cfg_.magnetic_field.z)) {
        throw std::invalid_argument("3D magnetic_field components must be finite");
    }
    if (cfg_.mode == RunMode::SteadyState) {
        if (cfg_.max_steps == 0) throw std::invalid_argument("3D max_steps must be positive for steady-state mode");
        if (cfg_.steady_window == 0) throw std::invalid_argument("3D steady_window must be positive");
        if (!std::isfinite(cfg_.steady_tolerance) || cfg_.steady_tolerance <= 0.0) {
            throw std::invalid_argument("3D steady_tolerance must be positive and finite");
        }
    }
    validate_runtime_policy(cfg_.runtime);
    if (cfg_.checkpoint_output && cfg_.checkpoint_interval == 0) {
        throw std::invalid_argument("3D checkpoint_interval must be positive when checkpoint_output is enabled");
    }
    cfg_.particle_boundary_config.left = resolve_particle_boundary(cfg_.particle_boundary_config.left, cfg_.boundary);
    cfg_.particle_boundary_config.right = resolve_particle_boundary(cfg_.particle_boundary_config.right, cfg_.boundary);
    cfg_.particle_boundary_config.bottom = resolve_particle_boundary(cfg_.particle_boundary_config.bottom, cfg_.boundary);
    cfg_.particle_boundary_config.top = resolve_particle_boundary(cfg_.particle_boundary_config.top, cfg_.boundary);
    cfg_.particle_boundary_config.back = resolve_particle_boundary(cfg_.particle_boundary_config.back, cfg_.boundary);
    cfg_.particle_boundary_config.front = resolve_particle_boundary(cfg_.particle_boundary_config.front, cfg_.boundary);
    for (const auto& sc : cfg_.species) species_.emplace_back(sc);
    if (species_.empty()) species_.emplace_back(Species3DConfig{});
}

void Simulation3D::initialize() {
    time_ = 0.0;
    step_ = 0;
    boundary_losses_ = {};
    for (auto& sp : species_) sp.initialize(mesh_, rng_);
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

void Simulation3D::deposit_and_solve() {
    mesh_.clear_charge();
    for (const auto& sp : species_) sp.deposit_charge(mesh_);
    solver_.solve(mesh_);
}

void Simulation3D::apply_particle_boundaries(Particle3D& particle) {
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
    if (!apply_lower_boundary(particle.position.z, particle.velocity_half.z, mesh_.length_z(),
                              cfg_.particle_boundary_config.back, boundary_losses_.absorbed_back)) {
        particle.alive = false;
        return;
    }
    if (!apply_upper_boundary(particle.position.z, particle.velocity_half.z, mesh_.length_z(),
                              cfg_.particle_boundary_config.front, boundary_losses_.absorbed_front)) {
        particle.alive = false;
        return;
    }
}

void Simulation3D::step() {
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

void Simulation3D::save_checkpoint(const std::filesystem::path& path) const {
    ensure_parent_directory(path);
    std::ofstream out(path);
    if (!out) throw std::runtime_error("cannot open 3D checkpoint for writing: " + path.string());
    out << std::setprecision(17);
    out << kCheckpointMagicV2 << '\n';
    out << "dimension 3\n";
    out << "units " << to_string(cfg_.units.system) << ' '
        << cfg_.units.relative_permittivity << ' '
        << cfg_.units.permittivity() << "\n";
    out << "step " << step_ << "\n";
    out << "time " << time_ << "\n";
    out << "boundary_losses " << boundary_losses_.absorbed_left << ' ' << boundary_losses_.absorbed_right << ' '
        << boundary_losses_.absorbed_bottom << ' ' << boundary_losses_.absorbed_top << ' '
        << boundary_losses_.absorbed_back << ' ' << boundary_losses_.absorbed_front << "\n";
    out << "species_count " << species_.size() << "\n";
    out << "rng " << rng_ << "\n";
    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        const auto& sp = species_[species_id];
        out << "species " << species_id << ' ' << sp.name() << ' ' << sp.particles().size() << "\n";
        for (const auto& p : sp.particles()) {
            out << p.position.x << ' ' << p.position.y << ' ' << p.position.z << ' '
                << p.velocity.x << ' ' << p.velocity.y << ' ' << p.velocity.z << ' '
                << p.velocity_half.x << ' ' << p.velocity_half.y << ' ' << p.velocity_half.z << ' '
                << (p.alive ? 1 : 0) << "\n";
        }
    }
    require_stream(out, "failed while writing 3D checkpoint: " + path.string());
}

void Simulation3D::load_checkpoint(const std::filesystem::path& path) {
    std::ifstream in(path);
    if (!in) throw std::runtime_error("cannot open 3D checkpoint for reading: " + path.string());
    std::string magic;
    std::getline(in, magic);
    const bool checkpoint_v1 = magic == kCheckpointMagicV1;
    const bool checkpoint_v2 = magic == kCheckpointMagicV2;
    if (!checkpoint_v1 && !checkpoint_v2) {
        throw std::runtime_error("invalid checkpoint magic in: " + path.string());
    }

    std::string key;
    unsigned dimension = 0;
    in >> key >> dimension;
    if (key != "dimension" || dimension != 3) throw std::runtime_error("checkpoint dimension does not match 3D simulation");
    in >> key;
    if (key == "units") {
        std::string unit_system;
        double relative_permittivity = 0.0;
        double permittivity = 0.0;
        in >> unit_system >> relative_permittivity >> permittivity;
        if (!checkpoint_v2 ||
            unit_system != to_string(cfg_.units.system) ||
            relative_permittivity != cfg_.units.relative_permittivity ||
            permittivity != cfg_.units.permittivity()) {
            throw std::runtime_error(
                "checkpoint unit system does not match 3D config");
        }
        in >> key;
    } else if (checkpoint_v2 ||
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
       >> boundary_losses_.absorbed_bottom >> boundary_losses_.absorbed_top
       >> boundary_losses_.absorbed_back >> boundary_losses_.absorbed_front;
    if (key != "boundary_losses") throw std::runtime_error("checkpoint missing 3D boundary loss counters");
    std::size_t species_count = 0;
    in >> key >> species_count;
    if (key != "species_count" || species_count != species_.size()) throw std::runtime_error("checkpoint species count does not match 3D config");
    in >> key;
    if (key != "rng") throw std::runtime_error("checkpoint missing rng state");
    in >> rng_;

    for (std::size_t species_id = 0; species_id < species_.size(); ++species_id) {
        std::size_t stored_species_id = 0;
        std::string stored_name;
        std::size_t particle_count = 0;
        in >> key >> stored_species_id >> stored_name >> particle_count;
        if (key != "species" || stored_species_id != species_id || stored_name != species_[species_id].name()) {
            throw std::runtime_error("checkpoint species metadata does not match 3D config");
        }
        auto& particles = species_[species_id].particles();
        particles.resize(particle_count);
        for (auto& p : particles) {
            int alive = 0;
            in >> p.position.x >> p.position.y >> p.position.z
               >> p.velocity.x >> p.velocity.y >> p.velocity.z
               >> p.velocity_half.x >> p.velocity_half.y >> p.velocity_half.z
               >> alive;
            p.alive = alive != 0;
        }
    }
    require_stream(in, "failed while reading 3D checkpoint: " + path.string());
    deposit_and_solve();
    initialized_ = true;
}

RunSummary3D Simulation3D::run() {
    if (!cfg_.restart_path.empty()) load_checkpoint(cfg_.restart_path);
    else initialize();

    write_unit_metadata(cfg_.output_dir, cfg_.units, 3);
    std::vector<InitializationSpeciesMoments> initialization_moments;
    initialization_moments.reserve(species_.size());
    for (const auto& species : species_) {
        initialization_moments.push_back(
            summarize_initialization(species));
    }
    write_initialization_report(
        cfg_.output_dir / "initialization.csv", 3,
        cfg_.restart_path.empty() ? "generated" : "restart",
        initialization_moments);
    Diagnostics3D diag(
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
    RunSummary3D summary;
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

DiagnosticSample3D Simulation3D::sample() const {
    DiagnosticSample3D s;
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
    for (std::size_t k = 0; k < mesh_.nz(); ++k) {
        for (std::size_t j = 0; j < mesh_.ny(); ++j) {
            for (std::size_t i = 0; i < mesh_.nx(); ++i) {
                const auto idx = mesh_.index(i, j, k);
                const double e2 = mesh_.electric_x()[idx] * mesh_.electric_x()[idx]
                                + mesh_.electric_y()[idx] * mesh_.electric_y()[idx]
                                + mesh_.electric_z()[idx] * mesh_.electric_z()[idx];
                const double volume = mesh_.node_volume(i, j, k);
                s.field_energy +=
                    0.5 * cfg_.units.permittivity() * e2 * volume;
                s.charge_l1 += std::abs(mesh_.rho()[idx]) * volume;
            }
        }
    }
    s.total_energy = s.kinetic_energy + s.field_energy;
    return s;
}
}
