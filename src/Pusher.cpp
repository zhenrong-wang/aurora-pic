#include "pic/Pusher.hpp"

namespace pic {
namespace {

Vec3 add(Vec3 a, Vec3 b) {
    return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}

Vec3 scale(Vec3 v, double factor) {
    return Vec3{factor * v.x, factor * v.y, factor * v.z};
}

Vec3 cross(Vec3 a, Vec3 b) {
    return Vec3{a.y * b.z - a.z * b.y,
                a.z * b.x - a.x * b.z,
                a.x * b.y - a.y * b.x};
}

Vec3 boris_advance(Vec3 velocity, Vec3 electric, Vec3 magnetic, double charge_to_mass, double dt) {
    const double half_qm_dt = 0.5 * charge_to_mass * dt;
    const Vec3 v_minus = add(velocity, scale(electric, half_qm_dt));
    const Vec3 t = scale(magnetic, half_qm_dt);
    const double t2 = t.x * t.x + t.y * t.y + t.z * t.z;
    const Vec3 s = scale(t, 2.0 / (1.0 + t2));
    const Vec3 v_prime = add(v_minus, cross(v_minus, t));
    const Vec3 v_plus = add(v_minus, cross(v_prime, s));
    return add(v_plus, scale(electric, half_qm_dt));
}

Vec3 to_vec3(Vec2 v) {
    return Vec3{v.x, v.y, 0.0};
}

Vec2 to_vec2(Vec3 v) {
    return Vec2{v.x, v.y};
}

Vec3 particle_velocity(const Particle2D& particle) {
    return {
        particle.velocity.x,
        particle.velocity.y,
        particle.velocity_z};
}

Vec3 particle_velocity_half(const Particle2D& particle) {
    return {
        particle.velocity_half.x,
        particle.velocity_half.y,
        particle.velocity_half_z};
}

void set_particle_velocity(Particle2D& particle, Vec3 velocity) {
    particle.velocity = to_vec2(velocity);
    particle.velocity_z = velocity.z;
}

void set_particle_velocity_half(
    Particle2D& particle, Vec3 velocity) {
    particle.velocity_half = to_vec2(velocity);
    particle.velocity_half_z = velocity.z;
}

} // namespace

void initialize_leapfrog_half_step(Particle& particle, double electric, double charge_to_mass, double dt) {
    particle.v_half = particle.v - 0.5 * charge_to_mass * electric * dt;
}

void kick_leapfrog(Particle& particle, double electric, double charge_to_mass, double dt) {
    particle.v_half += charge_to_mass * electric * dt;
}

void drift_leapfrog(Particle& particle, double dt) {
    particle.x += particle.v_half * dt;
}

void synchronize_leapfrog(Particle& particle, double electric, double charge_to_mass, double dt) {
    particle.v = particle.v_half + 0.5 * charge_to_mass * electric * dt;
}

void initialize_leapfrog_half_step(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt) {
    particle.velocity_half.x = particle.velocity.x - 0.5 * charge_to_mass * electric.x * dt;
    particle.velocity_half.y = particle.velocity.y - 0.5 * charge_to_mass * electric.y * dt;
    particle.velocity_half_z = particle.velocity_z;
}

void initialize_boris_half_step(
    Particle2D& particle, Vec2 electric, Vec3 magnetic,
    double charge_to_mass, double dt) {
    set_particle_velocity_half(
        particle,
        boris_advance(
            particle_velocity(particle), to_vec3(electric),
            magnetic, charge_to_mass, -0.5 * dt));
}

void initialize_boris_half_step(Particle2D& particle, Vec2 electric, double magnetic_z, double charge_to_mass, double dt) {
    initialize_boris_half_step(
        particle, electric, Vec3{0.0, 0.0, magnetic_z},
        charge_to_mass, dt);
}

void kick_leapfrog(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt) {
    particle.velocity_half.x += charge_to_mass * electric.x * dt;
    particle.velocity_half.y += charge_to_mass * electric.y * dt;
}

void kick_boris(
    Particle2D& particle, Vec2 electric, Vec3 magnetic,
    double charge_to_mass, double dt) {
    set_particle_velocity_half(
        particle,
        boris_advance(
            particle_velocity_half(particle), to_vec3(electric),
            magnetic, charge_to_mass, dt));
}

void kick_boris(Particle2D& particle, Vec2 electric, double magnetic_z, double charge_to_mass, double dt) {
    kick_boris(
        particle, electric, Vec3{0.0, 0.0, magnetic_z},
        charge_to_mass, dt);
}

void drift_leapfrog(Particle2D& particle, double dt) {
    particle.position.x += particle.velocity_half.x * dt;
    particle.position.y += particle.velocity_half.y * dt;
}

void synchronize_leapfrog(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt) {
    particle.velocity.x = particle.velocity_half.x + 0.5 * charge_to_mass * electric.x * dt;
    particle.velocity.y = particle.velocity_half.y + 0.5 * charge_to_mass * electric.y * dt;
    particle.velocity_z = particle.velocity_half_z;
}

void synchronize_boris(
    Particle2D& particle, Vec2 electric, Vec3 magnetic,
    double charge_to_mass, double dt) {
    set_particle_velocity(
        particle,
        boris_advance(
            particle_velocity_half(particle), to_vec3(electric),
            magnetic, charge_to_mass, 0.5 * dt));
}

void synchronize_boris(Particle2D& particle, Vec2 electric, double magnetic_z, double charge_to_mass, double dt) {
    synchronize_boris(
        particle, electric, Vec3{0.0, 0.0, magnetic_z},
        charge_to_mass, dt);
}

void initialize_leapfrog_half_step(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt) {
    particle.velocity_half.x = particle.velocity.x - 0.5 * charge_to_mass * electric.x * dt;
    particle.velocity_half.y = particle.velocity.y - 0.5 * charge_to_mass * electric.y * dt;
    particle.velocity_half.z = particle.velocity.z - 0.5 * charge_to_mass * electric.z * dt;
}

void initialize_boris_half_step(Particle3D& particle, Vec3 electric, Vec3 magnetic, double charge_to_mass, double dt) {
    particle.velocity_half = boris_advance(particle.velocity, electric, magnetic, charge_to_mass, -0.5 * dt);
}

void kick_leapfrog(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt) {
    particle.velocity_half.x += charge_to_mass * electric.x * dt;
    particle.velocity_half.y += charge_to_mass * electric.y * dt;
    particle.velocity_half.z += charge_to_mass * electric.z * dt;
}

void kick_boris(Particle3D& particle, Vec3 electric, Vec3 magnetic, double charge_to_mass, double dt) {
    particle.velocity_half = boris_advance(particle.velocity_half, electric, magnetic, charge_to_mass, dt);
}

void drift_leapfrog(Particle3D& particle, double dt) {
    particle.position.x += particle.velocity_half.x * dt;
    particle.position.y += particle.velocity_half.y * dt;
    particle.position.z += particle.velocity_half.z * dt;
}

void synchronize_leapfrog(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt) {
    particle.velocity.x = particle.velocity_half.x + 0.5 * charge_to_mass * electric.x * dt;
    particle.velocity.y = particle.velocity_half.y + 0.5 * charge_to_mass * electric.y * dt;
    particle.velocity.z = particle.velocity_half.z + 0.5 * charge_to_mass * electric.z * dt;
}

void synchronize_boris(Particle3D& particle, Vec3 electric, Vec3 magnetic, double charge_to_mass, double dt) {
    particle.velocity = boris_advance(particle.velocity_half, electric, magnetic, charge_to_mass, 0.5 * dt);
}

} // namespace pic
