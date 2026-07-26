#include "pic/Pusher.hpp"

namespace pic {

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
}

void kick_leapfrog(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt) {
    particle.velocity_half.x += charge_to_mass * electric.x * dt;
    particle.velocity_half.y += charge_to_mass * electric.y * dt;
}

void drift_leapfrog(Particle2D& particle, double dt) {
    particle.position.x += particle.velocity_half.x * dt;
    particle.position.y += particle.velocity_half.y * dt;
}

void synchronize_leapfrog(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt) {
    particle.velocity.x = particle.velocity_half.x + 0.5 * charge_to_mass * electric.x * dt;
    particle.velocity.y = particle.velocity_half.y + 0.5 * charge_to_mass * electric.y * dt;
}

void initialize_leapfrog_half_step(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt) {
    particle.velocity_half.x = particle.velocity.x - 0.5 * charge_to_mass * electric.x * dt;
    particle.velocity_half.y = particle.velocity.y - 0.5 * charge_to_mass * electric.y * dt;
    particle.velocity_half.z = particle.velocity.z - 0.5 * charge_to_mass * electric.z * dt;
}

void kick_leapfrog(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt) {
    particle.velocity_half.x += charge_to_mass * electric.x * dt;
    particle.velocity_half.y += charge_to_mass * electric.y * dt;
    particle.velocity_half.z += charge_to_mass * electric.z * dt;
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

} // namespace pic
