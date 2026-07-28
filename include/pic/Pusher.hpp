#pragma once
#include "pic/Types.hpp"

namespace pic {
// Electrostatic leapfrog helper routines. Particle::v / Particle2D::velocity /
// Particle3D::velocity are retained as time-centered velocities for diagnostics
// and output, while the *_half fields hold the adjacent half-step velocity used
// for leapfrog kicks and position drift. Initialization stores the previous
// half-step velocity.
void initialize_leapfrog_half_step(Particle& particle, double electric, double charge_to_mass, double dt);
void kick_leapfrog(Particle& particle, double electric, double charge_to_mass, double dt);
void drift_leapfrog(Particle& particle, double dt);
void synchronize_leapfrog(Particle& particle, double electric, double charge_to_mass, double dt);
void initialize_leapfrog_half_step(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt);
void initialize_boris_half_step(Particle2D& particle, Vec2 electric, Vec3 magnetic, double charge_to_mass, double dt);
void initialize_boris_half_step(Particle2D& particle, Vec2 electric, double magnetic_z, double charge_to_mass, double dt);
void kick_leapfrog(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt);
void kick_boris(Particle2D& particle, Vec2 electric, Vec3 magnetic, double charge_to_mass, double dt);
void kick_boris(Particle2D& particle, Vec2 electric, double magnetic_z, double charge_to_mass, double dt);
void drift_leapfrog(Particle2D& particle, double dt);
void synchronize_leapfrog(Particle2D& particle, Vec2 electric, double charge_to_mass, double dt);
void synchronize_boris(Particle2D& particle, Vec2 electric, Vec3 magnetic, double charge_to_mass, double dt);
void synchronize_boris(Particle2D& particle, Vec2 electric, double magnetic_z, double charge_to_mass, double dt);

void initialize_leapfrog_half_step(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt);
void initialize_boris_half_step(Particle3D& particle, Vec3 electric, Vec3 magnetic, double charge_to_mass, double dt);
void kick_leapfrog(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt);
void kick_boris(Particle3D& particle, Vec3 electric, Vec3 magnetic, double charge_to_mass, double dt);
void drift_leapfrog(Particle3D& particle, double dt);
void synchronize_leapfrog(Particle3D& particle, Vec3 electric, double charge_to_mass, double dt);
void synchronize_boris(Particle3D& particle, Vec3 electric, Vec3 magnetic, double charge_to_mass, double dt);
}
