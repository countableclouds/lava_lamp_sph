use crate::utility::*;
use rayon::prelude::*;

use std::{cmp::Eq, collections::HashMap, convert::TryInto, hash::Hash};

pub struct Map<T: Coords + Copy> {
    pub particles: [Particle<T>; 500],
    pub dim: T,
    pub radius: f64,
    pub particle_map: HashMap<T::Key, Vec<Particle<T>>>,
    gravity: f64,
    gas_constant: f64,
}

impl<T> Map<T>
where
    T: Copy
        + Coords
        + Poly6Kernel
        + SpikyKernel
        + ViscosityKernel
        + std::ops::Mul<f64, Output = T>
        + std::ops::AddAssign
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + Sync
        + Send
        + Default,
    T::Key: Hash + Eq + Send + Sync,
{
    pub fn update_hashmap(&mut self) {
        self.particle_map.clear();
        for particle in self.particles.iter() {
            (*self
                .particle_map
                .entry(particle.position.bin(self.radius))
                .or_insert(Vec::new()))
            .push(particle.clone());
        }
    }

    pub fn update_density(&self, particle: &Particle<T>) -> Particle<T> {
        let mut density = 0.;

        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for other_particle in particles {
                    density += other_particle.properties.mass
                        * particle
                            .position
                            .poly6(self.radius, other_particle.position);
                }
            }
        }
        particle.with_density(density)
    }

    pub fn update_velocity(&self, particle: &Particle<T>, delta_t: f64) -> Particle<T> {
        let mut acceleration = T::default();

        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for other_particle in particles {
                    acceleration += particle
                        .position
                        .grad_spiky(self.radius, other_particle.position)
                        * (-other_particle.volume() / 2.
                            * self.gas_constant
                            * (other_particle.density_disparity() + particle.density_disparity()));

                    acceleration += (other_particle.velocity - particle.velocity)
                        * other_particle.volume()
                        * (other_particle.properties.viscosity + particle.properties.viscosity)
                        * (particle
                            .position
                            .laplace_visc(self.radius, other_particle.position)
                            / 2.);
                }
            }
        }
        acceleration =
            (acceleration * particle.density.recip() + T::height(-self.gravity)) * delta_t;
        particle.with_velocity(particle.velocity + acceleration)
    }

    pub fn update(&mut self, delta_t: f64) {
        self.particles = match self
            .particles
            .par_iter()
            .map(|particle| self.update_density(particle))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
        {
            Ok(arr) => arr,
            Err(_) => panic!("Expected a Vec of a different length",),
        };
        self.particles = match self
            .particles
            .par_iter()
            .map(|particle| self.update_velocity(particle, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
        {
            Ok(arr) => arr,
            Err(_) => panic!("Expected a Vec of a different length",),
        };
    }
}
