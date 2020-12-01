use crate::utility::*;
use rayon::prelude::*;

use std::{cmp::Eq, collections::HashMap, convert::TryInto, hash::Hash};

pub struct Map<T: Coords + Copy> {
    pub particles: [Particle<T>; 10693],
    pub dim: T,
    pub radius: f64,
    pub particle_map: HashMap<T::Key, Vec<Particle<T>>>,
    pub gravity: f64,
    pub gas_constant: f64,
}

impl<T> Map<T>
where
    T: Copy
        + Coords
        + PartialEq
        + Poly6Kernel
        + SpikyKernel
        + ViscosityKernel
        + std::ops::Mul<f64, Output = T>
        + std::ops::AddAssign
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + Sync
        + Send
        + Default
        + std::fmt::Debug,
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
                    // if (other_particle.position - particle.position).mag() < 0.000000001
                    //     && other_particle != particle
                    // {
                    //     println!("{}", particle
                    //     .position
                    //     .poly6(self.radius, other_particle.position));
                    // }
                    density += other_particle.mass
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
        let mut water_benzyl_normal = T::default();
        let mut water_benzyl_curvature = 0.;
        let mut temperature = 0.;
        // println!("{:?}", particle.density_disparity());
        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for other_particle in particles {
                    acceleration += particle
                        .position
                        .grad_spiky(self.radius, other_particle.position)
                        * (-other_particle.volume() / 2.
                            * self.gas_constant
                            * (other_particle.gas_coefficient() + particle.gas_coefficient()));

                    // acceleration += (other_particle.velocity - particle.velocity)
                    //     * other_particle.volume()
                    //     * (other_particle.fluid_type.viscosity() + particle.fluid_type.viscosity())
                    //     * (particle
                    //         .position
                    //         .laplace_visc(self.radius, other_particle.position)
                    //         / 2.);

                    // water_benzyl_curvature += other_particle.volume()
                    //     * particle.fluid_type.color(other_particle.fluid_type)
                    //     * particle
                    //         .position
                    //         .laplace_poly6(self.radius, other_particle.position);

                    // water_benzyl_normal += particle
                    //     .position
                    //     .grad_poly6(self.radius, other_particle.position)
                    //     * (other_particle.volume()
                    //         * particle.fluid_type.color(other_particle.fluid_type));

                    // temperature += particle.fluid_type.diffusivity()
                    //     * (other_particle.temperature - particle.temperature)
                    //     * other_particle.volume()
                    //     * particle
                    //         .position
                    //         .laplace_visc(self.radius, other_particle.position);
                }
            }
        }
        // if water_benzyl_normal != T::default() {
        //     acceleration += water_benzyl_normal.normalize()
        //         * (Fluid::Saltwater.interfacial_tension(Fluid::BenzylAlcohol)
        //             * water_benzyl_curvature);
        // }

        acceleration =
            (acceleration * particle.density.recip() + T::height(self.gravity)) * delta_t;
        temperature = particle.temperature + temperature * delta_t;
        particle
            .with_velocity(particle.velocity + acceleration)
            .with_temperature(temperature)
            .control_update(self.dim, delta_t)
    }

    pub fn update(&mut self, delta_t: f64) {
        self.update_hashmap();

        self.particles = self
            .particles
            .par_iter()
            .map(|particle| self.update_density(particle))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");

        self.update_hashmap();

        self.particles = self
            .particles
            .par_iter()
            .map(|particle| self.update_velocity(particle, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length")
    }
}
