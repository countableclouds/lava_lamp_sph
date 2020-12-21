use crate::utility::*;
use crate::NUM_PARTICLES;
use rayon::prelude::*;
use std::{cmp::Eq, collections::HashMap, convert::TryInto, hash::Hash};

pub struct Map<T: Coords + Copy> {
    pub particles: [Particle<T>; NUM_PARTICLES],
    pub dim: T,
    pub radius: f64,
    pub particle_map: HashMap<T::Key, Vec<(usize, Particle<T>)>>,
    pub gravity: f64,
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
        for (i, particle) in self.particles.iter().enumerate() {
            (*self
                .particle_map
                .entry(particle.position.bin(self.radius))
                .or_insert(Vec::new()))
            .push((i, particle.clone()));
        }
    }

    pub fn update_density(&self, particle: &Particle<T>) -> Particle<T> {
        let mut density = 0.;
        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for (_, other_particle) in particles {
                    density += other_particle.mass
                        * particle
                            .position
                            .poly6(self.radius, other_particle.position);
                }
            }
        }
        particle.with_density(density)
    }

    pub fn update_nonpressure_velocity(&self, particle: &Particle<T>, delta_t: f64) -> Particle<T> {
        let mut acceleration = T::default();
        let mut water_benzyl_normal = T::default();
        let mut water_benzyl_curvature = 0.;
        let mut temperature = 0.;
        // println!("{:?}", particle.density_disparity());
        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for (_, other_particle) in particles {

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

        acceleration = acceleration * particle.density.recip() + T::height(self.gravity);
        temperature = particle.temperature + temperature * delta_t;
        particle
            .with_velocity(particle.velocity + acceleration * delta_t)
            .with_temperature(temperature)
    }
    pub fn get_diagonal(&self, delta_t: f64) -> [f64; NUM_PARTICLES] {
        let mut inverse_densities: [T; NUM_PARTICLES] = [T::default(); NUM_PARTICLES];
        let mut diagonal: [f64; NUM_PARTICLES] = [0.; NUM_PARTICLES];
        for (i, particle) in self.particles.iter().enumerate() {
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    for (_, other_particle) in particles {
                        inverse_densities[i] += particle
                            .position
                            .grad_poly6(self.radius, other_particle.position)
                            * (other_particle.volume() / other_particle.density);
                    }
                }
            }
        }
        for (i, particle) in self.particles.iter().enumerate() {
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    for (_, other_particle) in particles {
                        diagonal[i] -= particle
                            .position
                            .grad_poly6(self.radius, other_particle.position)
                            .dot(
                                inverse_densities[i]
                                    + particle
                                        .position
                                        .grad_poly6(self.radius, other_particle.position)
                                        * (particle.mass * particle.density.powi(-2)),
                            )
                            * (other_particle.mass * delta_t.powi(2));
                    }
                }
            }
        }
        diagonal
    }

    pub fn get_density_differences(&self, delta_t: f64) -> [f64; NUM_PARTICLES] {
        let mut density_diff: [f64; NUM_PARTICLES] = [0.; NUM_PARTICLES];
        for (i, particle) in self.particles.iter().enumerate() {
            density_diff[i] = -particle.delta_density();
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    for (_, other_particle) in particles {
                        density_diff[i] -= particle
                            .position
                            .grad_visc(self.radius, other_particle.position)
                            .dot(particle.velocity - other_particle.velocity)
                            * (other_particle.mass * delta_t);
                    }
                }
            }
        }
        density_diff
    }

    pub fn get_pressure_accelerations(
        &self,
        pressures: [f64; NUM_PARTICLES],
    ) -> [T; NUM_PARTICLES] {
        let mut accelerations: [T; NUM_PARTICLES] = [T::default(); NUM_PARTICLES];
        for (i, particle) in self.particles.iter().enumerate() {
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    for (j, other_particle) in particles {
                        accelerations[i] += particle
                            .position
                            .grad_spiky(self.radius, other_particle.position)
                            * (-other_particle.mass
                                * (pressures[i] * particle.density.powi(-2)
                                    + pressures[*j] * other_particle.density.powi(-2)));
                    }
                }
            }
        }
        accelerations
    }

    pub fn update_pressure_velocity(&mut self, delta_t: f64) {
        let num_iter: u64 = 100;
        let relaxation_coeff = 0.5;
        let mut pressures: [f64; NUM_PARTICLES] = [0.; NUM_PARTICLES];

        let mut pressure_accelerations: [T; NUM_PARTICLES] = [T::default(); NUM_PARTICLES];
        let diagonal = self.get_diagonal(delta_t);
        let density_diff = self.get_density_differences(delta_t);
        let mut error;
        for _ in 0..num_iter {
            let mut image: [f64; NUM_PARTICLES] = [0.; NUM_PARTICLES];
            error = 0.;
            for (i, particle) in self.particles.iter().enumerate() {
                for point in particle.position.along_axes(self.radius) {
                    if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                        for (j, other_particle) in particles {
                            image[i] += particle
                                .position
                                .grad_visc(self.radius, other_particle.position)
                                .dot(pressure_accelerations[i] - pressure_accelerations[*j])
                                * (other_particle.mass * delta_t.powi(2));
                        }
                    }
                }
                pressures[i] += relaxation_coeff * (density_diff[i] - image[i]) / diagonal[i];
                pressures[i] = pressures[i].max(0.)
            }
            pressure_accelerations = self.get_pressure_accelerations(pressures);

            // for i in 0..image.len() {
            //     error += (density_diff[i] - image[i]).abs();
            // }
            // println!("{:?}", error);
            // println!(
            //     "count: {}",
            //     pressures.len() - pressures.iter().filter(|&n| *n == 0.).count()
            // )
        }

        self.particles = self
            .particles
            .iter()
            .enumerate()
            .map(|(i, particle)| {
                particle.with_velocity(particle.velocity + pressure_accelerations[i] * delta_t)
            })
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");
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
            .map(|particle| self.update_nonpressure_velocity(particle, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");

        self.update_pressure_velocity(delta_t);
        self.particles = self
            .particles
            .iter()
            .map(|particle| particle.control_update(self.dim, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");
    }
}
