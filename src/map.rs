use crate::utility::*;
use crate::{NUM_BOUNDARY_PARTICLES, NUM_PARTICLES};
use rayon::prelude::*;
use std::fs;
use std::{cmp::Eq, collections::HashMap, convert::TryInto, fmt::Display, hash::Hash};
pub struct Map<T: Coords + Copy> {
    pub particles: [Particle<T>; NUM_PARTICLES],
    pub boundary_particles: [BoundaryParticle<T>; NUM_BOUNDARY_PARTICLES],
    pub dim: T,
    pub radius: f64,
    pub particle_map: HashMap<T::Key, Vec<(usize, Particle<T>)>>,
    pub boundary_particle_map: HashMap<T::Key, Vec<BoundaryParticle<T>>>,
    pub gravity: f64,
    pub pressures: [f64; NUM_PARTICLES],
    pub cfl: f64,
}

pub const TEST_NUM: usize = 0; //2223; // ;
pub const VOLUME_FACTOR: f64 = 1.;

impl<T> Map<T>
where
    T: Display
        + Copy
        + Coords
        + PartialEq
        + Poly6Kernel
        + SpikyKernel
        + ViscosityKernel
        + CubicSplineKernel
        + QuinticSplineKernel
        + std::ops::Mul<f64, Output = T>
        + std::ops::AddAssign
        + std::ops::Add<Output = T>
        + std::ops::Sub<Output = T>
        + Sync
        + Send
        + Default
        + std::fmt::Display,
    T::Key: Hash + Eq + Send + Sync,
{
    pub fn new(
        particles: [Particle<T>; NUM_PARTICLES],
        boundary_particles: [BoundaryParticle<T>; NUM_BOUNDARY_PARTICLES],
        radius: f64,
        dim: T,
        gravity: f64,
        cfl: f64,
    ) -> Self {
        Map {
            particles,
            boundary_particles,
            dim,
            radius,
            particle_map: HashMap::new(),
            boundary_particle_map: HashMap::new(),
            gravity,
            pressures: [0.; NUM_PARTICLES],
            cfl,
        }
    }
    pub fn update_hashmap(&mut self) {
        self.particle_map.clear();
        for (i, particle) in self.particles.iter().enumerate() {
            (*self
                .particle_map
                .entry(particle.position.bin(self.radius))
                .or_insert(Vec::new()))
            .push((i, particle.clone()));
        }
        self.boundary_particle_map.clear();
        for boundary_particle in self.boundary_particles.iter() {
            (*self
                .boundary_particle_map
                .entry(boundary_particle.position.bin(self.radius))
                .or_insert(Vec::new()))
            .push(boundary_particle.clone());
        }
    }

    pub fn update_density(&self, i: usize, particle: &Particle<T>) -> Particle<T> {
        let mut density = 0.;
        let mut point_count = 0;

        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for (_, fluid_particle) in particles {
                    if fluid_particle.mass
                        * particle
                            .position
                            .cubic(self.radius, fluid_particle.position)
                        > 0.
                    {
                        point_count += 1
                    }
                    density += fluid_particle.mass
                        * particle
                            .position
                            .cubic(self.radius, fluid_particle.position);
                }
            }

            if let Some(boundary_particles) =
                self.boundary_particle_map.get(&point.bin(self.radius))
            {
                for boundary_particle in boundary_particles {
                    if particle.mass
                        * particle
                            .position
                            .cubic(self.radius, boundary_particle.position)
                        > 0.
                    {
                        point_count += 1
                    }
                    density += boundary_particle.mass(particle.rest_density())
                        * particle
                            .position
                            .cubic(self.radius, boundary_particle.position);
                }
            }
        }

        if i == TEST_NUM {
            println!("Point count: {}", point_count);
        }

        particle.with_density(density)
    }

    pub fn update_boundary_volume(&self, particle: &BoundaryParticle<T>) -> BoundaryParticle<T> {
        let mut volume = 0.;
        for point in particle.position.along_axes(self.radius) {
            if let Some(boundary_particles) =
                self.boundary_particle_map.get(&point.bin(self.radius))
            {
                for boundary_particle in boundary_particles {
                    volume += VOLUME_FACTOR
                        * particle
                            .position
                            .cubic(self.radius, boundary_particle.position);
                }
            }
        }

        particle.with_volume(volume.recip())
    }

    pub fn update_nonpressure_velocity(
        &self,
        i: usize,
        particle: &Particle<T>,
        delta_t: f64,
    ) -> Particle<T> {
        let mut acceleration = T::default();
        let mut water_benzyl_normal = T::default();
        let mut water_benzyl_curvature = 0.;
        let mut temperature = 0.;
        // if (self.dim.height() - particle.position.height()) < 0.05 {
        //     factor = 100.;
        // }
        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for (_, fluid_particle) in particles {
                    acceleration += (fluid_particle.velocity - particle.velocity)
                        * (fluid_particle.volume()
                            * (fluid_particle.fluid_type.viscosity()
                                + particle.fluid_type.viscosity())
                            * (particle
                                .position
                                .laplace_cubic(self.radius, fluid_particle.position)
                                / 2.));

                    water_benzyl_curvature += fluid_particle.volume()
                        * particle.fluid_type.color(fluid_particle.fluid_type)
                        * particle
                            .position
                            .laplace_cubic(self.radius, fluid_particle.position);

                    water_benzyl_normal += particle
                        .position
                        .grad_cubic(self.radius, fluid_particle.position)
                        * (-fluid_particle.volume()
                            * particle.fluid_type.color(fluid_particle.fluid_type));

                    temperature += if particle.fluid_type == fluid_particle.fluid_type {
                        particle.fluid_type.diffusivity()
                    } else {
                        (particle.fluid_type.diffusivity()
                            + fluid_particle.fluid_type.diffusivity())
                            / 2.
                    } * (fluid_particle.temperature - particle.temperature)
                        * fluid_particle.volume()
                        * particle
                            .position
                            .laplace_cubic(self.radius, fluid_particle.position);
                }
            }
            if let Some(boundary_particles) =
                self.boundary_particle_map.get(&point.bin(self.radius))
            {
                for boundary_particle in boundary_particles {
                    temperature += particle.fluid_type.diffusivity()
                        * (boundary_particle.temperature - particle.temperature)
                        * boundary_particle.volume
                        * particle
                            .position
                            .laplace_cubic(self.radius, boundary_particle.position);
                }
            }
        }
        // if i == TEST_NUM {
        //     println!("{}", buoyancy);
        // }

        // if water_benzyl_normal != T::default() {
        //     acceleration += water_benzyl_normal.normalize()
        //         * (Fluid::Saltwater.interfacial_tension(Fluid::BenzylAlcohol)
        //             * water_benzyl_curvature);
        // }

        acceleration =
            acceleration * particle.density.recip() + T::default().with_height(self.gravity);
        // temperature = particle.temperature + temperature * delta_t;
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
                    for (_, fluid_particle) in particles {
                        inverse_densities[i] += particle
                            .position
                            .grad_cubic(self.radius, fluid_particle.position)
                            * (fluid_particle.mass * particle.density.powi(-2));
                    }
                }
                if let Some(boundary_particles) =
                    self.boundary_particle_map.get(&point.bin(self.radius))
                {
                    for boundary_particle in boundary_particles {
                        inverse_densities[i] += particle
                            .position
                            .grad_cubic(self.radius, boundary_particle.position)
                            * (boundary_particle.mass(particle.rest_density())
                                * particle.density.powi(-2));
                    }
                }
            }
        }
        for (i, particle) in self.particles.iter().enumerate() {
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    for (_, fluid_particle) in particles {
                        diagonal[i] -= particle
                            .position
                            .grad_cubic(self.radius, fluid_particle.position)
                            .dot(
                                inverse_densities[i]
                                    + particle
                                        .position
                                        .grad_cubic(self.radius, fluid_particle.position)
                                        * (particle.mass * particle.density.powi(-2)),
                            )
                            * fluid_particle.mass
                            * delta_t.powi(2);
                    }
                }

                if let Some(boundary_particles) =
                    self.boundary_particle_map.get(&point.bin(self.radius))
                {
                    for boundary_particle in boundary_particles {
                        diagonal[i] -= particle
                            .position
                            .grad_cubic(self.radius, boundary_particle.position)
                            .dot(inverse_densities[i])
                            * boundary_particle.mass(particle.rest_density())
                            * delta_t.powi(2);
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
            if i == TEST_NUM {
                println!("Boundary density diff thing: {}", density_diff[i]);
            }
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    for (_, fluid_particle) in particles {
                        density_diff[i] -= particle
                            .position
                            .grad_cubic(self.radius, fluid_particle.position)
                            .dot(particle.velocity - fluid_particle.velocity)
                            * fluid_particle.mass
                            * delta_t;
                    }
                }
                if let Some(boundary_particles) =
                    self.boundary_particle_map.get(&point.bin(self.radius))
                {
                    for boundary_particle in boundary_particles {
                        density_diff[i] -= particle
                            .position
                            .grad_cubic(self.radius, boundary_particle.position)
                            .dot(particle.velocity - boundary_particle.velocity)
                            * boundary_particle.mass(particle.rest_density())
                            * delta_t;
                    }
                }
            }
        }
        density_diff
    }

    pub fn get_pressure_accelerations(&self) -> [T; NUM_PARTICLES] {
        let mut accelerations: [T; NUM_PARTICLES] = [T::default(); NUM_PARTICLES];
        for (i, particle) in self.particles.iter().enumerate() {
            if particle.density == 0. {
                continue;
            }
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    if i == TEST_NUM {
                    } else {
                        for (j, fluid_particle) in particles {
                            if fluid_particle.density < 0.000000001 {
                                continue;
                            }
                            accelerations[i] += particle
                                .position
                                .grad_cubic(self.radius, fluid_particle.position)
                                * (-fluid_particle.mass
                                    * (self.pressures[i] * particle.density.powi(-2)
                                        + self.pressures[*j] * fluid_particle.density.powi(-2)));
                        }
                    }
                }
                if let Some(boundary_particles) =
                    self.boundary_particle_map.get(&point.bin(self.radius))
                {
                    for boundary_particle in boundary_particles {
                        accelerations[i] += particle
                            .position
                            .grad_cubic(self.radius, boundary_particle.position)
                            * (-boundary_particle.mass(particle.rest_density())
                                * self.pressures[i]
                                * particle.density.powi(-2));
                    }
                }
            }
        }

        accelerations
    }

    pub fn update_pressure_velocity(&mut self, delta_t: f64) {
        let num_iter: u64 = 20;
        let relaxation_coeff = 0.5;
        let mut pressure_accelerations: [T; NUM_PARTICLES] = [T::default(); NUM_PARTICLES];
        let diagonal = self.get_diagonal(delta_t);
        let density_diff = self.get_density_differences(delta_t);
        let mut error;
        let mut h: Vec<f64> = density_diff.iter().map(|e| e.abs()).collect();
        h.sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "Quantiles: {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8}",
            h[0],
            h[1 * h.len() / 10],
            h[2 * h.len() / 10],
            h[3 * h.len() / 10],
            h[4 * h.len() / 10],
            h[5 * h.len() / 10],
            h[6 * h.len() / 10],
            h[7 * h.len() / 10],
            h[8 * h.len() / 10],
            h[9 * h.len() / 10],
            h[h.len() - 1],
        );
        let mut h: Vec<(usize, f64)> = self
            .particles
            .iter()
            .enumerate()
            .map(|(i, p)| (i, p.density))
            .collect();
        h.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        println!(
            "Quantiles: {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8} {:.8}",
            h[0].1,
            h[1 * h.len() / 10].1,
            h[2 * h.len() / 10].1,
            h[3 * h.len() / 10].1,
            h[4 * h.len() / 10].1,
            h[5 * h.len() / 10].1,
            h[6 * h.len() / 10].1,
            h[7 * h.len() / 10].1,
            h[8 * h.len() / 10].1,
            h[9 * h.len() / 10].1,
            h[h.len() - 1].1,
        );
        println!(
            "mean: {}",
            h.iter().map(|(_, j)| j).sum::<f64>() / h.len() as f64
        );
        println!(
            "Location of min particle: {}, index: {}",
            self.particles[h[0].0].position, h[0].0
        );
        println!(
            "Location of max particle: {}, index: {}",
            self.particles[h[h.len() - 1].0].position,
            h[h.len() - 1].0
        );
        println!(
            "Offset (positive means we want denser): {:.8}",
            density_diff.iter().map(|e| *e).sum::<f64>() / density_diff.len() as f64
        );

        println!(
            "Proportion of too light particles: {}",
            density_diff
                .iter()
                .map(|e| (e > &0.) as i64 as f64)
                .sum::<f64>()
                / density_diff.len() as f64
        );
        println!(
            "Diagonal sum average: {:?}",
            diagonal.iter().map(|e| e).sum::<f64>() / density_diff.len() as f64
        );
        for j in 0..num_iter {
            pressure_accelerations = self.get_pressure_accelerations();
            let mut image: [f64; NUM_PARTICLES] = [0.; NUM_PARTICLES];
            error = 0.;
            let mut count = 0;
            for (i, particle) in self.particles.iter().enumerate() {
                if diagonal[i] == 0. {
                    continue;
                }
                for point in particle.position.along_axes(self.radius) {
                    if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                        for (j, fluid_particle) in particles {
                            image[i] += particle
                                .position
                                .grad_cubic(self.radius, fluid_particle.position)
                                .dot(
                                    pressure_accelerations[i]
                                        - if *j < NUM_PARTICLES {
                                            pressure_accelerations[*j]
                                        } else {
                                            T::default()
                                        },
                                )
                                * fluid_particle.mass
                                * delta_t.powi(2);
                        }
                    }

                    if let Some(boundary_particles) =
                        self.boundary_particle_map.get(&point.bin(self.radius))
                    {
                        for boundary_particle in boundary_particles {
                            image[i] += particle
                                .position
                                .grad_cubic(self.radius, boundary_particle.position)
                                .dot(pressure_accelerations[i])
                                * boundary_particle.mass(particle.rest_density())
                                * delta_t.powi(2);
                        }
                    }
                }
                // if i == TEST_NUM && CHECKER {
                //     println!(
                //         "Test Particles Velocity: {}",
                //         self.particles[TEST_NUM].velocity
                //     );
                //     println!(
                //         "Test Particles Pressure Acceleration: {}",
                //         pressure_accelerations[TEST_NUM]
                //     );
                //     println!("Test Particles Pressure : {}", self.pressures[TEST_NUM]);
                // }

                self.pressures[i] +=
                    relaxation_coeff * (density_diff[i] - image[i]) / (diagonal[i]);

                if i == TEST_NUM {
                    println!(
                        "Density Difference: {}, Image: {}, Diagonal: {}",
                        density_diff[i], image[i], diagonal[i]
                    );
                }

                if self.pressures[i] < 0. {
                    count += 1;
                }

                // pressures[i] = (0.29845 - self.particles[i].position.height()) * 100.;
                self.pressures[i] = self.pressures[i].max(0.);

                // if i == TEST_NUM {
                //     println!("Pressure: {}", pressures[i]);
                //     println!("Image: {}", image[i]);
                //     println!("Height: {}", self.particles[i].position.height());
                // }
                // if ERROR_CHECKER {
                //     if i == TEST_NUM {
                //         self.pressures[i] = 1.;
                //     } else {
                //         self.pressures[i] = 1.;
                //     }
                // }
                // if i != TEST_NUM {
                //     pressures[i] = 0.;
                // }
            }

            println!("TEST PRESSURE!: {}", self.pressures[TEST_NUM]);
            println!(
                "TEST PRESSURE ACCELERATIONS!: {}",
                pressure_accelerations[TEST_NUM]
            );
            // let max_pressure = pressures.iter().cloned().fold(0. / 0., f64::max);
            // for elem in &mut pressures {
            //     *elem /= max_pressure / 10.;
            // }

            // println!(
            //     "New velocity: {:?}",
            //     pressure_accelerations[TEST_NUM] * delta_t
            // );
            for i in 0..image.len() {
                error += (density_diff[i] - image[i]).abs() / self.particles[i].rest_density();
            }
            let mut image_average = 0.;
            for i in 0..image.len() {
                image_average += image[i].abs() / self.particles[i].rest_density();
            }
            // if ERROR_CHECKER {
            //     println!(
            //         "Density Differences: {:?}, {:?}, {:?}",
            //         image[TEST_NUM] + self.particles[TEST_NUM].density,
            //         self.pressures[TEST_NUM],
            //         image[TEST_NUM]
            //     );
            // }

            if j == 0 || j == num_iter - 1 || true {
                println!("ERROR Percentage: {:?}, {}", error / image.len() as f64, j);
            }

            // println!("Image Average: {:?}", image_average / image.len() as f64);

            // if j == 1 {
            //     break;
            // }
        }
        println!(
            "count: {}, pressure average: {}",
            self.pressures.len() - self.pressures.iter().filter(|&n| *n <= 0.).count(),
            self.pressures.iter().fold(0., |a, b| a + b) / self.pressures.len() as f64
        );
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
        println!("Velocity: {}", self.particles[TEST_NUM].velocity);
        for pressure in self.pressures.iter_mut() {
            *pressure *= 0.6;
        }
    }

    pub fn to_file(&mut self, path: String) {
        let mut data: String = "".to_owned();
        for &particle in &self.particles {
            data.push_str(&particle.position.to_string());
            data.push_str(", ");
        }
        data = data[..data.len() - 1].to_owned();

        fs::write(path, data).expect("Unable to write file");
    }

    pub fn update(&mut self, delta_t: f64) -> f64 {
        self.update_hashmap();

        self.boundary_particles = self
            .boundary_particles
            .par_iter()
            .map(|boundary_particle| self.update_boundary_volume(boundary_particle))
            .collect::<Vec<BoundaryParticle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");
        self.update_hashmap();
        self.particles = self
            .particles
            .par_iter()
            .enumerate()
            .map(|(i, particle)| self.update_density(i, particle))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");

        println!("Resultant density: {:?}", self.particles[TEST_NUM].density);
        println!("Particle position: {}", self.particles[TEST_NUM].position);

        self.update_hashmap();

        self.particles = self
            .particles
            .par_iter()
            .enumerate()
            .map(|(i, particle)| self.update_nonpressure_velocity(i, particle, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");

        self.update_hashmap();

        self.update_pressure_velocity(delta_t);

        let total_temperature: f64 = self
            .particles
            .iter()
            .map(|particle| particle.temperature)
            .fold(0., |x: f64, y: f64| x.max(y));

        println!("TOTAL TEMPERATURE: {}", total_temperature);

        self.particles = self
            .particles
            .iter()
            .map(|particle| particle.control_update(self.dim, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");

        let max_speed = self
            .particles
            .iter()
            .map(|&p| p.velocity.mag())
            .fold(0., f64::max);
        // self.particles = self
        //     .particles
        //     .iter()
        //     .map(|particle| particle.with_velocity(T::default()))
        //     .collect::<Vec<Particle<T>>>()
        //     .as_slice()
        //     .try_into()
        //     .expect("Expected a Vec of a different length");

        // let max_delta_t = -self.max_cfl * self.radius / self.gravity;
        if max_speed > 0.5 {
            return self.cfl * self.radius / max_speed;
        } else {
            return (-self.cfl * self.radius / self.gravity).sqrt();
        }
    }
}
