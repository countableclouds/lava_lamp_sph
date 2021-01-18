use crate::utility::*;
use crate::NUM_PARTICLES;
use rayon::prelude::*;
use std::{cmp::Eq, collections::HashMap, convert::TryInto, fmt::Display, hash::Hash};

pub struct Map<T: Coords + Copy> {
    pub particles: [Particle<T>; NUM_PARTICLES],
    pub dim: T,
    pub radius: f64,
    pub boundary_mass: f64,
    pub particle_map: HashMap<T::Key, Vec<(usize, Particle<T>)>>,
    pub gravity: f64,
}

pub const TEST_NUM: usize = 8782; //2223; // ;
const PRESSURE_FACTOR: f64 = 1.;
const MASS_FACTOR: f64 = 1.;
const checker: bool = false;
const error_checker: bool = false;

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
        boundary_mass: f64,
        radius: f64,
        dim: T,
        gravity: f64,
    ) -> Self {
        Map {
            particles,
            dim,
            radius,
            boundary_mass,
            particle_map: HashMap::new(),
            gravity,
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
    }

    pub fn update_density(&self, i: usize, particle: &Particle<T>) -> Particle<T> {
        let mut density = 0.;
        let mut point_count = 0;
        // if particle.position.height() < 0. {
        //     println!("{}", i);
        //     assert!(1 == 0);
        // }
        for point in particle.position.along_axes(self.radius) {
            if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                for (_, other_particle) in particles {
                    if if let Some(mass) = other_particle.mass {
                        mass
                    } else {
                        particle.mass.unwrap()
                    } * particle
                        .position
                        .cubic(self.radius, other_particle.position)
                        > 0.
                    {
                        // if i == TEST_NUM {
                        //     println!("{}", other_particle.position);
                        // };
                        point_count += 1
                    }
                    density += if let Some(mass) = other_particle.mass {
                        mass
                    } else {
                        particle.mass.unwrap()
                    } * particle
                        .position
                        .cubic(self.radius, other_particle.position);
                }
            }
        }
        if i == TEST_NUM {
            println!("Point count: {}", point_count);
        }
        let mut boundary_density = 0.;
        for boundary_particle in particle.position.proj() {
            boundary_density += self.boundary_mass
                * MASS_FACTOR
                * particle.position.cubic(self.radius, boundary_particle);
        }
        for boundary_particle in (self.dim - particle.position).proj() {
            boundary_density += self.boundary_mass
                * MASS_FACTOR
                * particle
                    .position
                    .cubic(self.radius, self.dim - boundary_particle);
        }
        if i == TEST_NUM {
            println!("Boundary density: {}", boundary_density);
        }
        // if density + boundary_density < 900. {
        //     println!(
        //         "Position: {:?}, Density: {}, Boundary Density: {}, Number: {}",
        //         particle.position, density, boundary_density, i
        //     )
        // }
        particle.with_density(density + boundary_density)
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

                    water_benzyl_curvature += other_particle.volume()
                        * particle.fluid_type.unwrap().color(
                            if let Some(fluid) = other_particle.fluid_type {
                                fluid
                            } else {
                                particle.fluid_type.unwrap()
                            },
                        )
                        * particle
                            .position
                            .laplace_quintic(self.radius, other_particle.position);

                    water_benzyl_normal += particle
                        .position
                        .grad_cubic(self.radius, other_particle.position)
                        * (-other_particle.volume()
                            * particle.fluid_type.unwrap().color(
                                if let Some(fluid) = other_particle.fluid_type {
                                    fluid
                                } else {
                                    particle.fluid_type.unwrap()
                                },
                            ));

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

        acceleration = acceleration * particle.density.unwrap().recip()
            + T::default().with_height(self.gravity);
        temperature = particle.temperature.unwrap() + temperature * delta_t;
        particle
            .with_velocity(particle.velocity.unwrap() + acceleration * delta_t)
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
                            .grad_cubic(self.radius, other_particle.position)
                            * (if let Some(mass) = other_particle.mass {
                                mass
                            } else {
                                particle.mass.unwrap()
                            } * particle.density.unwrap().powi(-2));
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
                            .grad_cubic(self.radius, other_particle.position)
                            .dot(
                                inverse_densities[i]
                                    + if !other_particle.mass.is_none() {
                                        particle
                                            .position
                                            .grad_cubic(self.radius, other_particle.position)
                                            * (particle.mass.unwrap()
                                                * particle.density.unwrap().powi(-2))
                                    } else {
                                        T::default()
                                    },
                            )
                            * (if let Some(mass) = other_particle.mass {
                                mass
                            } else {
                                particle.mass.unwrap()
                            } * delta_t.powi(2));
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
                    for (_, other_particle) in particles {
                        density_diff[i] -= particle
                            .position
                            .grad_cubic(self.radius, other_particle.position)
                            .dot(
                                particle.velocity.unwrap()
                                    - if let Some(velocity) = other_particle.velocity {
                                        velocity
                                    } else {
                                        T::default()
                                    },
                            )
                            * (if let Some(mass) = other_particle.mass {
                                mass
                            } else {
                                particle.mass.unwrap()
                            } * delta_t);
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
            if particle.density.unwrap() == 0. {
                continue;
            }
            for point in particle.position.along_axes(self.radius) {
                if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                    for (j, other_particle) in particles {
                        if other_particle.density.unwrap() < 0.000000001 {
                            continue;
                        }
                        accelerations[i] += particle
                            .position
                            .grad_cubic(self.radius, other_particle.position)
                            * (-if let Some(mass) = other_particle.mass {
                                mass
                            } else {
                                particle.mass.unwrap()
                            } * (pressures[i] * particle.density.unwrap().powi(-2)
                                + if let Some(density) = other_particle.density {
                                    pressures[*j] * density.powi(-2)
                                } else {
                                    0.
                                }));
                        if i == TEST_NUM
                            && checker
                            && false
                            && particle
                                .position
                                .grad_cubic(self.radius, other_particle.position)
                                != T::default()
                            && pressures[*j] != 0.
                        {
                            println!(
                                "Bordering Pressure: {}, Position: {}, Index: {}, Density: {}",
                                pressures[*j],
                                other_particle.position,
                                j,
                                self.particles[*j].density.unwrap()
                            );
                            println!("{}", accelerations[i]);
                        }
                    }
                }
            }
        }

        accelerations
    }

    pub fn update_pressure_velocity(&mut self, delta_t: f64) {
        let num_iter: u64 = 15;
        let relaxation_coeff = 0.5;
        let mut pressures: [f64; NUM_PARTICLES] = [0.; NUM_PARTICLES];

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
            .map(|(i, p)| (i, p.density.unwrap()))
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
            let mut image: [f64; NUM_PARTICLES] = [0.; NUM_PARTICLES];
            error = 0.;
            let mut count = 0;
            for (i, particle) in self.particles.iter().enumerate() {
                if diagonal[i] == 0. {
                    continue;
                }
                for point in particle.position.along_axes(self.radius) {
                    if let Some(particles) = self.particle_map.get(&point.bin(self.radius)) {
                        for (j, other_particle) in particles {
                            image[i] += particle
                                .position
                                .grad_cubic(self.radius, other_particle.position)
                                .dot(
                                    pressure_accelerations[i]
                                        - if *j < NUM_PARTICLES {
                                            pressure_accelerations[*j]
                                        } else {
                                            T::default()
                                        },
                                )
                                * (if let Some(mass) = other_particle.mass {
                                    mass
                                } else {
                                    particle.mass.unwrap()
                                } * delta_t.powi(2));
                        }
                    }
                }
                if i == TEST_NUM && checker {
                    println!(
                        "Test Particles Velocity: {}",
                        self.particles[TEST_NUM].velocity.unwrap()
                    );
                    println!(
                        "Test Particles Pressure Acceleration: {}",
                        pressure_accelerations[TEST_NUM]
                    );
                    println!("Test Particles Pressure : {}", pressures[TEST_NUM]);
                }

                pressures[i] += relaxation_coeff * (density_diff[i] - image[i]) / diagonal[i];
                if i == TEST_NUM {
                    println!(
                        "Density Difference: {}, Image: {}, Diagonal: {}",
                        density_diff[i], image[i], diagonal[i]
                    );
                }
                // if image[TEST_NUM] > density_diff[TEST_NUM] {
                //     println!("{}", i);
                //     assert!(1 == 0);
                // }
                // if image[i] > density_diff[i].max(0.) + 10000. {
                //     println!("{}", i);
                //     assert!(1 == 0);
                // }

                if pressures[i] < 0. {
                    count += 1;
                }

                // pressures[i] = (0.29845 - self.particles[i].position.height()) * 100.;
                pressures[i] = pressures[i].max(0.);

                // if i == TEST_NUM {
                //     println!("Pressure: {}", pressures[i]);
                //     println!("Image: {}", image[i]);x
                //     println!("Height: {}", self.particles[i].position.height());
                // }
                if error_checker {
                    if i == TEST_NUM {
                        pressures[i] = 1.;
                    } else {
                        pressures[i] = 1.;
                    }
                }
                // if i != TEST_NUM {
                //     pressures[i] = 0.;
                // }
            }
            println!("COUNT!: {}", count);
            // let max_pressure = pressures.iter().cloned().fold(0. / 0., f64::max);
            // for elem in &mut pressures {
            //     *elem /= max_pressure / 10.;
            // }
            pressure_accelerations = self.get_pressure_accelerations(pressures);
            // println!(
            //     "New velocity: {:?}",
            //     pressure_accelerations[TEST_NUM] * delta_t
            // );
            for i in 0..image.len() {
                error += (density_diff[i] - image[i]).abs() / self.particles[i].expected_density();
            }
            let mut image_average = 0.;
            for i in 0..image.len() {
                image_average += image[i].abs() / self.particles[i].expected_density();
            }
            if error_checker {
                println!(
                    "Density Differences: {:?}, {:?}, {:?}",
                    image[TEST_NUM] + self.particles[TEST_NUM].density.unwrap(),
                    pressures[TEST_NUM],
                    image[TEST_NUM]
                );
            }

            if j == 0 || j == num_iter - 1 || true {
                println!("Error Percentage: {:?}, {}", error / image.len() as f64, j);
            }

            // println!("Image Average: {:?}", image_average / image.len() as f64);

            // if j == 1 {
            //     break;
            // }
            // let bob: Vec<(usize, &f64)> = pressures
            //     .iter()
            //     .enumerate()
            //     .filter(|(i, &n)| n == 0.)
            //     .take(10)
            //     .collect();
            // let bob2: Vec<T> = bob
            //     .into_iter()
            //     .map(|(i, h)| self.particles[i].position)
            //     .collect();
            // println!("{:?}", bob2);  break;
            // }
        }
        println!(
            "count: {}, pressure average: {}",
            pressures.len() - pressures.iter().filter(|&n| *n <= 0.).count(),
            pressures.iter().fold(0., |a, b| a + b) / pressures.len() as f64
        );
        self.particles = self
            .particles
            .iter()
            .enumerate()
            .map(|(i, particle)| {
                particle
                    .with_velocity(particle.velocity.unwrap() + pressure_accelerations[i] * delta_t)
            })
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");
        println!("Velocity: {}", self.particles[TEST_NUM].velocity.unwrap());
    }

    pub fn update(&mut self, delta_t: f64) {
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
            .map(|particle| self.update_nonpressure_velocity(particle, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");
        self.update_hashmap();
        self.update_pressure_velocity(delta_t);
        self.particles = self
            .particles
            .iter()
            .map(|particle| particle.control_update(self.dim, delta_t))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");
        self.particles = self
            .particles
            .iter()
            .map(|particle| particle.with_velocity(T::default()))
            .collect::<Vec<Particle<T>>>()
            .as_slice()
            .try_into()
            .expect("Expected a Vec of a different length");
    }
}
