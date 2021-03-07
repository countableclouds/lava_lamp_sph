extern crate itertools;
extern crate kiss3d;
extern crate nalgebra as na;
pub mod map;
pub mod utility;
use itertools::Itertools;
use map::Map;
use rand::prelude::*;
use rand::rngs::StdRng;
use std::collections::HashSet;
use std::fs;
use std::str::FromStr;

use std::{convert::TryInto, time::Instant};

use utility::{BoundaryParticle, Coords, Fluid, Particle, Point3D};

use kiss3d::{camera::FirstPerson, light::Light, scene::SceneNode, window::Window};
use na::Point3;

const ROOM_TEMPERATURE: f64 = 293.15;
const MAP_SCALE: f32 = 1500.;
const LAYERS: usize = 6;
const DIM: Point3D = Point3D {
    x: 0.15,
    y: 0.15,
    z: 0.3 / 36. * LAYERS as f64,
};
const BOUNDARY_DIM: Point3D = Point3D {
    x: DIM.x,
    y: DIM.y,
    z: DIM.z * 3. / 2.,
};
const PARTICLES_DIM: [usize; 3] = [19, 19, LAYERS + 2];
const BOUNDARY_PARTICLES_DIM: [usize; 3] =
    [PARTICLES_DIM[0], PARTICLES_DIM[1], PARTICLES_DIM[2] * 3 / 2];
const NUM_PARTICLES: usize =
    (PARTICLES_DIM[0] - 1) * (PARTICLES_DIM[1] - 1) * (PARTICLES_DIM[2] - 1);
const NUM_BOUNDARY_PARTICLES: usize = (BOUNDARY_PARTICLES_DIM[0] + 1)
    * (BOUNDARY_PARTICLES_DIM[1] + 1)
    * (BOUNDARY_PARTICLES_DIM[2] + 1)
    - (BOUNDARY_PARTICLES_DIM[0] - 1)
        * (BOUNDARY_PARTICLES_DIM[1] - 1)
        * (BOUNDARY_PARTICLES_DIM[2] - 1)
    - (BOUNDARY_PARTICLES_DIM[0] - 1) * (BOUNDARY_PARTICLES_DIM[1] - 1);
const WRITE_TO_FILE: bool = false;

fn gen_points_random(
    dim_lower: &Point3D,
    dim_upper: &Point3D,
    num_particles: usize,
    rng: &mut StdRng,
) -> Vec<Point3D> {
    (0..num_particles)
        .into_iter()
        .map(|_| {
            Point3D::new(
                rng.gen_range(dim_lower.x, dim_upper.x),
                rng.gen_range(dim_lower.y, dim_upper.y),
                rng.gen_range(dim_lower.z, dim_upper.z),
            )
        })
        .collect()
}

fn gen_boundary_points(
    &start: &Point3D,
    &width: &Point3D,
    &height: &Point3D,
    width_particles: usize,
    height_particles: usize,
) -> Vec<Point3D> {
    let mut boundary_points: Vec<Point3D> = Vec::new();
    println!("start: {}, width: {}, height: {}", start, width, height);
    for i in 0..width_particles + 1 {
        for j in 0..height_particles + 1 {
            boundary_points.push(
                start
                    + width * ((i as f64) / (width_particles as f64))
                    + height * ((j as f64) / (height_particles as f64)),
            );
        }
    }

    boundary_points
}

fn gen_fluid_points(
    &start: &Point3D,
    &width: &Point3D,
    &length: &Point3D,
    &height: &Point3D,
    width_particles: usize,
    length_particles: usize,
    height_particles: usize,
) -> Vec<Point3D> {
    let mut fluid_points: Vec<Point3D> = Vec::new();
    for k in 1..height_particles {
        for i in 1..width_particles {
            for j in 1..length_particles {
                fluid_points.push(
                    start
                        + width * ((i as f64) / (width_particles as f64))
                        + length * ((j as f64) / (length_particles as f64))
                        + height * ((k as f64) / (height_particles as f64)),
                );
            }
        }
    }

    fluid_points
}

fn main() {
    let mut rng = StdRng::seed_from_u64(2);

    let mut points = gen_points_random(
        &Point3D::new(0., 0., DIM.z / 4.),
        &Point3D::new(DIM.x, DIM.y, DIM.z),
        NUM_PARTICLES * 3 / 4,
        &mut rng,
    );
    let mut benzyl_points = gen_points_random(
        &Point3D::default(),
        &Point3D::new(DIM.x, DIM.y, DIM.z / 4.),
        NUM_PARTICLES - NUM_PARTICLES * 3 / 4,
        &mut rng,
    );
    points.append(&mut benzyl_points);

    let mut points = Vec::new();
    let contents = fs::read_to_string("init_distros/dist_lava_lamp")
        .expect("Something went wrong reading the file");
    let particle_positions = contents.split("\n");

    for elem in particle_positions {
        let mut coord = elem.split(", ");
        points.push(Point3D {
            x: f64::from_str(coord.next().unwrap()).unwrap(),
            y: f64::from_str(coord.next().unwrap()).unwrap(),
            z: f64::from_str(coord.next().unwrap()).unwrap(),
        });
    }
    let points = gen_fluid_points(
        &Point3D::default(),
        &Point3D::default().with_x(DIM.x),
        &Point3D::default().with_y(DIM.y),
        &Point3D::default().with_z(DIM.z),
        PARTICLES_DIM[0],
        PARTICLES_DIM[1],
        PARTICLES_DIM[2],
    );
    // points.shuffle(&mut rng);
    let mut boundary_points = Vec::new();
    for elem in BOUNDARY_DIM
        .proj()
        .iter()
        .zip(BOUNDARY_PARTICLES_DIM.iter())
        .combinations(2)
    {
        //combinations(2) {
        boundary_points.append(&mut gen_boundary_points(
            &Point3D::default(),
            &(BOUNDARY_DIM - *elem[0].0),
            &(BOUNDARY_DIM - *elem[1].0),
            *elem[0].1,
            *elem[1].1,
        ));
    }

    for elem in BOUNDARY_DIM
        .proj()
        .iter()
        .zip(BOUNDARY_PARTICLES_DIM.iter())
        .combinations(2)
    {
        if elem[0].0.z != 0. && elem[1].0.z != 0. {
            continue;
        }
        boundary_points.append(&mut gen_boundary_points(
            &BOUNDARY_DIM,
            &(*elem[0].0 - BOUNDARY_DIM),
            &(*elem[1].0 - BOUNDARY_DIM),
            *elem[0].1,
            *elem[1].1,
        ));
    }
    let mut indices = HashSet::new();
    for (i, &elem) in boundary_points.iter().enumerate() {
        for (j, &elem2) in boundary_points.iter().enumerate() {
            if (elem - elem2).mag() < 0.0000001 && i != j {
                indices.insert(i.min(j));
            }
        }
    }
    boundary_points = boundary_points
        .iter()
        .enumerate()
        .filter(|(i, _)| !indices.contains(i))
        .map(|(_, point)| *point)
        .collect();
    println!("{}", boundary_points.len());
    let boundary_particles: [BoundaryParticle<Point3D>; NUM_BOUNDARY_PARTICLES] = boundary_points
        .into_iter()
        .map(|position| {
            BoundaryParticle::<Point3D>::new(position, Point3D::default(), ROOM_TEMPERATURE)
        })
        .collect::<Vec<BoundaryParticle<Point3D>>>()
        .try_into()
        .map_err(|v: Vec<BoundaryParticle<Point3D>>| {
            format!(
                "Expected a Vec of length {} instead of {}",
                NUM_BOUNDARY_PARTICLES,
                v.len()
            )
        })
        .unwrap();
    let mut window = Window::new("Simulation! 😎");
    // let particle_density = DIM.volume() / NUM_PARTICLES as f64;
    let particle_density = DIM.volume()
        / ((PARTICLES_DIM[0] + 1) * (PARTICLES_DIM[1] + 1) * (PARTICLES_DIM[2] + 1)) as f64;
    let particles: [Particle<Point3D>; NUM_PARTICLES] = points
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            Particle::<Point3D>::new(
                position,
                particle_density
                    * if i > (NUM_PARTICLES * 3 / 4) {
                        Fluid::BenzylAlcohol.density(ROOM_TEMPERATURE)
                    } else {
                        Fluid::Saltwater.density(ROOM_TEMPERATURE)
                    },
                Point3D::default(),
                ROOM_TEMPERATURE,
                if i > (NUM_PARTICLES * 3 / 4) {
                    Fluid::BenzylAlcohol
                } else {
                    Fluid::Saltwater
                },
            )
        })
        .collect::<Vec<Particle<Point3D>>>()
        .try_into()
        .map_err(|v: Vec<Particle<Point3D>>| {
            format!(
                "Expected a Vec of length {} instead of {}",
                NUM_PARTICLES,
                v.len()
            )
        })
        .unwrap();

    let mut map = Map::new(
        particles,
        boundary_particles,
        particle_density.cbrt() * 2.,
        DIM,
        -9.8,
        0.05,
    );
    println!("RADIUS: {:.8}", particle_density.cbrt() * 2.);
    let mut last_time: Instant = Instant::now();
    let mut time_elapsed = 0.;
    let mut real_time_elapsed = 0.;
    let mut i = 0;

    let vertices: [Point3<f32>; 8] = (DIM * MAP_SCALE as f64).cube(Point3D::default());

    window.set_light(Light::StickToCamera);

    window.set_background_color(1.0, 1.0, 1.0);
    let eye = (BOUNDARY_DIM * MAP_SCALE as f64 / 2.
        + Point3D::new(-0.35 * MAP_SCALE as f64, 0., 0.))
    .na_point();
    let mut first_person = FirstPerson::new(eye, (BOUNDARY_DIM * MAP_SCALE as f64 / 2.).na_point());

    let mut cubes: Vec<SceneNode> = vec![];
    first_person.set_move_step(10.);
    let mut delta_t = 0.01; //(10. as f64).recip();

    while window.render_with_camera(&mut first_person) {
        for cube in &mut cubes {
            window.remove_node(cube);
        }
        cubes = vec![];
        for (j, particle) in (&map.particles).iter().enumerate() {
            cubes.push(window.add_cube(0.001 * MAP_SCALE, 0.001 * MAP_SCALE, 0.001 * MAP_SCALE));

            let cube = cubes.get_mut(j).unwrap();

            let color = if particle.fluid_type == Fluid::BenzylAlcohol {
                particle
                    .fluid_type
                    .simulation_color(particle.temperature as f32, 300.15, 307.15)
            } else {
                particle
                    .fluid_type
                    .simulation_color(particle.temperature as f32, 293.15, 307.15)
            };
            cube.set_color(color.x, color.y, color.z);

            cube.append_translation(&(particle.position * MAP_SCALE as f64).translation())
        }

        for (j, particle) in (&map.boundary_particles).iter().enumerate() {
            cubes.push(window.add_cube(0.001 * MAP_SCALE, 0.001 * MAP_SCALE, 0.001 * MAP_SCALE));

            let cube = cubes.get_mut(j + map.particles.len()).unwrap();

            cube.set_color(1., 0., 0.);

            cube.append_translation(&(particle.position * MAP_SCALE as f64).translation())
        }

        map.update(delta_t);
        // delta_t = 0.026449456893745664;
        println!("The change in time is: {}", delta_t);
        i += 1;
        if i % 10 == 0 {
            println!(
                "{} seconds have passed, real time {}.",
                time_elapsed, real_time_elapsed
            );
            if WRITE_TO_FILE {
                let mut data: String = "".to_owned();
                for &particle in &map.particles {
                    data.push_str(&particle.position.to_string());
                    data.push_str("\n");
                }
                data = data[..data.len() - 1].to_owned();

                fs::write(format!("init_distros/dist_{}", i / 10), data)
                    .expect("Unable to write file");
            }
        }

        let real_delta_t = last_time.elapsed().as_secs_f64();
        real_time_elapsed += real_delta_t;
        time_elapsed += delta_t;

        // thread::sleep(Duration::from_secs_f64((delta_t - real_delta_t).max(0.)));
        last_time = Instant::now();
    }
}
