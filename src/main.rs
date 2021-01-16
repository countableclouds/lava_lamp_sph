extern crate kiss3d;
extern crate nalgebra as na;
pub mod map;
pub mod utility;
use map::{Map, TEST_NUM};
use rand::prelude::*;
use rand::rngs::StdRng;
use std::fs;
use std::str::FromStr;

use std::{
    collections::HashMap,
    convert::TryInto,
    time::{Duration, Instant},
};

use utility::{Coords, Fluid, Particle, Point, Point3D};

use kiss3d::{camera::FirstPerson, light::Light, scene::SceneNode, window::Window};
use na::Point3;

const ROOM_TEMPERATURE: f64 = 293.15;
const MAP_SCALE: f32 = 1000.;
const PARTICLES_UPPER_BOUND: usize = 12000;
const NUM_PARTICLES: usize = 10115;
const DIM: Point3D = Point3D {
    x: 0.14605,
    y: 0.14605,
    z: 0.29845,
};
fn gen_points(dim: &Point, num_particles: u64) -> Vec<Point> {
    let length = (-(dim.squared_mag() + dim.area() * (4. * (num_particles as f64) - 2.)).sqrt()
        + dim.x
        + dim.y)
        / (2. - (2 * num_particles) as f64);
    let width_particles = (dim.x / length - 1.) as u64;
    let height_particles = (dim.y / length - 1.) as u64;
    let width_offset = (dim.x - ((width_particles - 1) as f64 * length)) / 2.;
    let height_offset = (dim.y - ((height_particles - 1) as f64 * length)) / 2.;
    println!(
        "Using {} particles, with {} as the initial maximum. Radius {}.",
        width_particles * height_particles,
        num_particles,
        length
    );
    (0..(width_particles * height_particles))
        .into_iter()
        .map(|p| {
            Point::new(
                (p % width_particles) as f64 * length + width_offset,
                (p / width_particles) as f64 * length + height_offset,
            )
        })
        .collect()
}

fn gen_points_grid(dim: &Point3D, num_particles: u64, rng: &mut Option<StdRng>) -> Vec<Point3D> {
    let length = (dim.volume()
        / (num_particles as f64 - dim.x * dim.y - dim.y * dim.z - dim.z * dim.x))
        .cbrt();
    let width_particles = (dim.x / length - 1.) as u64;
    let length_particles = (dim.y / length - 1.) as u64;
    let height_particles = (dim.z / length - 1.) as u64;
    let width_offset = (dim.x - ((width_particles - 1) as f64 * length)) / 2.;
    let length_offset = (dim.y - ((length_particles - 1) as f64 * length)) / 2.;
    let height_offset = (dim.z - ((height_particles - 1) as f64 * length)) / 2.;
    println!(
        "Using {} particles, with {} as the initial maximum. Radius {}. x,y,z offset: {}, {}, {}",
        width_particles * height_particles * length_particles,
        num_particles,
        length,
        width_offset,
        length_offset,
        height_offset
    );
    (0..(width_particles * height_particles * length_particles))
        .into_iter()
        .map(|p| {
            Point3D::new(
                (p % width_particles) as f64 * length
                    + width_offset
                    + match rng {
                        Some(rng_gen) => (rng_gen.gen::<f64>() - 0.5) * length / 4.,
                        None => 0.,
                    },
                ((p % (width_particles * length_particles)) / width_particles as u64) as f64
                    * length
                    + length_offset
                    + match rng {
                        Some(rng_gen) => (rng_gen.gen::<f64>() - 0.5) * length / 4.,
                        None => 0.,
                    },
                (p / (width_particles * length_particles)) as f64 * length
                    + height_offset
                    + match rng {
                        Some(rng_gen) => (rng_gen.gen::<f64>() - 0.5) * length / 4.,
                        None => 0.,
                    },
            )
        })
        .collect()
}
fn gen_points_random(
    dim_lower: &Point3D,
    dim_upper: &Point3D,
    num_particles: u64,
    rng: &mut StdRng,
) -> Vec<Point3D> {
    let dim = dim_upper.clone() - dim_lower.clone();
    (0..num_particles)
        .into_iter()
        .map(|p| {
            Point3D::new(
                rng.gen_range(dim_lower.x, dim_upper.x),
                rng.gen_range(dim_lower.y, dim_upper.y),
                rng.gen_range(dim_lower.z, dim_upper.z),
            )
        })
        .collect()
}

fn disp_cube(window: &mut Window, vertices: [Point3<f32>; 8], color: Point3<f32>) {
    window.draw_line(&vertices[0], &vertices[1], &color);
    window.draw_line(&vertices[0], &vertices[2], &color);
    window.draw_line(&vertices[0], &vertices[4], &color);
    window.draw_line(&vertices[2], &vertices[3], &color);
    window.draw_line(&vertices[4], &vertices[5], &color);
    window.draw_line(&vertices[2], &vertices[6], &color);
    window.draw_line(&vertices[1], &vertices[5], &color);
    window.draw_line(&vertices[4], &vertices[6], &color);
    window.draw_line(&vertices[1], &vertices[3], &color);
    window.draw_line(&vertices[7], &vertices[3], &color);
    window.draw_line(&vertices[7], &vertices[5], &color);
    window.draw_line(&vertices[7], &vertices[6], &color);
}
fn main() {
    let mut rng = StdRng::seed_from_u64(2);
    // let mut rng = Some(StdRng::seed_from_u64(2));

    // let mut rng: Option<ThreadRng> = None;

    let mut points = gen_points_random(
        &Point3D::new(0., 0., DIM.z / 4.),
        &Point3D::new(DIM.x, DIM.y, DIM.z),
        (NUM_PARTICLES * 3 / 4) as u64,
        &mut rng,
    );
    let mut benzyl_points = gen_points_random(
        &Point3D::default(),
        &Point3D::new(DIM.x, DIM.y, DIM.z / 4.),
        NUM_PARTICLES as u64 - (NUM_PARTICLES * 3 / 4) as u64,
        &mut rng,
    );
    points.append(&mut benzyl_points);

    // let points = gen_points_grid(&DIM, PARTICLES_UPPER_BOUND as u64, &mut rng);

    // let mut points = Vec::new();
    // let contents = fs::read_to_string("particles/particle_2000")
    //     .expect("Something went wrong reading the file");
    // let particle_positions = contents.split("\n");

    // for elem in particle_positions {
    //     let mut coord = elem.split(", ");
    //     points.push(Point3D {
    //         x: f64::from_str(coord.next().unwrap()).unwrap(),
    //         y: f64::from_str(coord.next().unwrap()).unwrap(),
    //         z: f64::from_str(coord.next().unwrap()).unwrap(),
    //     })
    // }

    let mut window = Window::new("Simulation! 😎");
    let particle_density = DIM.volume() / NUM_PARTICLES as f64;

    let particles: [Particle<Point3D>; NUM_PARTICLES] = points
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            Particle::<Point3D>::new(
                position,
                particle_density
                    * if i < (NUM_PARTICLES * 3 / 4) {
                        Fluid::Saltwater.density(ROOM_TEMPERATURE)
                    } else {
                        Fluid::BenzylAlcohol.density(ROOM_TEMPERATURE)
                    },
                ROOM_TEMPERATURE,
                if i < (NUM_PARTICLES * 3 / 4) {
                    Fluid::Saltwater
                } else {
                    Fluid::BenzylAlcohol
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
        particle_density * 1000.,
        particle_density.cbrt() * 2.,
        DIM,
        -0.,
    );
    println!("RADIUS: {:.8}", particle_density.cbrt() * 2.);
    println!("BOUNDARY MASS: {:.8}", particle_density * 1200.);
    let mut last_time: Instant = Instant::now();
    let mut time_elapsed = 0.;
    let mut real_time_elapsed = 0.;
    let mut i = 0;

    let vertices: [Point3<f32>; 8] = (DIM * MAP_SCALE as f64).cube(Point3D::default());

    window.set_light(Light::StickToCamera);

    window.set_background_color(1.0, 1.0, 1.0);
    let eye =
        (DIM * MAP_SCALE as f64 / 2. + Point3D::new(0., -0.35 * MAP_SCALE as f64, 0.)).na_point();
    let mut first_person = FirstPerson::new(eye, (DIM * MAP_SCALE as f64 / 2.).na_point());
    let mut cubes: Vec<SceneNode> = vec![];
    let delta_t = (10. as f64).recip();
    while window.render_with_camera(&mut first_person) {
        for cube in &mut cubes {
            window.remove_node(cube);
        }
        cubes = vec![];
        for (j, particle) in (&map.particles).iter().enumerate() {
            cubes.push(window.add_cube(0.001 * MAP_SCALE, 0.001 * MAP_SCALE, 0.001 * MAP_SCALE));
            let cube = cubes.get_mut(j).unwrap();

            let color = particle.fluid_type.simulation_color();
            cube.set_color(color.x, color.y, color.z);
            if j == TEST_NUM {
                cube.set_color(0., 0., 0.);
            }
            cube.append_translation(&(particle.position * MAP_SCALE as f64).translation())
        }

        i += 1;
        if i % 1 == 0 {
            println!(
                "{} seconds have passed, real time {}.",
                time_elapsed, real_time_elapsed
            );
        }
        if i % 10 == 0 {
            let mut data: String = "".to_owned();
            for &particle in &map.particles {
                data.push_str(&particle.position.to_string());
                data.push_str("\n");
            }
            data = data[..data.len() - 1].to_owned();

            fs::write(format!("particles/particle_{}", i / 10), data)
                .expect("Unable to write file");
        }
        disp_cube(&mut window, vertices, Point3::new(0.0, 0.0, 0.0));

        if i > 0 {
            map.update(delta_t)
        };

        // assert!(1 == 0);
        let real_delta_t = last_time.elapsed().as_secs_f64();
        real_time_elapsed += real_delta_t;
        time_elapsed += delta_t;
        // thread::sleep(Duration::from_secs_f64((delta_t - real_delta_t).max(0.)));
        last_time = Instant::now();
    }
}
