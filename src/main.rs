extern crate kiss3d;
extern crate nalgebra as na;
pub mod map;
pub mod utility;
use map::Map;
use rand::prelude::*;
use rand::rngs::StdRng;
use std::fs;
use std::str::FromStr;

use std::{convert::TryInto, time::Instant};

use utility::{Fluid, Particle, Point3D};

use kiss3d::{camera::FirstPerson, light::Light, scene::SceneNode, window::Window};
use na::Point3;

const ROOM_TEMPERATURE: f64 = 293.15;
const MAP_SCALE: f32 = 1000.;
const NUM_PARTICLES: usize = 10115;
const DIM: Point3D = Point3D {
    x: 0.14605,
    y: 0.14605,
    z: 0.29845,
};

fn gen_points_random(
    dim_lower: &Point3D,
    dim_upper: &Point3D,
    num_particles: u64,
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

    let mut points = Vec::new();
    let contents = fs::read_to_string("particles/particle_updated_densities") //for particle without boundary, make the radius 2. and the boundary mass 1300. With large boundary needs radius 2. and boundary mass 2000., and 0.95 factor on the densities
        .expect("Something went wrong reading the file");
    let particle_positions = contents.split("\n");

    for elem in particle_positions {
        let mut coord = elem.split(", ");
        points.push(Point3D {
            x: f64::from_str(coord.next().unwrap()).unwrap(),
            y: f64::from_str(coord.next().unwrap()).unwrap(),
            z: f64::from_str(coord.next().unwrap()).unwrap(),
        })
    }
    // points.shuffle(&mut rng);

    let mut window = Window::new("Simulation! 😎");
    let particle_density = DIM.volume() / NUM_PARTICLES as f64;
    let particles: [Particle<Point3D>; NUM_PARTICLES] = points
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            Particle::<Point3D>::new(
                position,
                particle_density
                    * 1.
                    * if i > (NUM_PARTICLES * 3 / 4) {
                        Fluid::BenzylAlcohol.density(ROOM_TEMPERATURE)
                    } else {
                        Fluid::Saltwater.density(ROOM_TEMPERATURE)
                    },
                Point3D::default(),
                273.15,
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
        particle_density * 1200.,
        particle_density.cbrt() * 2.,
        DIM,
        -9.8,
    );
    println!("RADIUS: {:.8}", particle_density.cbrt() * 2.);
    println!("BOUNDARY MASS: {:.8}", particle_density * 5000.);
    let mut last_time: Instant = Instant::now();
    let mut time_elapsed = 0.;
    let mut real_time_elapsed = 0.;
    let mut i = 0;

    let vertices: [Point3<f32>; 8] = (DIM * MAP_SCALE as f64).cube(Point3D::default());

    window.set_light(Light::StickToCamera);

    window.set_background_color(1.0, 1.0, 1.0);
    let eye = (DIM * MAP_SCALE as f64 / 2.
        + Point3D::new(-0.35 * MAP_SCALE as f64, 0., -DIM.z * MAP_SCALE as f64 / 4.))
    .na_point();
    let mut first_person = FirstPerson::new(
        eye,
        Point3::new(
            DIM.x as f32 * MAP_SCALE / 2.,
            DIM.y as f32 * MAP_SCALE / 2.,
            DIM.z as f32 * MAP_SCALE / 4.,
        ),
    );

    let mut cubes: Vec<SceneNode> = vec![];
    first_person.set_move_step(10.);
    let delta_t = 0.012; //(10. as f64).recip();
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
        map.update(delta_t);
        i += 1;
        if i % 1 == 0 {
            println!(
                "{} seconds have passed, real time {}.",
                time_elapsed, real_time_elapsed
            );
        }
        disp_cube(&mut window, vertices, Point3::new(0.0, 0.0, 0.0));

        let real_delta_t = last_time.elapsed().as_secs_f64();
        real_time_elapsed += real_delta_t;
        time_elapsed += delta_t;

        // thread::sleep(Duration::from_secs_f64((delta_t - real_delta_t).max(0.)));
        last_time = Instant::now();
    }
}
