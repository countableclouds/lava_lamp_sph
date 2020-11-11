pub mod map;
pub mod utility;
use map::Map;
use piston_window::*;
use rand::prelude::*;
use std::{
    collections::HashMap,
    convert::TryInto,
    time::{Duration, Instant},
};

use utility::{Fluid, Particle, Point};
const ROOM_TEMPERATURE: f64 = 293.15;
const MAP_SCALE: f64 = 2000.;

fn gen_points(dim: &Point, num_particles: u64, rng: &mut ThreadRng) -> Vec<Point> {
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
                (p % width_particles as u64) as f64 * length
                    + width_offset
                    + (rng.gen::<f64>() - 0.5) * length / 2.,
                (p / width_particles as u64) as f64 * length
                    + height_offset
                    + (rng.gen::<f64>() - 0.5) * length / 2.,
            )
        })
        .collect()
}
fn main() {
    let mut window: PistonWindow = WindowSettings::new("Lava Lamp Simulation", [1000, 1000])
        .exit_on_esc(true)
        .build()
        .unwrap();
    let num_particles: u64 = 1000;
    let dim = Point::new(0.14605, 0.29845);
    let mut rng = rand::thread_rng();
    let points = gen_points(&dim, num_particles, &mut rng);
    let map_visual_margins: Point = Point::new(25., 25.);

    let particles: [Particle<Point>; 945] = points
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            Particle::<Point>::new(
                position,
                dim.area() / num_particles as f64
                    * if i < (num_particles * 3 / 4) as usize {
                        Fluid::Saltwater.density(ROOM_TEMPERATURE)
                    } else {
                        Fluid::BenzylAlcohol.density(ROOM_TEMPERATURE)
                    },
                ROOM_TEMPERATURE,
                if i < (num_particles * 3 / 4) as usize {
                    Fluid::Saltwater
                } else {
                    Fluid::BenzylAlcohol
                },
            )
        })
        .collect::<Vec<Particle<Point>>>()
        .as_slice()
        .try_into()
        .expect("Expected a Vec of a different length");

    let mut map = Map {
        particles,
        dim,
        radius: 0.01,
        particle_map: HashMap::new(),
        gravity: 9.8,
        gas_constant: 8.314,
    };
    let mut last_time: Instant = Instant::now();
    let mut time_elapsed = 0.;
    let mut real_time_elapsed = 0.;
    let mut i = 0;

    // for i in 0..40000 {
    //     if i == 1 {
    //         println!("{:?}", map.particles);
    //     }
    //     if i % 100 == 0 {
    //         println!("{} iterations have passed.", i);
    //     }
    //     map.update((11 as f64).recip());
    // }
    while let Some(e) = window.next() {
        if i % 1 == 0 {
            println!(
                "{} seconds have passed, real time {}.",
                time_elapsed, real_time_elapsed
            );
        }
        let delta_t = 250000f64.recip();

        map.update(delta_t);

        i += 1;

        window.draw_2d(&e, |c, g, _device| {
            clear([1.0; 4], g);
            line_from_to(
                [0., 0., 0., 1.],
                1.,
                Point::new(0., 0.) * MAP_SCALE + map_visual_margins,
                Point::new(0., dim.y) * MAP_SCALE + map_visual_margins,
                c.transform,
                g,
            );

            line_from_to(
                [0., 0., 0., 1.],
                1.,
                Point::new(dim.x, 0.) * MAP_SCALE + map_visual_margins,
                Point::new(dim.x, dim.y) * MAP_SCALE + map_visual_margins,
                c.transform,
                g,
            );

            line_from_to(
                [0., 0., 0., 1.],
                1.,
                Point::new(0., dim.y) * MAP_SCALE + map_visual_margins,
                Point::new(dim.x, dim.y) * MAP_SCALE + map_visual_margins,
                c.transform,
                g,
            );

            line_from_to(
                [0., 0., 0., 1.],
                1.,
                Point::new(0., 0.) * MAP_SCALE + map_visual_margins,
                Point::new(dim.x, 0.) * MAP_SCALE
                    // + Point::new(real_time_elapsed, 0.)
                    + map_visual_margins,
                c.transform,
                g,
            );
            // println!(
            //     "{:?}",
            //     map.particles
            //         .iter()
            //         .map(|particle| particle.position)
            //         .collect::<Vec<Point>>()
            // );
            // println!("{:?}", (&map.particles[269..271]));
            for particle in &map.particles {
                rectangle_from_to(
                    match particle.fluid_type {
                        Fluid::Saltwater => [0., 0., 0., 1.],
                        Fluid::BenzylAlcohol => [1., 0., 0., 1.],
                    },
                    particle.position * MAP_SCALE + Point::new(1., 1.) + map_visual_margins,
                    particle.position * MAP_SCALE - Point::new(1., 1.) + map_visual_margins,
                    c.transform,
                    g,
                );
            }
        });
        let real_delta_t = last_time.elapsed().as_secs_f64();
        real_time_elapsed += real_delta_t;
        time_elapsed += delta_t;
        // thread::sleep(Duration::from_secs_f64((delta_t - real_delta_t).max(0.)));
        last_time = Instant::now();
    }
}
