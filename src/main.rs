extern crate kiss3d;
extern crate nalgebra as na;
pub mod map;
pub mod utility;
use map::Map;
use rand::prelude::*;
use std::{
    collections::HashMap,
    convert::TryInto,
    time::{Duration, Instant},
};

use utility::{Fluid, Particle, Point, Point3D};

use kiss3d::context::Context;
use kiss3d::light::Light;
use kiss3d::window::Window;
use na::{Point3, UnitQuaternion};

const ROOM_TEMPERATURE: f64 = 293.15;

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
                (p % width_particles) as f64 * length
                    + width_offset
                    + (rng.gen::<f64>() - 0.5) * length / 2.,
                (p / width_particles) as f64 * length
                    + height_offset
                    + (rng.gen::<f64>() - 0.5) * length / 2.,
            )
        })
        .collect()
}

fn gen_points_3d(dim: &Point3D, num_particles: u64, rng: &mut ThreadRng) -> Vec<Point3D> {
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
        "Using {} particles, with {} as the initial maximum. Radius {}.",
        width_particles * height_particles * length_particles,
        num_particles,
        length
    );
    (0..(width_particles * height_particles * length_particles))
        .into_iter()
        .map(|p| {
            Point3D::new(
                (p % width_particles) as f64 * length
                    + width_offset
                    + (rng.gen::<f64>() - 0.5) * length / 2.,
                ((p % (width_particles * height_particles)) / width_particles as u64) as f64
                    * length
                    + length_offset
                    + (rng.gen::<f64>() - 0.5) * length / 2.,
                (p / width_particles * height_particles) as f64 * length
                    + height_offset
                    + (rng.gen::<f64>() - 0.5) * length / 2.,
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
    let num_particles: u64 = 14000;
    let dim = Point3D::new(0.14605, 0.14605, 0.29845);
    let mut rng = rand::thread_rng();
    let points = gen_points_3d(&dim, num_particles, &mut rng);
    let map_visual_margins: Point = Point::new(25., 25.);

    let particles: [Particle<Point3D>; 10693] = points
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            Particle::<Point3D>::new(
                position,
                dim.volume() / num_particles as f64
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
        .collect::<Vec<Particle<Point3D>>>()
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
    let mut window = Window::new("Simulation! 😎");
    // window.camera

    // let mut c = window.add_line(0.5, 0.5, 0.5);

    // c.set_color(1.0, 0.0, 0.0);
    let vertices: [Point3<f32>; 8] = (dim * 10.).cube(Point3D::default());

    window.set_light(Light::StickToCamera);

    window.set_background_color(1.0, 1.0, 1.0);

    // let rot = UnitQuaternion::from_axis_angle(&Vector3::y_axis(), 0.014);

    while window.render() {
        i += 1;
        if i % 1 == 0 {
            println!("{} iterations have passed.", i);
        }
        disp_cube(&mut window, vertices, Point3::new(0.0, 0.0, 0.0));
        window.set_point_size(10.);
        Context::get().point_size(10.);
        for particle in &map.particles {
            window.draw_point(
                &(particle.position * 10.).na_point(),
                &particle.fluid_type.simulation_color(),
            );
        }
        map.update((1100000000 as f64).recip());
    }
    // if i % 1 == 0 {
    //     println!(
    //         "{} seconds have passed, real time {}.",
    //         time_elapsed, real_time_elapsed
    //     );

    //     let real_delta_t = last_time.elapsed().as_secs_f64();
    //     real_time_elapsed += real_delta_t;
    //     time_elapsed += delta_t;
    //     // thread::sleep(Duration::from_secs_f64((delta_t - real_delta_t).max(0.)));
    //     last_time = Instant::now();
    // }
}
