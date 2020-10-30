pub mod map;
pub mod utility;
use map::Map;
use std::{collections::HashMap, convert::TryInto};
use utility::{Fluid, Particle, Point};
const ROOM_TEMPERATURE: f64 = 293.15;

fn gen_points(dim: &Point, num_particles: u64) -> Vec<Point> {
    let length = (-(dim.squared_mag() + dim.area() * (4. * (num_particles as f64) - 2.)).sqrt()
        + dim.x
        + dim.y)
        / (2. - (2 * num_particles) as f64);
    println!("{}", length);
    let width_particles = (dim.x / length - 1.) as u64;
    let height_particles = (dim.y / length - 1.) as u64;
    let width_offset = (dim.x - (width_particles as f64 * length)) / 2.;
    let height_offset = (dim.y - (height_particles as f64 * length)) / 2.;
    println!(
        "Using {} particles, with {} as the initial maximum.",
        width_particles * height_particles,
        num_particles
    );
    (0..(width_particles * height_particles))
        .into_iter()
        .map(|p| {
            Point::new(
                ((p + 1) % width_particles as u64) as f64 * length + width_offset,
                ((p + 1) / width_particles as u64) as f64 * length + height_offset,
            )
        })
        .collect()
}
fn main() {
    let num_particles: u64 = 11000;
    let dim = Point::new(0.14605, 0.29845);
    let points = gen_points(&dim, num_particles);
    let particles: [Particle<Point>; 10950] = match points
        .into_iter()
        .enumerate()
        .map(|(i, position)| {
            Particle::<Point>::new(
                position,
                dim.area() / num_particles as f64
                    * if i < (num_particles * 3 / 4) as usize {
                        Fluid::Saltwater.density(ROOM_TEMPERATURE) * 3. / 4.
                    } else {
                        Fluid::BenzylAlcohol.density(ROOM_TEMPERATURE) / 4.
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
    {
        Ok(arr) => arr,
        Err(_) => panic!("Expected a Vec of a different length",),
    };
    let mut map = Map {
        particles,
        dim,
        radius: 0.,
        particle_map: HashMap::new(),
        gravity: 9.8,
        gas_constant: 8.,
    };
    for i in 0..40000 {
        if i % 100 == 0 {
            println!("{} iterations have passed.", i);
        }
        map.update((11 as f64).recip());
    }
}
