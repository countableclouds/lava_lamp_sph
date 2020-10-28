use std::{cmp::Eq, collections::HashMap, f64::consts::PI, hash::Hash};

#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub struct Particle<T: Copy> {
    position: T,
    velocity: T,
    density: f64,
    temperature: f64,
    properties: ParticleProperties,
}

#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub struct ParticleProperties {
    /// Dynamic viscosity
    viscosity: f64,
    /// thermal conductivity divided by specific heat, no density yet
    diffusivity: f64,
    surface_tension: f64,
    thermal_expansion: f64,
}

impl ParticleProperties {
    // /// 5% salt water
    // fn salt_water() -> Self {
    //     Self {viscosity: , diffusivity}
    // }
}
/// Kernel for general smoothing, h is the radius of the support
pub trait Poly6Kernel {
    fn poly6(&self, h: f64, point: Self) -> f64;
    fn grad_poly6(&self, h: f64, point: Self) -> Self;
    fn laplace_poly6(&self, h: f64, point: Self) -> f64;
}
/// Kernel for smoothing pressure, h is the radius of the support
pub trait SpikyKernel {
    fn spiky(&self, h: f64, point: Self) -> f64;
    fn grad_spiky(&self, h: f64, point: Self) -> Self;
    fn laplace_spiky(&self, h: f64, point: Self) -> f64;
}
/// Kernel for smoothing viscosity, h is the radius of the support
pub trait ViscosityKernel {
    fn visc(&self, h: f64, point: Self) -> f64;
    fn grad_visc(&self, h: f64, point: Self) -> Self;
    fn laplace_visc(&self, h: f64, point: Self) -> f64;
}

pub trait Coords
where
    Self: std::marker::Sized,
{
    type Key;
    fn bin(&self, r: f64) -> Self::Key;
    fn along_axes(&self, r: f64) -> Vec<Self>;
}

#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}
impl Point {
    fn new(x: f64, y: f64) -> Point {
        Point { x, y }
    }
    fn squared_mag(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2)
    }
    fn mag(&self) -> f64 {
        self.squared_mag().sqrt()
    }
    fn with_x(&self, x: f64) -> Point {
        Point { x, y: self.y }
    }
    fn with_y(&self, y: f64) -> Point {
        Point { x: self.x, y: y }
    }
}

impl Coords for Point {
    type Key = (f64, f64);
    fn bin(&self, r: f64) -> Self::Key {
        (self.clone() * r.recip()).into()
    }
    fn along_axes(&self, r: f64) -> Vec<Point> {
        vec![
            self.with_x(self.x + r),
            self.with_x(self.x - r),
            self.with_y(self.y + r),
            self.with_y(self.y - r),
        ]
    }
}

impl std::ops::Add for Point {
    type Output = Point;

    fn add(self, other: Point) -> Point {
        Point {
            x: self.x + other.x,
            y: self.y + other.y,
        }
    }
}

impl std::ops::Sub for Point {
    type Output = Point;

    fn sub(self, other: Point) -> Point {
        Point {
            x: self.x - other.x,
            y: self.y - other.y,
        }
    }
}

impl Into<(f64, f64)> for Point {
    fn into(self) -> (f64, f64) {
        (self.x, self.y)
    }
}

impl std::ops::Mul<f64> for Point {
    type Output = Point;

    fn mul(self, other: f64) -> Point {
        Point {
            x: self.x * other,
            y: self.y * other,
        }
    }
}

impl Poly6Kernel for Point {
    fn poly6(&self, h: f64, point: Point) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return 0.;
        }
        h.powi(-9) * squared_diff.powi(3) * 315. * (64. * PI).recip()
    }

    fn grad_poly6(&self, h: f64, point: Point) -> Point {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return Point::default();
        }
        dist * (-h.powi(-9) * squared_diff.powi(2) * 6. * 315. * (64. * PI).recip())
    }

    fn laplace_poly6(&self, h: f64, point: Point) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return 0.;
        }
        6. * squared_diff * (-squared_diff * 2. + 4. * dist.squared_mag())
    }
}

impl SpikyKernel for Point {
    fn spiky(&self, h: f64, point: Point) -> f64 {
        let dist = (self.clone() - point).mag();
        if dist > h {
            return 0.;
        }
        15. * PI.recip() * h.powi(-6) * dist.powi(3)
    }
    fn grad_spiky(&self, h: f64, point: Point) -> Point {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h {
            return Point::default();
        }
        dist_point * ((-3.) * (h - dist).powi(2) * dist.recip() * 15. * PI.recip() * h.powi(-6))
    }
    fn laplace_spiky(&self, h: f64, point: Point) -> f64 {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h {
            return 0.;
        }
        (-3.) * (h - dist).powi(2) * dist.recip() + 6. * (h - dist) * 15. * PI.recip() * h.powi(-6)
    }
}

impl ViscosityKernel for Point {
    fn visc(&self, h: f64, point: Self) -> f64 {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h {
            return 0.;
        }
        15. * (2. * PI).recip()
            * h.powi(-3)
            * (-dist.powi(3) * h.powi(-3) / 2.
                + dist.powi(2) * h.powi(-2)
                + h * (2. * dist).recip()
                - 1.)
    }
    fn grad_visc(&self, h: f64, point: Self) -> Point {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        dist_point * (-3. / 2. * dist * h.powi(-3) + 2. * h.powi(-2) - h * dist.powi(-3) / 2.)
    }
    fn laplace_visc(&self, h: f64, point: Self) -> f64 {
        let dist = (self.clone() - point).mag();
        if dist > h {
            return 0.;
        }
        45. * PI.recip() * h.powi(-6) * (h - dist)
    }
}

pub struct Map<T: Coords + Copy> {
    pub particles: [Particle<T>; 500],
    pub dim: T,
    pub radius: f64,
    pub particle_map: HashMap<T::Key, Vec<Particle<T>>>,
}

impl<T> Map<T>
where
    T: Copy + Coords + Poly6Kernel + SpikyKernel + ViscosityKernel + std::ops::Mul<f64>,
    T::Key: Hash + Eq,
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
}
