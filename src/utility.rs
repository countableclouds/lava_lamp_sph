use std::f64::consts::PI;
const WATER_BENZYL_TENSION: f64 = 0.03;
const WATER_VISCOSITY: f64 = 1.12189;
const WATER_DIFFUSIVITY: f64 = 1.69e-7;
const BENZYL_VISCOSITY: f64 = 5.7684;
const BENZYL_DIFFUSIVITY: f64 = 6.912e-7;

pub trait Coords
where
    Self: std::marker::Sized,
{
    type Key;
    fn bin(&self, r: f64) -> Self::Key;
    fn along_axes(&self, r: f64) -> Vec<Self>;
    fn height(height: f64) -> Self;
    fn normalize(self) -> Self;
    fn control_update(
        &mut self,
        velocity: &mut Self,
        temperature: f64,
        dim: &Self,
        delta_t: f64,
    ) -> f64;
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
#[derive(Debug, PartialEq, Clone, Copy)]
pub struct Particle<T: Copy> {
    pub position: T,
    pub velocity: T,
    pub density: f64,
    pub mass: f64,
    pub temperature: f64,
    pub properties: Fluid,
}

impl<T: Copy + Coords + Default> Particle<T> {
    pub fn new(position: T, mass: f64, temperature: f64, properties: Fluid) -> Self {
        Particle {
            position,
            velocity: T::default(),
            density: 0.,
            mass,
            temperature,
            properties,
        }
    }
    pub fn with_density(mut self, density: f64) -> Self {
        self.density = density;
        self
    }
    pub fn with_temperature(mut self, temperature: f64) -> Self {
        self.temperature = temperature;
        self
    }
    pub fn with_velocity(mut self, velocity: T) -> Self {
        self.velocity = velocity;
        self
    }

    pub fn density_disparity(&self) -> f64 {
        self.density - self.properties.density(self.temperature)
    }
    pub fn volume(&self) -> f64 {
        self.mass * self.density.recip()
    }
    pub fn control_update(mut self, dim: T, delta_t: f64) -> Self {
        self.temperature =
            self.position
                .control_update(&mut self.velocity, self.temperature, &dim, delta_t);
        self
    }
}

#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}
impl Point {
    pub fn new(x: f64, y: f64) -> Point {
        Point { x, y }
    }
    pub fn squared_mag(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2)
    }
    pub fn mag(&self) -> f64 {
        self.squared_mag().sqrt()
    }
    pub fn area(&self) -> f64 {
        self.x * self.y
    }
    pub fn with_x(&self, x: f64) -> Point {
        Point { x, y: self.y }
    }
    pub fn with_y(&self, y: f64) -> Point {
        Point { x: self.x, y: y }
    }
}

impl Coords for Point {
    type Key = (u64, u64);
    fn bin(&self, r: f64) -> Self::Key {
        (self.clone() * r.recip()).into()
    }
    fn along_axes(&self, r: f64) -> Vec<Point> {
        vec![
            self.clone(),
            self.with_x(self.x + r),
            self.with_x(self.x - r),
            self.with_y(self.y + r),
            self.with_y(self.y - r),
        ]
    }
    fn height(height: f64) -> Point {
        Point::new(0., height)
    }

    fn normalize(self) -> Self {
        self / self.mag()
    }

    fn control_update(
        &mut self,
        velocity: &mut Point,
        temperature: f64,
        dim: &Self,
        delta_t: f64,
    ) -> f64 {
        self.x = self.x + velocity.x * delta_t;
        self.y = self.y + velocity.y * delta_t;
        if self.x > dim.x || self.x < 0. {
            self.x = self.x.min(dim.x).max(0.);
            velocity.x = -velocity.x;
        }
        if self.y > dim.y {
            self.y = dim.y;
            velocity.y = -velocity.y;
            return 333.15;
        }
        if self.y < 0. {
            self.y = 0.;
            velocity.y = 0.;
            return 293.15;
        }
        temperature
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

impl std::ops::AddAssign for Point {
    fn add_assign(&mut self, other: Point) {
        self.x += other.x;
        self.y += other.y;
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

impl Into<(u64, u64)> for Point {
    fn into(self) -> (u64, u64) {
        (self.x as u64, self.y as u64)
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

impl std::ops::Div<f64> for Point {
    type Output = Point;

    fn div(self, other: f64) -> Point {
        Point {
            x: self.x / other,
            y: self.y / other,
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
#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Fluid {
    Saltwater,
    BenzylAlcohol,
}

impl Fluid {
    /// Dynamic viscosity
    pub fn viscosity(&self) -> f64 {
        match self {
            Fluid::Saltwater => WATER_VISCOSITY,
            Fluid::BenzylAlcohol => BENZYL_VISCOSITY,
        }
    }
    pub fn diffusivity(&self) -> f64 {
        match self {
            Fluid::Saltwater => WATER_DIFFUSIVITY,
            Fluid::BenzylAlcohol => BENZYL_DIFFUSIVITY,
        }
    }
    pub fn interfacial_tension(&self, fluid: Fluid) -> f64 {
        match self {
            Fluid::Saltwater => match fluid {
                Fluid::Saltwater => 0.,
                Fluid::BenzylAlcohol => WATER_BENZYL_TENSION,
            },
            Fluid::BenzylAlcohol => match fluid {
                Fluid::Saltwater => WATER_BENZYL_TENSION,
                Fluid::BenzylAlcohol => 0.,
            },
        }
    }

    pub fn color(&self, fluid: Fluid) -> f64 {
        match self {
            Fluid::Saltwater => match fluid {
                Fluid::Saltwater => 0.,
                Fluid::BenzylAlcohol => 0.5,
            },
            Fluid::BenzylAlcohol => match fluid {
                Fluid::Saltwater => -0.5,
                Fluid::BenzylAlcohol => 0.,
            },
        }
    }

    pub fn density(&self, temperature: f64) -> f64 {
        match self {
            Fluid::Saltwater => 1034. / ((293.15 - temperature) * 3.09e-4 + 1.),
            Fluid::BenzylAlcohol => 1030. / ((313.15 - temperature) * 7.63e-4 + 1.),
        }
    }
}
