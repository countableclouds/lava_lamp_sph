use na::{Point3, Translation3, UnitQuaternion, Vector3};
use std::f64::consts::PI;

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
    fn mag(&self) -> f64;
}
/// Kernel for general smoothing, h is the radius of the support
pub trait Poly6Kernel {
    fn poly6_coeff(h: f64) -> f64;
    fn poly6(&self, h: f64, point: Self) -> f64;
    fn grad_poly6(&self, h: f64, point: Self) -> Self;
    fn laplace_poly6(&self, h: f64, point: Self) -> f64;
}
/// Kernel for smoothing pressure, h is the radius of the support
pub trait SpikyKernel {
    fn spiky_coeff(h: f64) -> f64;
    fn spiky(&self, h: f64, point: Self) -> f64;
    fn grad_spiky(&self, h: f64, point: Self) -> Self;
    fn laplace_spiky(&self, h: f64, point: Self) -> f64;
}
/// Kernel for smoothing viscosity, h is the radius of the support
pub trait ViscosityKernel {
    fn visc_coeff(h: f64) -> f64;
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
    pub fluid_type: Fluid,
}

impl<T: Copy + Coords + Default> Particle<T> {
    pub fn new(position: T, mass: f64, temperature: f64, fluid_type: Fluid) -> Self {
        Particle {
            position,
            velocity: T::default(),
            density: 0.,
            mass,
            temperature,
            fluid_type,
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

    pub fn gas_coefficient(&self) -> f64 {
        self.fluid_type.molar_mass().recip()
            * self.temperature
            * (self.density - self.fluid_type.density(self.temperature))
    }

    pub fn volume(&self) -> f64 {
        self.mass * self.density.recip()
    }
    pub fn control_update(mut self, dim: T, delta_t: f64) -> Self {
        // self.temperature =
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
    fn mag(&self) -> f64 {
        self.squared_mag().sqrt()
    }

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
            return 293.15;
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

impl Into<[f64; 2]> for Point {
    fn into(self) -> [f64; 2] {
        [self.x, self.y]
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
    fn poly6_coeff(h: f64) -> f64 {
        4. / PI * h.powi(-8)
    }
    fn poly6(&self, h: f64, point: Point) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return 0.;
        }
        squared_diff.powi(3) * Self::poly6_coeff(h)
    }

    fn grad_poly6(&self, h: f64, point: Point) -> Point {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return Point::default();
        }
        dist * (squared_diff.powi(2) * 6. * Self::poly6_coeff(h))
    }

    fn laplace_poly6(&self, h: f64, point: Point) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return 0.;
        }
        -6. * squared_diff * (squared_diff * 2. - 4. * dist.squared_mag()) * Self::poly6_coeff(h)
    }
}

impl SpikyKernel for Point {
    fn spiky_coeff(h: f64) -> f64 {
        10. / PI * h.powi(-5)
    }
    fn spiky(&self, h: f64, point: Point) -> f64 {
        let dist = (self.clone() - point).mag();
        if dist > h {
            return 0.;
        }
        dist.powi(3) * Self::spiky_coeff(h)
    }
    fn grad_spiky(&self, h: f64, point: Point) -> Point {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h || dist == 0. {
            return Point::default();
        }
        dist_point * ((-3.) * (h - dist).powi(2) * dist.recip() * Self::spiky_coeff(h))
    }
    fn laplace_spiky(&self, h: f64, point: Point) -> f64 {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h || dist == 0. {
            return 0.;
        }
        (-3. * h.powi(2) * dist.recip() + 12. * h - 9. * dist) * Self::spiky_coeff(h)
    }
}

impl ViscosityKernel for Point {
    fn visc_coeff(h: f64) -> f64 {
        10. * (3. * PI * h.powi(2)).recip()
    }
    fn visc(&self, h: f64, point: Self) -> f64 {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h {
            return 0.;
        }

        (-dist.powi(3) * h.powi(-3) / 2. + dist.powi(2) * h.powi(-2) + h * (2. * dist).recip() - 1.)
            * Self::visc_coeff(h)
    }
    fn grad_visc(&self, h: f64, point: Self) -> Point {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        assert!(
            dist != 0.,
            "Gradient of viscosity function produced infinity."
        );
        dist_point
            * ((-3. / 2. * dist * h.powi(-3) + 2. * h.powi(-2) - h * dist.powi(-3) / 2.)
                * Self::visc_coeff(h))
    }
    fn laplace_visc(&self, h: f64, point: Self) -> f64 {
        let dist = (self.clone() - point).mag();
        if dist > h {
            return 0.;
        }
        (h - dist) * h.powi(-3) * (6. * Self::visc_coeff(h))
    }
}

#[derive(Default, Debug, PartialEq, Clone, Copy)]
pub struct Point3D {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}
impl Point3D {
    pub fn new(x: f64, y: f64, z: f64) -> Point3D {
        Point3D { x, y, z }
    }
    pub fn squared_mag(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }

    pub fn volume(&self) -> f64 {
        self.x * self.y * self.z
    }
    pub fn with_x(&self, x: f64) -> Point3D {
        Point3D {
            x,
            y: self.y,
            z: self.z,
        }
    }
    pub fn with_y(&self, y: f64) -> Point3D {
        Point3D {
            x: self.x,
            y,
            z: self.z,
        }
    }
    pub fn with_z(&self, z: f64) -> Point3D {
        Point3D {
            x: self.x,
            y: self.y,
            z,
        }
    }
    pub fn cube(&self, other: Point3D) -> [Point3<f32>; 8] {
        [
            self.na_point(),
            Point3D::new(self.x, self.y, other.z).na_point(),
            Point3D::new(self.x, other.y, self.z).na_point(),
            Point3D::new(self.x, other.y, other.z).na_point(),
            Point3D::new(other.x, self.y, self.z).na_point(),
            Point3D::new(other.x, self.y, other.z).na_point(),
            Point3D::new(other.x, other.y, self.z).na_point(),
            Point3D::new(other.x, other.y, other.z).na_point(),
        ]
    }
    pub fn na_point(&self) -> Point3<f32> {
        Point3::<f32>::new(self.x as f32, self.y as f32, self.z as f32)
    }
    pub fn translation(&self) -> Translation3<f32> {
        Translation3::new(self.x as f32, self.y as f32, self.z as f32)
    }
}

impl Coords for Point3D {
    type Key = (u64, u64, u64);
    fn mag(&self) -> f64 {
        self.squared_mag().sqrt()
    }

    fn bin(&self, r: f64) -> Self::Key {
        (self.clone() * r.recip()).into()
    }
    fn along_axes(&self, r: f64) -> Vec<Point3D> {
        vec![
            self.clone(),
            self.with_x(self.x + r),
            self.with_x(self.x - r),
            self.with_y(self.y + r),
            self.with_y(self.y - r),
            self.with_z(self.z + r),
            self.with_z(self.z - r),
        ]
    }
    fn height(height: f64) -> Self {
        Self::new(0., 0., height)
    }

    fn normalize(self) -> Self {
        self / self.mag()
    }

    fn control_update(
        &mut self,
        velocity: &mut Point3D,
        temperature: f64,
        dim: &Self,
        delta_t: f64,
    ) -> f64 {
        self.x = self.x + velocity.x * delta_t;
        self.y = self.y + velocity.y * delta_t;
        self.z = self.z + velocity.z * delta_t;
        if self.x > dim.x || self.x < 0. {
            self.x = self.x.min(dim.x).max(0.);
            velocity.x = -velocity.x;
        }
        if self.y > dim.y || self.y < 0. {
            self.y = self.y.min(dim.y).max(0.);
            velocity.y = -velocity.y;
        }
        if self.z > dim.z {
            self.z = dim.z;
            velocity.z = -velocity.z;
            return 293.15;
        }
        if self.z < 0. {
            self.z = 0.;
            velocity.z = 0.;
            return 293.15;
        }
        temperature
    }
}

impl std::ops::Add for Point3D {
    type Output = Point3D;

    fn add(self, other: Point3D) -> Point3D {
        Point3D {
            x: self.x + other.x,
            y: self.y + other.y,
            z: self.z + other.z,
        }
    }
}

impl std::ops::AddAssign for Point3D {
    fn add_assign(&mut self, other: Point3D) {
        self.x += other.x;
        self.y += other.y;
        self.z += other.z;
    }
}

impl std::ops::Sub for Point3D {
    type Output = Point3D;

    fn sub(self, other: Point3D) -> Point3D {
        Point3D {
            x: self.x - other.x,
            y: self.y - other.y,
            z: self.z - other.z,
        }
    }
}

impl Into<(u64, u64, u64)> for Point3D {
    fn into(self) -> (u64, u64, u64) {
        (self.x as u64, self.y as u64, self.z as u64)
    }
}

impl Into<[f64; 3]> for Point3D {
    fn into(self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

impl From<(f64, f64, f64)> for Point3D {
    fn from(other: (f64, f64, f64)) -> Self {
        Point3D::new(other.0, other.1, other.2)
    }
}

impl std::ops::Mul<f64> for Point3D {
    type Output = Point3D;

    fn mul(self, other: f64) -> Point3D {
        Point3D {
            x: self.x * other,
            y: self.y * other,
            z: self.z * other,
        }
    }
}

impl std::ops::Div<f64> for Point3D {
    type Output = Point3D;

    fn div(self, other: f64) -> Point3D {
        Point3D {
            x: self.x / other,
            y: self.y / other,
            z: self.z / other,
        }
    }
}

impl Poly6Kernel for Point3D {
    fn poly6_coeff(h: f64) -> f64 {
        315. / 64. / PI * h.powi(-9)
    }
    fn poly6(&self, h: f64, point: Point3D) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return 0.;
        }
        squared_diff.powi(3) * Self::poly6_coeff(h)
    }

    fn grad_poly6(&self, h: f64, point: Point3D) -> Point3D {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return Point3D::default();
        }
        dist * (squared_diff.powi(2) * 6. * Self::poly6_coeff(h))
    }

    fn laplace_poly6(&self, h: f64, point: Point3D) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return 0.;
        }
        -6. * squared_diff * (squared_diff * 2. - 4. * dist.squared_mag()) * Self::poly6_coeff(h)
    }
}

impl SpikyKernel for Point3D {
    fn spiky_coeff(h: f64) -> f64 {
        15. / PI * h.powi(-6)
    }
    fn spiky(&self, h: f64, point: Point3D) -> f64 {
        let dist = (self.clone() - point).mag();
        if dist > h {
            return 0.;
        }
        dist.powi(3) * Self::spiky_coeff(h)
    }
    fn grad_spiky(&self, h: f64, point: Point3D) -> Point3D {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h || dist == 0. {
            return Point3D::default();
        }
        dist_point * ((-3.) * (h - dist).powi(2) * dist.recip() * Self::spiky_coeff(h))
    }
    fn laplace_spiky(&self, h: f64, point: Point3D) -> f64 {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h || dist == 0. {
            return 0.;
        }
        (-3. * h.powi(2) * dist.recip() + 12. * h - 9. * dist) * Self::spiky_coeff(h)
    }
}

impl ViscosityKernel for Point3D {
    fn visc_coeff(h: f64) -> f64 {
        15. / 2. / PI * h.powi(-3)
    }
    fn visc(&self, h: f64, point: Self) -> f64 {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h {
            return 0.;
        }

        (-dist.powi(3) * h.powi(-3) / 2. + dist.powi(2) * h.powi(-2) + h * (2. * dist).recip() - 1.)
            * Self::visc_coeff(h)
    }
    fn grad_visc(&self, h: f64, point: Self) -> Point3D {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        assert!(
            dist != 0.,
            "Gradient of viscosity function produced infinity."
        );
        dist_point
            * ((-3. / 2. * dist * h.powi(-3) + 2. * h.powi(-2) - h * dist.powi(-3) / 2.)
                * Self::visc_coeff(h))
    }
    fn laplace_visc(&self, h: f64, point: Self) -> f64 {
        let dist = (self.clone() - point).mag();
        if dist > h {
            return 0.;
        }
        (h - dist) * h.powi(-3) * (6. * Self::visc_coeff(h))
    }
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Fluid {
    Saltwater,
    BenzylAlcohol,
}

impl Fluid {
    const WATER_BENZYL_TENSION: f64 = 0.03297;
    const WATER_VISCOSITY: f64 = 1.12189;
    const WATER_DIFFUSIVITY: f64 = 1.69e-7;
    const WATER_DENSITY: f64 = 1034.;
    const WATER_THERMAL_EXPANSION: f64 = 3.09e-4;
    const WATER_MOLAR_MASS: f64 = 1.8015e-2;
    const BENZYL_VISCOSITY: f64 = 5.7684;
    const BENZYL_DIFFUSIVITY: f64 = 6.912e-7;
    const BENZYL_DENSITY: f64 = 1030.;
    const BENZYL_THERMAL_EXPANSION: f64 = 7.63e-4;
    const BENZYL_MOLAR_MASS: f64 = 1.0814e-1;
    /// Dynamic viscosity
    pub fn viscosity(&self) -> f64 {
        match self {
            Fluid::Saltwater => Fluid::WATER_VISCOSITY,
            Fluid::BenzylAlcohol => Fluid::BENZYL_VISCOSITY,
        }
    }
    pub fn diffusivity(&self) -> f64 {
        match self {
            Fluid::Saltwater => Fluid::WATER_DIFFUSIVITY,
            Fluid::BenzylAlcohol => Fluid::BENZYL_DIFFUSIVITY,
        }
    }
    // molar mass in kg/mol
    pub fn molar_mass(&self) -> f64 {
        match self {
            Fluid::Saltwater => Fluid::WATER_MOLAR_MASS,
            Fluid::BenzylAlcohol => Fluid::BENZYL_MOLAR_MASS,
        }
    }
    pub fn interfacial_tension(&self, fluid: Fluid) -> f64 {
        match self {
            Fluid::Saltwater => match fluid {
                Fluid::Saltwater => 0.,
                Fluid::BenzylAlcohol => Fluid::WATER_BENZYL_TENSION,
            },
            Fluid::BenzylAlcohol => match fluid {
                Fluid::Saltwater => Fluid::WATER_BENZYL_TENSION,
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

    pub fn simulation_color(&self) -> Point3<f32> {
        match self {
            Fluid::Saltwater => Point3::new(1., 0., 0.),
            Fluid::BenzylAlcohol => Point3::new(0., 1., 1.),
        }
    }

    pub fn density(&self, temperature: f64) -> f64 {
        match self {
            Fluid::Saltwater => {
                Fluid::WATER_DENSITY
                    / ((temperature - 293.15) * Fluid::WATER_THERMAL_EXPANSION + 1.)
            }
            Fluid::BenzylAlcohol => {
                Fluid::BENZYL_DENSITY
                    / ((temperature - 313.15) * Fluid::BENZYL_THERMAL_EXPANSION + 1.)
            }
        }
    }
}
