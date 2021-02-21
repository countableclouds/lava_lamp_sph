use na::{Point3, Translation3, UnitQuaternion, Vector3};
use std::f64::consts::PI;

pub trait Coords
where
    Self: std::marker::Sized,
{
    type Key;
    fn bin(&self, r: f64) -> Self::Key;
    fn along_axes(&self, r: f64) -> Vec<Self>;
    fn with_height(self, height: f64) -> Self;
    fn height(&self) -> f64;
    fn normalize(self) -> Self;
    fn control_update(&mut self, velocity: &mut Self, dim: &Self, delta_t: f64);
    fn mag(&self) -> f64;
    fn squared_mag(&self) -> f64;
    fn dot(&self, other: Self) -> f64;
    fn proj(&self) -> Vec<Self>;
    fn apply(&self, func: Box<dyn Fn(f64) -> f64>) -> Self;
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

pub trait QuinticSplineKernel {
    fn quintic_coeff(h: f64) -> f64;
    fn quintic(&self, h: f64, point: Self) -> f64;
    fn grad_quintic(&self, h: f64, point: Self) -> Self;
    fn laplace_quintic(&self, h: f64, point: Self) -> f64;
}

/// Kernel for smoothing viscosity, h is the radius of the support
pub trait CubicSplineKernel {
    fn cubic_coeff(h: f64) -> f64;
    fn cubic(&self, h: f64, point: Self) -> f64;
    fn grad_cubic(&self, h: f64, point: Self) -> Self;
    fn laplace_cubic(&self, h: f64, point: Self) -> f64;
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
    pub fn new(position: T, mass: f64, velocity: T, temperature: f64, fluid_type: Fluid) -> Self {
        Particle {
            position,
            velocity,
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
    pub fn expected_density(&self) -> f64 {
        self.fluid_type.density(self.temperature)
    }

    pub fn delta_density(&self) -> f64 {
        self.density - self.expected_density()
    }

    pub fn volume(&self) -> f64 {
        self.mass * self.density.recip()
    }
    pub fn control_update(mut self, dim: T, delta_t: f64) -> Self {
        // self.temperature =
        self.position
            .control_update(&mut self.velocity, &dim, delta_t);
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
    type Key = (i64, i64);
    fn squared_mag(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2)
    }
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
            self.with_x(self.x + r).with_y(self.y + r),
            self.with_x(self.x + r).with_y(self.y - r),
            self.with_x(self.x - r).with_y(self.y + r),
            self.with_x(self.x - r).with_y(self.y - r),
        ]
    }
    fn dot(&self, other: Point) -> f64 {
        other.x * self.x + other.y * self.y
    }
    fn with_height(mut self, height: f64) -> Point {
        self.y = height;
        self
    }

    fn height(&self) -> f64 {
        self.y
    }

    fn normalize(self) -> Self {
        self / self.mag()
    }

    fn control_update(&mut self, velocity: &mut Point, dim: &Self, delta_t: f64) {
        self.x = self.x + velocity.x * delta_t;
        self.y = self.y + velocity.y * delta_t;
        if self.x > dim.x || self.x < 0. {
            self.x = self.x.min(dim.x).max(0.);
            velocity.x = -velocity.x;
        }
        if self.y > dim.y {
            self.y = dim.y;
            velocity.y = -velocity.y;
        }
        if self.y < 0. {
            self.y = 0.;
            velocity.y = 0.;
        }
    }
    fn proj(&self) -> Vec<Self> {
        vec![Point::new(0., self.y), Point::new(self.x, 0.)]
    }
    fn apply(&self, func: Box<dyn Fn(f64) -> f64>) -> Self {
        Point::new(func(self.x), func(self.y))
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

impl Into<(i64, i64)> for Point {
    fn into(self) -> (i64, i64) {
        (self.x.floor() as i64, self.y.floor() as i64)
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
        if dist > h || dist == 0. {
            return 0.;
        }

        (-dist.powi(3) * h.powi(-3) / 2. + dist.powi(2) * h.powi(-2) + h * (2. * dist).recip() - 1.)
            * Self::visc_coeff(h)
    }
    fn grad_visc(&self, h: f64, point: Self) -> Point {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h || dist == 0. {
            return Point::default();
        }

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
    type Key = (i64, i64, i64);
    fn squared_mag(&self) -> f64 {
        self.x.powi(2) + self.y.powi(2) + self.z.powi(2)
    }
    fn mag(&self) -> f64 {
        self.squared_mag().sqrt()
    }

    fn bin(&self, r: f64) -> Self::Key {
        (self.clone() * r.recip()).into()
    }
    fn along_axes(&self, r: f64) -> Vec<Point3D> {
        let mut points: Vec<Point3D> = Vec::new();
        for i in -1..2 {
            for j in -1..2 {
                for k in -1..2 {
                    points.push(Point3D::new(
                        self.x + (i as f64) * r,
                        self.y + (j as f64) * r,
                        self.z + (k as f64) * r,
                    ))
                }
            }
        }
        points
    }
    fn dot(&self, other: Point3D) -> f64 {
        other.x * self.x + other.y * self.y + other.z * self.z
    }
    fn with_height(mut self, height: f64) -> Point3D {
        self.z = height;
        self
    }
    fn height(&self) -> f64 {
        self.z
    }

    fn normalize(self) -> Self {
        self / self.mag()
    }

    fn control_update(&mut self, velocity: &mut Point3D, dim: &Self, delta_t: f64) {
        let margin: f64 = 0.001;
        self.x = self.x + velocity.x * delta_t;
        self.y = self.y + velocity.y * delta_t;
        self.z = self.z + velocity.z * delta_t;
        // if self.x > dim.x - margin || self.x < margin {
        //     self.x = self.x.min(dim.x - margin).max(margin);
        //     velocity.x = 0.;
        // }
        // if self.y > dim.y - margin || self.y < margin {
        //     self.y = self.y.min(dim.y - margin).max(margin);
        //     velocity.y = 0.;
        // }
        // if self.z > dim.z - margin {
        //     self.z = dim.z - margin;
        //     velocity.z = 0.;
        // }
        // if self.z < margin {
        //     self.z = margin;
        //     velocity.z = 0.;
        // };
    }
    fn proj(&self) -> Vec<Self> {
        vec![
            Point3D::new(0., self.y, self.z),
            Point3D::new(self.x, 0., self.z),
            Point3D::new(self.x, self.y, 0.),
        ]
    }
    fn apply(&self, func: Box<dyn Fn(f64) -> f64>) -> Self {
        Point3D::new(func(self.x), func(self.y), func(self.z))
    }
}

impl std::fmt::Display for Point3D {
    // This trait requires `fmt` with this exact signature.
    fn fmt(&self, f: &mut std::fmt::Formatter) -> std::fmt::Result {
        // Write strictly the first element into the supplied output
        // stream: `f`. Returns `fmt::Result` which indicates whether the
        // operation succeeded or failed. Note that `write!` uses syntax which
        // is very similar to `println!`.
        write!(f, "{}, {}, {}", self.x, self.y, self.z)
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

impl Into<(i64, i64, i64)> for Point3D {
    fn into(self) -> (i64, i64, i64) {
        (
            self.x.floor() as i64,
            self.y.floor() as i64,
            self.z.floor() as i64,
        )
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

impl CubicSplineKernel for Point3D {
    fn cubic_coeff(h: f64) -> f64 {
        PI.recip() * h.powi(-3)
    }
    fn cubic(&self, h: f64, point: Point3D) -> f64 {
        let support = h / 2.;
        let dist = self.clone() - point;
        let q = dist.mag() / support;
        Self::cubic_coeff(support)
            * if q > 2. {
                0.
            } else if q > 1. {
                0.25 * (2. - q).powi(3)
            } else {
                1. - 1.5 * q.powi(2) + 0.75 * q.powi(3)
            }
    }

    fn grad_cubic(&self, h: f64, point: Point3D) -> Point3D {
        let support = h / 2.;
        let dist = self.clone() - point;
        let diff = dist.mag();
        let q = diff / support;
        if diff == 0. {
            return Point3D::default();
        }
        dist * Self::cubic_coeff(support)
            * (support * diff).recip()
            * if q > 2. {
                0.
            } else if q > 1. {
                -0.75 * (2. - q).powi(2)
            } else {
                -3. * q + 2.25 * q.powi(2)
            }
    }

    fn laplace_cubic(&self, h: f64, point: Point3D) -> f64 {
        let support = h / 2.;
        let dist = self.clone() - point;
        let diff = dist.mag();
        if diff == 0. {
            return 0.;
        }
        let q = diff / support;
        Self::cubic_coeff(support)
            * if q > 2. {
                0.
            } else if q > 1. {
                3. * (-2. * (diff * support).recip() + 3. * support.powi(-2)
                    - diff * support.powi(-3))
            } else {
                -9. * support.powi(-2) + 9. * diff * support.powi(-3)
            }
    }
}

impl QuinticSplineKernel for Point3D {
    fn quintic_coeff(h: f64) -> f64 {
        (120. * PI).recip() * h.powi(-3)
    }
    fn quintic(&self, h: f64, point: Point3D) -> f64 {
        let support = h / 3.;
        let dist = self.clone() - point;
        let q = dist.mag() / support;
        Self::quintic_coeff(support)
            * (if q < 3. { (3. - q).powi(5) } else { 0. }
                + if q < 2. { -6. * (2. - q).powi(5) } else { 0. }
                + if q < 1. { 15. * (1. - q).powi(5) } else { 0. })
    }

    fn grad_quintic(&self, h: f64, point: Point3D) -> Point3D {
        let support = h / 3.;
        let dist = self.clone() - point;
        let diff = dist.mag();
        if diff == 0. {
            return Point3D::default();
        }
        let q = diff / support;
        dist * (Self::quintic_coeff(support)
            * (support * diff).recip()
            * (if q < 3. { -5. * (3. - q).powi(4) } else { 0. }
                + if q < 2. { 30. * (2. - q).powi(4) } else { 0. }
                + if q < 1. { -75. * (1. - q).powi(4) } else { 0. }))
    }

    fn laplace_quintic(&self, h: f64, point: Point3D) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        if squared_diff < 0. {
            return 0.;
        }
        -6. * squared_diff * (squared_diff * 2. - 4. * dist.squared_mag()) * Self::quintic_coeff(h)
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
        if dist > h || dist == 0. {
            return 0.;
        }

        (-dist.powi(3) * h.powi(-3) / 2. + dist.powi(2) * h.powi(-2) + h * (2. * dist).recip() - 1.)
            * Self::visc_coeff(h)
    }
    fn grad_visc(&self, h: f64, point: Self) -> Point3D {
        let dist_point = self.clone() - point;
        let dist = dist_point.mag();
        if dist > h || dist == 0. {
            return Point3D::default();
        }
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
    const WATER_VISCOSITY: f64 = 1.12189e-3;
    const WATER_DIFFUSIVITY: f64 = 1.69e-7;
    const WATER_DENSITY: f64 = 1017.7;
    const WATER_THERMAL_EXPANSION: f64 = 3.09e-4;
    const BENZYL_VISCOSITY: f64 = 5.7684e-3;
    const BENZYL_DIFFUSIVITY: f64 = 6.912e-7;
    const BENZYL_DENSITY: f64 = 1023.75;
    const BENZYL_THERMAL_EXPANSION: f64 = 7.63e-4;
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

    pub fn disp(&self) -> &str {
        match self {
            Fluid::Saltwater => "Salt Water",
            Fluid::BenzylAlcohol => "Benzyl Alcohol",
        }
    }

    pub fn simulation_color(
        &self,
        temperature: f32,
        min_temperature: f32,
        max_temperature: f32,
    ) -> Point3<f32> {
        match self {
            Fluid::Saltwater => Point3::new(
                (temperature - min_temperature) / (max_temperature - min_temperature),
                1.,
                0.,
            ),
            Fluid::BenzylAlcohol => Point3::new(
                ((temperature - min_temperature) / (max_temperature - min_temperature))
                    * (255. / 255.),
                ((temperature - min_temperature) / (max_temperature - min_temperature))
                    * (0. / 255.),
                1. + ((temperature - min_temperature) / (max_temperature - min_temperature))
                    * (0. / 255. - 1.),
            ),
        }
    }

    pub fn density(&self, temperature: f64) -> f64 {
        let temperature = 293.15;
        match self {
            Fluid::Saltwater => {
                Fluid::WATER_DENSITY
                    / ((temperature - 288.55) * Fluid::WATER_THERMAL_EXPANSION + 1.)
            }
            Fluid::BenzylAlcohol => {
                Fluid::BENZYL_DENSITY
                    / ((temperature - 289.85) * Fluid::BENZYL_THERMAL_EXPANSION + 1.)
            }
        }
    }
}
impl std::fmt::Display for Fluid {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Fluid::Saltwater => write!(f, "Salt Water"),
            Fluid::BenzylAlcohol => write!(f, "Benzyl Alcohol"),
        }
    }
}
