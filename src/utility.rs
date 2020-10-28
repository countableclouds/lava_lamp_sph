use std::f64::consts::PI;
/// Kernel for general smoothing
pub trait Poly6Kernel {
    fn poly6(&self, h: f64, point: Self) -> f64;
    fn grad_poly6(&self, h: f64, point: Self) -> Self;
    fn laplace_poly6(&self, h: f64, point: Self) -> f64;
}
/// Kernel for smoothing pressure
pub trait SpikyKernel {
    fn spiky(&self, h: f64, point: Self) -> f64;
    fn grad_spiky(&self, h: f64, point: Self) -> Self;
    fn laplace_spiky(&self, h: f64, point: Self) -> f64;
}
/// Kernel for smoothing viscosity
pub trait ViscosityKernel {
    fn visc(&self, h: f64, point: Self) -> f64;
    fn grad_visc(&self, h: f64, point: Self) -> Self;
    fn laplace_visc(&self, h: f64, point: Self) -> f64;
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
        h.powi(-9) * (h.powi(2) - dist.squared_mag()).powi(3) * 315. / (64. * PI)
    }

    fn grad_poly6(&self, h: f64, point: Point) -> Point {
        let dist = self.clone() - point;
        dist * (-h.powi(-9) * (h.powi(2) - dist.squared_mag()).powi(2) * 6. * 315. / (64. * PI))
    }

    fn laplace_poly6(&self, h: f64, point: Point) -> f64 {
        let dist = self.clone() - point;
        let squared_diff = h.powi(2) - dist.squared_mag();
        6. * squared_diff * (-squared_diff * 2. + 4. * dist.squared_mag())
    }
}

pub struct Map2D {
    dim: Point,
    points: Vec<Point>,
}

impl Map2D {}
