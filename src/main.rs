// extern crate kiss3d;
// extern crate nalgebra as na;
// pub mod map;
// pub mod utility;
// use map::Map;
// use std::{
//     collections::HashMap,
//     convert::TryInto,
//     time::{Duration, Instant},
// };

// use utility::{Coords, Fluid, Particle, Point, Point3D};

// use kiss3d::{camera::FirstPerson, light::Light, scene::SceneNode, window::Window};
// use na::Point3;

// const ROOM_TEMPERATURE: f64 = 293.15;
// const MAP_SCALE: f32 = 1000.;
// const PARTICLE_UPPER_BOUND: u64 = 15000;
// const NUM_PARTICLES: usize = 12312;
// fn gen_points(dim: &Point, num_particles: u64) -> Vec<Point> {
//     let length = (-(dim.squared_mag() + dim.area() * (4. * (num_particles as f64) - 2.)).sqrt()
//         + dim.x
//         + dim.y)
//         / (2. - (2 * num_particles) as f64);
//     let width_particles = (dim.x / length - 1.) as u64;
//     let height_particles = (dim.y / length - 1.) as u64;
//     let width_offset = (dim.x - ((width_particles - 1) as f64 * length)) / 2.;
//     let height_offset = (dim.y - ((height_particles - 1) as f64 * length)) / 2.;
//     println!(
//         "Using {} particles, with {} as the initial maximum. Radius {}.",
//         width_particles * height_particles,
//         num_particles,
//         length
//     );
//     (0..(width_particles * height_particles))
//         .into_iter()
//         .map(|p| {
//             Point::new(
//                 (p % width_particles) as f64 * length + width_offset,
//                 (p / width_particles) as f64 * length + height_offset,
//             )
//         })
//         .collect()
// }

// fn gen_points_3d(dim: &Point3D, num_particles: u64) -> Vec<Point3D> {
//     let length = (dim.volume()
//         / (num_particles as f64 - dim.x * dim.y - dim.y * dim.z - dim.z * dim.x))
//         .cbrt();
//     let width_particles = (dim.x / length - 1.) as u64;
//     let length_particles = (dim.y / length - 1.) as u64;
//     let height_particles = (dim.z / length - 1.) as u64;
//     let width_offset = (dim.x - ((width_particles - 1) as f64 * length)) / 2.;
//     let length_offset = (dim.y - ((length_particles - 1) as f64 * length)) / 2.;
//     let height_offset = (dim.z - ((height_particles - 1) as f64 * length)) / 2.;
//     println!(
//         "Using {} particles, with {} as the initial maximum. Radius {}.",
//         width_particles * height_particles * length_particles,
//         num_particles,
//         length
//     );
//     (0..(width_particles * height_particles * length_particles))
//         .into_iter()
//         .map(|p| {
//             Point3D::new(
//                 (p % width_particles) as f64 * length + width_offset,
//                 // + (rng.gen::<f64>() - 0.5) * length / 2.,
//                 ((p % (width_particles * length_particles)) / width_particles as u64) as f64
//                     * length
//                     + length_offset,
//                 // + (rng.gen::<f64>() - 0.5) * length / 2.,
//                 (p / (width_particles * length_particles)) as f64 * length + height_offset, // + (rng.gen::<f64>() - 0.5) * length / 2.,
//             )
//         })
//         .collect()
// }

// fn disp_cube(window: &mut Window, vertices: [Point3<f32>; 8], color: Point3<f32>) {
//     window.draw_line(&vertices[0], &vertices[1], &color);
//     window.draw_line(&vertices[0], &vertices[2], &color);
//     window.draw_line(&vertices[0], &vertices[4], &color);
//     window.draw_line(&vertices[2], &vertices[3], &color);
//     window.draw_line(&vertices[4], &vertices[5], &color);
//     window.draw_line(&vertices[2], &vertices[6], &color);
//     window.draw_line(&vertices[1], &vertices[5], &color);
//     window.draw_line(&vertices[4], &vertices[6], &color);
//     window.draw_line(&vertices[1], &vertices[3], &color);
//     window.draw_line(&vertices[7], &vertices[3], &color);
//     window.draw_line(&vertices[7], &vertices[5], &color);
//     window.draw_line(&vertices[7], &vertices[6], &color);
// }
// fn main() {
//     let num_particles: u64 = PARTICLE_UPPER_BOUND;
//     let dim = Point3D::new(0.14605, 0.14605, 0.29845);
//     let points = gen_points_3d(&dim, num_particles);
//     let mut window = Window::new("Simulation! 😎");
//     let particles: [Particle<Point3D>; NUM_PARTICLES] = points
//         .into_iter()
//         .enumerate()
//         .map(|(i, position)| {
//             Particle::<Point3D>::new(
//                 position,
//                 dim.volume() / NUM_PARTICLES as f64
//                     * if i < (NUM_PARTICLES * 3 / 4) {
//                         Fluid::Saltwater.density(ROOM_TEMPERATURE)
//                     } else {
//                         Fluid::BenzylAlcohol.density(ROOM_TEMPERATURE)
//                     },
//                 ROOM_TEMPERATURE,
//                 if i < (NUM_PARTICLES * 3 / 4) {
//                     Fluid::Saltwater
//                 } else {
//                     Fluid::BenzylAlcohol
//                 },
//             )
//         })
//         .collect::<Vec<Particle<Point3D>>>()
//         .try_into()
//         .map_err(|v: Vec<Particle<Point3D>>| {
//             format!(
//                 "Expected a Vec of length {} instead of {}",
//                 NUM_PARTICLES,
//                 v.len()
//             )
//         })
//         .unwrap();

//     let mut map = Map {
//         particles,
//         dim,
//         radius: 0.015,
//         particle_map: HashMap::new(),
//         gravity: 0.,
//     };
//     let mut last_time: Instant = Instant::now();
//     let mut time_elapsed = 0.;
//     let mut real_time_elapsed = 0.;
//     let mut i = 0;

//     let vertices: [Point3<f32>; 8] = (dim * MAP_SCALE as f64).cube(Point3D::default());

//     window.set_light(Light::StickToCamera);

//     window.set_background_color(1.0, 1.0, 1.0);
//     let eye =
//         (dim * MAP_SCALE as f64 / 2. + Point3D::new(0.0, 1.0 * MAP_SCALE as f64, 0.0)).na_point();
//     let mut first_person = FirstPerson::new(eye, (dim * MAP_SCALE as f64 / 2.).na_point());
//     let mut cubes: Vec<SceneNode> = vec![];
//     let delta_t = (11 as f64).recip();
//     while window.render_with_camera(&mut first_person) {
//         // window.scene_mut().unlink();
//         for cube in &mut cubes {
//             window.remove_node(cube);
//         }
//         cubes = vec![];
//         for (j, particle) in (&map.particles).iter().enumerate() {
//             cubes.push(window.add_cube(0.001 * MAP_SCALE, 0.001 * MAP_SCALE, 0.001 * MAP_SCALE));
//             let cube = cubes.get_mut(j).unwrap();
//             let color = particle.fluid_type.simulation_color();
//             cube.set_color(color.x, color.y, color.z);
//             cube.append_translation(&(particle.position * MAP_SCALE as f64).translation())
//         }

//         i += 1;
//         if i % 1 == 0 {
//             println!(
//                 "{} seconds have passed, real time {}.",
//                 time_elapsed, real_time_elapsed
//             );
//         }
//         disp_cube(&mut window, vertices, Point3::new(0.0, 0.0, 0.0));

//         if i == 1 || false {
//             map.update(delta_t)
//         };

//         // assert!(1 == 0);
//         let real_delta_t = last_time.elapsed().as_secs_f64();
//         real_time_elapsed += real_delta_t;
//         time_elapsed += delta_t;
//         // thread::sleep(Duration::from_secs_f64((delta_t - real_delta_t).max(0.)));
//         last_time = Instant::now();
//     }
// }

extern crate kiss3d;
extern crate nalgebra as na;

use kiss3d::camera::Camera;
use kiss3d::context::Context;
use kiss3d::planar_camera::PlanarCamera;
use kiss3d::post_processing::PostProcessingEffect;
use kiss3d::renderer::Renderer;
use kiss3d::resource::{
    AllocationType, BufferType, Effect, GPUVec, ShaderAttribute, ShaderUniform,
};
use kiss3d::text::Font;
use kiss3d::window::{State, Window};
use na::{Matrix4, Point2, Point3, Vector3};

// Custom renderers are used to allow rendering objects that are not necessarily
// represented as meshes. In this example, we will render a large, growing, point cloud
// with a color associated to each point.

// Writing a custom renderer requires the main loop to be
// handled by the `State` trait instead of a `while window.render()`
// like other examples.

struct AppState {
    point_cloud_renderer: PointCloudRenderer,
}

impl State for AppState {
    // Return the custom renderer that will be called at each
    // render loop.
    fn cameras_and_effect_and_renderer(
        &mut self,
    ) -> (
        Option<&mut dyn Camera>,
        Option<&mut dyn PlanarCamera>,
        Option<&mut dyn Renderer>,
        Option<&mut dyn PostProcessingEffect>,
    ) {
        (None, None, Some(&mut self.point_cloud_renderer), None)
    }

    fn step(&mut self, window: &mut Window) {
        if self.point_cloud_renderer.num_points() < 1_000_000 {
            // Add some random points to the point cloud.
            for _ in 0..1_000 {
                let random: Point3<f32> = rand::random();
                self.point_cloud_renderer
                    .push((random - Vector3::repeat(0.5)) * 0.5, rand::random());
            }
        }

        let num_points_text = format!(
            "Number of points: {}",
            self.point_cloud_renderer.num_points()
        );
        window.draw_text(
            &num_points_text,
            &Point2::new(0.0, 20.0),
            60.0,
            &Font::default(),
            &Point3::new(1.0, 1.0, 1.0),
        );
    }
}

fn main() {
    let window = Window::new("Kiss3d: persistent_point_cloud");
    let app = AppState {
        point_cloud_renderer: PointCloudRenderer::new(4000.0),
    };

    window.render_loop(app)
}

/// Structure which manages the display of long-living points.
struct PointCloudRenderer {
    shader: Effect,
    pos: ShaderAttribute<Point3<f32>>,
    color: ShaderAttribute<Point3<f32>>,
    proj: ShaderUniform<Matrix4<f32>>,
    view: ShaderUniform<Matrix4<f32>>,
    colored_points: GPUVec<Point3<f32>>,
    point_size: f32,
}

impl PointCloudRenderer {
    /// Creates a new points renderer.
    fn new(point_size: f32) -> PointCloudRenderer {
        let mut shader = Effect::new_from_str(VERTEX_SHADER_SRC, FRAGMENT_SHADER_SRC);

        shader.use_program();

        PointCloudRenderer {
            colored_points: GPUVec::new(Vec::new(), BufferType::Array, AllocationType::StreamDraw),
            pos: shader.get_attrib::<Point3<f32>>("position").unwrap(),
            color: shader.get_attrib::<Point3<f32>>("color").unwrap(),
            proj: shader.get_uniform::<Matrix4<f32>>("proj").unwrap(),
            view: shader.get_uniform::<Matrix4<f32>>("view").unwrap(),
            shader,
            point_size,
        }
    }

    fn push(&mut self, point: Point3<f32>, color: Point3<f32>) {
        if let Some(colored_points) = self.colored_points.data_mut() {
            colored_points.push(point);
            colored_points.push(color);
        }
    }

    fn num_points(&self) -> usize {
        self.colored_points.len() / 2
    }
}

impl Renderer for PointCloudRenderer {
    /// Actually draws the points.
    fn render(&mut self, pass: usize, camera: &mut dyn Camera) {
        if self.colored_points.len() == 0 {
            return;
        }

        self.shader.use_program();
        self.pos.enable();
        self.color.enable();

        camera.upload(pass, &mut self.proj, &mut self.view);

        self.color.bind_sub_buffer(&mut self.colored_points, 1, 1);
        self.pos.bind_sub_buffer(&mut self.colored_points, 1, 0);

        let ctxt = Context::get();
        ctxt.point_size(self.point_size);
        ctxt.draw_arrays(Context::POINTS, 0, (self.colored_points.len() / 2) as i32);

        self.pos.disable();
        self.color.disable();
    }
}

const VERTEX_SHADER_SRC: &'static str = "#version 100
    attribute vec3 position;
    attribute vec3 color;
    varying   vec3 Color;
    uniform   mat4 proj;
    uniform   mat4 view;
    void main() {
        gl_Position = proj * view * vec4(position, 1.0);
        Color = color;
    }";

const FRAGMENT_SHADER_SRC: &'static str = "#version 100
#ifdef GL_FRAGMENT_PRECISION_HIGH
   precision highp float;
#else
   precision mediump float;
#endif
    varying vec3 Color;
    void main() {
        gl_FragColor = vec4(Color, 1.0);
    }";
