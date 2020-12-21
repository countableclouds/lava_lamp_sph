use std::collections::{HashMap, HashSet};

use conrod::{
    widget::{Canvas, Slider},
    Scalar,
};
use kiss3d::widget_ids;

use kiss3d::conrod::color::{Color, Colorable};
use kiss3d::conrod::position::{Positionable, Sizeable};
use kiss3d::conrod::widget::{button::Style, range_slider::RangeSlider, Button, Text, Widget};
use kiss3d::conrod::Labelable;
use petgraph::{
    graph::{EdgeIndex, NodeIndex},
    visit::{EdgeRef, IntoEdgeReferences, IntoNodeReferences},
};
use wasm_bindgen::prelude::*;

use kiss3d::conrod;
use kiss3d::light::Light;
use kiss3d::scene::SceneNode;
use kiss3d::window::{State, Window};
use nalgebra::{distance, Point3, Rotation3, Translation3, UnitQuaternion, Vector3, VectorN};

use connectome_model::{
    sim::Simulation,
    simplex::{faces, SimplicialComplex},
};
use rand::{distributions::Uniform, rngs::ThreadRng, seq::IteratorRandom, Rng};

struct AppState {
    ids: Ids,
    app: DemoApp,
}

impl AppState {
    fn new(window: &mut Window, ids: Ids, app: DemoApp) -> Self {
        Self { ids, app }
    }
}

impl State for AppState {
    fn step(&mut self, window: &mut Window) {
        self.app.step(window);

        let mut ui = window.conrod_ui_mut().set_widgets();
        self.app.gui(&mut ui, &self.ids);
    }
}

#[wasm_bindgen]
pub fn init_panic_hook() {
    console_error_panic_hook::set_once();
}

#[wasm_bindgen(start)]
pub fn main() -> Result<(), JsValue> {
    init_panic_hook();

    let mut window = Window::new("Synaptogenesis model");
    window.set_background_color(1.0, 1.0, 1.0);

    window.set_light(Light::StickToCamera);

    // Generate the widget identifiers.
    let ids = Ids::new(window.conrod_ui_mut().widget_id_generator());
    let app = DemoApp::new();
    window.conrod_ui_mut().theme = theme();

    let state = AppState::new(&mut window, ids, app);

    window.render_loop(state);
    Ok(())
}

pub fn theme() -> conrod::Theme {
    use conrod::position::{Align, Direction, Padding, Position, Relative};
    conrod::Theme {
        name: "Demo Theme".to_string(),
        padding: Padding::none(),
        x_position: Position::Relative(Relative::Align(Align::Start), None),
        y_position: Position::Relative(Relative::Direction(Direction::Backwards, 20.0), None),
        background_color: conrod::color::DARK_CHARCOAL,
        shape_color: conrod::color::LIGHT_CHARCOAL,
        border_color: conrod::color::BLACK,
        border_width: 0.0,
        label_color: conrod::color::WHITE,
        font_id: None,
        font_size_large: 26,
        font_size_medium: 18,
        font_size_small: 12,
        widget_styling: conrod::theme::StyleMap::default(),
        mouse_drag_threshold: 0.0,
        double_click_threshold: std::time::Duration::from_millis(500),
    }
}

// Generate a unique `WidgetId` for each widget.
widget_ids! {
    pub struct Ids {
        canvas,
        connectivity_text,
        connectivity_slider,
        myelination_text,
        myelination_slider,
        decay_text,
        decay_slider,
        max_myelination_text,
        max_myelination_slider,
        distance_exp_text,
        distance_exp_slider,
        refractory_period_text,
        refractory_period_slider,
        output_text,
    }
}

pub const WIN_W: u32 = 600;
pub const WIN_H: u32 = 420;

pub struct DemoApp {
    sim_state: Option<SimState>,
    sim_params: SimParams,
    sim_invalidated: bool,
}

impl DemoApp {
    pub fn new() -> Self {
        DemoApp {
            sim_state: Default::default(),
            sim_params: Default::default(),
            sim_invalidated: Default::default(),
        }
    }

    fn gui(&mut self, ui: &mut conrod::UiCell, ids: &Ids) {
        const PADDING: Scalar = 20.0;
        const VERTICAL_SPACING: Scalar = 20.0;

        Canvas::new()
            .pad(PADDING)
            .align_top()
            .align_left()
            .w(ui.win_w * 0.2)
            .h(ui.win_h * 0.5)
            .scroll_kids_vertically()
            .color(Color::Rgba(0.0, 0.0, 0.0, 0.5))
            .set(ids.canvas, ui);

        Text::new("connectivity rate")
            .parent(ids.canvas)
            .align_top()
            .align_middle_x()
            .padded_w_of(ids.canvas, PADDING)
            .set(ids.connectivity_text, ui);

        for value in Slider::new(self.sim_params.connectivity_rate, 0.0, 1.0)
            .h(20.0)
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.connectivity_text, VERTICAL_SPACING)
            .set(ids.connectivity_slider, ui)
        {
            self.sim_params.connectivity_rate = value;
            self.sim_invalidated = true;
        }

        Text::new("myelination rate")
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.connectivity_slider, VERTICAL_SPACING)
            .set(ids.myelination_text, ui);

        for value in Slider::new(self.sim_params.myelination_rate, 0.0, 1.0)
            .h(20.0)
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.myelination_text, VERTICAL_SPACING)
            .set(ids.myelination_slider, ui)
        {
            self.sim_params.myelination_rate = value;
            self.sim_invalidated = true;
        }

        Text::new("myelination decay rate")
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.myelination_slider, VERTICAL_SPACING)
            .set(ids.decay_text, ui);

        for value in Slider::new(self.sim_params.decay_rate, 0.0, 0.1)
            .h(20.0)
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.decay_text, VERTICAL_SPACING)
            .set(ids.decay_slider, ui)
        {
            self.sim_params.decay_rate = value;
            self.sim_invalidated = true;
        }

        Text::new("max myelination")
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.decay_slider, VERTICAL_SPACING)
            .set(ids.max_myelination_text, ui);

        for value in Slider::new(self.sim_params.max_myelination as f32, 0.0, 20.0)
            .h(20.0)
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.max_myelination_text, VERTICAL_SPACING)
            .set(ids.max_myelination_slider, ui)
        {
            self.sim_params.max_myelination = value as usize;
            self.sim_invalidated = true;
        }

        Text::new("distance exponent")
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.max_myelination_slider, VERTICAL_SPACING)
            .set(ids.distance_exp_text, ui);

        for value in Slider::new(self.sim_params.distance_exp as f32, 0.0, 10.0)
            .h(20.0)
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.distance_exp_text, VERTICAL_SPACING)
            .set(ids.distance_exp_slider, ui)
        {
            self.sim_params.distance_exp = value as i32;
            self.sim_invalidated = true;
        }

        Text::new("refractory period")
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.distance_exp_slider, VERTICAL_SPACING)
            .set(ids.refractory_period_text, ui);

        for value in Slider::new(self.sim_params.refractory_period as f32, 0.0, 10.0)
            .h(20.0)
            .parent(ids.canvas)
            .padded_w_of(ids.canvas, PADDING)
            .down_from(ids.refractory_period_text, VERTICAL_SPACING)
            .set(ids.refractory_period_slider, ui)
        {
            self.sim_params.refractory_period = value as usize;
            self.sim_invalidated = true;
        }

        Text::new(&format!(
            "simplex sizes: {:?}\nbetti numbers: {:?}",
            self.sim_state.as_ref().unwrap().cached_outputs.0,
            self.sim_state.as_ref().unwrap().cached_outputs.1,
        ))
        .parent(ids.canvas)
        .padded_w_of(ids.canvas, PADDING)
        .down_from(ids.refractory_period_slider, VERTICAL_SPACING)
        .set(ids.output_text, ui);
    }
}

impl State for DemoApp {
    fn step(&mut self, window: &mut Window) {
        if self.sim_invalidated {
            self.sim_state.take().map(|mut state| state.cleanup(window));
            self.sim_invalidated = false;
        }

        let params = self.sim_params;

        self.sim_state
            .get_or_insert_with(|| SimState::new(params, window))
            .step(window);
    }
}

#[derive(Clone, Copy)]
struct SimParams {
    connectivity_rate: f64,
    myelination_rate: f64,
    decay_rate: f64,
    max_myelination: usize,
    distance_exp: i32,
    refractory_period: usize,
}

impl Default for SimParams {
    fn default() -> Self {
        Self {
            connectivity_rate: 0.5,
            myelination_rate: 0.5,
            decay_rate: 0.01,
            max_myelination: 10,
            distance_exp: 4,
            refractory_period: 3,
        }
    }
}

struct SimState {
    rng: ThreadRng,
    sim: Simulation<ThreadRng>,
    simplicial_complex: SimplicialComplex,
    cached_outputs: (Vec<usize>, Vec<i64>),
    neuron_nodes: HashMap<usize, SceneNode>,
    synapse_nodes: HashMap<(usize, usize), SceneNode>,
    activity_nodes: Vec<SceneNode>,
}

impl SimState {
    const NUM_NODES: u32 = 6;
    const VIS_SCALE: f32 = 0.2;
    const SYNAPSE_MAX_WIDTH: f32 = 0.01;

    fn new(params: SimParams, window: &mut Window) -> Self {
        let rng = rand::thread_rng();
        let mut sim = Simulation::new(
            params.connectivity_rate,
            params.myelination_rate,
            params.decay_rate,
            params.max_myelination,
            params.distance_exp,
            params.refractory_period,
            rand::thread_rng(),
        );

        sim.init_uniform(3, Self::NUM_NODES);

        let simplicial_complex =
            SimplicialComplex::new((0..Self::NUM_NODES.pow(3) as usize).collect());

        let mut neuron_nodes = HashMap::with_capacity(sim.graph.node_count());

        for (id, node) in sim.graph.node_references() {
            let mut scene_node = window.add_sphere(0.05);
            scene_node.append_translation(&Translation3::from(
                nalgebra::convert::<_, Vector3<f32>>(node.position.coords) * Self::VIS_SCALE,
            ));

            neuron_nodes.insert(id.index(), scene_node);
        }

        let synapse_nodes = HashMap::new();
        let activity_nodes = Vec::new();

        Self {
            rng,
            sim,
            simplicial_complex,
            cached_outputs: Default::default(),
            neuron_nodes,
            synapse_nodes,
            activity_nodes,
        }
    }

    fn cleanup(&mut self, window: &mut Window) {
        for scene_node in self
            .neuron_nodes
            .values_mut()
            .chain(self.synapse_nodes.values_mut())
            .chain(self.activity_nodes.iter_mut())
        {
            window.remove_node(scene_node);
        }
    }
}

impl State for SimState {
    fn step(&mut self, window: &mut Window) {
        let id_range = Uniform::new(0, Self::NUM_NODES.pow(3) as usize);

        let result = self.sim.step(
            &self
                .rng
                .sample_iter(id_range)
                .take(if self.sim.timestep < 4000 { 10 } else { 2 })
                .collect::<Vec<_>>(),
        );

        for (id, node) in self.sim.graph.node_references() {
            let scene_node = self.neuron_nodes.get_mut(&id.index()).unwrap();

            if node.is_active(self.sim.timestep) {
                scene_node.set_color(1.0, 0.9, 0.0);
            } else {
                scene_node.set_color(0.75, 0.75, 0.75)
            }
        }

        for pair in result.removed_edges {
            self.simplicial_complex.remove(vec![pair.0, pair.1]);

            window.remove_node(&mut self.synapse_nodes.remove(&pair).unwrap());
        }

        for (source_id, target_id) in result.added_edges {
            self.simplicial_complex.add(vec![source_id, target_id]);

            let source = self
                .sim
                .graph
                .node_weight(NodeIndex::new(source_id))
                .unwrap();
            let target = self
                .sim
                .graph
                .node_weight(NodeIndex::new(target_id))
                .unwrap();

            let source_vec =
                nalgebra::convert::<_, Vector3<f32>>(source.position.coords) * Self::VIS_SCALE;
            let target_vec =
                nalgebra::convert::<_, Vector3<f32>>(target.position.coords) * Self::VIS_SCALE;

            let dist = distance(&source.position, &target.position) as f32 * Self::VIS_SCALE;
            let trans = Translation3::from((source_vec + target_vec) * 0.5);
            let rot = UnitQuaternion::face_towards(&(target_vec - source_vec), &Vector3::y());

            let mut scene_node = window.add_cylinder(1.0, dist);
            scene_node.append_translation(&trans);

            if ((target_vec - source_vec).normalize().abs() - Vector3::y()).magnitude() != 0.0 {
                scene_node.append_rotation_wrt_center(&UnitQuaternion::face_towards(
                    &Vector3::y(),
                    &Vector3::z(),
                ));

                scene_node.append_rotation_wrt_center(&rot);
            }

            self.synapse_nodes
                .insert((source_id, target_id), scene_node);
        }

        for node in self.activity_nodes.iter_mut() {
            window.remove_node(node);
        }

        self.activity_nodes.clear();

        for edge in self.sim.graph.edge_references() {
            let source = self.sim.graph.node_weight(edge.source()).unwrap();
            let target = self.sim.graph.node_weight(edge.target()).unwrap();

            let weight = edge.weight();

            let source_vec =
                nalgebra::convert::<_, Vector3<f32>>(source.position.coords) * Self::VIS_SCALE;
            let target_vec =
                nalgebra::convert::<_, Vector3<f32>>(target.position.coords) * Self::VIS_SCALE;

            let dist = distance(&source.position, &target.position) as f32 * Self::VIS_SCALE;
            let radius = (weight.myelination + 1) as f32
                * (Self::SYNAPSE_MAX_WIDTH / (self.sim.max_myelination + 1) as f32);

            let scene_node = self
                .synapse_nodes
                .get_mut(&(edge.source().index(), edge.target().index()))
                .unwrap();

            scene_node.set_local_scale(radius * 2.0, dist, radius * 2.0);

            for activation in weight.activation_queue.iter() {
                let mut scene_node = window.add_sphere(0.025);

                scene_node.set_color(1.0, 0.9, 0.0);
                scene_node.append_translation(&Translation3::from(
                    source_vec
                        + ((target_vec - source_vec).normalize()
                            * (dist / (1 + (activation.at - activation.queued_at)) as f32)
                            * (self.sim.timestep - activation.queued_at + 1) as f32),
                ));

                self.activity_nodes.push(scene_node);
            }
        }

        if self.sim.timestep % 10 == 0 {
            let lengths = self
                .simplicial_complex
                .simplices
                .iter()
                .map(|simplex| simplex.len())
                .collect();

            let betti_numbers = self.simplicial_complex.betti_numbers();

            self.cached_outputs = (lengths, betti_numbers);
        }
    }
}
