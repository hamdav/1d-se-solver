//#![deny(warnings)]
//

use std::mem;

use speedy2d::color::Color;
use speedy2d::dimen::Vector2;
use speedy2d::window::{
    KeyScancode,
    ModifiersState,
    MouseButton,
    MouseScrollDistance,
    VirtualKeyCode,
    WindowHandler,
    WindowHelper,
    WindowStartupInfo
};
use speedy2d::{Graphics2D, Window};

mod numerov;
mod exactdiag;


fn main()
{

    let window =
        Window::new_centered("1D SE solver", (1280, 960)).unwrap();

    window.run_loop(MyWindowHandler::new(Vector2{x: 1280, y:960}));
}

fn convert_curve(curve: Vec<Vector2<f64>>) -> Vec<Vector2<f64>>{
    /*
     * Converts a curve to be a valid function over x
     * i.e. curve has only one y value per x value
     * also, reverses curve if it goes from right to left
     */

    // First, decide which way you should monotonize, 
    // is the curve mostly drawn from right to left or from left to right.
    let left_to_right = curve[0].x < curve[curve.len()-1].x;

    // Only add new points if they are to the right (or left) 
    // of the old ones
    let mut rv = vec![curve[0]];
    for v in curve {
        if left_to_right && v.x > rv[rv.len()-1].x {
            rv.push(v);
        }
        if !left_to_right && v.x < rv[rv.len()-1].x {
            rv.push(v);
        }
    }

    // Reverse if needed
    if !left_to_right {
        rv.reverse();
    }

    rv
}

fn resample_curve(curve: Vec<Vector2<f64>>, n_samples: usize) -> Vec<f64> {
    /*
     * Returns the same curve but linearly interpolated at n_samples equally
     * spaced points in the same span as the original curve. 
     */
    let mut resampled_curve = vec![];
    let x_lo = curve[0].x;
    let x_hi = curve[curve.len()-1].x;
    let dx = (x_hi - x_lo) / (n_samples as f64 - 1.);

    // xs is the x coordinates of the original curve
    let xs = curve.iter().map(|v| v.x).collect::<Vec<f64>>();

    for i in 0..n_samples {

        let x = x_lo + dx * i as f64;
        // Check if the x is exactly at one of the original curves
        // points or if not, which two it is in between
        // If it's between two, interpolate and add that to the samples
        let res = xs.binary_search_by(|a| a.partial_cmp(&x).unwrap());

        match res {
            Ok(ind) => resampled_curve.push(curve[ind].y),
            Err(ind) => {
                if ind == curve.len() {
                    resampled_curve.push(curve[curve.len()-1].y);
                } else {
                    let t = (x - xs[ind-1]) / (xs[ind] - xs[ind-1]);
                    let y = curve[ind-1].y * (1.-t) + curve[ind].y * t;
                    resampled_curve.push(y);
                }
            }
        }
    }

    resampled_curve
}



struct MyWindowHandler
{
    drawn_curve: Vec<Vector2<f32>>,
    potential: Vec<f64>,
    support: (f64, f64),
    energy_scale: f64,
    wf_scale: f64,
    wfs: Vec<numerov::State>,
    mouse_button_down: bool,
    window_size: Vector2<u32>,
    marked_wf: Option<usize>,
}

impl MyWindowHandler {
    fn new(window_size: Vector2<u32>) -> Self {
        MyWindowHandler{
            energy_scale: 100.,
            wf_scale: 3.,
            window_size,
            drawn_curve: Vec::new(),
            potential: Vec::new(),
            support: (0., 0.),
            wfs: Vec::new(),
            mouse_button_down: false,
            marked_wf: None,
        }
    }
    fn update_potential_and_wf(&mut self) {
        /*
         * Recalculates the potential from the currently drawn curve
         * and recalculates the bound states of the potential,
         * updating self.wfs
         * If the recalculation invalidates the marked wf, the marking is
         * cleared
         */

        // Convert from screen coordinates to physical coordinates (-1, 1)
        let phys_curve: Vec<Vector2<f64>> = self.drawn_curve.iter()
            .map(|v| self.px2phys(*v))
            .collect();

        // Calculate the support of the potential
        self.support = (phys_curve[0].x, phys_curve[phys_curve.len()-1].x);

        // Calculate the potential by resampling and converting from
        // physical coordinates to energy coordinates
        self.potential = resample_curve(convert_curve(phys_curve), 300)
            .iter()
            .map(|v| self.energy_scale * v)
            .collect();
        // Extend the drawn potential by inserting zeros at the edges
        // Note that the support also needs to be altered
        let dx = (self.support.1 - self.support.0) 
            / (self.potential.len() as f64 - 1.);
        self.support.0 -= dx;
        self.support.1 += dx;
        self.potential.insert(0, 0.);
        self.potential.push(0.);

        // find the bound states
        //self.wfs = numerov::find_bound_states((-1., 1.), self.support, &self.potential);
        self.wfs = exactdiag::find_bound_states((-1., 1.), self.support, &self.potential);

        // If this recalculation made our marked state not exist,
        // stop marking it
        if let Some(idx) = self.marked_wf {
            if idx >= self.wfs.len() {
                self.marked_wf = None;
            }
        }
    }
    fn px2phys(&self, pxpos: Vector2<f32>) -> Vector2<f64> {
        /*
         * physical coordinates: 
         * 0,0 in middle of screen,
         * 1,1 on top right corner
         */
        Vector2{
            x: pxpos.x as f64 / (self.window_size.x as f64 / 2.)  - 1.,
            y: -pxpos.y as f64 / (self.window_size.y as f64 / 2.)  + 1.,
        }
    }
    fn phys2px(&self, physpos: Vector2<f64>) -> Vector2<f32> {
        /*
         * physical coordinates:
         * 0,0 in middle of screen,
         * 1,1 on top right corner
         */
        Vector2{
            x: (physpos.x as f32 + 1.) * (self.window_size.x as f32 / 2.),
            y: (1. - physpos.y as f32) * (self.window_size.y as f32 / 2.),

        }
    }
}

impl WindowHandler for MyWindowHandler
{
    fn on_start(&mut self, helper: &mut WindowHelper, info: WindowStartupInfo)
    {
        println!("Got on_start callback: {:?}", info);
        helper.set_resizable(true);

    }

    fn on_draw(&mut self, _helper: &mut WindowHelper, graphics: &mut Graphics2D)
    {
        // Clear the screen
        graphics.clear_screen(Color::from_rgb(0.1, 0.1, 0.1));

        // Draw coordinate axes
        graphics.draw_line(self.phys2px(Vector2{x: -1., y: 0.}), 
                           self.phys2px(Vector2{x:  1., y: 0.}), 1.0, Color::GRAY);
        graphics.draw_line(self.phys2px(Vector2{x: 0., y: -1.}), 
                           self.phys2px(Vector2{x: 0., y:  1.}), 1.0, Color::GRAY);

        // Draw the curve
        for (v1, v2) in self.drawn_curve.iter().zip(self.drawn_curve.iter().skip(1)) {
            graphics.draw_line(*v1, *v2, 2.0, Color::GRAY);
        }

        //Draw the potential
        let dx = (self.support.1 - self.support.0) / (self.potential.len() as f64 - 1.);
        for i in 1..self.potential.len() {
            let x = self.support.0 + i as f64 * dx;
            graphics.draw_line(
                self.phys2px(Vector2{
                    x: x - dx, y: self.potential[i-1] / self.energy_scale}), 
                self.phys2px(Vector2{
                    x: x, y: self.potential[i] / self.energy_scale}),
                    3.0, Color::RED);
        }


        // Draw the wf energy lines
        for (idx, wf) in self.wfs.iter().enumerate() {

            let dx = (wf.x_right - wf.x_left) / (wf.wf.len() -1) as f64;
            for i in 0..wf.wf.len()-1 {
                let color = if let Some(marked_idx) = self.marked_wf {
                    if marked_idx == idx {Color::GREEN} 
                    else {Color::GRAY}
                } else { Color::GRAY };
                graphics.draw_line(
                    self.phys2px(Vector2{
                        x: wf.x_left + i as f64 * dx, y: wf.energy / self.energy_scale}), 
                    self.phys2px(Vector2{
                        x: wf.x_left + (i+1) as f64 * dx, y: wf.energy / self.energy_scale}),
                    3.0, Color::from_rgba(color.r(), color.g(), color.b(), wf.wf[i].powi(2) as f32))
            }
        }
        // Draw the wavefunction
        if let Some(idx) = self.marked_wf {
            let wf = &self.wfs[idx];
            let dx = (wf.x_right - wf.x_left) / (wf.wf.len() -1) as f64;
            for i in 0..wf.wf.len()-1 {
                graphics.draw_line(
                    self.phys2px(Vector2{
                        x: wf.x_left + i as f64 * dx, y: wf.wf[i] / self.wf_scale}), 
                    self.phys2px(Vector2{
                        x: wf.x_left + (i+1) as f64 * dx, y: wf.wf[i+1] / self.wf_scale}),
                    3.0, Color::GREEN);
            }
        }
    }

    fn on_mouse_move(&mut self, helper: &mut WindowHelper, position: Vector2<f32>)
    {
        if self.mouse_button_down {
            self.drawn_curve.push(position);
            helper.request_redraw();
        }
    }

    fn on_mouse_button_down(&mut self, _helper: &mut WindowHelper, button: MouseButton)
    {
        if button == MouseButton::Left {
            self.mouse_button_down = true;
        }
    }

    fn on_mouse_button_up(&mut self, _helper: &mut WindowHelper, button: MouseButton)
    {
        if button == MouseButton::Left {
            self.mouse_button_down = false;
        }
    }

    fn on_key_down(
        &mut self,
        helper: &mut WindowHelper,
        virtual_key_code: Option<VirtualKeyCode>,
        scancode: KeyScancode
    )
    {
        println!(
            "Got on_key_down callback: {:?}, scancode {}",
            virtual_key_code,
            scancode
        );
        // Press Q to (re)calculate the bound states
        if virtual_key_code == Some(VirtualKeyCode::Q) {
            self.update_potential_and_wf();
            helper.request_redraw();
        }
        // Press K to decrease the energy scaling,
        // effectively decreasing the depth of the potential
        if virtual_key_code == Some(VirtualKeyCode::K) {
            self.energy_scale /= 1.03;
            self.update_potential_and_wf();
            helper.request_redraw();
        }
        // Press J to increase the energy scaling,
        // effectively increasing the depth of the potential
        if virtual_key_code == Some(VirtualKeyCode::J) {
            self.energy_scale *= 1.03;
            self.update_potential_and_wf();
            helper.request_redraw();
        }
        
        // Press C to clear everything, potential, drawn curve and wfs
        if virtual_key_code == Some(VirtualKeyCode::C) {
            self.potential = Vec::new();
            self.drawn_curve = Vec::new();
            self.wfs = Vec::new();
            self.marked_wf = None;
            helper.request_redraw();
        }
        // Press the up arrow to make the wf above marked,
        // or mark the bottom one if nothing is marked
        if virtual_key_code == Some(VirtualKeyCode::Up) {
            if let Some(idx) = self.marked_wf {
                if idx+1 < self.wfs.len() {
                    self.marked_wf = Some(idx+1);
                }
            } else {
                self.marked_wf = Some(0);
            }
            helper.request_redraw();
        }
        // Press the down arrow to make the wf below marked,
        // remove the mark if the bottom one is marked
        if virtual_key_code == Some(VirtualKeyCode::Down) {
            if let Some(idx) = self.marked_wf {
                if idx != 0 {
                    self.marked_wf = Some(idx-1);
                } else {
                    self.marked_wf = None;
                }
            }
            helper.request_redraw();
        }
        // Press H to decrease the wf scaling,
        if virtual_key_code == Some(VirtualKeyCode::H) {
            self.wf_scale /= 1.03;
            helper.request_redraw();
        }
        // Press G to increase the wf scaling,
        if virtual_key_code == Some(VirtualKeyCode::G) {
            self.wf_scale *= 1.03;
            helper.request_redraw();
        }

        // TEMPORARY FOR BUGS
        // Press D and F to change the energy which the wavefunction is calculated
        // from.
        if virtual_key_code == Some(VirtualKeyCode::F) {
            if let Some(i) = self.marked_wf {
                let new_energy = self.wfs[i].energy * 1.0001;
                let (new_psi, _) = numerov::bidirectional_shooting(new_energy,
                                                                   (-1., 1.),
                                                                   self.support,
                                                                   &self.potential);
                self.wfs[i] = new_psi;
                helper.request_redraw();
            }

        }
        if virtual_key_code == Some(VirtualKeyCode::D) {
            if let Some(i) = self.marked_wf {
                let new_energy = self.wfs[i].energy / 1.0001;
                let (new_psi, _) = numerov::bidirectional_shooting(new_energy,
                                                                   (-1., 1.),
                                                                   self.support,
                                                                   &self.potential);
                self.wfs[i] = new_psi;
                helper.request_redraw();
            }

        }

    }
    fn on_resize(&mut self, _helper: &mut WindowHelper, size_pixels: Vector2<u32>)
    {
        self.window_size = size_pixels;
        println!("Got on_resize callback: {:?}", size_pixels);
    }

    fn on_key_up(
        &mut self,
        _helper: &mut WindowHelper,
        virtual_key_code: Option<VirtualKeyCode>,
        scancode: KeyScancode
    )
    {
        println!(
            "Got on_key_up callback: {:?}, scancode {}",
            virtual_key_code,
            scancode
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resample() {
        let curve = vec![Vector2{x: 0., y: 1.},
            Vector2{x: 1., y: 2.},
            Vector2{x: 2., y: 2.}];
        
        let resampled_curve = resample_curve(curve, 5);
        assert_eq!(resampled_curve, vec![1., 1.5, 2., 2., 2.]);
    }
}
