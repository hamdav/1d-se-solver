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

fn main()
{

    let window =
        Window::new_centered("Speedy2D: Input Callbacks Example", (1280, 960)).unwrap();

    window.run_loop(MyWindowHandler {
        drawn_curve: Vec::new(),
        potential: Vec::new(),
        support: (0., 0.),
        energy_scale: 100.,
        wf_scale: 2.,
        mouse_button_down: false,
        window_size: Vector2{x: 1280, y:960},
        wfs: vec![],
    })
}

fn convert_curve(curve: Vec<Vector2<f32>>) -> Vec<Vector2<f32>>{
    /*
     * Converts a curve to be a valid function over x
     * i.e. curve has only one y value per x value
     * also, reverses curve if it goes from right to left
     */

    // First, decide which way you should monotonize, 
    // is the curve mostly drawn from right to left or from left to right.
    let left_to_right = curve[0].x < curve[curve.len()-1].x;

    let mut rv = vec![curve[0]];
    for v in curve {
        if left_to_right && v.x > rv[rv.len()-1].x {
            rv.push(v);
        }
        if !left_to_right && v.x < rv[rv.len()-1].x {
            rv.push(v);
        }
    }
    if !left_to_right {
        rv.reverse();
    }
    rv
}

fn resample_curve(curve: Vec<Vector2<f32>>, n_samples: usize) -> Vec<f32> {
    /*
     * Returns the same curve but linearly interpolated at n_samples equally
     * spaced points in the same span as the original curve. 
     */
    let mut resampled_curve = vec![];
    let x_lo = curve[0].x;
    let x_hi = curve[curve.len()-1].x;
    let dx = (x_hi - x_lo) / (n_samples as f32 - 1.);
    let xs = curve.iter().map(|v| v.x).collect::<Vec<f32>>();

    for i in 0..n_samples {

        let x = x_lo + dx * i as f32;
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
    //println!("{:?}", curve);
    //println!("{:?}", resampled_curve);
    //
    resampled_curve
}



struct MyWindowHandler
{
    drawn_curve: Vec<Vector2<f32>>,
    potential: Vec<f32>,
    support: (f32, f32),
    energy_scale: f32,
    wf_scale: f32,
    wfs: Vec<numerov::State>,
    mouse_button_down: bool,
    window_size: Vector2<u32>,
}

impl MyWindowHandler {
    fn px2phys(&self, pxpos: Vector2<f32>) -> Vector2<f32> {
        /*
         * physical coordinates: 0,0 in middle of screen, 1,1 on top right corner
         */
        Vector2{
            x: pxpos.x / (self.window_size.x as f32 / 2.)  - 1.,
            y: -pxpos.y / (self.window_size.y as f32 / 2.)  + 1.,
        }
    }
    fn phys2px(&self, physpos: Vector2<f32>) -> Vector2<f32> {
        /*
         * physical coordinates: 0,0 in middle of screen, 1,1 on top right corner
         */
        Vector2{
            x: (physpos.x + 1.) * (self.window_size.x as f32 / 2.),
            y: (1. - physpos.y) * (self.window_size.y as f32 / 2.),
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
        let dx = (self.support.1 - self.support.0) / (self.potential.len() as f32 - 1.);
        for i in 1..self.potential.len() {
            let x = self.support.0 + i as f32 * dx;
            graphics.draw_line(
                self.phys2px(Vector2{
                    x: x - dx, y: self.potential[i-1] / self.energy_scale}), 
                self.phys2px(Vector2{
                    x: x, y: self.potential[i] / self.energy_scale}),
                    3.0, Color::RED);
        }


        // Draw the wf
        for wf in self.wfs.iter() {
            graphics.draw_line(self.phys2px(Vector2{
                x: -1., y: wf.energy / self.energy_scale}),
                self.phys2px(Vector2{
                x: 1., y: wf.energy / self.energy_scale}),
                3.0, Color::GREEN);

            let dx = (wf.x_right - wf.x_left) / (wf.wf.len() -1) as f32;
            for i in 0..wf.wf.len()-1 {
                graphics.draw_line(
                    self.phys2px(Vector2{
                        x: wf.x_left + i as f32 * dx, y: wf.wf[i] / self.wf_scale}), 
                    self.phys2px(Vector2{
                        x: wf.x_left + (i+1) as f32 * dx, y: wf.wf[i+1] / self.wf_scale}),
                    3.0, Color::GREEN);
            }
        }
    }

    fn on_mouse_move(&mut self, helper: &mut WindowHelper, position: Vector2<f32>)
    {
        // println!(
        //     "Got on_mouse_move callback: ({:.1}, {:.1})",
        //     position.x,
        //     position.y
        // );
        if self.mouse_button_down {
            self.drawn_curve.push(position);
            helper.request_redraw();
        }
    }

    fn on_mouse_button_down(&mut self, helper: &mut WindowHelper, button: MouseButton)
    {
        if button == MouseButton::Left {
            self.mouse_button_down = true;
        }
    }

    fn on_mouse_button_up(&mut self, helper: &mut WindowHelper, button: MouseButton)
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
        if virtual_key_code == Some(VirtualKeyCode::P) {
            let phys_curve: Vec<Vector2<f32>> = self.drawn_curve.iter()
                .map(|v| self.px2phys(*v))
                .collect();
            //println!("{:?}", phys_curve);
            self.support = (phys_curve[0].x, phys_curve[phys_curve.len()-1].x);
            self.potential = resample_curve(convert_curve(phys_curve), 100)
                .iter()
                .map(|v| self.energy_scale*v)
                .collect();
            //println!("{:?}", self.potential);
            helper.request_redraw();
        }
        if virtual_key_code == Some(VirtualKeyCode::Q) {
            self.wfs.extend(numerov::find_bound_states((-1., 1.), self.support, &self.potential));
            println!("{:?}", self.wfs);
            helper.request_redraw();
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
