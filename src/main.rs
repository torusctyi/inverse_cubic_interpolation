mod lib;
use lib::inv_cubic_solve;

// function for which we want to find the root and its derivative
fn func(x: f64)  -> f64{x.sin() + x.powi(3) }
fn dfunc(x: f64) -> f64{x.cos() + 3.0*x.powi(2) }

fn main(){
   let root =  inv_cubic_solve(-1.0, 0.5, 1e-14, func, dfunc);
   println!("\nroot = {}\nf(root) = {}\nlog10(|f(root)|) = {}", root, func(root), func(root).abs().log10());
}
           
