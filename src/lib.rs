use std::f64;

const BISECTION_FREQ : usize = 5;

#[derive(PartialEq)]
enum CoordinateChangeFlag {
   First,
   Second,
   Reset
}

/** Standard bisection method */
fn naive_bisection (x0: f64, x1: f64, f0: f64, f1: f64,  func: fn(f64) -> f64) -> (f64, f64){

   let x_new = (x0 + x1)/2.0; 
   let f_new = func(x_new);
   
   if f_new * f0 <= 0.0{
       (x0, x_new)
   }
   else if f_new * f1 <= 0.0 {
       (x_new, x1)
   }
   else {
        panic!("Bisection Failure");
   }
}

fn false_position (x0: f64, x1: f64, f0: f64, f1: f64,  func: fn(f64) -> f64) -> (f64, f64){

   let x_new = (x0*f1 - x1*f0)/(f1 - f0); 
   let f_new = func(x_new);
   
   if f_new * f0 <= 0.0{
       return (x0, x_new);
   }
   else if f_new * f1 <= 0.0 {
       return (x_new, x1);
   }
   else {
        panic!("Bisection Failure");
   }
}


/** Evaluate the  cubic that matches f(x0), f'(x0), f(x1),  and f'(x1) at the value x */
fn two_point_cubic(x: f64, x0: f64, x1: f64, f0: f64, df0: f64, f1: f64, df1: f64) -> f64{

     let t: f64 = (x - x0)/(x1 - x0);

     let h00 = 2.0*t*t*t - 3.0*t*t + 1.0;
     let h10 = t*t*t - 2.0*t*t + t;
     let h01 = -2.0*t*t*t + 3.0*t*t;
     let h11 =  t*t*t - t*t;

     f0*h00 + df0*h10 + f1*h01 + df1*h11
}

/** The same as two_point_cubic, except the inverse of f(x) is approximated, and evaluated at 0 */
fn two_point_cubic_inverse ( x0: f64, x1: f64, f0: f64, df0: f64, f1: f64, df1: f64) -> f64{
     let y0 = f0;
     let y1 = f1;
 
     let g0 = x0;
     let g1 = x1;
       
     let dg0 = 1.0/df0;
     let dg1 = 1.0/df1;

     two_point_cubic(0.0, y0, y1, g0, dg0, g1, dg1)
}

/** Given points x0 and x1 such that f(x0)*f(x1) < 0, use cubic interpolation to find a point between that is an approximate 
  * root of f(x).
  */
fn cubic_bisection (x0: f64, x1: f64, f0: f64, df0:f64,  f1: f64, df1: f64,  func: fn(f64) -> f64) -> Result<(f64, f64), f64>{
   let x_new = two_point_cubic_inverse(x0, x1, f0, df0, f1, df1);

   if x_new <= x0 || x_new >= x1 {
       Err(std::f64::NAN)
   } 
   else {
       let f_new = func(x_new);
   
       if f_new * f0 <= 0.0{
           Ok((x0, x_new))
       }
       else if f_new * f1 <= 0.0{
           Ok((x_new, x1))
       }
       else { 
           panic!("Bisection Failure"); // this condition should never happen in theory
       }
   }
}  

/** The actual solver */
pub fn inv_cubic_solve(x0: f64, x1: f64, tol: f64, func: fn(f64) -> f64, deriv: fn(f64) -> f64) ->  f64{

    let mut x_best;
    let mut x  = (x0, x1);
    let mut dx  = (x0 - x1).abs();

    // function values
    let mut f0 = func(x0);
    let mut f1 = func(x1);

    // derivatives
    let mut df0 = deriv(x0);
    let mut df1 = deriv(x1);

    let mut last_coord_changed =  CoordinateChangeFlag::Reset;
    let mut should_bisect  = false;

    // number of iterations
    let mut n_iters :usize = 1;

    /* 
       NB: The loop below will terminate so long as f(x0) and f(x1) are of opposite sign.  In the worst case, the width will halve every BISECTION_FREQ iterations
    */
    loop{

        // temporarily store the current values
        let x_old = x;

        let f_max = f0.abs().max(f1.abs());
        let f_min = f0.abs().min(f1.abs());

        // print data
        println!("{0:0<02} x1 = {1:0<022.19} x2= {2:0<022.19}  min(|f(x1)|, |f(x2)|) = {3:0<022.19}  max(|f(x1)|, |f(x2)|) = {4:0<024.19} log10(|x2 - x1|) = {5:0<+012.10}",
                        n_iters, x.0 , x.1, f_min, f_max, dx.abs().log10());

        // get the best point found so far
        x_best =  if f0.abs() < f1.abs() {x.0} else {x.1};

        // if the method has converged, return the best point
        if dx.abs() < tol || f_min < tol {
            return x_best;
        }
              
        // perform bisection every nth iteration, or if a point hasn't been changed in two iterations.
        if n_iters % BISECTION_FREQ  == 0 || should_bisect {
            x = naive_bisection(x.0, x.1, f0, f1, func);
            should_bisect = false;
            // reset the flag to tell if a point hasn't changed in two iterations
            last_coord_changed = CoordinateChangeFlag::Reset; 
        } else {
            let result = cubic_bisection(x.0, x.1, f0, df0,  f1, df1, func);
            x = match result{
               Ok(_x) => _x,
               Err(_) => false_position(x.0, x.1, f0, f1, func)
           };      
        }

        // update the function values and derivatives depending on which point was changed
        // and determine if the same point was changed twice in a row
        let second_last_coord_changed = last_coord_changed;
        if x_old.0 != x.0 {
              f0 = func(x.0);
             df0 = deriv(x.0);
             last_coord_changed = CoordinateChangeFlag::First;
        } else {
              f1 = func(x.1);
             df1 = deriv(x.1);
             last_coord_changed = CoordinateChangeFlag::Second;
        };
        
        if last_coord_changed == second_last_coord_changed{
            should_bisect = true;
        }

        // update various statistics
        n_iters = n_iters + 1;
        dx = (x.1 - x.0).abs();

    }   
}

// Let's test it out!!!




