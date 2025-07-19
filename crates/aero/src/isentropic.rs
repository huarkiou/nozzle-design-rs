use math::{Tolerance, rootfinding::toms748};

use crate::material::Cp;

pub fn cal_total_temperature(cp: &Cp, t_static: f64, velocity: f64) -> f64 {
    match cp {
        Cp::Constant(v) => t_static + velocity.powi(2) / (2. * v),
        Cp::Variable(f) => {
            let equation = |t_total| t_static + velocity.powi(2) / (2. * f(t_static)) - t_total;
            let res =
                toms748::solve_bracket(0., 6000., &equation, Tolerance::new(1e-12, 1e-12), 30)
                    .unwrap();
            res.root()
        }
    }
}

pub fn cal_static_temperature(cp: &Cp, t_total: f64, velocity: f64) -> f64 {
    match cp {
        Cp::Constant(v) => t_total - velocity.powi(2) / (2. * v),
        Cp::Variable(f) => {
            let equation = |t_static| t_total - velocity.powi(2) / (2. * f(t_static)) - t_static;
            let res =
                toms748::solve_bracket(0., t_total, &equation, Tolerance::new(1e-12, 1e-12), 30)
                    .unwrap();
            res.root()
        }
    }
}
