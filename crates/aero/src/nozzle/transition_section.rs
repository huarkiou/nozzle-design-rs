use std::rc::Rc;

use crate::{
    moc::{
        unitprocess::{Context, ExitLineFunc, UnitProcess},
        CharLine, CharLines, MocPoint,
    },
    nozzle::Section,
};

pub struct TransitionSection {
    pub char_lines: CharLines,
    pub line_init: CharLine,
    cal_exit: Rc<dyn Fn(f64, &MocPoint) -> (f64, f64)>,
    pub max_length: f64,
    calculated: bool,
}

impl TransitionSection {
    pub fn new(cal_exit: ExitLineFunc, max_length: f64) -> Self {
        Self {
            char_lines: CharLines::new(),
            line_init: CharLine::new(),
            cal_exit: Rc::from(cal_exit),
            max_length,
            calculated: false,
        }
    }

    pub fn set_line_init(&mut self, line_init: CharLine) {
        self.line_init = line_init;
        self.char_lines = CharLines::new();
        self.calculated = false;
    }

    fn cal_next_transition_line(
        unitprocess: &dyn UnitProcess,
        line_prev: &CharLine,
        cal_exit: &Rc<dyn Fn(f64, &MocPoint) -> (f64, f64)>,
        _max_length: f64,
        dist_exit: f64,
    ) -> CharLine {
        let n = line_prev.len();
        if n < 2 {
            return CharLine::new();
        }

        let mfr_in = f64::abs(CharLine::mass_flow_rate(line_prev, unitprocess.area_type()));

        let mut line_cur = CharLine::with_capacity(n);
        for p in line_prev.iter() {
            line_cur.push(p.clone());
        }

        let mut j = n as i64 - 1;
        while j >= 0 {
            let ju = j as usize;

            if ju == line_cur.len() - 1 {
                // Exit point: step from prev line's last point along LEFT characteristic
                if line_prev.len() >= 2 {
                    let context = Context {
                        prev: line_prev,
                        next: &line_cur,
                        idx_prev: line_prev.len() - 1,
                        idx_next: 0,
                    };
                    if let Some(point) = unitprocess.exit_characteristics_point_fixed_dist(
                        context,
                        Box::new({
                            let f = Rc::clone(cal_exit);
                            move |y, p| f(y, p)
                        }),
                        dist_exit,
                    ) {
                        if point.is_valid() && point.y >= 0.0 {
                            line_cur[ju] = point;
                        }
                    }
                }

                if line_cur.len() == 2 {
                    // C++: 质量流量达标后缩减为单出口点
                    let exit_point = &line_cur[ju];
                    let last_prev = line_prev.last().unwrap();
                    let mfr_slice = [last_prev.clone(), exit_point.clone()];
                    let mfr_exit = f64::abs(CharLine::mass_flow_rate(
                        &mfr_slice,
                        unitprocess.area_type(),
                    ));
                    if (mfr_exit - mfr_in) > -1e-7 {
                        let context = Context {
                            prev: line_prev,
                            next: &line_cur,
                            idx_prev: 0,
                            idx_next: 0,
                        };
                        if let Some(point) = unitprocess.last_point(
                            context,
                            Box::new({
                                let f = Rc::clone(cal_exit);
                                move |y, p| f(y, p)
                            }),
                            mfr_in,
                        ) {
                            if point.is_valid() && point.x >= -1e-9 {
                                line_cur = CharLine::new();
                                line_cur.push(point);
                                break;
                            }
                        }
                    }
                }
            } else if ju == 0 {
                // Wall point
                if line_cur.len() > 1 {
                    let old_wall_x = line_cur[0].x;
                    let max_k = 3_usize.min(line_prev.len().saturating_sub(1));
                    let mut k: usize = 0;
                    while k <= max_k && 1 + k < line_prev.len() {
                        let context = Context {
                            prev: line_prev,
                            next: &line_cur,
                            idx_prev: k,
                            idx_next: 1,
                        };
                        if let Some(point) = unitprocess.transition_wall_point(context) {
                            let px = point.x;
                            if point.is_valid() && px >= -1e-9 && point.y >= 0.0 {
                                if line_cur.len() > 2 && px >= line_cur[1].x {
                                    line_cur.remove(1);
                                    continue;
                                }
                                if px > old_wall_x + 1e-12 {
                                    line_cur[0] = point;
                                    break;
                                }
                            }
                        }
                        k += 1;
                    }
                    // C++ distance cleanup: 壁面点与第一内点过近则移除内点
                    if line_cur.len() >= 3 {
                        let d01 = line_cur[0].distance_to(&line_cur[1]);
                        let d12 = line_cur[1].distance_to(&line_cur[2]);
                        if d01 < 0.2 * d12 {
                            line_cur.remove(1);
                        }
                    }
                }
            } else {
                // Interior point
                let ju_cur = ju.min(line_cur.len().saturating_sub(2));
                if ju_cur + 1 < line_cur.len() && ju < line_prev.len() {
                    let context = Context {
                        prev: line_prev,
                        next: &line_cur,
                        idx_prev: ju,
                        idx_next: ju_cur + 1,
                    };
                    if let Some(point) = unitprocess.transition_interior_point(context) {
                        if point.is_valid() && point.x >= -1e-9 && point.y >= -1e-9 {
                            line_cur[ju_cur] = point;
                        }
                    }
                }
            }
            j -= 1;
        }

        line_cur.retain(|p| p.is_valid());
        line_cur
    }
}

impl Section for TransitionSection {
    fn run(&mut self, unitprocess: &dyn UnitProcess, _config: &super::NozzleConfig) {
        assert!(
            !self.line_init.is_empty(),
            "TransitionSection: line_init must be set"
        );
        if self.calculated {
            return;
        }
        if self.line_init.len() < 2 {
            self.calculated = true;
            return;
        }

        self.char_lines = CharLines::new();
        let avg_dist: f64 = {
            let mut sum = 0.0;
            for i in 1..self.line_init.len() {
                sum += self.line_init[i].distance_to(&self.line_init[i - 1]);
            }
            sum / (self.line_init.len() - 1) as f64
        };

        let mut line_cur = self.line_init.clone();
        let max_iter = 2000u32;
        let mut iter = 0u32;

        while line_cur.len() > 1 {
            iter += 1;
            if iter > max_iter {
                break;
            }
            let next = Self::cal_next_transition_line(
                unitprocess,
                &line_cur,
                &self.cal_exit,
                self.max_length,
                avg_dist,
            );
            if next.is_empty() {
                continue;
            }
            self.char_lines.push(next.clone());
            line_cur = next;
            // 只要特征线未缩减为单点且壁面未超出长度，就继续推进
            if line_cur[0].x > self.max_length + 1e-6 {
                break;
            }
        }
        self.calculated = true;
    }

    fn get_charlines(&self) -> &CharLines {
        assert!(self.calculated, "TransitionSection has not run");
        &self.char_lines
    }

    fn inherit_last_line(&mut self, line: &CharLine) {
        self.set_line_init(line.clone());
    }
}

/// 计算 Rao 最大推力喷管理论中当前喷管的最佳出口背压 (OTN).
///
/// # 参数
/// - `point`: 壁面出口参考点
///
/// # 返回值
/// 最佳出口背压 (Pa)
///
/// # 公式
/// ```text
/// p_a = p - 0.5 * ρ * V² * sin(2θ) / cot(asin(1/Ma))
///     = p - 0.5 * ρ * V² * sin(2θ) / √(Ma² - 1)
/// ```
pub fn cal_pb_otn(point: &MocPoint) -> f64 {
    let v_sq = point.velocity_squared();
    let theta = point.flow_direction();
    let ma = point.mach_number();
    let ma_sq = ma * ma;
    // cot(asin(1/Ma)) = √(Ma² - 1)
    let cot_mach_angle = (ma_sq - 1.0).sqrt();
    point.p - 0.5 * point.rho * v_sq * (2.0 * theta).sin() / cot_mach_angle
}

/// Rao 最大推力喷管出口边界条件 (OTN).
/// 求解 y 坐标处满足 Riemann 不变量 C1 和质量流量不变量 C2 的 (V, θ).
pub fn cal_exit_otn(y: f64, p_ref: &MocPoint, axisym: bool) -> (f64, f64) {
    let ma_ref = p_ref.mach_number();
    let mu_ref = (1.0 / ma_ref).asin();
    let v_ref = p_ref.velocity();
    let theta_ref = p_ref.flow_direction();
    let gamma_ref = p_ref.gamma(p_ref.t);
    let rg = p_ref.rg();
    let (tt, _pt, rt) = p_ref.total_temperature_pressure_density();
    let cp = gamma_ref * rg / (gamma_ref - 1.0);

    // 不变量
    let c1 = v_ref * (mu_ref - theta_ref).cos() / mu_ref.cos();
    let c2 = if axisym {
        p_ref.y * p_ref.rho * v_ref * v_ref * theta_ref.sin() * theta_ref.sin() * mu_ref.tan()
    } else {
        p_ref.rho * v_ref * v_ref * theta_ref.sin() * theta_ref.sin() * mu_ref.tan()
    };

    let mut x = [v_ref, theta_ref];

    let _ = math::rootfinding::newton_2d(
        &|x: &[f64; 2]| -> [f64; 2] {
            let v_sq = x[0] * x[0];
            let sin_t = x[1].sin();
            let t_static = (tt - v_sq / (2.0 * cp)).max(1.0);
            let soundspeed_sq = gamma_ref * rg * t_static;
            let ma = (v_sq / soundspeed_sq).sqrt().max(1.0 + 1e-12);
            let mu = (1.0 / ma).asin();
            let rho = rt / (1.0 + (gamma_ref - 1.0) / 2.0 * ma * ma).powf(1.0 / (gamma_ref - 1.0));
            let delta_y = if axisym { y } else { 1.0 };

            let f1 = x[0] * (mu - x[1]).cos() / mu.cos() - c1;
            let f2 = delta_y * rho * v_sq * sin_t * sin_t * mu.tan() - c2;
            [f1, f2]
        },
        &mut x,
        1e-8,
        50,
    );

    (x[0] * x[1].cos(), x[0] * x[1].sin())
}

pub fn make_exit_otn(axisym: bool) -> ExitLineFunc {
    Box::new(move |y: f64, p_ref: &MocPoint| cal_exit_otn(y, p_ref, axisym))
}
