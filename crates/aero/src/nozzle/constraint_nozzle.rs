use crate::{
    moc::{unitprocess::UnitProcess, CharLine, CharLines},
    nozzle::{
        initial_line::InitialLine,
        transition_section::{make_exit_otn, TransitionSection},
        ExpansionSection, InitialSection, NozzleConfig, Section,
    },
};
use math::Tolerance;

pub struct ConstraintNozzle {
    config: NozzleConfig,
    unitprocess: Box<dyn UnitProcess>,
    sections: Vec<Box<dyn Section>>,
    /// 转向段出口边界条件工厂
    cal_exit_factory: fn(bool) -> crate::moc::unitprocess::ExitLineFunc,
}

impl ConstraintNozzle {
    /// 构建最大推力喷管OTN
    pub fn new_otn(config: NozzleConfig) -> Self {
        let length = config.geometry.length;
        let r_t = config.throat.radius_throat;
        let theta_a = config.throat.theta_a;
        let axisym = config.control.axisymmetric;
        Self {
            unitprocess: config.to_unitprocess(),
            config,
            sections: vec![
                Box::new(InitialLine::new()),
                Box::new(InitialSection::new()),
                Box::new(ExpansionSection::new(r_t, theta_a, length)),
                Box::new(TransitionSection::new(make_exit_otn(axisym), length)),
            ],
            cal_exit_factory: make_exit_otn,
        }
    }

    /// 计算喷管各段内流场数据
    ///
    /// 按顺序对每个段执行特征线法计算步骤，
    /// 前一个截面段的最后一条特征线作为下一截面段的初始边界条件传入。
    pub fn run(&mut self) {
        // 前 3 段流水线
        for i in 0..3 {
            if i > 0 {
                let last_line = {
                    let prev_lines = self.sections[i - 1].get_charlines();
                    prev_lines.last().cloned()
                };
                if let Some(line) = last_line {
                    self.sections[i].inherit_last_line(&line);
                }
            }
            self.sections[i].run(self.unitprocess.as_ref(), &self.config);
        }

        // 第 4 段：TransitionSection — 先试全段，再二分搜索最优子集
        {
            let last_line = {
                let prev_lines = self.sections[2].get_charlines();
                prev_lines.last().cloned()
            };
            if let Some(line) = last_line {
                self.sections[3].inherit_last_line(&line);
            }
        }
        self.sections[3].run(self.unitprocess.as_ref(), &self.config);
        self.cal_transition_to_target_length();
    }

    /// 二分搜索确定转向段初始线子集，使出口 x 接近目标长度
    fn cal_transition_to_target_length(&mut self) {
        let exp_lines = self.sections[2].get_charlines();
        let exp_last = match exp_lines.last() {
            Some(l) => l.clone(),
            None => return,
        };
        let n_exp = exp_last.len();
        if n_exp < 2 {
            return;
        }

        let tol = Tolerance::new(1e-4, 1e-4);
        let axisym = self.config.control.axisymmetric;
        let length = self.config.geometry.length;
        let up = self.unitprocess.as_ref();
        let cfg = &self.config;
        let exit_factory = self.cal_exit_factory;

        // 辅助函数：用 expansion 线前 i+1 个点作为 line_init，运行转向段，返回 exit_x - length
        let trial = |i: usize| -> f64 {
            let mut sub_line = CharLine::with_capacity(i + 1);
            for pt in exp_last.iter().take(i + 1) {
                sub_line.push(pt.clone());
            }
            let cal_exit = exit_factory(axisym);
            let mut ts = TransitionSection::new(cal_exit, length);
            ts.inherit_last_line(&sub_line);
            ts.run(up, cfg);
            ts.get_charlines()
                .last()
                .and_then(|line| line.last())
                .map_or(f64::NAN, |p| p.x - length)
        };

        // 先试全段
        let f_full = trial(n_exp - 1);
        if f_full.is_nan() || f_full.abs() < tol.abs {
            return;
        }

        // 二分搜索：找合适的 line_init 子集长度
        let lo = 1_i64;
        let hi = (n_exp - 1) as i64;
        let flo = trial(lo as usize);

        if flo.is_nan() || flo * f_full > 0.0 {
            return;
        }

        let mut left = lo;
        let mut right = hi;
        #[allow(unused_assignments)]
        let mut f_left = flo;
        for _ in 0..50 {
            if (right - left) <= 1 {
                break;
            }
            let mid = (left + right) / 2;
            let f_mid = trial(mid as usize);
            if f_mid.is_nan() {
                break;
            }
            if f_mid.abs() < tol.abs {
                left = mid;
                right = mid;
                break;
            }
            if f_left * f_mid < 0.0 {
                right = mid;
            } else {
                left = mid;
                f_left = f_mid;
            }
        }

        let best_i = if (trial(left as usize).abs()) < (trial(right as usize).abs()) {
            left as usize
        } else {
            right as usize
        };

        let mut sub_line = CharLine::with_capacity(best_i + 1);
        for pt in exp_last.iter().take(best_i + 1) {
            sub_line.push(pt.clone());
        }
        let cal_exit = exit_factory(axisym);
        let mut ts = TransitionSection::new(cal_exit, length);
        ts.inherit_last_line(&sub_line);
        ts.run(up, cfg);
        self.sections[3] = Box::new(ts);
    }

    pub fn get_assembly_charlines(&self) -> CharLines {
        let mut lines: CharLines = CharLines::new();
        for section in &self.sections {
            lines.extend(section.get_charlines());
        }
        lines
    }
}

#[cfg(test)]
mod tests {
    use std::path::PathBuf;

    use crate::{nozzle::config::*, Material};

    use super::*;
    #[test]
    fn test_new_and_run() {
        let config = NozzleConfig {
            control: Control::default(),
            material: Material::air_piecewise_polynomial(), // 使用理想空气
            inlet: Inlet::default(),
            geometry: Geometry::default(),
            throat: Throat::default(),
            outlet: Outlet::default(),
            io: IO::default(),
        };
        let target_len = config.geometry.length;
        let mut n = ConstraintNozzle::new_otn(config);
        n.run();
        let lines = n.get_assembly_charlines();

        // 验证至少产生了特征线数据
        assert!(!lines.is_empty(), "应产生流场数据");

        // 验证所有点坐标不越界
        for (li, line) in lines.iter().enumerate() {
            for (pi, point) in line.iter().enumerate() {
                assert!(
                    point.is_valid(),
                    "存在无效点: line={li} pt={pi}: {:?}",
                    point
                );
                if point.x < 0.0 || point.y < 0.0 {
                    panic!(
                        "line={li} pt={pi}: x={}, y={}, u={}, v={}",
                        point.x, point.y, point.u, point.v
                    );
                }
            }
        }

        // 验证出口 x 不超出目标长度
        for line in lines.iter() {
            for point in line.iter() {
                assert!(
                    point.x < target_len + 1.0,
                    "x={} 超出长度 {target_len}",
                    point.x
                );
            }
        }
    }

    /// 使用非零膨胀角、非零喉部半径测试完整喷管（含膨胀段）
    #[test]
    fn test_with_expansion() {
        let config = NozzleConfig {
            control: Control::default(),
            material: Material::from_rgas_gamma(287.042, 1.4),
            inlet: Inlet::default(),
            geometry: Geometry::default(),
            throat: Throat {
                radius_throat: 0.0,             // 喉部过渡圆弧半径
                theta_a: 19.0_f64.to_radians(), // 初始膨胀角 19°（接近 C++ 自动优化值）
            },
            outlet: Outlet::default(),
            io: IO::default(),
        };
        let target_len = config.geometry.length;
        let mut n = ConstraintNozzle::new_otn(config);
        n.run();
        let lines = n.get_assembly_charlines();

        // 验证至少产生了特征线数据
        assert!(!lines.is_empty(), "应产生流场数据");

        // 验证所有点坐标不越界
        for (li, line) in lines.iter().enumerate() {
            for (pi, point) in line.iter().enumerate() {
                assert!(
                    point.is_valid(),
                    "存在无效点: line={li} pt={pi}: {:?}",
                    point
                );
                if point.x < 0.0 || point.y < 0.0 {
                    panic!(
                        "line={li} pt={pi}: x={}, y={}, u={}, v={}",
                        point.x, point.y, point.u, point.v
                    );
                }
            }
        }

        // 验证出口 x 不超出目标长度
        for line in lines.iter() {
            for point in line.iter() {
                assert!(
                    point.x < target_len + 1.0,
                    "x={} 超出长度 {target_len}",
                    point.x
                );
            }
        }

        // 写入文件以便可视化检查
        let mut output_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
        output_dir.push("../../target/tmp/fluid_field_expansion.txt");
        lines.write_to_file(output_dir, false).unwrap();
    }
}
