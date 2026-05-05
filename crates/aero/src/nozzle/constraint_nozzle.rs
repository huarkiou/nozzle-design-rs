use std::f64::consts::PI;

use crate::{
    moc::{unitprocess::UnitProcess, CharLine, CharLines, MocPoint},
    nozzle::{
        initial_line::InitialLine,
        transition_section::{cal_pb_otn, make_exit_otn, TransitionSection},
        uniform_section::UniformSection,
        ExpansionSection, InitialSection, NozzleConfig, Section,
    },
};
use math::{rootfinding::toms748, Tolerance};

pub struct ConstraintNozzle {
    config: NozzleConfig,
    unitprocess: Box<dyn UnitProcess>,
    sections: Vec<Box<dyn Section>>,
    /// 转向段出口边界条件工厂
    cal_exit_factory: fn(bool) -> crate::moc::unitprocess::ExitLineFunc,
    /// 均一区（填补转向段起始线子集与膨胀段完整最后特征线之间的缺口）
    uniform_section: UniformSection,
}

/// 给定膨胀段最后一条特征线，运行转向段并二分搜索使出口 x 逼近目标长度。
///
/// 分两阶段搜索（对应 C++ constraint-nozzle.hpp）：
/// 1. 整数索引二分 — 确定转向段起始线在膨胀段最后特征线中的网格区间
/// 2. 弧长插值二分 — 在区间内精确找到使出口 x 匹配 target_length 的起始点
///
/// 返回构建好的 `TransitionSection`；若无法收敛则返回 `None`。
fn run_transition_to_target_length(
    exp_last: &CharLine,
    unitprocess: &dyn UnitProcess,
    config: &NozzleConfig,
    cal_exit_factory: fn(bool) -> crate::moc::unitprocess::ExitLineFunc,
) -> Option<TransitionSection> {
    let n_exp = exp_last.len();
    if n_exp < 2 {
        return None;
    }

    let eps = config.control.eps;
    let axisym = config.control.axisymmetric;
    let length = config.geometry.length;
    let exit_factory = cal_exit_factory;

    // 辅助函数：用 expansion 线前 i+1 个点作为 line_init，运行转向段，返回 exit_x - length
    let trial = |i: usize| -> f64 {
        let mut sub_line = CharLine::with_capacity(i + 1);
        for pt in exp_last.iter().take(i + 1) {
            sub_line.push(pt.clone());
        }
        let cal_exit = exit_factory(axisym);
        let mut ts = TransitionSection::new(cal_exit, length);
        ts.inherit_last_line(&sub_line);
        ts.run(unitprocess, config);
        ts.get_charlines()
            .last()
            .and_then(|line| line.last())
            .map_or(f64::NAN, |p| p.x - length)
    };

    // ── 阶段 1：整数索引二分 ──
    // 先试全段
    let f_full = trial(n_exp - 1);
    if f_full.is_nan() {
        return None;
    }
    if f_full.abs() < eps {
        // 全段即收敛，直接用全段结果
        let cal_exit = exit_factory(axisym);
        let mut ts = TransitionSection::new(cal_exit, length);
        let mut full_line = CharLine::with_capacity(n_exp);
        for pt in exp_last.iter() {
            full_line.push(pt.clone());
        }
        ts.inherit_last_line(&full_line);
        ts.run(unitprocess, config);
        return Some(ts);
    }

    // ── 阶段 1：整数索引二分 — 用通用整数二分法确定网格区间 ──
    let lo = 1_i64;
    let hi = (n_exp - 1) as i64;
    let flo = trial(lo as usize);

    if flo.is_nan() || flo * f_full > 0.0 {
        return None;
    }

    let ib = match math::rootfinding::solve_bracket(
        lo,
        hi,
        &|i| trial(i as usize),
        |a, b| (a + b) / 2,
        |a, b| b - a <= 1,
        50,
    ) {
        Ok(b) => b,
        Err(_) => return None,
    };

    let i_min = ib.lo as usize;
    let i_max = ib.hi as usize;

    // ── 阶段 2：弧长插值二分 — 在网格区间内精确找目标点 ──
    let p_min = &exp_last[i_min];
    let p_max = &exp_last[i_max];
    let l_max = p_max.distance_to(p_min);
    let alpha = ((p_max.y - p_min.y) / (p_max.x - p_min.x)).atan();

    // 辅助函数：在弧长 L_cur 处插值一点，用它作为 line_init 最后一点，运行转向段
    let cal_transition_l = |l_cur: f64| -> f64 {
        let mut point_tmp = MocPoint {
            x: p_min.x + l_cur * alpha.cos(),
            y: p_min.y + l_cur * alpha.sin(),
            u: 0.0,
            v: 0.0,
            p: 0.0,
            t: 0.0,
            rho: 0.0,
            mat: p_min.mat.clone(),
        };
        point_tmp = point_tmp.interpolate_along(p_min, p_max);

        // 当插值点靠近 p_min 时用 i_min 个前导点，否则用 i_min+1 个
        let n_before = if l_cur / l_max < 0.5 {
            i_min.saturating_sub(1)
        } else {
            i_min
        };
        let mut sub_line = CharLine::with_capacity(n_before + 2);
        for pt in exp_last.iter().take(n_before + 1) {
            sub_line.push(pt.clone());
        }
        sub_line.push(point_tmp);

        let cal_exit = exit_factory(axisym);
        let mut ts = TransitionSection::new(cal_exit, length);
        ts.inherit_last_line(&sub_line);
        ts.run(unitprocess, config);
        ts.get_charlines()
            .last()
            .and_then(|line| line.last())
            .map_or(f64::NAN, |p| p.x - length)
    };

    // 检查阶段 2 边界是否有解
    let f0 = cal_transition_l(0.0);
    let f_max = cal_transition_l(l_max);
    if f0.is_nan() || f_max.is_nan() {
        // 插值不可用，回退到阶段 1 的结果
        return trial_best(
            exp_last,
            unitprocess,
            config,
            exit_factory,
            axisym,
            length,
            i_min,
            i_max,
            &trial,
        );
    }
    if f0.abs() < eps {
        return trial_best(
            exp_last,
            unitprocess,
            config,
            exit_factory,
            axisym,
            length,
            i_min,
            i_min,
            &trial,
        );
    }
    if f_max.abs() < eps {
        return trial_best(
            exp_last,
            unitprocess,
            config,
            exit_factory,
            axisym,
            length,
            i_max,
            i_max,
            &trial,
        );
    }
    if f0 * f_max > 0.0 {
        // 同号，回退到阶段 1 的最佳结果
        return trial_best(
            exp_last,
            unitprocess,
            config,
            exit_factory,
            axisym,
            length,
            i_min,
            i_max,
            &trial,
        );
    }

    // 阶段 2: TOMS748 寻根 — 在弧长区间内精确找使出口 x 匹配 target_length 的起始点
    let toms748_f = |l_cur: f64| -> f64 { cal_transition_l(l_cur) };

    let l_best = match toms748::solve_bracket(0.0, l_max, &toms748_f, Tolerance::new(eps, eps), 30)
    {
        Ok(bracket) => {
            // 检查是否收敛到边界
            if bracket.lo.abs() < eps || (bracket.hi - l_max).abs() < eps {
                // 根靠近区间端点，可能不在区间内
            }
            0.5 * (bracket.lo + bracket.hi)
        }
        Err(_) => {
            // TOMS748 失败，回退到 f0/f_max 中更接近零的端点
            if f0.abs() < f_max.abs() {
                0.0
            } else {
                l_max
            }
        }
    };
    let mut point_tmp = MocPoint {
        x: p_min.x + l_best * alpha.cos(),
        y: p_min.y + l_best * alpha.sin(),
        u: 0.0,
        v: 0.0,
        p: 0.0,
        t: 0.0,
        rho: 0.0,
        mat: p_min.mat.clone(),
    };
    point_tmp = point_tmp.interpolate_along(p_min, p_max);

    let n_before = if l_best / l_max < 0.5 {
        i_min.saturating_sub(1)
    } else {
        i_min
    };
    let mut sub_line = CharLine::with_capacity(n_before + 2);
    for pt in exp_last.iter().take(n_before + 1) {
        sub_line.push(pt.clone());
    }
    sub_line.push(point_tmp);

    let cal_exit = exit_factory(axisym);
    let mut ts = TransitionSection::new(cal_exit, length);
    ts.inherit_last_line(&sub_line);
    ts.run(unitprocess, config);
    Some(ts)
}

/// 回退函数：当阶段 2 插值不可用时，使用阶段 1 的最佳索引结果
fn trial_best(
    exp_last: &CharLine,
    unitprocess: &dyn UnitProcess,
    config: &NozzleConfig,
    cal_exit_factory: fn(bool) -> crate::moc::unitprocess::ExitLineFunc,
    axisym: bool,
    length: f64,
    i_min: usize,
    i_max: usize,
    trial: &dyn Fn(usize) -> f64,
) -> Option<TransitionSection> {
    let best_i = if (trial(i_min).abs()) < (trial(i_max).abs()) {
        i_min
    } else {
        i_max
    };
    let mut sub_line = CharLine::with_capacity(best_i + 1);
    for pt in exp_last.iter().take(best_i + 1) {
        sub_line.push(pt.clone());
    }
    let cal_exit = cal_exit_factory(axisym);
    let mut ts = TransitionSection::new(cal_exit, length);
    ts.inherit_last_line(&sub_line);
    ts.run(unitprocess, config);
    Some(ts)
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
            uniform_section: UniformSection::new(length),
        }
    }

    /// 计算喷管各段内流场数据
    ///
    /// 按顺序对每个段执行特征线法计算步骤，
    /// 前一个截面段的最后一条特征线作为下一截面段的初始边界条件传入。
    ///
    /// 当 `theta_a` 为 NaN、负值或 >= 90° 时，自动迭代选取满足背压约束的初始膨胀角。
    pub fn run(&mut self) {
        // ── 第 0～1 段：初值线 + 初值问题（与 theta_a 无关） ──
        for i in 0..2 {
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

        // ── 判断是否需要自动迭代 theta_a ──
        let theta_a = self.config.throat.theta_a;
        let theta_auto = !theta_a.is_finite() || theta_a <= 0.0 || theta_a >= PI / 2.0;

        if theta_auto {
            // 自动迭代：从初值段最后一条特征线出发，反复试算膨胀+转向，寻找最优 theta_a
            let initial_last_line = self.sections[1]
                .get_charlines()
                .last()
                .cloned()
                .expect("InitialSection 未产生任何特征线");
            let optimal_theta = self.find_optimal_theta_a(&initial_last_line);

            // 将最佳 theta_a 写回配置和膨胀段
            self.config.throat.theta_a = optimal_theta;
            self.sections[2].set_theta_a(optimal_theta);

            if !self.config.geometry.height_e.is_finite() {
                eprintln!(
                    "  →  selected theta_a = {:.6}° for p_ambient = {:.1} Pa",
                    optimal_theta.to_degrees(),
                    self.config.outlet.p_ambient
                );
            } else {
                eprintln!(
                    "  →  selected theta_a = {:.6}° for height_e = {:.4} m",
                    optimal_theta.to_degrees(),
                    self.config.geometry.height_e
                );
            }
        }

        // ── 第 2 段：膨胀段 ──
        {
            let last_line = {
                let prev_lines = self.sections[1].get_charlines();
                prev_lines.last().cloned()
            };
            if let Some(line) = last_line {
                self.sections[2].inherit_last_line(&line);
            }
        }
        self.sections[2].run(self.unitprocess.as_ref(), &self.config);

        // ── 第 3 段：转向段 — 先试全段，再二分搜索最优子集 ──
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

        // ── 均一区：填补转向段起始线子集与膨胀段完整最后特征线之间的缺口 ──
        self.run_uniform_section();
    }

    /// 构建并运行均一区：从膨胀段最后特征线的"未使用"部分出发，
    /// 对转向段每个出口点计算对应的右行特征线。
    fn run_uniform_section(&mut self) {
        let exp_last_line = match self.sections[2].get_charlines().last() {
            Some(l) => l,
            None => return,
        };
        let n_exp = exp_last_line.len();
        if n_exp < 2 {
            return;
        }

        // 转向段起始线使用了膨胀段最后特征线的子集（前 trans_init_size 个点），
        // 均一区 line_init = 剩余的点（靠近对称轴侧）
        let trans_init_size = self.sections[3].line_init_len();

        let uniform_line_init = if trans_init_size >= n_exp {
            // 转向段使用了全部点，均一区仅用对称轴点
            let mut line = CharLine::with_capacity(1);
            line.push(exp_last_line.last().unwrap().clone());
            line
        } else {
            let count = n_exp - trans_init_size;
            let mut line = CharLine::with_capacity(count);
            for pt in exp_last_line.iter().skip(trans_init_size) {
                line.push(pt.clone());
            }
            line
        };

        // 均一区 line_exit = 每条转向段特征线的最末点（出口点）
        let mut uniform_line_exit = CharLine::new();
        for line in self.sections[3].get_charlines().iter() {
            if let Some(pt) = line.last() {
                uniform_line_exit.push(pt.clone());
            }
        }

        self.uniform_section.line_init = uniform_line_init;
        self.uniform_section.line_exit = uniform_line_exit;
        self.uniform_section.length = self.config.geometry.length;
        self.uniform_section
            .run(self.unitprocess.as_ref(), &self.config);
    }

    /// 自动迭代寻找满足背压（或出口高度）约束的最优初始膨胀角 `theta_a`。
    ///
    /// 算法流程（参考 C++ constraint-nozzle.hpp）：
    /// 1. 从大到小试探上界（80°→10°），从小到大试探下界（0.1°→30°）
    /// 2. 用 TOMS748 在上下界之间寻根
    /// 3. 返回使 `p_a - p_ambient = 0`（或 `h_e - height_e = 0`）的 theta_a
    fn find_optimal_theta_a(&self, initial_last_line: &CharLine) -> f64 {
        eprintln!("  Finding best initial expansion angle theta_a ...");

        let p_ambient = self.config.outlet.p_ambient;
        let height_e_fixed = self.config.geometry.height_e.is_finite();

        // ── 辅助函数：给定 theta_a，计算约束差值 ──
        let cal_cond_diff = |theta_a: f64, verbose: bool| -> Option<f64> {
            let r_t = self.config.throat.radius_throat;
            let length = self.config.geometry.length;

            // 1. 创建并运行膨胀段
            let mut exp = ExpansionSection::new(r_t, theta_a, length);
            exp.inherit_last_line(initial_last_line);
            exp.run(self.unitprocess.as_ref(), &self.config);

            let exp_lines = exp.get_charlines();
            let exp_last = exp_lines.last()?;
            if exp_last.len() < 2 {
                if verbose {
                    eprintln!(
                        "    θ = {:.3}°  →  expansion produced < 2 points, skip",
                        theta_a.to_degrees()
                    );
                }
                return None;
            }

            // 2. 运行转向段至目标长度
            let ts = run_transition_to_target_length(
                exp_last,
                self.unitprocess.as_ref(),
                &self.config,
                self.cal_exit_factory,
            )?;

            // 3. 取出壁面出口点，计算约束差值
            let ts_lines = ts.get_charlines();
            let wall_point = ts_lines.last()?.first()?;
            if !wall_point.is_valid() {
                if verbose {
                    eprintln!(
                        "    θ = {:.3}°  →  exit wall point invalid, skip",
                        theta_a.to_degrees()
                    );
                }
                return None;
            }

            // 喷嘴长度 = 出口轴线点 x - 初值线壁面点 x
            let nozzle_len = ts_lines
                .last()
                .and_then(|line| line.last())
                .map(|p| p.x)
                .unwrap_or(f64::NAN);

            if height_e_fixed {
                // 固定出口高度模式：基于出口高度收敛
                let h_e = wall_point.y;
                let diff = h_e - self.config.geometry.height_e;
                let pct =
                    (h_e - self.config.geometry.height_e) / self.config.geometry.height_e * 100.0;
                if verbose {
                    eprintln!(
                        "    θ = {:>10.6}°  →  h_e = {:>10.6} m  ({:>+9.4}%)  L = {:.6} m",
                        theta_a.to_degrees(),
                        h_e,
                        pct,
                        nozzle_len
                    );
                }
                Some(diff)
            } else {
                // 自由出口高度模式：基于背压收敛
                let p_a = cal_pb_otn(wall_point);
                let diff = p_a - p_ambient;
                let pct = (p_a - p_ambient) / p_ambient * 100.0;
                if verbose {
                    eprintln!(
                        "    θ = {:>10.6}°  →  p_a = {:>10.3} Pa  ({:>+9.4}%)  L = {:.6} m",
                        theta_a.to_degrees(),
                        p_a,
                        pct,
                        nozzle_len
                    );
                }
                Some(diff)
            }
        };

        // ── 1. 上界试探：从大到小找能使计算不崩溃的角度 ──
        eprintln!("    searching upper bound (80° → 10°) ...");
        let upper_result = (|| -> Option<(f64, f64)> {
            for degree in [80.0_f64, 70.0, 60.0, 50.0, 40.0, 30.0, 20.0, 10.0] {
                let theta = degree.to_radians();
                if let Some(diff) = cal_cond_diff(theta, false) {
                    if diff.is_finite() {
                        // 打印成功的上界以便用户了解范围
                        cal_cond_diff(theta, true);
                        return Some((theta, diff));
                    }
                }
            }
            None
        })();

        // ── 2. 下界试探：从小到上找能使计算不崩溃的角度 ──
        eprintln!("    searching lower bound (0.1° → 30°) ...");
        let lower_result = (|| -> Option<(f64, f64)> {
            for degree in [0.1_f64, 1.0, 10.0, 20.0, 30.0] {
                let theta = degree.to_radians();
                if let Some(diff) = cal_cond_diff(theta, false) {
                    if diff.is_finite() {
                        cal_cond_diff(theta, true);
                        return Some((theta, diff));
                    }
                }
            }
            None
        })();

        // ── 3a. 两边界均失败 — 退回默认值 ──
        if upper_result.is_none() && lower_result.is_none() {
            eprintln!("    all bracket candidates failed, using default 20°");
            return 20.0_f64.to_radians();
        }

        // ── 3b. 仅有一边界成功 — 直接使用该角度 ──
        if upper_result.is_none() {
            let (theta, _) = lower_result.unwrap();
            eprintln!(
                "    only lower bound found, using θ = {:.3}°",
                theta.to_degrees()
            );
            return theta;
        }
        if lower_result.is_none() {
            let (theta, _) = upper_result.unwrap();
            eprintln!(
                "    only upper bound found, using θ = {:.3}°",
                theta.to_degrees()
            );
            return theta;
        }

        let (xa, fa) = lower_result.unwrap();
        let (xb, fb) = upper_result.unwrap();

        eprintln!(
            "    bracket: [{:.3}°, {:.3}°]  fa={:+.3e}  fb={:+.3e}",
            xa.to_degrees(),
            xb.to_degrees(),
            fa,
            fb
        );

        // ── 3c. 两端点函数值同号 — 区间内可能无根，取更接近零的一端 ──
        if fa * fb > 0.0 {
            let theta = if fa.abs() < fb.abs() { xa } else { xb };
            eprintln!(
                "    fa·fb > 0 (no sign change in bracket), using best endpoint θ = {:.3}°",
                theta.to_degrees()
            );
            return theta;
        }

        // ── 4. TOMS748 寻根 ──
        eprintln!("    solving with TOMS748 (tol = 0.01°, max_iter = 20) ...");
        let theta_tol = 0.01_f64.to_radians();
        let max_iter = 20;

        // 包装函数：每次求值输出进度
        use std::cell::Cell;
        let eval_count = Cell::new(0u32);

        let f = |theta: f64| -> f64 {
            let n = eval_count.get() + 1;
            eval_count.set(n);
            match cal_cond_diff(theta, true) {
                Some(diff) => diff,
                None => {
                    eprintln!(
                        "    TOMS748 iter {}: θ = {:.3}°  →  computation failed, using NaN",
                        n,
                        theta.to_degrees()
                    );
                    f64::NAN
                }
            }
        };

        match toms748::solve_bracket(xa, xb, &f, Tolerance::new(theta_tol, theta_tol), max_iter) {
            Ok(bracket) => {
                let theta_opt = 0.5 * (bracket.lo + bracket.hi);

                eprintln!(
                    "    TOMS748 converged in {} iterations, bracket = [{:.6}°, {:.6}°]",
                    bracket.iterations,
                    bracket.lo.to_degrees(),
                    bracket.hi.to_degrees()
                );

                // 检查是否收敛到边界（可能未真正找到根）
                if (theta_opt - xa).abs() < theta_tol || (theta_opt - xb).abs() < theta_tol {
                    eprintln!(
                        "    WARNING: solution hit bracket boundary [{:.3}°, {:.3}°], root may be outside",
                        xa.to_degrees(),
                        xb.to_degrees()
                    );
                }
                theta_opt
            }
            Err(_) => {
                let theta = if fa.abs() < fb.abs() { xa } else { xb };
                eprintln!(
                    "    TOMS748 failed to converge, using best endpoint θ = {:.3}°",
                    theta.to_degrees()
                );
                theta
            }
        }
    }

    /// 二分搜索确定转向段初始线子集，使出口 x 接近目标长度
    fn cal_transition_to_target_length(&mut self) {
        let exp_lines = self.sections[2].get_charlines();
        let exp_last = match exp_lines.last() {
            Some(l) => l.clone(),
            None => return,
        };
        if exp_last.len() < 2 {
            return;
        }

        match run_transition_to_target_length(
            &exp_last,
            self.unitprocess.as_ref(),
            &self.config,
            self.cal_exit_factory,
        ) {
            Some(ts) => {
                self.sections[3] = Box::new(ts);
            }
            None => {
                // 转向段无法收敛到目标长度，保留现有结果
            }
        }
    }

    pub fn get_assembly_charlines(&self) -> CharLines {
        let mut lines: CharLines = CharLines::new();

        // ── 段 0～2：初值线 + 初值问题 + 膨胀段 ──
        for i in 0..3 {
            lines.extend(self.sections[i].get_charlines());
        }

        // ── 段 3 + 均一区：逐对合并转向段特征线与均一区扩展线 ──
        let trans_lines = self.sections[3].get_charlines();
        let uniform_lines = self.uniform_section.get_charlines();

        if uniform_lines.is_empty() {
            // 均一区未产生数据时直接追加转向段结果
            lines.extend(trans_lines);
        } else {
            let n = trans_lines.len().min(uniform_lines.len());
            for i in 0..n {
                let mut combined =
                    CharLine::with_capacity(trans_lines[i].len() + uniform_lines[i].len());
                for pt in trans_lines[i].iter() {
                    combined.push(pt.clone());
                }
                for pt in uniform_lines[i].iter() {
                    combined.push(pt.clone());
                }
                lines.push(combined);
            }
            // 剩余转向段线（如有）
            for i in n..trans_lines.len() {
                lines.push(trans_lines[i].clone());
            }
        }

        // ── 出口边界线：收集所有截面段的截短线并合并 ──
        let mut all_cuts: Vec<CharLine> = Vec::new();
        // 初值线 + 初值问题 + 膨胀段 + 转向段
        for i in 0..4 {
            let cut = self.sections[i].get_line_cut();
            if !cut.is_empty() {
                all_cuts.push(cut);
            }
        }
        // 均一区
        let uniform_cut = self.uniform_section.get_line_cut();
        if !uniform_cut.is_empty() {
            all_cuts.push(uniform_cut);
        }

        let exit_line = cal_exit_line(&all_cuts, self.config.geometry.length);
        if !exit_line.is_empty() {
            lines.push(exit_line);
        }

        lines
    }
}

/// 合并所有截面段的截短线为统一的出口边界线。
///
/// 收集所有 section 的 `get_line_cut()` 输出，合并去重，
/// 补充 y=0（对称轴）端点，统一修正所有点 x = length。
fn cal_exit_line(all_cuts: &[CharLine], length: f64) -> CharLine {
    let mut exit_line = CharLine::new();

    for cut in all_cuts {
        for pt in cut.iter() {
            let mut p = pt.clone();
            p.x = length;
            exit_line.push(p);
        }
    }

    if exit_line.is_empty() {
        return exit_line;
    }

    // 按 y 降序排列（壁面在上，对称轴在下）
    exit_line.sort_by(|a, b| b.y.partial_cmp(&a.y).unwrap_or(std::cmp::Ordering::Equal));

    // 剔除 y 坐标过于接近的重复点
    exit_line.dedup_by(|a, b| (a.y - b.y).abs() < 1e-10);

    // 补充 y = 0 点（对称轴点）
    if exit_line.last().map(|p| p.y).unwrap_or(0.0) > 0.0 {
        let n = exit_line.len();
        if n >= 2 {
            let p1 = &exit_line[n - 2];
            let p2 = &exit_line[n - 1];
            let mut p = p2.clone();
            p.y = 0.0;
            p.x = length;
            p = p.interpolate_along(p1, p2);
            exit_line.push(p);
        } else {
            let mut p = exit_line[0].clone();
            p.y = 0.0;
            p.x = length;
            exit_line.push(p);
        }
    }

    // 统一强制修正所有点 x = length，确保出口为严格竖直直线
    for pt in exit_line.iter_mut() {
        pt.x = length;
    }

    exit_line
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
            material: Material::from_rgas_gamma(287.042, 1.4), // 使用理想空气
            inlet: Inlet::default(),
            geometry: Geometry::default(),
            throat: Throat {
                radius_throat: 0.0,
                theta_a: f64::NAN, // 让算法自动选择初始膨胀角
            },
            outlet: Outlet::default(),
            io: IO::default(),
        };
        let target_len = config.geometry.length;
        let mut n = ConstraintNozzle::new_otn(config);
        n.run();
        let lines = n.get_assembly_charlines();

        // 写入文件以便可视化检查
        let mut output_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
        output_dir.push("../../target/tmp/fluid_field.txt");
        lines.write_to_file(output_dir, false).unwrap();

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

        // 验证出口 x 不超出目标长度，且至少到达目标附近
        let mut max_x: f64 = 0.0;
        for line in lines.iter() {
            for point in line.iter() {
                max_x = max_x.max(point.x);
                assert!(
                    point.x < target_len + 1.0,
                    "x={} 超出长度 {target_len}",
                    point.x
                );
            }
        }
        assert!(
            max_x >= target_len - 1e-4,
            "出口未到达目标长度: max_x={:.6}, target={:.6}",
            max_x,
            target_len
        );
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

    /// 诊断测试：非零喉部半径下膨胀段-转向段壁面连续性检查
    #[test]
    fn test_wall_continuity_with_arc() {
        let config = NozzleConfig {
            control: Control::default(),
            material: Material::from_rgas_gamma(287.042, 1.4),
            inlet: Inlet::default(),
            geometry: Geometry::default(),
            throat: Throat {
                radius_throat: 0.5,             // 非零喉部过渡圆弧半径
                theta_a: 25.0_f64.to_radians(), // 25° 初始膨胀角
            },
            outlet: Outlet::default(),
            io: IO::default(),
        };
        let target_len = config.geometry.length;
        let mut n = ConstraintNozzle::new_otn(config);
        n.run();

        let exp_lines = n.sections[2].get_charlines();
        let trans_lines = n.sections[3].get_charlines();

        println!("\n=== 膨胀段-转向段壁面连续性诊断 ===");
        println!("膨胀段特征线条数: {}", exp_lines.len());
        println!("转向段特征线条数: {}", trans_lines.len());

        // 膨胀段最后一个壁面点
        if let Some(exp_last) = exp_lines.last() {
            let exp_wall = &exp_last[0];
            println!(
                "\n膨胀段最后壁面点: x={:.6}, y={:.6}, θ={:.4}°",
                exp_wall.x,
                exp_wall.y,
                exp_wall.flow_direction().to_degrees()
            );
        }

        // 转向段第一个壁面点
        if let Some(trans_first) = trans_lines.first() {
            if !trans_first.is_empty() {
                let trans_wall = &trans_first[0];
                println!(
                    "转向段第一个壁面点: x={:.6}, y={:.6}, θ={:.4}°",
                    trans_wall.x,
                    trans_wall.y,
                    trans_wall.flow_direction().to_degrees()
                );
            } else {
                println!("转向段第一条特征线为空！");
            }
        }

        // 计算两段壁面点之间的距离
        if let (Some(exp_last), Some(trans_first)) = (exp_lines.last(), trans_lines.first()) {
            if !exp_last.is_empty() && !trans_first.is_empty() {
                let exp_wall = &exp_last[0];
                let trans_wall = &trans_first[0];
                let gap = exp_wall.distance_to(trans_wall);
                let dx = trans_wall.x - exp_wall.x;
                let dy = trans_wall.y - exp_wall.y;
                println!("\n壁面点间距: {:.6e} m", gap);
                println!("  Δx = {:.6e} m", dx);
                println!("  Δy = {:.6e} m", dy);

                // 打印膨胀段最后5个壁面点
                println!("\n膨胀段最后5个壁面点:");
                let start = exp_lines.len().saturating_sub(5);
                for i in start..exp_lines.len() {
                    let p = &exp_lines[i][0];
                    println!(
                        "  [{:>3}] x={:>10.6}, y={:>10.6}, θ={:>6.2}°",
                        i,
                        p.x,
                        p.y,
                        p.flow_direction().to_degrees()
                    );
                }

                // 打印转向段前5个壁面点
                println!("\n转向段前5个壁面点:");
                let end = (trans_lines.len()).min(5);
                for i in 0..end {
                    let p = &trans_lines[i][0];
                    println!(
                        "  [{:>3}] x={:>10.6}, y={:>10.6}, θ={:>6.2}°",
                        i,
                        p.x,
                        p.y,
                        p.flow_direction().to_degrees()
                    );
                }

                if gap > 0.01 {
                    println!("\n⚠️  壁面点间距过大 ({:.4} m)，可能存在不连续！", gap);
                } else if gap > 1e-6 {
                    println!("\n⚠️  壁面有微小间隙 ({:.6e} m)，可能是求解器衔接问题", gap);
                } else {
                    println!("\n✅ 壁面连续，间距在容差范围内");
                }

                // 检查膨胀段壁面点是否都在圆弧上
                println!("\n膨胀段壁面点与理论圆弧的偏差:");
                let x0 = exp_lines[0][0].x;
                let y0 = exp_lines[0][0].y;
                let r_t = 0.5;
                for i in 0..exp_lines.len() {
                    let p = &exp_lines[i][0];
                    // 理论圆弧: x = x0 + r_t*sin(θ), y = y0 + r_t*(1-cos(θ))
                    // 反推 theta = asin((x - x0) / r_t)
                    let sin_theta = ((p.x - x0) / r_t).clamp(-1.0, 1.0);
                    let theta_arc = sin_theta.asin();
                    let x_theory = x0 + r_t * theta_arc.sin();
                    let y_theory = y0 + r_t * (1.0 - theta_arc.cos());
                    let err = ((p.x - x_theory).powi(2) + (p.y - y_theory).powi(2)).sqrt();
                    if err > 1e-6 {
                        println!(
                            "  [{:>3}] 偏差={:.3e} m  (实际 x={:.6} y={:.6}  理论 x={:.6} y={:.6})",
                            i, err, p.x, p.y, x_theory, y_theory
                        );
                    }
                }
            }
        }

        // 写入输出文件
        let lines = n.get_assembly_charlines();
        let mut output_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
        output_dir.push("../../target/tmp/fluid_field_continuity_check.txt");
        lines.write_to_file(output_dir, false).unwrap();

        // 验证所有点有效
        let mut max_x: f64 = 0.0;
        for (li, line) in lines.iter().enumerate() {
            for (pi, point) in line.iter().enumerate() {
                assert!(point.is_valid(), "无效点: line={li} pt={pi}");
                assert!(
                    point.x >= -1e-9 && point.y >= -1e-9,
                    "坐标越界: line={li} pt={pi}: x={}, y={}",
                    point.x,
                    point.y
                );
                max_x = max_x.max(point.x);
            }
        }
        assert!(
            max_x >= target_len - 1e-4,
            "出口未到达目标长度: max_x={:.6}, target={:.6}",
            max_x,
            target_len
        );
    }
}
