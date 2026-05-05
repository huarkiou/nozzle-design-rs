use crate::{
    moc::{
        CharLine, CharLines, MocPoint,
        unitprocess::{Context, UnitProcess},
    },
    nozzle::Section,
};

/// 超声速特征线法的初始膨胀段
///
/// 从给定的初始线（来自初值段的最后一条右行特征线）开始，
/// 沿喉部过渡圆弧逐条计算膨胀特征线。
///
/// 每条右行特征线按从壁面(下标0)到对称轴(末下标)的顺序排列。
pub struct ExpansionSection {
    /// 依次为初始线上各点发出的右行特征线（含 line_init 作为第一条）
    pub char_lines: CharLines,
    /// 开始计算的起始线（来自上一段的最后一条特征线）
    pub line_init: CharLine,
    /// 喉部过渡圆弧半径 (m)
    pub radius_throat: f64,
    /// 初始膨胀角 (rad)
    pub theta_a: f64,
    /// 最大允许 x 坐标 (m)，超过将被截短
    pub max_length: f64,
    /// 超出最大长度的截短线
    pub exit_boundary_segment: CharLine,
}

impl ExpansionSection {
    /// 创建膨胀段，调用 [`run`](Section::run) 前需设置 `line_init`
    pub fn new(radius_throat: f64, theta_a: f64, max_length: f64) -> Self {
        Self {
            char_lines: CharLines::new(),
            line_init: CharLine::new(),
            radius_throat,
            theta_a,
            max_length,
            exit_boundary_segment: CharLine::new(),
        }
    }

    /// 设置初始线，并清除已有的计算结果
    pub fn set_line_init(&mut self, line_init: CharLine) {
        self.line_init = line_init;
        self.char_lines = CharLines::new();
        self.exit_boundary_segment = CharLine::new();
    }

    /// 根据前一条右行特征线，和当前右行特征线壁面点坐标和切向，
    /// 计算当前右行特征线。
    ///
    /// # 参数
    /// * `unitprocess` - 特征线法基本过程
    /// * `line_prev` - 前一条右行特征线
    /// * `p_init` - 当前壁面点（已知 x, y, θ 的喉部圆弧点）
    /// * `max_length` - 最大允许 x 坐标
    ///
    /// # 返回值
    /// 计算得到的当前右行特征线，从壁面点(下标0)排到对称轴点(末下标)
    fn cal_next_expansion_line(
        unitprocess: &dyn UnitProcess,
        line_prev: &CharLine,
        p_init: &MocPoint,
        max_length: f64,
    ) -> CharLine {
        let n_prev = line_prev.len();
        let mut line_cur = CharLine::with_capacity(n_prev + 1);

        // ── 步骤 1：壁面点（逆推壁面点） ──
        // 从 p_init（喉部圆弧上的已知点）出发，
        // 寻找 line_prev 上合适的区间使左行特征线通过 p_init
        let wall_info = (p_init.x, p_init.y, p_init.flow_direction());
        let mut wall_point = p_init.clone();
        // start_idx 至少为 1，避免 interior_point 将壁面点与自身比较
        let mut start_idx: usize = 1;

        for idx_prev in 0..n_prev.saturating_sub(1) {
            let mut tmp_cur = CharLine::new();
            tmp_cur.push(wall_point.clone());

            let context = Context {
                prev: line_prev,
                next: &tmp_cur,
                idx_prev,
                idx_next: 0,
            };
            if let Some(refined) = unitprocess.inverse_wall_point(context, wall_info) {
                // 壁面点不应跑到上游 (x<0) 或超出最大长度
                if refined.is_valid() && refined.x >= -1e-9 && refined.x < max_length {
                    wall_point = refined;
                    start_idx = idx_prev + 1;
                    break;
                }
            }
        }

        line_cur.push(wall_point.clone());

        // ── 步骤 2：流场内点 ──
        // 从 start_idx 开始，对 line_prev 上剩余的点依次计算 interior_point
        for idx_prev in start_idx..n_prev {
            let context = Context {
                prev: line_prev,
                next: &line_cur,
                idx_prev,
                idx_next: line_cur.len() - 1,
            };
            if let Some(point) = unitprocess.interior_point(context) {
                if point.is_valid() && point.y > 0.0 {
                    // 首个内点若在壁面上游，跳过（对应 C++ cal_next_expansion_line L149）
                    if line_cur.len() == 1 && line_cur[0].x > point.x {
                        continue;
                    }
                    let x_check = point.x;
                    line_cur.push(point);

                    // 如果当前点 x 已超出 max_length，截短后停止
                    if x_check >= max_length {
                        // 截短最后一点到 max_length
                        let n = line_cur.len();
                        if n >= 2 {
                            let p_prev = &line_cur[n - 2];
                            let p_last = &line_cur[n - 1];
                            let ratio = (max_length - p_prev.x) / (p_last.x - p_prev.x);
                            let mut p_cut = p_last.clone();
                            p_cut.x = max_length;
                            p_cut.y = p_prev.y + ratio * (p_last.y - p_prev.y);
                            p_cut = p_cut.interpolate_along(p_prev, p_last);
                            line_cur.pop();
                            line_cur.push(p_cut);
                        }
                        break;
                    }
                }
            }
        }

        // ── 步骤 3：对称轴边界条件 ──
        // 始终尝试计算对称轴点（不因截短而跳过），确保特征线到达 y=0
        if !line_cur.is_empty() {
            let context = Context {
                prev: line_prev,
                next: &line_cur,
                idx_prev: n_prev.saturating_sub(1),
                idx_next: line_cur.len() - 1,
            };
            if let Some(point) = unitprocess.symmetry_axis_point(context) {
                if point.is_valid() {
                    if point.x > max_length {
                        let p_prev = line_cur.last().unwrap();
                        let ratio = (max_length - p_prev.x) / (point.x - p_prev.x);
                        let mut p_cut = point.clone();
                        p_cut.x = max_length;
                        p_cut.y = p_prev.y + ratio * (point.y - p_prev.y);
                        p_cut = p_cut.interpolate_along(p_prev, &point);
                        line_cur.push(p_cut);
                    } else {
                        line_cur.push(point);
                    }
                }
            }
        }

        // 确保最后一点满足约定：在对称轴上(y=0) 或 在右边界上(x=max_length)
        if let Some(last) = line_cur.last() {
            if last.y > 0.0 && (last.x - max_length).abs() > 1e-9 {
                // 既不在轴上也不在右边界上 → 强制截断到右边界
                if line_cur.len() >= 2 {
                    let p_prev = &line_cur[line_cur.len() - 2];
                    let ratio = (max_length - p_prev.x) / (last.x - p_prev.x);
                    let mut p_cut = last.clone();
                    p_cut.x = max_length;
                    p_cut.y = p_prev.y + ratio * (p_cut.x - p_prev.x);
                    p_cut = p_cut.interpolate_along(p_prev, last);
                    *line_cur.last_mut().unwrap() = p_cut;
                }
            }
        }

        line_cur
    }
}

impl Section for ExpansionSection {
    fn set_theta_a(&mut self, theta_a: f64) {
        self.theta_a = theta_a;
    }

    fn run(&mut self, unitprocess: &dyn UnitProcess, _config: &super::NozzleConfig) {
        assert!(
            !self.line_init.is_empty(),
            "ExpansionSection: line_init 必须在 run() 之前设置"
        );

        self.char_lines = CharLines::new();

        // 喉部圆弧起点 = 初始线的壁面点（line_init 反转后 [0]=壁面）
        let p_wall = self.line_init.first().unwrap();
        let x0 = p_wall.x;
        let y0 = p_wall.y;
        let _theta_0 = p_wall.flow_direction();

        let r_t = self.radius_throat;
        let theta_a = self.theta_a;

        // 无效膨胀角时直接跳过（NaN 或 ≤0 均无膨胀），
        // 但仍将 line_init 传递给下游段
        if !theta_a.is_finite() || theta_a <= 0.0 {
            self.char_lines.push(self.line_init.clone());
            return;
        }

        // 初始角度增量
        const DTHETA_MIN: f64 = 0.1_f64.to_radians();
        const DTHETA_MAX: f64 = 0.5_f64.to_radians();
        let mut dtheta = DTHETA_MIN;
        let mut need_dtheta = true;
        let mut is_ending = false;
        let mut iter_count = 0u32;
        const MAX_ITER: u32 = 500;

        // 从 line_init 开始逐条计算膨胀特征线
        let mut prev_line = self.line_init.clone();
        let mut prev_theta = prev_line[0].flow_direction();

        // 主循环：沿喉部圆弧逐条计算膨胀特征线
        while prev_theta < theta_a - 1e-12 {
            iter_count += 1;
            if iter_count > MAX_ITER {
                // 防止死循环：达到最大迭代次数后强制退出
                break;
            }
            let target_theta = if is_ending {
                // 结束阶段：直接以 theta_a 为目标，避免反复缩小步长导致特征线过密
                theta_a
            } else {
                prev_theta + dtheta
            };

            // 喉部圆弧上的壁面点（速度方向沿圆弧切线 target_theta）
            let velo = prev_line[0].velocity();
            let p_init = MocPoint {
                x: x0 + r_t * target_theta.sin(),
                y: y0 + r_t * (1.0 - target_theta.cos()),
                u: velo * target_theta.cos(),
                v: velo * target_theta.sin(),
                p: prev_line[0].p,
                rho: prev_line[0].rho,
                t: prev_line[0].t,
                mat: prev_line[0].mat.clone(),
            };

            let line_cur =
                Self::cal_next_expansion_line(unitprocess, &prev_line, &p_init, self.max_length);

            if line_cur.is_empty() {
                break;
            }

            // 若壁面角度超出目标值，进入结束阶段并重试
            if !is_ending && line_cur[0].flow_direction() > theta_a {
                is_ending = true;
                continue; // 放弃超出 theta_a 的线，下一轮直接用 theta_a
            }

            // dtheta 自适应调整 — 基于轴线侧末两点间距反馈
            // 仅在线未截断（最后一点在对称轴上）时执行，
            // 截断线的末两点不在轴线侧，其间距不能反映特征线密度
            if need_dtheta
                && line_cur.len() >= 3
                && line_cur.last().map(|p| p.y <= 0.0).unwrap_or(false)
            {
                let last_point = line_cur.last().unwrap();
                let l_prev = prev_line
                    .last()
                    .unwrap()
                    .distance_to(&prev_line[prev_line.len() - 2]);
                let l_cur = last_point.distance_to(&line_cur[line_cur.len() - 2]);
                let target_l = l_prev * line_cur[0].y / prev_line[0].y;

                if (l_cur - target_l).abs() > 1e-6 {
                    let ratio = target_l / l_cur;
                    // 阻尼修正：70% 保持 + 30% 朝目标方向，避免振荡
                    dtheta *= 0.7 + 0.3 * ratio;
                    dtheta = dtheta.clamp(DTHETA_MIN * 0.2, DTHETA_MAX);
                }

                if dtheta + line_cur[0].flow_direction() >= theta_a {
                    need_dtheta = false;
                }
            }

            let reached_theta_a = (line_cur[0].flow_direction() - theta_a).abs() < 1e-9;

            self.char_lines.push(line_cur.clone());
            prev_theta = line_cur[0].flow_direction();
            prev_line = line_cur;

            if reached_theta_a {
                break;
            }
        }

        // ── 去重叠：若最后两条线壁面侧间距不足前面平均间距的一半，移除倒数第二条 ──
        let n = self.char_lines.len();
        if n >= 3 {
            let last_spacing = self.char_lines[n - 1][0].distance_to(&self.char_lines[n - 2][0]);
            // 计算前面若干条线的平均壁面间距
            let mut sum = 0.0;
            let mut count = 0;
            for i in (0..n.saturating_sub(2)).rev() {
                let d = self.char_lines[i + 1][0].distance_to(&self.char_lines[i][0]);
                sum += d;
                count += 1;
                if count >= 5 {
                    break;
                }
            }
            if count > 0 {
                let avg = sum / count as f64;
                if last_spacing < 0.5 * avg {
                    self.char_lines.remove(n - 2);
                }
            }
        }

        // ── 计算出口边界贡献点 ──
        self.exit_boundary_segment = CharLine::new();
        for line in self.char_lines.iter().rev() {
            let p = line.last().unwrap();
            self.exit_boundary_segment.push(p.clone());
            if p.y <= 0.0 {
                break;
            }
        }
        // 反转使 exit_boundary_segment 从上壁面到轴线排列
        self.exit_boundary_segment.reverse();
    }

    fn get_charlines(&self) -> &CharLines {
        &self.char_lines
    }

    fn inherit_last_line(&mut self, line: &CharLine) {
        self.set_line_init(line.clone());
    }

    fn exit_boundary_segment(&self) -> CharLine {
        self.exit_boundary_segment.clone()
    }
}
