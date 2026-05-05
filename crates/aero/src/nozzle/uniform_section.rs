use crate::moc::{
    unitprocess::{Context, UnitProcess},
    CharLine, CharLines, MocPoint,
};

/// MOC 容差（与 Control.eps 默认值一致）
const EPS: f64 = 1e-7;

/// 均一区（Uniform Section）：填补转向段起始线子集与膨胀段完整最后一条特征线之间的"缺口"。
///
/// 对应 C++ `UniformSection`（uniform-section.hpp）：
/// - `line_init` = 膨胀段最后一条特征线中未被转向段使用的下游点（靠近对称轴的部分）
/// - `line_exit` = 转向段每条特征线的最末点（出口点），按从壁面到对称轴排列
/// - `run()` 从每个出口点出发，沿 line_init 方向逐条计算右行特征线
/// - `cal_exit_line()` 整合出口边界线
pub struct UniformSection {
    /// 依次计算得到的右行特征线
    pub char_lines: CharLines,
    /// 膨胀段最后特征线中未被转向段使用的下游点（从壁面侧到对称轴侧排列）
    pub line_init: CharLine,
    /// 转向段出口左行特征线（每条转向段特征线的最后一个点）
    pub line_exit: CharLine,
    /// 喷管目标长度
    pub length: f64,
}

impl UniformSection {
    pub fn new(length: f64) -> Self {
        Self {
            char_lines: CharLines::new(),
            line_init: CharLine::new(),
            line_exit: CharLine::new(),
            length,
        }
    }

    /// 执行均一区计算
    ///
    /// 对于 line_exit 中每个出口点 p_init，在前一条特征线 line_prev
    /// 的基础上计算一条新的右行特征线。
    ///
    /// 首条 line_prev = line_init（膨胀段下游点），
    /// 后续 line_prev = 前一步算出的 line_cur。
    ///
    /// 如果 `line_init` 或 `line_exit` 为空，则直接返回不做计算。
    pub fn run(&mut self, unitprocess: &dyn UnitProcess) {
        if self.line_init.is_empty() || self.line_exit.is_empty() {
            return;
        }

        self.char_lines = CharLines::new();
        let mut line_prev = self.line_init.clone();

        for p_init in self.line_exit.iter() {
            let line_cur =
                Self::cal_next_uniform_line(unitprocess, &line_prev, p_init, self.length);
            if !line_cur.is_empty() {
                self.char_lines.push(line_cur.clone());
                line_prev = line_cur;
            }
        }
    }

    /// 计算出口边界特征线。
    ///
    /// 将膨胀段截短线 `line_cut` 与均一区出口点合并，
    /// 补充 y=0（对称轴）点和 y=ymax（壁面）点，按 y 降序排列。
    ///
    /// 最后统一将所有点的 x 强制修正为 `self.length`，
    /// 确保出口为严格竖直直线（对应 C++ uniform-section.hpp:24-29 注释）。
    pub fn cal_exit_line(&self, _unitprocess: &dyn UnitProcess, line_cut: &CharLine) -> CharLine {
        let mut exit_line = line_cut.clone();
        // 移除 line_cut 末尾的对称轴点（后续重新补充）
        if exit_line.len() > 1 {
            exit_line.pop();
        }

        // 添加均一区出口点
        for line in self.char_lines.iter() {
            if let Some(p) = line.last() {
                if (p.x - self.length).abs() < EPS {
                    let mut pt = p.clone();
                    pt.x = self.length;
                    exit_line.push(pt);
                }
            }
        }

        // 按 y 降序排列（壁面在上，对称轴在下）
        exit_line.sort_by(|a, b| b.y.partial_cmp(&a.y).unwrap_or(std::cmp::Ordering::Equal));

        // 补充 y = y_max 点（壁面点）
        if let Some(last_line) = self.char_lines.last() {
            if let Some(wall_pt) = last_line.first() {
                if wall_pt.y > exit_line.first().map(|p| p.y).unwrap_or(0.0) {
                    let mut p = wall_pt.clone();
                    p.x = self.length;
                    exit_line.insert(0, p);
                }
            }
        }

        // 补充 y = 0 点（对称轴点）
        if let Some(last_pt) = exit_line.last() {
            if last_pt.y > 0.0 {
                let p1 = &exit_line[exit_line.len() - 2];
                let p2 = last_pt;
                let mut p = p2.clone();
                p.y = 0.0;
                p.x = self.length;
                p = p.interpolate_along(p1, p2);
                exit_line.push(p);
            }
        }

        // 统一强制修正所有点 x = length，确保出口为严格竖直直线
        for pt in exit_line.iter_mut() {
            pt.x = self.length;
        }

        exit_line
    }

    /// 已知前一条特征线 `line_prev` 和出口点 `p_init`，
    /// 计算一条从出口延伸到 line_prev 区域的右行特征线。
    ///
    /// # 参数
    /// - `line_prev`: 前一条右行特征线（首条为膨胀段下游点，后续为前一步结果）
    /// - `p_init`: 当前出口初始点（来自转向段出口左行特征线）
    /// - `max_length`: 喷管目标长度
    fn cal_next_uniform_line(
        unitprocess: &dyn UnitProcess,
        line_prev: &CharLine,
        p_init: &MocPoint,
        max_length: f64,
    ) -> CharLine {
        let n_prev = line_prev.len();
        let mut line_cur = CharLine::with_capacity(n_prev + 1);
        // j: 因无效点跳过的偏移量
        let mut j: usize = 0;

        for i in 0..n_prev {
            // 已知点：i==0 时用 p_init，否则用 line_cur 上一步的结果
            let known = if i == 0 {
                p_init.clone()
            } else {
                let idx = i.saturating_sub(1).saturating_sub(j);
                if idx >= line_cur.len() {
                    break;
                }
                line_cur[idx].clone()
            };

            let mut tmp_cur = CharLine::with_capacity(2);
            tmp_cur.push(known.clone());

            let context = Context {
                prev: line_prev,
                next: &tmp_cur,
                idx_prev: i,
                idx_next: 0,
            };

            let mut point = match unitprocess.interior_point(context) {
                Some(p) => p,
                None => {
                    j += 1;
                    continue;
                }
            };

            // 剔除比 p_init 还靠左的点 和 无效点
            if !point.is_valid() || point.x <= p_init.x {
                j += 1;
                continue;
            }

            // 对称轴边界条件
            if let Some(last_prev) = line_prev.last() {
                if last_prev.y == 0.0 && (point.y <= 0.0 || !point.is_valid()) {
                    if i == 0 {
                        break;
                    }
                    let sym_idx = i.saturating_sub(1).saturating_sub(j);
                    if sym_idx >= line_cur.len() {
                        break;
                    }
                    let sym_known = line_cur[sym_idx].clone();
                    let mut sym_tmp = CharLine::new();
                    sym_tmp.push(sym_known);
                    let sym_ctx = Context {
                        prev: line_prev,
                        next: &sym_tmp,
                        idx_prev: i,
                        idx_next: 0,
                    };
                    match unitprocess.symmetry_axis_point(sym_ctx) {
                        Some(p) => point = p,
                        None => {
                            j += 1;
                            continue;
                        }
                    }
                }
            }

            // 超出喷管长度则截短
            if point.x > max_length {
                if i == 0 {
                    // 首点超出：从 p_init 到 point 线性插值到 max_length
                    let k = (p_init.y - point.y) / (p_init.x - point.x);
                    let orig = point.clone();
                    point.x = max_length;
                    point.y = p_init.y + k * (max_length - p_init.x);
                    point = point.interpolate_along(p_init, &orig);
                } else {
                    // 内点超出：在 line_prev[i-1] 和 line_prev[i] 之间二分找交点
                    let p1 = &line_prev[i - 1];
                    let _p2 = &line_prev[i];
                    // 简单截短
                    let ratio = (max_length - p1.x) / (point.x - p1.x);
                    point.x = max_length;
                    point.y = p1.y + ratio * (point.y - p1.y);
                    point = point.interpolate_along(p1, &point.clone());
                }
            }

            // 剔除位置倒退的点
            if !line_cur.is_empty() && point.x <= line_cur.last().unwrap().x {
                j += 1;
                continue;
            }

            if !point.is_valid() {
                j += 1;
                continue;
            }

            line_cur.push(point.clone());

            // 到达对称轴或到达出口长度则停止
            if point.y <= 0.0 || (point.x - max_length).abs() < 1e-9 {
                break;
            }
        }

        // ── 后处理：确保最后一点到达对称轴或出口 ──
        if !line_cur.is_empty() {
            let last_y = line_cur.last().unwrap().y;
            let last_x = line_cur.last().unwrap().x;

            if last_y > 0.0 && last_x < max_length {
                // 尝试补对称轴点
                if let Some(lp_last) = line_prev.last() {
                    if lp_last.y == 0.0 {
                        let sym_known = line_cur.last().unwrap().clone();
                        let mut sym_tmp = CharLine::new();
                        sym_tmp.push(sym_known.clone());
                        let sym_ctx = Context {
                            prev: line_prev,
                            next: &sym_tmp,
                            idx_prev: n_prev.saturating_sub(1),
                            idx_next: 0,
                        };
                        if let Some(p) = unitprocess.symmetry_axis_point(sym_ctx) {
                            if p.is_valid() {
                                let p_orig = p.clone();
                                let mut p_final = p;
                                if p_final.x > max_length {
                                    let ratio =
                                        (max_length - sym_known.x) / (p_final.x - sym_known.x);
                                    p_final.x = max_length;
                                    p_final.y = sym_known.y + ratio * (p_final.y - sym_known.y);
                                    p_final = p_final.interpolate_along(&sym_known, &p_orig);
                                }
                                line_cur.push(p_final);
                            }
                        }
                    }
                }
            }
        }

        // 剔除过于接近的相邻点
        if line_cur.len() > 3 {
            let n = line_cur.len();
            let p_last = &line_cur[n - 1];
            let p_prev = &line_cur[n - 2];
            let p_pprev = &line_cur[n - 3];
            let d1 = p_last.distance_to(p_prev);
            let d2 = p_prev.distance_to(p_pprev);
            if d1 < 0.2 * d2 {
                line_cur.remove(n - 2);
            }
        }

        line_cur
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_uniform_section_empty_init() {
        let us = UniformSection::new(6.0);
        assert!(us.line_init.is_empty());
        assert!(us.line_exit.is_empty());
    }
}
