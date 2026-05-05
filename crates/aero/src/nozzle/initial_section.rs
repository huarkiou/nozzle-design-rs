use crate::{
    moc::{
        unitprocess::{Context, UnitProcess},
        CharLine, CharLines, MocPoint,
    },
    nozzle::Section,
};

/// 超声速特征线法的初值问题计算段
///
/// 从给定的初值线(`line_init`)开始，逐条计算其上的右行特征线，
/// 通过迭代推进的方式求解超声速流场。
///
/// 每条右行特征线按从壁面(下标0)到对称轴(末下标)的顺序排列。
pub struct InitialSection {
    /// 依次为初值线上各点发出的右行特征线
    pub char_lines: CharLines,
    /// 开始计算的起始线(初值线)，应在调用 [`run`](Section::run) 之前设置
    pub line_init: CharLine,
}

impl InitialSection {
    /// 创建一个空的初值段，调用 [`run`](Section::run) 前需设置 `line_init`
    pub fn new() -> Self {
        Self {
            char_lines: CharLines::new(),
            line_init: CharLine::new(),
        }
    }

    /// 创建并指定初值线
    pub fn with_line_init(line_init: CharLine) -> Self {
        Self {
            char_lines: CharLines::new(),
            line_init,
        }
    }

    /// 设置初值线，并清除已有的计算结果
    pub fn set_line_init(&mut self, line_init: CharLine) {
        self.line_init = line_init;
        self.char_lines = CharLines::new();
    }

    /// 根据前一条右行特征线和当前右行特征线的第一点（壁面点），
    /// 计算当前右行特征线。
    ///
    /// # 参数
    /// * `unitprocess` - 控制特征线法基本过程的 trait 对象
    /// * `line_prev` - 前一条右行特征线，最后一点为对称轴点(轴对称问题)
    ///                  或下壁面点(二维平面问题)
    /// * `p_init` - 当前右行特征线的已知初始壁面点
    ///
    /// # 返回值
    /// 计算得到的当前右行特征线，从壁面点(下标0)排到对称轴点(末下标)
    fn cal_next_ivp_line(
        unitprocess: &dyn UnitProcess,
        line_prev: &CharLine,
        p_init: &MocPoint,
    ) -> CharLine {
        let n_prev = line_prev.len();
        // 动态构建当前特征线，避免用 MocPoint::default() 预填 NaN
        let mut line_cur = CharLine::with_capacity(n_prev + 2);

        // j=0：起始边界条件 — 已知壁面点
        line_cur.push(p_init.clone());

        // j=1..n_prev：流场内点
        // 对前一线上的每个点，用右行(p1)+左行(p2)特征线相交求解
        for idx_prev in 0..n_prev {
            let context = Context {
                prev: line_prev,
                next: &line_cur,
                idx_prev,
                idx_next: line_cur.len() - 1, // 当前线上刚算出的前一点
            };
            if let Some(point) = unitprocess.interior_point(context) {
                if point.is_valid() {
                    line_cur.push(point);
                }
                // 结果无效时跳过此点，不阻塞后续计算
            }
        }

        // 最后一个点：对称轴边界条件
        if !line_cur.is_empty() {
            let context = Context {
                prev: line_prev,
                next: &line_cur,
                idx_prev: n_prev.saturating_sub(1), // line_prev 的最后一点
                idx_next: line_cur.len() - 1,
            };
            if let Some(point) = unitprocess.symmetry_axis_point(context) {
                if point.is_valid() {
                    line_cur.push(point);
                }
            }
        }

        line_cur
    }
}

impl Section for InitialSection {
    fn run(&mut self, unitprocess: &dyn UnitProcess, _config: &super::NozzleConfig) {
        assert!(
            !self.line_init.is_empty(),
            "InitialSection: line_init 必须在 run() 之前设置"
        );

        self.char_lines = CharLines::new();

        // 从对称轴点开始（初值线反转后，line_init[last] = 对称轴点 y=0）。
        // interior_point 的约定：p1（next 线）的 y 必须高于 p2（prev 线）。
        // 从轴→壁面迭代可保证 line_cur（新建线）始终在 line_prev 上方，
        // 特征线交点落在下游 x>0 且不会越过中心线 y<0。
        let mut line_cur = CharLine::new();
        line_cur.push(self.line_init.last().unwrap().clone());

        // 沿初值线从对称轴向壁面方向反向迭代
        for point in self.line_init.iter().rev().skip(1) {
            line_cur = Self::cal_next_ivp_line(unitprocess, &line_cur, point);
            if !line_cur.is_empty() {
                self.char_lines.push(line_cur.clone());
            }
        }
    }

    fn get_charlines(&self) -> &CharLines {
        &self.char_lines
    }

    fn inherit_last_line(&mut self, line: &CharLine) {
        self.set_line_init(line.clone());
    }
}
