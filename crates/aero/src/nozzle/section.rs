use crate::{
    moc::{CharLine, CharLines, unitprocess::UnitProcess},
    nozzle::NozzleConfig,
};

/// 表示喷管中的一个截面段，负责执行特征线计算
///
/// 每个 `Section` 实例代表喷管几何中的一个特定区域，
/// 通过特征线方法进行流场计算。
///
/// # 生命周期
///
/// Section 遵循 **先设置边界条件，再运行** 的模式：
/// 1. 构造或通过 `inherit_last_line` 传入上一段的出口线作为本段入口条件
/// 2. 调用 `run()` 执行特征线推进计算
/// 3. 调用 `get_charlines()` 获取计算结果
///
/// 如果 `run()` 尚未被调用，`get_charlines()` 返回**空特征线集合**（而非 panic）。
/// 调用者应通过 `char_lines.is_empty()` 判断是否已有可用结果。
pub trait Section {
    /// 执行当前截面段的计算步骤
    ///
    /// 使用给定的单元过程和喷管配置参数，
    /// 更新截面段的内部状态和特征线数据。
    ///
    /// # 参数
    /// * `unitprocess` - 单元过程对象，提供超音速流的特征线法计算方法
    /// * `config` - 喷管配置参数，包含几何和物理参数等
    ///
    /// # 幂等性
    ///
    /// 重复调用 `run()` 会**重新计算**并覆盖之前的结果。
    /// 不会跳过已计算状态（无内部 `calculated` 标记）。
    fn run(&mut self, unitprocess: &dyn UnitProcess, config: &NozzleConfig);

    /// 获取当前截面段的特征线数据
    ///
    /// 返回当前已有的特征线集合。如果 `run()` 尚未被调用，
    /// 返回空集合（不包含之前段的初始线）。
    ///
    /// # 返回值
    /// 指向当前截面段所有特征线的不可变引用
    fn get_charlines(&self) -> &CharLines;

    /// 将上一截面段最后一条特征线作为当前段的初始边界条件传入。
    ///
    /// 默认实现为空操作。只有需要从前一段接收数据的截面段
    /// （如 `InitialSection` 需要 `InitialLine` 的输出）才需覆盖此方法。
    fn inherit_last_line(&mut self, _line: &CharLine) {}

    /// 设置初始膨胀角（仅 ExpansionSection 使用）
    ///
    /// 用于自动迭代选取 theta_a 时重新运行膨胀段。
    /// 默认实现为空操作。
    fn set_theta_a(&mut self, _theta_a: f64) {}

    /// 获取该截面段在出口边界（`x = geometry.length`）上的贡献点。
    ///
    /// 返回该段计算出的、位于或截断于出口边界处的特征线点。
    /// 不同截面段的出口边界语义：
    /// - `ExpansionSection`：膨胀线被长度截断后的右端点集合（壁面→对称轴）
    /// - `TransitionSection`：特征线收敛为单点后的唯一点（最终出口点）
    /// - `UniformSection`：所有延伸到出口边界处的特征线末点
    /// - `InitialLine` / `InitialSection`：空（不到达出口边界）
    ///
    /// 默认实现返回空特征线。
    fn exit_boundary_segment(&self) -> CharLine {
        CharLine::new()
    }

    /// 获取起始线长度（仅 TransitionSection 使用）
    ///
    /// 返回 `line_init.len()`，用于均一区计算其起始线在膨胀段
    /// 最后特征线中的分割位置。默认返回 0。
    fn line_init_len(&self) -> usize {
        0
    }
}
