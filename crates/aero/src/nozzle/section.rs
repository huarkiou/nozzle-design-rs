use crate::{
    moc::{unitprocess::UnitProcess, CharLine, CharLines},
    nozzle::NozzleConfig,
};

/// 表示喷管中的一个截面段，负责执行特征线计算
///
/// 每个 `Section` 实例代表喷管几何中的一个特定区域，
/// 通过特征线方法进行流场计算。
pub trait Section {
    /// 执行当前截面段的计算步骤
    ///
    /// 使用给定的单元过程和喷管配置参数，
    /// 更新截面段的内部状态和特征线数据。
    ///
    /// # 参数
    /// * `unitprocess` - 单元过程对象，提供超音速流的特征线法计算方法
    /// * `config` - 喷管配置参数，包含几何和物理参数等
    fn run(&mut self, unitprocess: &dyn UnitProcess, config: &NozzleConfig);

    /// 获取当前截面段的特征线数据，不包括作为初始边界条件的第一条特征线，仅包括通过run计算出来的新结果
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

    /// 获取截短线（仅 ExpansionSection 使用）
    ///
    /// 膨胀段中因超出最大允许长度而被截断产生的右边界线，
    /// 从上壁面排到对称轴。
    /// 默认实现返回空特征线。
    fn get_line_cut(&self) -> CharLine {
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
