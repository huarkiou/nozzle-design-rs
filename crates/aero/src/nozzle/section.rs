use crate::{
    moc::{CharLines, unitprocess::UnitProcess},
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
}
