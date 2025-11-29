use crate::{
    moc::unitprocess::UnitProcess,
    nozzle::{NozzleConfig, Section, initial_line::InitialLine},
};

pub struct ConstraintNozzle {
    config: NozzleConfig,
    unitprocess: Box<dyn UnitProcess>,
    sections: Vec<Box<dyn Section>>,
}

impl ConstraintNozzle {
    /// 构建最大推力喷管OTN
    pub fn new_otn(config: NozzleConfig) -> Self {
        let mut sections: Vec<Box<dyn Section>> = Vec::new();
        sections.push(Box::new(InitialLine::new()));
        Self {
            unitprocess: config.to_unitprocess(),
            config,
            sections,
        }
    }

    /// 计算喷管各段内流场数据
    ///
    /// 按顺序对每个段执行特征线法计算步骤。
    ///
    /// 计算顺序为 `self.sections` 中存储的顺序，
    /// 前一个截面段的计算结果作为下一截面段的初始边界条件。
    pub fn run(&mut self) {
        for section in self.sections.iter_mut() {
            section.run(self.unitprocess.as_ref(), &self.config);
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::{Material, nozzle::config::*};

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
        let mut n = ConstraintNozzle::new_otn(config);
        n.run();
    }
}
