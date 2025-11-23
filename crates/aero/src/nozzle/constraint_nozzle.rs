use math::Tolerance;

use crate::{
    moc::{
        AreaType,
        unitprocess::{GeneralConfig, Irrotational, Rotational, UnitProcess},
    },
    nozzle::{NozzleConfig, Section},
};

pub struct ConstraintNozzle {
    config: NozzleConfig,
    unitprocess: Box<dyn UnitProcess>,
    sections: Vec<Box<dyn Section>>,
}

impl ConstraintNozzle {
    /// 构建最大推力喷管OTN
    pub fn new_otn(config: NozzleConfig) -> Self {
        let unitprocess_config = GeneralConfig {
            axisym: if config.control.axisymmetric {
                AreaType::Axisymmetric
            } else {
                AreaType::Planar(config.geometry.width)
            },
            tol: Tolerance::new(config.control.eps, config.control.eps),
            n_corr: config.control.n_correction_max,
        };

        let sections: Vec<Box<dyn Section>> = Vec::new();

        Self {
            unitprocess: if config.control.irrotational {
                Box::new(Irrotational {
                    conf: unitprocess_config,
                })
            } else {
                Box::new(Rotational {
                    conf: unitprocess_config,
                })
            },
            config: config,
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
