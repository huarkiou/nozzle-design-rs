use crate::{
    moc::{unitprocess::UnitProcess, CharLines},
    nozzle::{initial_line::InitialLine, ExpansionSection, InitialSection, NozzleConfig, Section},
};

pub struct ConstraintNozzle {
    config: NozzleConfig,
    unitprocess: Box<dyn UnitProcess>,
    sections: Vec<Box<dyn Section>>,
}

impl ConstraintNozzle {
    /// 构建最大推力喷管OTN
    pub fn new_otn(config: NozzleConfig) -> Self {
        let length = config.geometry.length;
        let r_t = config.throat.radius_throat;
        let theta_a = config.throat.theta_a;
        Self {
            unitprocess: config.to_unitprocess(),
            config,
            sections: vec![
                Box::new(InitialLine::new()),
                Box::new(InitialSection::new()),
                Box::new(ExpansionSection::new(r_t, theta_a, length)),
            ],
        }
    }

    /// 计算喷管各段内流场数据
    ///
    /// 按顺序对每个段执行特征线法计算步骤，
    /// 前一个截面段的最后一条特征线作为下一截面段的初始边界条件传入。
    pub fn run(&mut self) {
        for i in 0..self.sections.len() {
            if i > 0 {
                // 先克隆上一段最后一条特征线（释放不可变借用），再传给当前段
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
        let mut n = ConstraintNozzle::new_otn(config);
        n.run();
        let lines = n.get_assembly_charlines();
        let mut output_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
        output_dir.push("../../target/tmp/fluid_field.txt");
        lines.write_to_file(output_dir, false).unwrap();
    }

    /// 使用非零膨胀角、非零喉部半径测试完整喷管（含膨胀段）
    #[test]
    fn test_with_expansion() {
        let config = NozzleConfig {
            control: Control::default(),
            material: Material::air_piecewise_polynomial(),
            inlet: Inlet {
                ma: 1.2,
                ..Inlet::default()
            },
            geometry: Geometry {
                height_i: 0.1, // 缩小进口高度使喷管更紧凑
                length: 2.0,   // 缩短长度
                ..Geometry::default()
            },
            throat: Throat {
                radius_throat: 0.05,           // 喉部过渡圆弧半径
                theta_a: 5.0_f64.to_radians(), // 初始膨胀角 5°
            },
            outlet: Outlet::default(),
            io: IO::default(),
        };
        let mut n = ConstraintNozzle::new_otn(config);
        n.run();
        let lines = n.get_assembly_charlines();

        // 验证至少产生了特征线数据
        assert!(!lines.is_empty(), "应产生流场数据");

        // 验证所有点坐标不越界
        for line in lines.iter() {
            for point in line.iter() {
                assert!(point.is_valid(), "存在无效点: {:?}", point);
                assert!(point.x >= 0.0, "x 不应为负: {}", point.x);
                assert!(point.y >= 0.0, "y 不应为负: {}", point.y);
            }
        }

        // 验证出口（最后一条线）的 x 接近 target length
        let last_line = lines.last().unwrap();
        let last_point = last_line.last().unwrap();
        assert!(last_point.x <= 2.0 + 1e-6, "出口 x 不应超过目标长度");

        // 写入文件以便可视化检查
        let mut output_dir = PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap());
        output_dir.push("../../target/tmp/fluid_field_expansion.txt");
        lines.write_to_file(output_dir, false).unwrap();
    }
}
