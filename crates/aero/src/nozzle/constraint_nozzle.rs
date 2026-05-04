use crate::{
    moc::{unitprocess::UnitProcess, CharLines},
    nozzle::{initial_line::InitialLine, InitialSection, NozzleConfig, Section},
};

pub struct ConstraintNozzle {
    config: NozzleConfig,
    unitprocess: Box<dyn UnitProcess>,
    sections: Vec<Box<dyn Section>>,
}

impl ConstraintNozzle {
    /// 构建最大推力喷管OTN
    pub fn new_otn(config: NozzleConfig) -> Self {
        Self {
            unitprocess: config.to_unitprocess(),
            config,
            sections: vec![
                Box::new(InitialLine::new()),
                Box::new(InitialSection::new()),
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
}
