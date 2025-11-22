use crate::{
    moc::unitprocess::UnitProcess,
    nozzle::{NozzleConfig, Section},
};

pub struct ConstraintNozzle {
    config: NozzleConfig,
    unitprocess: Box<dyn UnitProcess>,
    sections: Vec<Box<dyn Section>>,
}

impl ConstraintNozzle {
    pub fn run(&mut self) {
        for section in self.sections.iter_mut() {
            section.run(&self.unitprocess, &self.config);
        }
    }
}
