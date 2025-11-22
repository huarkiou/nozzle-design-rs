use crate::{moc::unitprocess::UnitProcess, nozzle::NozzleConfig};

pub trait Section {
    fn run(&mut self, unitprocess: &Box<dyn UnitProcess>, config: &NozzleConfig);
}
