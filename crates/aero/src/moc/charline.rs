use std::ops::{Deref, DerefMut};

use crate::moc::MocPoint;

#[derive(Clone)]
pub struct CharLine {
    data: Vec<MocPoint>,
}

impl Deref for CharLine {
    type Target = Vec<MocPoint>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for CharLine {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}
