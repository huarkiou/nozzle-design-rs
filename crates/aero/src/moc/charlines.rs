use std::ops::{Deref, DerefMut};

use crate::moc::CharLine;

#[derive(Clone)]
pub struct CharLines {
    data: Vec<CharLine>,
}

impl Deref for CharLines {
    type Target = Vec<CharLine>;

    fn deref(&self) -> &Self::Target {
        &self.data
    }
}

impl DerefMut for CharLines {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.data
    }
}

impl CharLines {
    pub fn new() -> Self {
        Self { data: Vec::new() }
    }
}
