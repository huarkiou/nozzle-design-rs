use crate::{
    moc::{CharLine, CharLines, MocPoint, unitprocess::Context},
    nozzle::Section,
};

struct InitialSection {
    char_lines: CharLines,
    line_init: CharLine,
    calculated: bool,
}

impl InitialSection {
    pub fn new(line_init: CharLine) -> Self {
        Self {
            char_lines: CharLines::new(),
            line_init,
            calculated: false,
        }
    }

    fn cal_next_ivp_line(
        unitprocess: &dyn crate::moc::unitprocess::UnitProcess,
        line_prev: &CharLine,
        p_init: &MocPoint,
    ) -> CharLine {
        let mut line_cur = Vec::with_capacity(line_prev.len() + 2);
        let context: Context;

        let mut j2 = line_prev.len() as i64 - 1;
        let mut n = 0i32; // Counter for inserted points

        for j in 0..=j2 {
            if j == 0 {
                j2 = j2 + 2; // Update final index of current characteristic line
                points[3] = *p_init; // Known initial point
            } else if j == j2 {
                // Boundary condition for symmetry axis
                points[0] = line_cur[j as usize - 1];
                points[2] = *line_prev.last().unwrap();
                unitprocess.symmetry_axis_point(&mut points);
            } else {
                points[0] = line_cur[j as usize - 1];
                let prev_idx = (j as i32 - 1 - n) as usize;
                if prev_idx < line_prev.len() {
                    points[1] = line_prev[prev_idx];

                    let offset = unitprocess.interior_point(&mut points);
                    if offset == -999 {
                        j2 -= 1;
                        continue;
                    }
                }
            }

            // Ensure line_cur has enough capacity
            if j >= line_cur.len() as i64 {
                line_cur.resize((j + 1) as usize, MocPoint::default());
            }
            line_cur[j as usize] = points[3];
        }

        // Final resize to match expected length
        if line_cur.len() != (j2 + 1) as usize {
            line_cur.truncate((j2 + 1) as usize);
        }

        line_cur
    }
}

impl Section for InitialSection {
    fn run(
        &mut self,
        unitprocess: &dyn crate::moc::unitprocess::UnitProcess,
        config: &super::NozzleConfig,
    ) {
        assert!(self.char_lines.first().is_some_and(|line| !line.is_empty()));
        if self.calculated {
            return;
        }
        let mut line_cur = CharLine::new();
        line_cur.push(self.line_init.last().unwrap().clone());

        let init_iter = self.line_init.iter().rev().skip(1);

        for point in init_iter {
            line_cur = Self::cal_next_ivp_line(unitprocess, &line_cur, point);
            if !line_cur.is_empty() {
                self.char_lines.push(line_cur);
            }
        }
        self.calculated = true;
    }

    fn get_charlines(&self) -> &crate::moc::CharLines {
        if !self.calculated {
            panic!("InitialSection has not run, cannot get result!");
        }
        &self.char_lines
    }
}
