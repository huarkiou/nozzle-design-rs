use crate::{
    moc::{CharLine, CharLines, MocPoint},
    nozzle::{NozzleConfig, Section},
};

pub struct InitialLine {
    char_lines: CharLines,
    calculated: bool,
}

impl InitialLine {
    pub fn new() -> Self {
        Self {
            char_lines: CharLines::new(),
            calculated: false,
        }
    }
}

impl Section for InitialLine {
    fn run(
        &mut self,
        _unitprocess: &dyn crate::moc::unitprocess::UnitProcess,
        config: &super::NozzleConfig,
    ) {
        // 根据入口马赫数选择计算方法
        let initial_line = if config.inlet.ma < 1.0 {
            panic!("the Mach number of inlet flow cannot be lower than 1.0");
        } else if config.inlet.ma < 1.01 {
            cal_initial_value_line_sauers(&config)
        } else {
            cal_initial_value_line_direct(&config)
        };

        self.char_lines = initial_line;
        self.calculated = true;
    }

    fn get_charlines(&self) -> &crate::moc::CharLines {
        &self.char_lines
    }
}

fn cal_initial_value_line_direct(config: &NozzleConfig) -> CharLines {
    let n_inlet = config.control.n_inlet as usize;
    let y_t = config.geometry.height_i;
    let p_total = config.inlet.p_total;
    let t_total = config.inlet.temperature_total;
    let ma_i = config.inlet.ma;
    let theta_i = config.inlet.theta;
    let gamma = config.material.gamma(config.inlet.temperature_total);
    let rg = config.material.rgas();
    let dy = if n_inlet > 1 {
        y_t / (n_inlet as f64 - 1.0)
    } else {
        0.0
    };

    // 入口静温
    let t = t_total / (1.0 + (gamma - 1.0) / 2.0 * ma_i * ma_i);
    // γ/(γ-1)
    let tmp = gamma / (gamma - 1.0);
    // 气流定压比热容
    let cp = tmp * rg;
    // 入口气流速度
    let v = ((t_total - t) * 2.0 * cp).sqrt();
    let p = p_total * (t / t_total).powf(tmp);
    let rho = p / rg / t;

    let mut line_iv = CharLine::with_capacity(n_inlet);
    for i in 0..n_inlet {
        let x = 0.0;
        let y = if i > 0 { line_iv[i - 1].y + dy } else { 0.0 };

        // 创建 MocPoint
        let point = MocPoint {
            x,
            y,
            u: v * theta_i.cos(),
            v: v * theta_i.sin(),
            p,
            rho,
            t,
            mat: config.material.clone(),
        };

        if point.mach_number() < 1.0 {
            panic!("the Mach number of inlet flow cannot be lower than 1.0")
        }

        line_iv.push(point);
    }

    if !line_iv.is_empty() {
        let p = line_iv.last_mut().unwrap();
        p.x = 0.0;
        p.y = y_t;
    }

    line_iv.reverse();

    let mut char_lines = CharLines::new();
    char_lines.push(line_iv);
    char_lines
}

fn cal_initial_value_line_sauers(config: &NozzleConfig) -> CharLines {
    let n_inlet = config.control.n_inlet as usize;
    let y_t = config.geometry.height_i; // 假设使用入口高度作为喉部高度
    let axisymmetric = config.control.axisymmetric;
    let rtu = 2.01 * y_t;

    let mut line_iv = CharLine::with_capacity(n_inlet);

    let gamma = config.material.gamma(config.inlet.temperature_total);
    let rg = config.material.rgas();
    let p_total = config.inlet.p_total;
    let t_total = config.inlet.temperature_total;

    let dy = if n_inlet > 1 {
        y_t / (n_inlet as f64 - 1.0)
    } else {
        0.0
    };

    // Sauers法 v=0 -> x = -(γ+1)·α/[2·(3+δ)]·y^2 = C1·y^2
    if (rtu / y_t) < 2.0 {
        panic!("In Sauers method, RtU/y_t < 2.0 possibily results in low accuracy!");
    }

    let axisymmetric_val = if axisymmetric { 1.0 } else { 0.0 };

    // 计算参考参数
    let alpha = ((1.0 + axisymmetric_val) / ((gamma + 1.0) * rtu * y_t)).sqrt();
    let epsilon = -(gamma + 1.0) * alpha * y_t * y_t / (2.0 * (3.0 + axisymmetric_val));
    // x = C1·y ^ 2;
    let c1 = -(gamma + 1.0) * alpha / (2.0 * (3.0 + axisymmetric_val));
    // du / dy = C2·y - > u = 0.5·C2·y ^ 2
    let c2 = (gamma + 1.0) * alpha * alpha / (2.0 * (1.0 + axisymmetric_val));

    // 临界声速
    let c_crit = (2.0 * gamma * rg * t_total / (gamma + 1.0)).sqrt();

    for i in 0..n_inlet {
        let mut x;
        let y;
        if i > 0 {
            y = line_iv[i - 1].y + dy;
            x = c1 * y * y;
        } else {
            y = 0.0;
            x = 0.0;
        }

        let v = c_crit * (1.0 + alpha * x + c2 * y * y);
        x = x - epsilon;

        let t = t_total / (1.0 + (gamma - 1.0) / 2.0 * (v * v) / (gamma * rg));
        let p = p_total * (t / t_total).powf(gamma / (gamma - 1.0));
        let rho = p / rg / t;

        // 创建 MocPoint
        let point = MocPoint {
            x,
            y,
            u: v,   // Sauers法中theta=0，所以u=v
            v: 0.0, // Sauers法中theta=0，所以v=0
            p,
            rho,
            t,
            mat: config.material.clone(),
        };

        if point.mach_number() < 1.0 {
            panic!("the Mach number of inlet flow cannot be lower than 1.0")
        }

        line_iv.push(point);
    }

    if !line_iv.is_empty() {
        let p = line_iv.last_mut().unwrap();
        p.x = 0.0;
        p.y = y_t;
    }

    line_iv.reverse();

    let mut char_lines = CharLines::new();
    char_lines.push(line_iv);
    char_lines
}
