#[derive(Debug, Clone, PartialEq)]
pub struct Polynomial {
    /// 常数项
    constant: f64,

    /// 正次幂系数数组，coefficients_positive[i] 对应 x^(i+1) 的系数
    coefficients_positive: Vec<f64>,

    /// 负次幂系数数组，coefficients_negative[i] 对应 x^-(i+1) 的系数
    coefficients_negative: Vec<f64>,
}

impl From<f64> for Polynomial {
    fn from(value: f64) -> Self {
        Polynomial::from_constant(value)
    }
}

impl Default for Polynomial {
    fn default() -> Self {
        Self {
            constant: Default::default(),
            coefficients_positive: Default::default(),
            coefficients_negative: Default::default(),
        }
    }
}

impl Polynomial {
    /// 从各部分创建多项式
    pub fn from_constant(constant: f64) -> Self {
        Self {
            constant,
            coefficients_positive: Default::default(),
            coefficients_negative: Default::default(),
        }
    }

    /// 从各部分创建多项式
    pub fn from_parts(constant: f64, positive: Vec<f64>, negative: Vec<f64>) -> Self {
        Self {
            constant,
            coefficients_positive: positive,
            coefficients_negative: negative,
        }
    }

    /// 获取指定次幂的系数
    pub fn coefficient(&self, power: i32) -> f64 {
        if power == 0 {
            self.constant
        } else if power > 0 {
            let idx = (power - 1) as usize;
            if idx < self.coefficients_positive.len() {
                self.coefficients_positive[idx]
            } else {
                0.0
            }
        } else {
            // power < 0
            let idx = (-power - 1) as usize;
            if idx < self.coefficients_negative.len() {
                self.coefficients_negative[idx]
            } else {
                0.0
            }
        }
    }

    /// 设置指定次幂的系数
    pub fn set_coefficient(&mut self, power: i32, coeff: f64) {
        if power == 0 {
            self.constant = coeff;
        } else if power > 0 {
            let idx = (power - 1) as usize;
            if idx >= self.coefficients_positive.len() {
                self.coefficients_positive.resize(idx + 1, 0.0);
            }
            self.coefficients_positive[idx] = coeff;
        } else {
            // power < 0
            let idx = (-power - 1) as usize;
            if idx >= self.coefficients_negative.len() {
                self.coefficients_negative.resize(idx + 1, 0.0);
            }
            self.coefficients_negative[idx] = coeff;
        }
    }

    /// 多项式求值（x不能为0）
    pub fn evaluate(&self, x: f64) -> f64 {
        if !self.coefficients_negative.is_empty() && x == 0.0 {
            return f64::NAN; // x=0时负次幂未定义
        }

        let mut result = self.constant;

        // 计算正次幂部分
        let mut x_power = x;
        for &coeff in &self.coefficients_positive {
            result += coeff * x_power;
            x_power *= x;
        }

        // 计算负次幂部分
        let mut x_power = 1.0 / x; // x^-1
        for &coeff in &self.coefficients_negative {
            result += coeff * x_power;
            x_power /= x;
        }

        result
    }

    /// 获取常数项
    pub fn constant(&self) -> f64 {
        self.constant
    }

    /// 获取正次幂系数数组的引用
    pub fn positive_coefficients(&self) -> &Vec<f64> {
        &self.coefficients_positive
    }

    /// 获取负次幂系数数组的引用
    pub fn negative_coefficients(&self) -> &Vec<f64> {
        &self.coefficients_negative
    }

    /// 判断是否为零多项式
    pub fn is_zero(&self) -> bool {
        self.constant.abs() < f64::EPSILON
            && self
                .coefficients_positive
                .iter()
                .all(|&c| c.abs() < f64::EPSILON)
            && self
                .coefficients_negative
                .iter()
                .all(|&c| c.abs() < f64::EPSILON)
    }

    // 判断是否为常多项式
    pub fn is_constant(&self) -> bool {
        self.coefficients_negative.is_empty() && self.coefficients_positive.is_empty()
    }
}

impl std::fmt::Display for Polynomial {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.is_zero() {
            return write!(f, "0");
        }

        let mut parts = Vec::new();

        // 添加常数项
        if self.constant.abs() >= f64::EPSILON {
            parts.push(format!("{}", self.constant));
        }

        // 添加正次幂项
        for (i, &coeff) in self.coefficients_positive.iter().enumerate() {
            if coeff.abs() >= f64::EPSILON {
                let power = i + 1;
                let abs_coeff = coeff.abs();
                let sign = if coeff > 0.0 { "+" } else { "-" };

                if (abs_coeff - 1.0).abs() < f64::EPSILON {
                    parts.push(format!("{} x^{}", sign, power));
                } else {
                    parts.push(format!("{} {}*x^{}", sign, abs_coeff, power));
                }
            }
        }

        // 添加负次幂项
        for (i, &coeff) in self.coefficients_negative.iter().enumerate() {
            if coeff.abs() >= f64::EPSILON {
                let power = -(i as i32 + 1);
                let abs_coeff = coeff.abs();
                let sign = if coeff > 0.0 { "+" } else { "-" };

                if (abs_coeff - 1.0).abs() < f64::EPSILON {
                    parts.push(format!("{} x^{}", sign, power));
                } else {
                    parts.push(format!("{} {}*x^{}", sign, abs_coeff, power));
                }
            }
        }

        write!(f, "{}", parts.join(" "))?;
        Ok(())
    }
}

impl std::ops::Add for Polynomial {
    type Output = Polynomial;

    fn add(self, other: Polynomial) -> Polynomial {
        let mut result = self.clone();
        result += other;
        result
    }
}

impl std::ops::AddAssign for Polynomial {
    fn add_assign(&mut self, other: Polynomial) {
        self.constant += other.constant;

        // 处理正次幂
        let max_len = self
            .coefficients_positive
            .len()
            .max(other.coefficients_positive.len());
        self.coefficients_positive.resize(max_len, 0.0);

        for (i, &coeff) in other.coefficients_positive.iter().enumerate() {
            self.coefficients_positive[i] += coeff;
        }

        // 处理负次幂
        let max_len = self
            .coefficients_negative
            .len()
            .max(other.coefficients_negative.len());
        self.coefficients_negative.resize(max_len, 0.0);

        for (i, &coeff) in other.coefficients_negative.iter().enumerate() {
            self.coefficients_negative[i] += coeff;
        }
    }
}

impl std::ops::Sub for Polynomial {
    type Output = Polynomial;

    fn sub(self, other: Polynomial) -> Polynomial {
        let mut result = self.clone();
        result -= other;
        result
    }
}

impl std::ops::SubAssign for Polynomial {
    fn sub_assign(&mut self, other: Polynomial) {
        self.constant -= other.constant;

        // 处理正次幂
        let max_len = self
            .coefficients_positive
            .len()
            .max(other.coefficients_positive.len());
        self.coefficients_positive.resize(max_len, 0.0);

        for (i, &coeff) in other.coefficients_positive.iter().enumerate() {
            self.coefficients_positive[i] -= coeff;
        }

        // 处理负次幂
        let max_len = self
            .coefficients_negative
            .len()
            .max(other.coefficients_negative.len());
        self.coefficients_negative.resize(max_len, 0.0);

        for (i, &coeff) in other.coefficients_negative.iter().enumerate() {
            self.coefficients_negative[i] -= coeff;
        }
    }
}
