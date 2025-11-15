use std::fmt;

use crate::polynomial::Polynomial;

#[derive(Debug, Clone, PartialEq)]
pub struct PiecewisePolynomial {
    /// 分段点，升序排列
    breakpoints: Vec<f64>,

    /// 每个区间上的多项式
    /// 区间 [breakpoints[i], breakpoints[i+1]) 对应 polynomials[i]
    polynomials: Vec<Polynomial>,
}

impl From<f64> for PiecewisePolynomial {
    fn from(value: f64) -> Self {
        PiecewisePolynomial::from_constant(value)
    }
}

impl PiecewisePolynomial {
    /// 创建分段多项式
    /// breakpoints: 分段点，需要升序排列
    /// polynomials: 每个区间上的多项式，数量应比分段点数量多1
    pub fn new(breakpoints: Vec<f64>, polynomials: Vec<Polynomial>) -> Result<Self, &'static str> {
        if breakpoints.is_empty() && polynomials.len() == 1 {
            // 单个多项式的情况
            return Ok(Self {
                breakpoints: Vec::new(),
                polynomials,
            });
        }

        if breakpoints.len() + 1 != polynomials.len() {
            return Err("分段点数量必须比多项式数量少1");
        }

        // 检查分段点是否升序
        for i in 1..breakpoints.len() {
            if breakpoints[i] <= breakpoints[i - 1] {
                return Err("分段点必须严格升序");
            }
        }

        Ok(Self {
            breakpoints,
            polynomials,
        })
    }

    pub fn from_constant(c0: f64) -> Self {
        Self {
            breakpoints: Default::default(),
            polynomials: vec![c0.into()],
        }
    }

    /// 计算给定温度下的值
    pub fn evaluate(&self, temperature: f64) -> f64 {
        if self.breakpoints.is_empty() {
            // 只有一个多项式的情况
            return self.polynomials[0].evaluate(temperature);
        }

        // 找到对应的区间
        let mut idx = 0;
        for (i, &bp) in self.breakpoints.iter().enumerate() {
            if temperature >= bp {
                idx = i + 1;
            } else {
                break;
            }
        }

        // idx 对应的多项式
        self.polynomials[idx].evaluate(temperature)
    }

    /// 获取分段点
    pub fn breakpoints(&self) -> &[f64] {
        &self.breakpoints
    }

    /// 获取多项式段
    pub fn polynomials(&self) -> &[Polynomial] {
        &self.polynomials
    }

    /// 判断是否为常多项式
    pub fn is_constant(&self) -> bool {
        for window in self.polynomials.windows(2) {
            if !window[0].is_constant() || (window[0].constant() != window[1].constant()) {
                return false;
            }
        }
        self.polynomials.last().unwrap().is_constant()
    }
}

impl fmt::Display for PiecewisePolynomial {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.breakpoints.is_empty() {
            // 只有一个多项式的情况
            if let Some(poly) = self.polynomials.first() {
                return write!(f, "{}", poly);
            } else {
                return write!(f, "0");
            }
        }

        // 分段多项式的情况
        write!(f, "Piecewise[")?;

        for (i, (bp, poly)) in self
            .breakpoints
            .iter()
            .zip(self.polynomials.iter())
            .enumerate()
        {
            if i > 0 {
                write!(f, ", ")?;
            }
            write!(f, "{} for T < {}", poly, bp)?;
        }

        // 最后一个区间
        if let Some(last_poly) = self.polynomials.last() {
            write!(
                f,
                ", {} for T ≥ {}",
                last_poly,
                self.breakpoints.last().unwrap()
            )?;
        }

        write!(f, "]")?;
        Ok(())
    }
}
