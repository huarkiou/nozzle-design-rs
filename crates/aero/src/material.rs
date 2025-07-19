use std::sync::Arc;

use math::quadrature;

// 热容类型枚举
#[derive(Clone)]
pub enum Cp {
    Constant(f64),
    Variable(Arc<dyn Fn(f64) -> f64>),
}

impl PartialEq for Cp {
    fn eq(&self, other: &Self) -> bool {
        match (self, other) {
            (Self::Constant(l0), Self::Constant(r0)) => l0 == r0,
            (Self::Variable(l0), Self::Variable(r0)) => Arc::ptr_eq(l0, r0),
            _ => false,
        }
    }
}

impl Cp {
    pub fn new(cp: impl Fn(f64) -> f64 + 'static) -> Self {
        Self::Variable(Arc::new(cp))
    }

    pub fn eval(&self, temperature: f64) -> f64 {
        match self {
            Cp::Constant(v) => *v,
            Cp::Variable(f) => f(temperature),
        }
    }
}

impl Into<Cp> for f64 {
    fn into(self) -> Cp {
        Cp::Constant(self)
    }
}

// MaterialProperty 结构体定义
#[derive(Clone)]
pub struct Material {
    /// 摩尔质量 (kg/kmol)
    molecular_weight: f64,

    /// 定压比热容 Cp(J/(kg·K))，输入温度 T(K)
    cp: Cp,
}

impl PartialEq for Material {
    fn eq(&self, other: &Self) -> bool {
        self.molecular_weight == other.molecular_weight && self.cp == other.cp
    }
}

impl Material {
    pub fn borrow_cp(&self) -> &Cp {
        &self.cp
    }

    pub fn borrom_mw(&self) -> &f64 {
        &self.molecular_weight
    }
}

impl Material {
    // 常量定义为 associated constants
    pub const UNIVERSAL_GAS_CONSTANT: f64 = 8.31446261815324; // J/(mol·K) R = N_A * k
    pub const AVOGADRO_CONSTANT: f64 = 6.02214076e23; // mol⁻¹
    pub const BOLTZMANN_CONSTANT: f64 = 1.380649e-23; // J/K
    pub const SPEED_OF_LIGHT: f64 = 299792458.0; // m/s

    /// 用摩尔质量kg/kmol和Cp(T)函数构造
    /// # 参数
    /// - `molecular_weight`: 摩尔质量(单位kg/kmol) 注意此处不是国际单位制
    /// - `cp`: 定压比热容(单位J/(kg·K))关于温度T(单位K)的函数
    pub fn new(molecular_weight: f64, cp: impl Fn(f64) -> f64 + 'static) -> Self {
        Self {
            molecular_weight,
            cp: (Cp::new(cp)),
        }
    }

    /// 用摩尔质量kg/kmol和Cp(T)函数构造
    /// # 参数
    /// - `molecular_weight`: 摩尔质量(单位kg/kmol) 注意此处不是国际单位制
    /// - `cp`: 定压比热容(单位J/(kg·K)) 常数
    pub fn from_mw_cp(molecular_weight: f64, cp: f64) -> Self {
        Self {
            molecular_weight,
            cp: (Cp::Constant(cp)),
        }
    }

    /// 用气体常数和比热比构造
    /// # 参数
    /// - `r_gas`: 气体常数
    /// - `gamma`: 比热比
    pub fn from_rgas_gamma(r_gas: f64, gamma: f64) -> Self {
        Self {
            molecular_weight: Self::UNIVERSAL_GAS_CONSTANT / r_gas * 1e3, // kg/kmol
            cp: Cp::Constant(gamma / (gamma - 1.0) * r_gas),
        }
    }

    // 获取 Cp 值（在给定温度下）
    pub fn cp(&self, temperature: f64) -> f64 {
        self.cp.eval(temperature)
    }

    // 比热比 gamma，默认使用 Cp 计算 Cv = Cp - Rgas
    pub fn gamma(&self, temperature: f64) -> f64 {
        let cp = self.cp(temperature);
        let cv = cp - self.rgas();
        cp / cv
    }

    // 气体常数 Rgas = R / Mw （单位：J/(kg·K)）
    pub fn rgas(&self) -> f64 {
        Self::UNIVERSAL_GAS_CONSTANT / self.molecular_weight * 1e3
    }

    // 焓值 h(T)，参考点默认为 T_ref = 0.0
    pub fn enthalpy(&self, temperature: f64, t_ref: f64) -> f64 {
        if temperature == t_ref {
            return 0.0;
        }
        quadrature::gauss_legendre::integrate::<20, _>(|t| self.cp(t), t_ref, temperature)
    }
}

// 静态方法模拟：工厂函数
impl Material {
    // 空气（常数 Cp）
    pub fn air_constant() -> Self {
        Self::from_mw_cp(28.968, 1004.675)
    }

    // 空气（分段多项式 Cp）
    pub fn air_piecewise_polynomial() -> Self {
        Self::new(28.968, move |t| {
            let poly1 = |x: f64| {
                1161.482
                    + (-2.368819
                        + (0.01485511
                            + (-5.034909e-5
                                + (9.928570e-8
                                    + (-1.111097e-10 + (6.540196e-14 - 1.573588e-17 * x) * x)
                                        * x)
                                    * x)
                                * x)
                            * x)
                        * x
            };
            let poly2 = |x: f64| {
                -7069.814
                    + (33.70605
                        + (-0.05812760
                            + (5.421615e-5
                                + (-2.936679e-8
                                    + (9.237533e-12 + (-1.565553e-15 + 1.112335e-19 * x) * x)
                                        * x)
                                    * x)
                                * x)
                            * x)
                        * x
            };
            if 100. <= t && t <= 1000. {
                poly1(t)
            } else if 1000. < t && t <= 3000. {
                poly2(t)
            } else if t < 100. {
                poly1(100.)
            } else if t > 3000. {
                poly2(3000.)
            } else {
                panic!("MaterialProperty::air_piecewise_polynomial: Temperature is lower than 0!");
            }
        })
    }

    // 空气(NASA 9系数模型）
    pub fn air_nasa9piecewise_polynomial() -> Self {
        Self::new(28.968, move |t| {
            let poly1 = |x: f64| {
                (2898903. / x - 56496.26) / x
                    + 1437.799
                    + (-1.653609 + (0.003062254 + (-2.279138e-06 + 6.272365e-10 * x) * x) * x) * x
            };
            let poly2 = |x: f64| {
                (6.932494e+07 / x - 361053.2) / x
                    + 1476.665
                    + (-0.06138349 + (2.027963e-05 + (-3.075525e-09 + 1.888054e-13 * x) * x) * x)
                        * x
            };
            if 200. <= t && t <= 1000. {
                poly1(t)
            } else if 1000. < t && t <= 6000. {
                poly2(t)
            } else if t < 200. {
                poly1(200.)
            } else if t > 6000. {
                poly2(6000.)
            } else {
                panic!("MaterialProperty::air_piecewise_polynomial: Temperature is lower than 0!");
            }
        })
    }
}

#[cfg(test)]
mod test {
    use super::*;
    #[test]
    fn test_material() {
        let _1 = Material::new(25.0, |_| 1006.1);
        assert_eq!(_1.cp(100.), 1006.1);
        let _2 = Material::air_constant();
        let _3 = Material::air_piecewise_polynomial();
        let _3 = Material::air_nasa9piecewise_polynomial();
    }
}
