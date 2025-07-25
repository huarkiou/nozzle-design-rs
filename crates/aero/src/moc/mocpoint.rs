use std::{
    fmt::Display,
    ops::{Add, Div, Mul, Sub},
};

use math::Tolerance;

use crate::{Material, isentropic};

/// 描述超声速流场中的任一点的气流参数
#[derive(Clone)]
pub struct MocPoint {
    /// x坐标 单位：m
    pub x: f64,
    /// y坐标 单位：m
    pub y: f64,
    /// x方向速度 单位：m/s
    pub u: f64,
    /// y方向速度 单位：m/s
    pub v: f64,
    /// 静压 单位：Pa
    pub p: f64,
    /// 静温 单位：K
    pub t: f64,
    /// 静密度 单位：kg/m^3
    pub rho: f64,
    /// 物性参数
    pub mat: Material,
}

// 构造方法
impl MocPoint {
    pub fn new(x: f64, y: f64, u: f64, v: f64, p: f64, t: f64, rho: f64, mat: Material) -> Self {
        Self {
            x,
            y,
            u,
            v,
            p,
            t,
            rho,
            mat,
        }
    }

    pub fn from_compatible(
        x: f64,
        y: f64,
        u: f64,
        v: f64,
        p: f64,
        t_total: f64,
        rho: f64,
        mat: Material,
    ) -> Self {
        let t = isentropic::cal_static_temperature(
            mat.borrow_cp(),
            t_total,
            (u.powi(2) + v.powi(2)).sqrt(),
        );
        Self {
            x,
            y,
            u,
            v,
            p,
            t,
            rho,
            mat,
        }
    }
}

// 设置方法
impl MocPoint {
    pub fn set_temperature_pressure_density(&mut self, t_total: f64, p_total: f64, rho_total: f64) {
        let t_static =
            isentropic::cal_static_temperature(self.mat.borrow_cp(), t_total, self.velocity());
        self.t = t_static;

        let ma = self.mach_number();
        let sub_exp1 = 1.0 + (self.gamma(self.t) - 1.0) / 2.0 * ma.powi(2);
        let sub_exp2 = 1. / (self.gamma(self.t) - 1.);

        self.p = p_total / sub_exp1.powf(self.gamma(self.t) * sub_exp2);
        self.rho = rho_total / sub_exp1.powf(sub_exp2);
    }
}

// 计算和获取气流参数
impl MocPoint {
    /// 计算比热比
    pub fn gamma(&self, t: f64) -> f64 {
        self.mat.gamma(t)
    }

    /// 气体常数
    pub fn rg(&self) -> f64 {
        self.mat.rgas()
    }

    /// 计算定压比热容Cp 单位：J/(kg·K)
    pub fn cp(&self, t: f64) -> f64 {
        self.mat.cp(t)
    }

    /// 该点到另一点的距离的平方 单位：m^2
    pub fn distance_squared_to(&self, other: Self) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2)
    }

    /// 该点到另一点的距离 单位：m
    pub fn distance_to(&self, other: Self) -> f64 {
        self.distance_squared_to(other).sqrt()
    }

    /// 气流速度的平方 单位：(m/s)^2
    pub fn velocity_squared(&self) -> f64 {
        self.u * self.u + self.v * self.v
    }

    /// 气流速度 单位：m/s
    pub fn velocity(&self) -> f64 {
        self.velocity_squared().sqrt()
    }

    /// 气流方向 单位：rad
    pub fn flow_direction(&self) -> f64 {
        (self.v / self.u).atan()
    }

    /// 计算声速平方c^2 单位：(m/s)^2
    pub fn sound_speed_squared(&self) -> f64 {
        self.gamma(self.t) * self.rg() * self.t
    }

    /// 计算声速c 单位：m/s
    pub fn sound_speed(&self) -> f64 {
        self.sound_speed_squared().sqrt()
    }

    /// 计算马赫数Ma 单位：1
    pub fn mach_number(&self) -> f64 {
        let ma = self.velocity_squared() / self.sound_speed_squared();
        if ma < 1. { 1. } else { ma.sqrt() }
    }

    /// 同时计算声速和马赫数
    pub fn sound_speed_and_mach_number(&self) -> (f64, f64) {
        let c2 = self.sound_speed_squared();
        let v2 = self.velocity_squared();
        let ma2 = v2 / c2;
        let ma = if ma2 < 1.0 { 1.0 } else { ma2.sqrt() };
        (c2.sqrt(), ma)
    }

    /// 计算总温T* 单位：K
    pub fn total_temperature(&self) -> f64 {
        return self.t + self.velocity_squared() / (2. * self.cp(self.t));
    }

    /// 计算总压p* 单位：Pa
    pub fn total_pressure(&self) -> f64 {
        let ma = self.mach_number();
        let sub_exp1 = self.gamma(self.t) - 1.;
        let sub_exp2 = 1. + sub_exp1 / 2. * ma.powi(2);
        return self.p * sub_exp2.powf(self.gamma(self.t) / sub_exp1);
    }

    /// 计算总密度ρ* 单位：kg/m^3
    pub fn total_density(&self) -> f64 {
        let ma = self.mach_number();
        let sub_exp1 = self.gamma(self.t) - 1.;
        let sub_exp2 = 1. + sub_exp1 / 2. * ma.powi(2);
        return self.rho * sub_exp2.powf(1. / sub_exp1);
    }

    /// 计算总参数 (总温、总压、总密度)
    pub fn total_temperature_pressure_density(&self) -> (f64, f64, f64) {
        let ma = self.mach_number();
        let sub_exp1 = self.gamma(self.t) - 1.;
        let sub_exp2 = 1. + sub_exp1 / 2. * ma.powi(2);
        let tt = self.t + self.velocity_squared() / (2. * self.cp(self.t));
        let tp = self.p * sub_exp2.powf(self.gamma(self.t) / sub_exp1);
        let td = self.rho * sub_exp2.powf(1. / sub_exp1);
        (tt, tp, td)
    }
}

// 坐标判断
impl MocPoint {
    /// 判断当前点是否在 point1 和 point2 构成的轴对齐矩形包围盒内（考虑浮点误差）
    pub fn is_between(&self, point1: &Self, point2: &Self, eps: f64) -> bool {
        // 判断是否在包围盒内（包括边界）
        let in_x =
            (self.x >= point1.x.min(point2.x) - eps) && (self.x <= point1.x.max(point2.x) + eps);
        let in_y =
            (self.y >= point1.y.min(point2.y) - eps) && (self.y <= point1.y.max(point2.y) + eps);
        in_x && in_y
    }

    /// 判断当前点是否与point1和point2共线（考虑浮点误差）
    pub fn is_collinear_with(&self, point1: &Self, point2: &Self, eps: f64) -> bool {
        let cross_product = (point2.y - point1.y) * (self.x - point1.x)
            - (point2.x - point1.x) * (self.y - point1.y);
        cross_product.abs() < eps
    }

    /// 判断当前点是否位于线段 point1 -- point2 上（考虑浮点误差）
    pub fn is_on_segment(&self, point1: &Self, point2: &Self, eps: f64) -> bool {
        // 1. 叉积判断是否共线
        if !self.is_collinear_with(point1, point2, eps) {
            return false;
        }
        // 2. 判断是否在线段范围内
        self.is_between(point1, point2, eps)
    }
}

// 收敛判断
impl MocPoint {
    // 注意：is_*_converged_with()系列函数依赖 IEEE-754 浮点运算语义。
    // 若启用了 fast-math 编译选项（如 RUSTFLAGS="-C fast-math"），
    // 则可能导致 NaN 比较行为异常，从而影响收敛判断。

    /// 坐标(x,y)是否收敛
    pub fn is_position_converged_with(&self, other: &Self, tol: Tolerance) -> bool {
        tol.approx_eq(self.x, other.x) && tol.approx_eq(self.y, other.y)
    }

    /// 速度(u,v)是否收敛
    pub fn is_velocity_converged_with(&self, other: &Self, tol: Tolerance) -> bool {
        tol.approx_eq(self.u, other.u) && tol.approx_eq(self.v, self.v)
    }

    /// 所有参数是否收敛
    pub fn is_converged_with(&self, other: &Self, tol: Tolerance) -> bool {
        // 注意：Rust 中 NaN != NaN，所以 (a - b).abs() < eps 会正确返回 false 如果任意是 NaN
        self.is_position_converged_with(other, tol)
            && self.is_velocity_converged_with(other, tol)
            && tol.approx_eq(self.p, other.p)
            && tol.approx_eq(self.t, other.t)
            && tol.approx_eq(self.rho, other.rho)
    }

    /// 坐标(x,y)是否有效
    pub fn is_position_valid(&self) -> bool {
        self.x.is_finite() && self.y.is_finite()
    }

    /// 速度(u,v)是否有效
    pub fn is_velocity_valid(&self) -> bool {
        self.u.is_finite() && self.v.is_finite()
    }

    /// 材料参数比热比和气体常数是否有效
    pub fn is_material_valid(&self) -> bool {
        self.gamma(self.t) > 1. && self.rg() > Material::UNIVERSAL_GAS_CONSTANT
    }

    /// 所有参数是否有效
    pub fn is_valid(&self) -> bool {
        self.is_position_valid()
            && self.is_velocity_valid() & self.is_material_valid()
            && self.p > 0.
            && self.t > 0.
            && self.rho > 0.
    }
}

impl Default for MocPoint {
    fn default() -> Self {
        Self {
            x: f64::NAN,
            y: f64::NAN,
            u: f64::NAN,
            v: f64::NAN,
            p: f64::NAN,
            t: f64::NAN,
            rho: f64::NAN,
            mat: Material::air_constant(),
        }
    }
}

impl Display for MocPoint {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{0}, {1}, {2}, {3}, {4}, {5}, {6}, {7}, {8}, {9}",
            self.x,
            self.y,
            self.velocity(),
            self.flow_direction(),
            self.p,
            self.rho,
            self.t,
            self.rg(),
            self.gamma(self.t),
            self.mach_number()
        )
    }
}

impl PartialEq for MocPoint {
    fn eq(&self, other: &Self) -> bool {
        self.x == other.x
            && self.y == other.y
            && self.u == other.u
            && self.v == other.v
            && self.p == other.p
            && self.t == other.t
            && self.rho == other.rho
            && self.mat == other.mat
    }
}

impl Add for &MocPoint {
    type Output = MocPoint;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            u: self.u + rhs.u,
            v: self.v + rhs.v,
            p: self.p + rhs.p,
            t: self.t + rhs.t,
            rho: self.rho + rhs.rho,
            mat: self.mat.clone(),
        }
    }
}

impl Add for MocPoint {
    type Output = MocPoint;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
            u: self.u + rhs.u,
            v: self.v + rhs.v,
            p: self.p + rhs.p,
            t: self.t + rhs.t,
            rho: self.rho + rhs.rho,
            mat: self.mat.clone(),
        }
    }
}

impl Sub for &MocPoint {
    type Output = MocPoint;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            u: self.u - rhs.u,
            v: self.v - rhs.v,
            p: self.p - rhs.p,
            t: self.t - rhs.t,
            rho: self.rho - rhs.rho,
            mat: self.mat.clone(),
        }
    }
}

impl Sub for MocPoint {
    type Output = MocPoint;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
            u: self.u - rhs.u,
            v: self.v - rhs.v,
            p: self.p - rhs.p,
            t: self.t - rhs.t,
            rho: self.rho - rhs.rho,
            mat: self.mat.clone(),
        }
    }
}

impl Mul<f64> for &MocPoint {
    type Output = MocPoint;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x * rhs,
            y: self.y * rhs,
            u: self.u * rhs,
            v: self.v * rhs,
            p: self.p * rhs,
            t: self.t * rhs,
            rho: self.rho * rhs,
            mat: self.mat.clone(),
        }
    }
}

impl Mul<f64> for MocPoint {
    type Output = MocPoint;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x * rhs,
            y: self.y * rhs,
            u: self.u * rhs,
            v: self.v * rhs,
            p: self.p * rhs,
            t: self.t * rhs,
            rho: self.rho * rhs,
            mat: self.mat.clone(),
        }
    }
}

impl Div<f64> for &MocPoint {
    type Output = MocPoint;

    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x / rhs,
            y: self.y / rhs,
            u: self.u / rhs,
            v: self.v / rhs,
            p: self.p / rhs,
            t: self.t / rhs,
            rho: self.rho / rhs,
            mat: self.mat.clone(),
        }
    }
}

impl Div<f64> for MocPoint {
    type Output = MocPoint;

    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x / rhs,
            y: self.y / rhs,
            u: self.u / rhs,
            v: self.v / rhs,
            p: self.p / rhs,
            t: self.t / rhs,
            rho: self.rho / rhs,
            mat: self.mat.clone(),
        }
    }
}
