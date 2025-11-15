use std::ops::{Add, Div, Mul, Sub};

#[derive(Copy, Clone, PartialEq)]
pub struct Coord2d {
    pub x: f64,
    pub y: f64,
}

impl Coord2d {
    pub fn new(x: f64, y: f64) -> Self {
        Self { x, y }
    }

    /// 计算两点的距离的平方
    ///
    /// # 输入参数
    /// - `other:Coord2d`: 另一点
    /// # 返回值
    /// - `f64`: 当前坐标到other的距离的平方
    #[inline]
    pub fn distance_squared_to(&self, other: &Self) -> f64 {
        (self.x - other.x).powi(2) + (self.y - other.y).powi(2)
    }

    /// 计算两点的距离
    ///
    /// # 输入参数
    /// - `other:Coord2d`: 另一点
    /// # 返回值
    /// - `f64`: 当前坐标到other的距离
    #[inline]
    pub fn distance_to(&self, other: &Self) -> f64 {
        self.distance_squared_to(other).sqrt()
    }

    /// 计算当前坐标点相对点origin的方位角
    ///
    /// # 输入参数
    /// - `origin:Coord2d`: 参考点origin，以其为原点
    /// # 返回值
    /// - `f64`: 当前坐标相对origin的方位角、极角
    #[inline]
    pub fn bearing_to(&self, origin: &Self) -> f64 {
        (self.y - origin.y).atan2(self.x - origin.x)
    }

    /// 绕指定点逆时针旋转
    ///
    /// # 参数
    /// - `origin`: 旋转中心点
    /// - `alpha`: 旋转角度（弧度），正数表示逆时针旋转
    ///
    /// # 示例
    /// ```
    /// # use math::geometry::Coord2d;
    /// let mut point = Coord2d::new(1.0, 0.0);
    /// let origin = Coord2d::new(0.0, 0.0);
    /// point.rotate_around(&origin, std::f64::consts::PI / 2.0); // 逆时针旋转90度
    /// // point 现在大约是 (-0.0, 1.0)
    /// ```
    #[inline]
    pub fn rotate_around(&mut self, origin: &Self, alpha: f64) {
        let theta = self.bearing_to(origin) + alpha;
        let length = self.distance_to(origin);
        self.x = origin.x + length * theta.cos();
        self.y = origin.y + length * theta.sin();
    }
}

impl Default for Coord2d {
    fn default() -> Self {
        Self { x: 0., y: 0. }
    }
}

impl Add for Coord2d {
    type Output = Coord2d;

    fn add(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x + rhs.x,
            y: self.y + rhs.y,
        }
    }
}

impl Sub for Coord2d {
    type Output = Coord2d;

    fn sub(self, rhs: Self) -> Self::Output {
        Self::Output {
            x: self.x - rhs.x,
            y: self.y - rhs.y,
        }
    }
}

impl Mul<f64> for Coord2d {
    type Output = Coord2d;

    fn mul(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x * rhs,
            y: self.y * rhs,
        }
    }
}

impl Div<f64> for Coord2d {
    type Output = Coord2d;

    fn div(self, rhs: f64) -> Self::Output {
        Self::Output {
            x: self.x / rhs,
            y: self.y / rhs,
        }
    }
}
