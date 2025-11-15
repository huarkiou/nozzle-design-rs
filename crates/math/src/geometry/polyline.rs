use crate::geometry::Coord2d;

#[derive(Clone, PartialEq)]
pub struct Polyline2d {
    pub points: Vec<Coord2d>,
}

impl Polyline2d {
    pub fn new() -> Self {
        Self { points: Vec::new() }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            points: Vec::with_capacity(capacity),
        }
    }

    pub fn from_points(points: Vec<Coord2d>) -> Self {
        Self { points }
    }

    /// 添加一个点到折线末尾
    pub fn push(&mut self, point: Coord2d) {
        self.points.push(point);
    }

    /// 获取折线中点的数量
    pub fn len(&self) -> usize {
        self.points.len()
    }

    /// 折线是否为空
    pub fn is_empty(&self) -> bool {
        self.points.is_empty()
    }

    /// 获取指定索引的点的引用
    pub fn get(&self, index: usize) -> Option<&Coord2d> {
        self.points.get(index)
    }

    /// 获取指定索引的点的可变引用
    pub fn get_mut(&mut self, index: usize) -> Option<&mut Coord2d> {
        self.points.get_mut(index)
    }

    /// 计算折线的总长度
    pub fn length(&self) -> f64 {
        if self.points.len() < 2 {
            return 0.0;
        }

        self.points
            .windows(2)
            .map(|w| w[0].distance_to(&w[1]))
            .sum()
    }

    /// 绕指定点逆时针旋转整个折线
    pub fn rotate_around(&mut self, origin: &Coord2d, alpha: f64) {
        for point in &mut self.points {
            point.rotate_around(origin, alpha);
        }
    }

    /// 平移整个折线
    pub fn translate(&mut self, dx: f64, dy: f64) {
        for point in &mut self.points {
            point.x += dx;
            point.y += dy;
        }
    }

    /// 获取折线的边界框
    pub fn bounding_box(&self) -> Option<(Coord2d, Coord2d)> {
        if self.points.is_empty() {
            return None;
        }

        let mut min_x = self.points[0].x;
        let mut min_y = self.points[0].y;
        let mut max_x = self.points[0].x;
        let mut max_y = self.points[0].y;

        for point in &self.points {
            min_x = min_x.min(point.x);
            min_y = min_y.min(point.y);
            max_x = max_x.max(point.x);
            max_y = max_y.max(point.y);
        }

        Some((Coord2d::new(min_x, min_y), Coord2d::new(max_x, max_y)))
    }

    /// 获取折线的质心（所有点的平均位置）
    pub fn centroid(&self) -> Option<Coord2d> {
        if self.points.is_empty() {
            return None;
        }

        let sum = self.points.iter().fold(Coord2d::new(0.0, 0.0), |acc, p| {
            Coord2d::new(acc.x + p.x, acc.y + p.y)
        });
        let n = self.points.len() as f64;

        Some(Coord2d::new(sum.x / n, sum.y / n))
    }
}

impl Default for Polyline2d {
    fn default() -> Self {
        Self::new()
    }
}
