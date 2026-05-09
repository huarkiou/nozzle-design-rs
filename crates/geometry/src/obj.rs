use crate::point::Point3d;
use std::io::{BufWriter, Write};

/// A triangle mesh ready for Wavefront OBJ export.
#[derive(Debug, Clone)]
pub struct ObjModel {
    /// Vertex coordinates: `v x y z`
    pub vertices: Vec<[f64; 3]>,
    /// Vertex normals: `vn x y z` (one per vertex)
    pub normals: Vec<[f64; 3]>,
    /// Triangular faces: 1-based indices into `vertices` / `normals`.
    pub faces: Vec<[usize; 3]>,
}

impl ObjModel {
    /// Create an empty model.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            normals: Vec::new(),
            faces: Vec::new(),
        }
    }

    /// Build a triangle mesh by connecting adjacent profile rings.
    ///
    /// Each profile is a closed loop of `n_theta` points.  All profiles must
    /// have the same number of points.  Consecutive profiles are connected by
    /// quad strips, each quad being split into two triangles.
    ///
    /// Vertex normals are automatically computed via [`auto_normal`] after
    /// the faces have been created.
    pub fn from_wallpoints(profiles: &[Vec<Point3d>]) -> Self {
        // ---- edge cases ----
        if profiles.is_empty() {
            return Self::new();
        }
        let n_theta = profiles[0].len();
        if n_theta < 2 {
            return Self::new();
        }
        // Validate all profiles have the same point count.
        for p in profiles {
            if p.len() != n_theta {
                return Self::new();
            }
        }
        let n_profiles = profiles.len();

        // ---- collect vertices ----
        let mut vertices: Vec<[f64; 3]> = Vec::with_capacity(n_profiles * n_theta);
        for profile in profiles {
            for pt in profile {
                vertices.push([pt.x, pt.y, pt.z]);
            }
        }

        // ---- build faces (triangulate quad strips) ----
        let mut faces = Vec::new();
        let degenerate_eps: f64 = 1e-15;
        // Only create faces when there are at least 2 profiles.
        for i in 0..n_profiles.saturating_sub(1) {
            let base_i = i * n_theta;
            let base_ip1 = (i + 1) * n_theta;
            for j in 0..n_theta {
                let j_next = (j + 1) % n_theta;

                // 0-based vertex indices
                let ai = base_i + j;
                let bi = base_i + j_next;
                let ci = base_ip1 + j_next;
                let di = base_ip1 + j;

                // Triangle 1: a, b, c
                if triangle_area(&vertices[ai], &vertices[bi], &vertices[ci]) > degenerate_eps {
                    faces.push([ai + 1, bi + 1, ci + 1]);
                }
                // Triangle 2: a, c, d
                if triangle_area(&vertices[ai], &vertices[ci], &vertices[di]) > degenerate_eps {
                    faces.push([ai + 1, ci + 1, di + 1]);
                }
            }
        }

        // ---- normals placeholder ----
        let normals = vec![[0.0_f64; 3]; vertices.len()];

        let mut model = Self {
            vertices,
            normals,
            faces,
        };
        model.auto_normal();
        model
    }

    /// Write the model as a Wavefront OBJ file.
    ///
    /// Coordinates are written with 6 decimal places.
    /// Face lines use the format `f v1//vn1 v2//vn2 v3//vn3`.
    pub fn write_obj(&self, path: &str) -> std::io::Result<()> {
        let file = std::fs::File::create(path)?;
        let mut w = BufWriter::new(file);

        writeln!(w, "# OBJ file")?;

        // vertices
        for v in &self.vertices {
            writeln!(w, "v {:.6} {:.6} {:.6}", v[0], v[1], v[2])?;
        }

        // normals
        for n in &self.normals {
            writeln!(w, "vn {:.6} {:.6} {:.6}", n[0], n[1], n[2])?;
        }

        // faces
        for face in &self.faces {
            writeln!(
                w,
                "f {}//{} {}//{} {}//{}",
                face[0], face[0], face[1], face[1], face[2], face[2],
            )?;
        }

        w.flush()?;
        Ok(())
    }

    /// Compute per-vertex normals as the area-weighted average of adjacent
    /// face normals.
    ///
    /// The `normals` vector is resized to match `vertices` before computation
    /// and each entry is initialised to zero.  Zero-length accumulated
    /// normals are left as `[0, 0, 0]` instead of producing NaN.
    pub fn auto_normal(&mut self) {
        // Ensure normals vector matches vertices length
        self.normals.resize(self.vertices.len(), [0.0; 3]);
        for n in &mut self.normals {
            *n = [0.0; 3];
        }

        for face in &self.faces {
            let i1 = face[0].saturating_sub(1);
            let i2 = face[1].saturating_sub(1);
            let i3 = face[2].saturating_sub(1);

            if i1 >= self.vertices.len() || i2 >= self.vertices.len() || i3 >= self.vertices.len() {
                continue;
            }

            let v1 = &self.vertices[i1];
            let v2 = &self.vertices[i2];
            let v3 = &self.vertices[i3];

            let n = face_normal(v1, v2, v3);

            // Area = 0.5 * |cross product|
            let area = 0.5 * vector_length(&n);

            // Normalize the face normal
            let n_unit = normalize(&n);

            // Accumulate area-weighted normal for each vertex
            for idx in &[i1, i2, i3] {
                self.normals[*idx][0] += area * n_unit[0];
                self.normals[*idx][1] += area * n_unit[1];
                self.normals[*idx][2] += area * n_unit[2];
            }
        }

        // Normalize accumulated vertex normals
        for n in &mut self.normals {
            *n = normalize(n);
        }
    }
}

impl Default for ObjModel {
    fn default() -> Self {
        Self::new()
    }
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Compute the unsigned area of triangle (v1, v2, v3).
/// Returns 0.0 for degenerate (collinear) triangles.
fn triangle_area(v1: &[f64; 3], v2: &[f64; 3], v3: &[f64; 3]) -> f64 {
    let u = [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]];
    let v = [v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]];
    let cr = cross(&u, &v);
    0.5 * vector_length(&cr)
}

/// Compute the unit face normal (cross product of two edge vectors).
///
/// Returns a unit vector pointing outward according to the right-hand rule
/// `(v2 - v1) × (v3 - v1)`.  If the cross product is zero the zero vector
/// is returned.
pub fn face_normal(v1: &[f64; 3], v2: &[f64; 3], v3: &[f64; 3]) -> [f64; 3] {
    let u = [v2[0] - v1[0], v2[1] - v1[1], v2[2] - v1[2]];
    let v = [v3[0] - v1[0], v3[1] - v1[1], v3[2] - v1[2]];
    cross(&u, &v)
}

/// Cross product of two 3D vectors.
fn cross(a: &[f64; 3], b: &[f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Euclidean length of a 3D vector.
fn vector_length(v: &[f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Normalize a 3D vector in place.  Returns the zero vector unchanged.
fn normalize(v: &[f64; 3]) -> [f64; 3] {
    let len = vector_length(v);
    if len < f64::EPSILON {
        [0.0; 3]
    } else {
        [v[0] / len, v[1] / len, v[2] / len]
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_empty() {
        let m = ObjModel::new();
        assert!(m.vertices.is_empty());
        assert!(m.normals.is_empty());
        assert!(m.faces.is_empty());
    }

    #[test]
    fn test_default_is_empty() {
        let m = ObjModel::default();
        assert!(m.vertices.is_empty());
    }

    #[test]
    fn test_empty_profiles() {
        let m = ObjModel::from_wallpoints(&[]);
        assert!(m.vertices.is_empty());
        assert!(m.faces.is_empty());
    }

    #[test]
    fn test_single_profile_no_faces() {
        // Single profile with 4 points should produce vertices but no faces.
        let profile = vec![vec![
            Point3d::new(0.0, 0.0, 0.0),
            Point3d::new(1.0, 0.0, 0.0),
            Point3d::new(1.0, 1.0, 0.0),
            Point3d::new(0.0, 1.0, 0.0),
        ]];
        let m = ObjModel::from_wallpoints(&profile);
        assert_eq!(m.vertices.len(), 4);
        assert_eq!(m.faces.len(), 0);
        assert_eq!(m.normals.len(), 4);
    }

    #[test]
    fn test_zero_point_profile() {
        let profiles = vec![vec![]];
        let m = ObjModel::from_wallpoints(&profiles);
        assert!(m.vertices.is_empty());
        assert!(m.faces.is_empty());
    }

    #[test]
    fn test_two_profiles_quad() {
        // Two profiles with 4 points each → 4 quads → 8 triangles.
        let profile0 = vec![
            Point3d::new(0.0, 0.0, 0.0),
            Point3d::new(1.0, 0.0, 0.0),
            Point3d::new(1.0, 1.0, 0.0),
            Point3d::new(0.0, 1.0, 0.0),
        ];
        let profile1 = vec![
            Point3d::new(0.0, 0.0, 1.0),
            Point3d::new(1.0, 0.0, 1.0),
            Point3d::new(1.0, 1.0, 1.0),
            Point3d::new(0.0, 1.0, 1.0),
        ];
        let m = ObjModel::from_wallpoints(&[profile0, profile1]);
        assert_eq!(m.vertices.len(), 8);
        assert_eq!(m.faces.len(), 8); // 4 quads × 2 triangles
        assert_eq!(m.normals.len(), 8);

        // Check 1-based indexing on first face
        let first_face = m.faces[0];
        assert!(first_face[0] >= 1 && first_face[0] <= 8);
        assert!(first_face[1] >= 1 && first_face[1] <= 8);
        assert!(first_face[2] >= 1 && first_face[2] <= 8);

        // Verify the specific triangle pattern for j=0
        // p0[0]=idx1, p0[1]=idx2, p1[1]=idx6 → triangle [1,2,6]
        // p0[0]=idx1, p1[1]=idx6, p1[0]=idx5 → triangle [1,6,5]
        assert_eq!(m.faces[0], [1, 2, 6]);
        assert_eq!(m.faces[1], [1, 6, 5]);
    }

    #[test]
    fn test_face_normal_unit_square() {
        // A triangle in the XY plane → normal should be [0, 0, 1]
        let v1 = [0.0, 0.0, 0.0];
        let v2 = [1.0, 0.0, 0.0];
        let v3 = [0.0, 1.0, 0.0];
        let n = face_normal(&v1, &v2, &v3);
        // (v2-v1) × (v3-v1) = (1,0,0) × (0,1,0) = (0,0,1)
        assert!((n[0] - 0.0).abs() < 1e-12);
        assert!((n[1] - 0.0).abs() < 1e-12);
        assert!((n[2] - 1.0).abs() < 1e-12);
    }

    #[test]
    fn test_face_normal_degenerate() {
        // Collinear points → zero normal
        let v1 = [0.0, 0.0, 0.0];
        let v2 = [1.0, 0.0, 0.0];
        let v3 = [2.0, 0.0, 0.0];
        let n = face_normal(&v1, &v2, &v3);
        assert!((n[0]).abs() < 1e-12);
        assert!((n[1]).abs() < 1e-12);
        assert!((n[2]).abs() < 1e-12);
    }

    #[test]
    fn test_auto_normal_flat_square() {
        // Two triangles forming a unit square in the XY plane
        let mut m = ObjModel {
            vertices: vec![
                [0.0, 0.0, 0.0], // 1
                [1.0, 0.0, 0.0], // 2
                [1.0, 1.0, 0.0], // 3
                [0.0, 1.0, 0.0], // 4
            ],
            normals: vec![[0.0; 3]; 4],
            faces: vec![[1, 2, 3], [1, 3, 4]],
        };
        m.auto_normal();
        // All normals should point in +z direction
        for n in &m.normals {
            let len = (n[0] * n[0] + n[1] * n[1] + n[2] * n[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-12 || len < 1e-12);
            if len > 0.5 {
                assert!((n[2].abs() - 1.0).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn test_write_obj_output() {
        let m = ObjModel {
            vertices: vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            normals: vec![[0.0, 0.0, 1.0], [0.0, 0.0, 1.0], [0.0, 0.0, 1.0]],
            faces: vec![[1, 2, 3]],
        };

        let path = std::env::temp_dir().join("test_nozzle.obj");
        let path_str = path.to_string_lossy();
        m.write_obj(&path_str).expect("write should succeed");

        let content = std::fs::read_to_string(&path).expect("read should succeed");
        let lines: Vec<&str> = content.lines().collect();

        assert!(lines[0].starts_with("# OBJ file"));
        assert_eq!(lines[1], "v 0.000000 0.000000 0.000000");
        assert_eq!(lines[2], "v 1.000000 0.000000 0.000000");
        assert_eq!(lines[3], "v 0.000000 1.000000 0.000000");
        assert_eq!(lines[4], "vn 0.000000 0.000000 1.000000");
        assert_eq!(lines[5], "vn 0.000000 0.000000 1.000000");
        assert_eq!(lines[6], "vn 0.000000 0.000000 1.000000");
        assert_eq!(lines[7], "f 1//1 2//2 3//3");

        let _ = std::fs::remove_file(&path);
    }

    #[test]
    fn test_inconsistent_profile_lengths() {
        let p0 = vec![Point3d::new(0.0, 0.0, 0.0), Point3d::new(1.0, 0.0, 0.0)];
        let p1 = vec![
            Point3d::new(0.0, 0.0, 1.0),
            Point3d::new(1.0, 0.0, 1.0),
            Point3d::new(2.0, 0.0, 1.0),
        ];
        let m = ObjModel::from_wallpoints(&[p0, p1]);
        assert!(m.vertices.is_empty());
    }
}
