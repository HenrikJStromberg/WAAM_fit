use std::ops::{Add, Mul};

#[derive(Clone, Copy, PartialEq, Debug)]
struct Vector3D {
    i: f64,
    j: f64,
    k: f64,
}

impl Vector3D {
    pub fn new(i: impl Into<f64>, j: impl Into<f64>, k: impl Into<f64>) -> Self {
        Vector3D {
            i: i.into(),
            j: j.into(),
            k: k.into(),
        }
    }

    pub fn length(&self) -> f64 {
        (self.i * self.i + self.j * self.j + self.k * self.k).sqrt()
    }

    pub fn to_array(self) -> [f64; 3] {
        [self.i, self.j, self.k]
    }

    pub fn dot(&self, other: &Vector3D) -> f64 {
        self.i * other.i + self.j * other.j + self.k * other.k
    }
}

impl Add for Vector3D {
    type Output = Vector3D;

    fn add(self, rhs: Self) -> Self::Output {
        Self {
            i: self.i + rhs.i,
            j: self.j + rhs.j,
            k: self.k + rhs.k,
        }
    }
}

impl<T> Mul<T> for Vector3D
where
    T: Into<f64>,
{
    type Output = Vector3D;

    fn mul(self, scalar: T) -> Vector3D {
        let scalar_f64: f64 = scalar.into();
        Vector3D {
            i: self.i * scalar_f64,
            j: self.j * scalar_f64,
            k: self.k * scalar_f64,
        }
    }
}
#[cfg(test)]
mod brep_element_tests {
    use super::*;

    #[test]
    fn test_create_vector_from_non_f64() {
        assert_eq!(
            Vector3D::new(1, 2, 3),
            Vector3D {
                i: 1.0,
                j: 2.0,
                k: 3.0
            }
        )
    }
    #[test]
    fn test_add_vector3d_to_vector3d() {
        let v1 = Vector3D::new(1.0, 2.0, 3.0);
        let v2 = Vector3D::new(2.0, 3.0, 4.0);

        assert_eq!(v1 + v2, Vector3D::new(3.0, 5.0, 7.0))
    }

    #[test]
    fn test_mul_vector3d_scalar() {
        let v = Vector3D::new(1.0, 2.0, 3.0);
        assert_eq!(v * 3, Vector3D::new(3.0, 6.0, 9.0))
    }

    #[test]
    fn test_vector_length() {
        let v = Vector3D::new(2.0, 3.0, 4.0);
        assert_eq!(v.length(), 29.0_f64.sqrt())
    }

    #[test]
    fn test_vector_to_array() {
        let v = Vector3D::new(2.0, 3.0, 4.0);
        assert_eq!(v.to_array(), [2.0, 3.0, 4.0])
    }

    #[test]
    fn test_dot_product() {
        let v1 = Vector3D::new(1, 0, 0);
        let v2 = Vector3D::new(2, 2, 2);

        assert_eq!(v1.dot(&v2), 2.0);
    }
}
