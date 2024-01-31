use std::{
    borrow::Borrow,
    ops::{Add, Mul, Sub},
};

#[derive(Clone, Copy, PartialEq, Debug)]
pub struct Vector3D {
    pub i: f64,
    pub j: f64,
    pub k: f64,
}

impl Vector3D {
    pub fn new(i: impl Into<f64>, j: impl Into<f64>, k: impl Into<f64>) -> Self {
        Self {
            i: i.into(),
            j: j.into(),
            k: k.into(),
        }
    }

    pub fn length(&self) -> f64 {
        (self.dot(self)).sqrt()
    }

    pub const fn to_array(self) -> [f64; 3] {
        [self.i, self.j, self.k]
    }

    pub fn dot(&self, other: &Self) -> f64 {
        self.i
            .mul_add(other.i, self.j.mul_add(other.j, self.k * other.k))
    }

    pub fn distance_to(&self, b: &Self) -> f64 {
        (b - self).length()
    }
}

impl<'a, B> Add<B> for &'a Vector3D
where
    B: Borrow<Vector3D>,
{
    type Output = Vector3D;

    fn add(self, rhs: B) -> Self::Output {
        Vector3D {
            i: self.i + rhs.borrow().i,
            j: self.j + rhs.borrow().j,
            k: self.k + rhs.borrow().k,
        }
    }
}

impl<B> Add<B> for Vector3D
where
    B: Borrow<Self>,
{
    type Output = Self;

    fn add(self, rhs: B) -> Self::Output {
        &self + rhs
    }
}

impl<'a, B> Sub<B> for &'a Vector3D
where
    B: Borrow<Vector3D>,
{
    type Output = Vector3D;

    fn sub(self, rhs: B) -> Self::Output {
        Vector3D {
            i: self.i - rhs.borrow().i,
            j: self.j - rhs.borrow().j,
            k: self.k - rhs.borrow().k,
        }
    }
}

impl<B> Sub<B> for Vector3D
where
    B: Borrow<Self>,
{
    type Output = Self;
    fn sub(self, rhs: B) -> Self::Output {
        &self - rhs
    }
}

impl<T> Mul<T> for Vector3D
where
    T: Into<f64>,
{
    type Output = Self;

    fn mul(self, scalar: T) -> Self {
        let scalar_f64: f64 = scalar.into();
        Self {
            i: self.i * scalar_f64,
            j: self.j * scalar_f64,
            k: self.k * scalar_f64,
        }
    }
}
#[cfg(test)]
mod linear_algebra_tests {
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

        assert_eq!(&v1 + &v2, Vector3D::new(3.0, 5.0, 7.0));
        assert_eq!(v1 + &v2, Vector3D::new(3.0, 5.0, 7.0));
        assert_eq!(&v1 + v2, Vector3D::new(3.0, 5.0, 7.0));
        assert_eq!(v1 + v2, Vector3D::new(3.0, 5.0, 7.0));
    }

    #[test]
    fn test_sub_vector3d_to_vector3d() {
        let v1 = Vector3D::new(1.0, 2.0, 3.0);
        let v2 = Vector3D::new(2.0, 3.0, 4.0);

        assert_eq!(&v2 - &v1, Vector3D::new(1.0, 1.0, 1.0));
        assert_eq!(v2 - &v1, Vector3D::new(1.0, 1.0, 1.0));
        assert_eq!(&v2 - v1, Vector3D::new(1.0, 1.0, 1.0));
        assert_eq!(v2 - v1, Vector3D::new(1.0, 1.0, 1.0));
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

    #[test]
    fn test_vector_distance() {
        let v0 = Vector3D::new(0, 0, 0);
        let v1 = Vector3D::new(1, 1, 1);

        assert_eq!(Vector3D::distance_to(&v0, &v1), 3.0_f64.sqrt());
        assert_eq!(v1.distance_to(&v0), 3.0_f64.sqrt());
    }
}
