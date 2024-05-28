# OpenQuaternion
It is a tiny c++ codes consisting solely of header files.
## Brief guide for functions
###CSMQuaternion Rotate(const CSMQuaternion P)
It returns a rotated quaternion.
### Arithmetic operations
Addition(+, +=), Multiplication(%, %=)
### Conjugate()
It returns a conjugate quaternion.
### Partial derivative of quaternion rotation for 3D point(x, y, z)
Respect to w, x, y, and z, respectively
Rotate_dw(x, y, z)
Rotate_dx(x, y, z)
Rotate_dy(x, y, z)
Rotate_dz(x, y, z)
### Norm2()
It returns the square sum of w, x, y, and z
### Norm()
It returns the size of a quaternion.
### Normalize()
It normalize itself.
### Quaternion2Euler()
It returns a 3by3 rotation matrix equivalent to the quaternion.
### Euler2Quaternion(R)
It returns a quaternion equivalent to the rotation matrix ,R.
### Interpolation(qa, qb, ta, tb, t, th)
It returns a spherical interpolated quaternion.
### InterpolationLerp(qa, qb, u)
It returns a lerp interpolated quaternion.
### Rotation about an arbitrary axis
inline CSMQuaternion rotateArbitraryAxis(
    const double x, const double y, const double z,
    const double rx, const double ry, const double rz,
    const double theta)
It rotates a point(x, y, z) about the vector(rx, ry, rz) by a rotation angle, theta.
