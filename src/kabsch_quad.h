#pragma once

#include <cmath>
#include <algorithm>
#include <cassert>

struct Vec3
{
    double x, y, z;
};

struct RigidTransform
{
    double R[3][3]; // row-major rotation
    Vec3 t;         // translation
    double rmsd;
};

static inline Vec3 operator+(Vec3 a, Vec3 b)
{
    return { a.x + b.x, a.y + b.y, a.z + b.z };
}

static inline Vec3 operator-(Vec3 a, Vec3 b)
{
    return { a.x - b.x, a.y - b.y, a.z - b.z };
}

static inline Vec3 operator*(double s, Vec3 a)
{
    return { s*a.x, s*a.y, s*a.z };
}

static inline double Dot(Vec3 a, Vec3 b)
{
    return a.x*b.x + a.y*b.y + a.z*b.z;
}

static inline double Norm2(Vec3 a)
{
    return Dot(a, a);
}

static inline Vec3 ApplyRotation(const double R[3][3], Vec3 p)
{
    return {
        R[0][0]*p.x + R[0][1]*p.y + R[0][2]*p.z,
        R[1][0]*p.x + R[1][1]*p.y + R[1][2]*p.z,
        R[2][0]*p.x + R[2][1]*p.y + R[2][2]*p.z
    };
}

static inline Vec3 ApplyTransform(const RigidTransform& T, Vec3 p)
{
    Vec3 r = ApplyRotation(T.R, p);
    return r + T.t;
}

static inline void QuaternionToMatrix(
    double q0, double q1, double q2, double q3,
    double R[3][3])
{
    // q0 = scalar part, q1/q2/q3 = vector part
    const double n = std::sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    assert(n > 0);

    q0 /= n;
    q1 /= n;
    q2 /= n;
    q3 /= n;

    R[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R[0][1] = 2.0*(q1*q2 - q0*q3);
    R[0][2] = 2.0*(q1*q3 + q0*q2);

    R[1][0] = 2.0*(q2*q1 + q0*q3);
    R[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    R[1][2] = 2.0*(q2*q3 - q0*q1);

    R[2][0] = 2.0*(q3*q1 - q0*q2);
    R[2][1] = 2.0*(q3*q2 + q0*q1);
    R[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}

static inline void LargestEigenvector4x4_PowerIteration(
    const double A[4][4],
    double q[4])
{
    // Horn's N has trace 0, so eigenvalues are mixed sign. Plain power
    // iteration converges to an eigenvector for max |lambda|, not the
    // algebraically largest eigenvalue required for the optimal quaternion.
    // Shifting by alpha > ||A||_2 makes all eigenvalues of (A + alpha*I)
    // positive (lambda + alpha > 0), and the dominant mode is then the
    // eigenvector for max lambda(A), unchanged by the shift.
    double nf2 = 0.0;
    for (unsigned i = 0; i < 4; ++i)
        for (unsigned j = 0; j < 4; ++j)
        {
            const double a = A[i][j];
            nf2 += a * a;
        }
    const double nf = std::sqrt(nf2);
    const double shift = nf + 1.0;

    q[0] = 1.0;
    q[1] = 0.0;
    q[2] = 0.0;
    q[3] = 0.0;

    for (unsigned iter = 0; iter < 128; ++iter)
    {
        double v[4] = {};

        for (unsigned i = 0; i < 4; ++i)
        {
            double s = shift * q[i];
            for (unsigned j = 0; j < 4; ++j)
                s += A[i][j] * q[j];
            v[i] = s;
        }

        const double n = std::sqrt(
            v[0]*v[0] + v[1]*v[1] + v[2]*v[2] + v[3]*v[3]);

        if (n <= 0.0)
            break;

        for (unsigned i = 0; i < 4; ++i)
            q[i] = v[i] / n;
    }
}

static inline RigidTransform FitTriangleRigid(
    Vec3 a0, Vec3 a1, Vec3 a2,
    Vec3 b0, Vec3 b1, Vec3 b2)
{
    RigidTransform T;
    const Vec3 A[3] = { a0, a1, a2 };
    const Vec3 B[3] = { b0, b1, b2 };

    Vec3 ca = (1.0/3.0)*(a0 + a1 + a2);
    Vec3 cb = (1.0/3.0)*(b0 + b1 + b2);

    // Cross-covariance S = sum_i x_i y_i^T,
    // where x_i = A_i - ca, y_i = B_i - cb.
    double Sxx = 0, Sxy = 0, Sxz = 0;
    double Syx = 0, Syy = 0, Syz = 0;
    double Szx = 0, Szy = 0, Szz = 0;

    for (unsigned i = 0; i < 3; ++i)
    {
        Vec3 x = A[i] - ca;
        Vec3 y = B[i] - cb;

        Sxx += x.x*y.x; Sxy += x.x*y.y; Sxz += x.x*y.z;
        Syx += x.y*y.x; Syy += x.y*y.y; Syz += x.y*y.z;
        Szx += x.z*y.x; Szy += x.z*y.y; Szz += x.z*y.z;
    }

    // Horn's 4x4 symmetric matrix.
    double N[4][4];

    N[0][0] = Sxx + Syy + Szz;
    N[0][1] = Syz - Szy;
    N[0][2] = Szx - Sxz;
    N[0][3] = Sxy - Syx;

    N[1][0] = N[0][1];
    N[1][1] = Sxx - Syy - Szz;
    N[1][2] = Sxy + Syx;
    N[1][3] = Szx + Sxz;

    N[2][0] = N[0][2];
    N[2][1] = N[1][2];
    N[2][2] = -Sxx + Syy - Szz;
    N[2][3] = Syz + Szy;

    N[3][0] = N[0][3];
    N[3][1] = N[1][3];
    N[3][2] = N[2][3];
    N[3][3] = -Sxx - Syy + Szz;

    double q[4];
    LargestEigenvector4x4_PowerIteration(N, q);

    QuaternionToMatrix(q[0], q[1], q[2], q[3], T.R);

    // t = centroid_B - R * centroid_A
    Vec3 Rca = ApplyRotation(T.R, ca);
    T.t = cb - Rca;

    double sum2 = 0.0;
    for (unsigned i = 0; i < 3; ++i)
    {
        Vec3 p = ApplyTransform(T, A[i]);
        Vec3 d = p - B[i];
        sum2 += Norm2(d);
    }

    T.rmsd = std::sqrt(sum2 / 3.0);
    return T;
}

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <random>
#include <cassert>

static double Dot3(const double a[3], const double b[3])
{
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

static void MatVec(const double R[3][3], const double p[3], double q[3])
{
    q[0] = R[0][0]*p[0] + R[0][1]*p[1] + R[0][2]*p[2];
    q[1] = R[1][0]*p[0] + R[1][1]*p[1] + R[1][2]*p[2];
    q[2] = R[2][0]*p[0] + R[2][1]*p[1] + R[2][2]*p[2];
}

static inline void Transpose3x3(double R[3][3])
{
    double tmp;

    tmp = R[0][1]; R[0][1] = R[1][0]; R[1][0] = tmp;
    tmp = R[0][2]; R[0][2] = R[2][0]; R[2][0] = tmp;
    tmp = R[1][2]; R[1][2] = R[2][1]; R[2][1] = tmp;
}

static void ApplyTransform(
    const double R[3][3],
    const double t[3],
    const double p[3],
    double q[3])
{
    MatVec(R, p, q);
    q[0] += t[0];
    q[1] += t[1];
    q[2] += t[2];
}

static void MatMulATB(
    const double A[3][3],
    const double B[3][3],
    double C[3][3])
{
    // C = A^T * B
    for (unsigned i = 0; i < 3; ++i)
        for (unsigned j = 0; j < 3; ++j)
            C[i][j] =
                A[0][i]*B[0][j] +
                A[1][i]*B[1][j] +
                A[2][i]*B[2][j];
}

static double RotationAngleBetween(
    const double R1[3][3],
    const double R2[3][3])
{
    // Relative rotation: D = R1^T * R2
    double D[3][3];
    MatMulATB(R1, R2, D);

    const double tr = D[0][0] + D[1][1] + D[2][2];
    double c = 0.5*(tr - 1.0);

    if (c < -1.0) c = -1.0;
    if (c >  1.0) c =  1.0;

    return std::acos(c); // radians
}

static void RandomRotation(std::mt19937_64& rng, double R[3][3])
{
    std::normal_distribution<double> normal(0.0, 1.0);

    double q0 = normal(rng);
    double q1 = normal(rng);
    double q2 = normal(rng);
    double q3 = normal(rng);

    double n = std::sqrt(q0*q0 + q1*q1 + q2*q2 + q3*q3);
    q0 /= n;
    q1 /= n;
    q2 /= n;
    q3 /= n;

    R[0][0] = q0*q0 + q1*q1 - q2*q2 - q3*q3;
    R[0][1] = 2.0*(q1*q2 - q0*q3);
    R[0][2] = 2.0*(q1*q3 + q0*q2);

    R[1][0] = 2.0*(q2*q1 + q0*q3);
    R[1][1] = q0*q0 - q1*q1 + q2*q2 - q3*q3;
    R[1][2] = 2.0*(q2*q3 - q0*q1);

    R[2][0] = 2.0*(q3*q1 - q0*q2);
    R[2][1] = 2.0*(q3*q2 + q0*q1);
    R[2][2] = q0*q0 - q1*q1 - q2*q2 + q3*q3;
}

static double Dist3(const double a[3], const double b[3])
{
    const double dx = a[0] - b[0];
    const double dy = a[1] - b[1];
    const double dz = a[2] - b[2];
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

/** Four-point rigid fit via Kabsch (closed-form); prefer for generic point-cloud seeds. */
RigidTransform FitQuadRigid(
    Vec3 a0, Vec3 a1, Vec3 a2, Vec3 a3,
    Vec3 b0, Vec3 b1, Vec3 b2, Vec3 b3);