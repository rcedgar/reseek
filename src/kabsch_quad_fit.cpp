#include "myutils.h"
#include "kabsch.h"
#include "kabsch_quad.h"

RigidTransform FitQuadRigid(
	Vec3 a0, Vec3 a1, Vec3 a2, Vec3 a3,
	Vec3 b0, Vec3 b1, Vec3 b2, Vec3 b3)
	{
	RigidTransform T;
	const Vec3 A[4] = { a0, a1, a2, a3 };
	const Vec3 B[4] = { b0, b1, b2, b3 };

	double x[4][3];
	double y[4][3];
	for (unsigned i = 0; i < 4; ++i)
		{
		x[i][0] = A[i].x;
		x[i][1] = A[i].y;
		x[i][2] = A[i].z;
		y[i][0] = B[i].x;
		y[i][1] = B[i].y;
		y[i][2] = B[i].z;
		}

	const double *px[4] = { x[0], x[1], x[2], x[3] };
	const double *py[4] = { y[0], y[1], y[2], y[3] };

	double t[3];
	double u[3][3];
	Kabsch(px, py, 4, t, u);

	for (unsigned i = 0; i < 3; ++i)
		for (unsigned j = 0; j < 3; ++j)
			T.R[i][j] = u[i][j];

	T.t.x = t[0];
	T.t.y = t[1];
	T.t.z = t[2];

	double sum2 = 0.0;
	for (unsigned i = 0; i < 4; ++i)
		{
		Vec3 p = ApplyTransform(T, A[i]);
		Vec3 d = p - B[i];
		sum2 += Norm2(d);
		}

	T.rmsd = std::sqrt(sum2 / 4.0);
	return T;
	}
