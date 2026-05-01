#include "myutils.h"
#include "kabsch_quad.h"

static double AbsTetraVolume6(const double A[4][3])
	{
	const double ax = A[1][0] - A[0][0];
	const double ay = A[1][1] - A[0][1];
	const double az = A[1][2] - A[0][2];
	const double bx = A[2][0] - A[0][0];
	const double by = A[2][1] - A[0][1];
	const double bz = A[2][2] - A[0][2];
	const double cx = A[3][0] - A[0][0];
	const double cy = A[3][1] - A[0][1];
	const double cz = A[3][2] - A[0][2];
	const double rx = by * cz - bz * cy;
	const double ry = bz * cx - bx * cz;
	const double rz = bx * cy - by * cx;
	return std::fabs(ax * rx + ay * ry + az * rz);
	}

static void TestFitQuadRigid()
	{
	std::mt19937_64 rng(1);

	std::uniform_real_distribution<double> coord(-10.0, 10.0);
	std::uniform_real_distribution<double> trans(-100.0, 100.0);
	std::normal_distribution<double> noise(0.0, 0.01);

	for (unsigned test = 0; test < 100; ++test)
		{
		double A[4][3];
		double Ap[4][3];
		double B[4][3];

		// 1. Generate a non-degenerate tetrahedron (|scalar triple product|).
		unsigned gen = 0;
		for (;;)
			{
			for (unsigned i = 0; i < 4; ++i)
				for (unsigned k = 0; k < 3; ++k)
					A[i][k] = coord(rng);
			if (AbsTetraVolume6(A) > 0.5)
				break;
			asserta(++gen < 1000);
			}

		// 2. Perturbed copy.
		for (unsigned i = 0; i < 4; ++i)
			for (unsigned k = 0; k < 3; ++k)
				Ap[i][k] = A[i][k] + noise(rng);

		// 3. Random transform.
		double Rtrue[3][3];
		double ttrue[3];

		RandomRotation(rng, Rtrue);

		ttrue[0] = trans(rng);
		ttrue[1] = trans(rng);
		ttrue[2] = trans(rng);

		// 4. Apply transform to perturbed copy.
		for (unsigned i = 0; i < 4; ++i)
			ApplyTransform(Rtrue, ttrue, Ap[i], B[i]);

		Vec3 a0{ A[0][0], A[0][1], A[0][2] };
		Vec3 a1{ A[1][0], A[1][1], A[1][2] };
		Vec3 a2{ A[2][0], A[2][1], A[2][2] };
		Vec3 a3{ A[3][0], A[3][1], A[3][2] };

		Vec3 b0{ B[0][0], B[0][1], B[0][2] };
		Vec3 b1{ B[1][0], B[1][1], B[1][2] };
		Vec3 b2{ B[2][0], B[2][1], B[2][2] };
		Vec3 b3{ B[3][0], B[3][1], B[3][2] };

		const RigidTransform T = FitQuadRigid(a0, a1, a2, a3, b0, b1, b2, b3);

		double sum2 = 0.0;
		for (unsigned i = 0; i < 4; ++i)
			{
			double q[3];
			const double tfit[3] = { T.t.x, T.t.y, T.t.z };

			ApplyTransform(T.R, tfit, A[i], q);

			const double d = Dist3(q, B[i]);
			sum2 += d * d;
			}

		const double rmsd_check = std::sqrt(sum2 / 4.0);

		const double rot_err_rad = RotationAngleBetween(Rtrue, T.R);
		const double rot_err_deg = rot_err_rad * 180.0 / 3.14159265358979323846;

		const double tfit[3] = { T.t.x, T.t.y, T.t.z };
		const double trans_err = Dist3(ttrue, tfit);

		ProgressLog(
			"test %u  rmsd %.3g  rmsd_check %.3g  rot_err_deg %.3g  trans_err %.3g\n",
			test,
			T.rmsd,
			rmsd_check,
			rot_err_deg,
			trans_err);

		assert(std::fabs(T.rmsd - rmsd_check) < 1e-9);
		}
	}

void cmd_kabsch_quad()
	{
	TestFitQuadRigid();
	}
