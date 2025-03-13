#include "myutils.h"
#include "coords.h"

void coords::assert_unit_basis(const coords &X,
						const coords &Y,
						const coords &Z)
	{
	asserta(feq(X.norm(), 1));
	asserta(feq(Y.norm(), 1));
	asserta(feq(Z.norm(), 1));

	float angle = coords::angle_radians(X, Y);
	asserta(feq(angle, PI/2));
	}

// https://en.wikipedia.org/wiki/Rodrigues%27_rotation_formula
void coords::rotate(const coords &v, const coords &axis, float theta,
					coords &v_rotated)
	{
	coords k;
	axis.normalize(k);

	float k_dot_v = dot(k, v);

	coords k_cross_v;
	coords::cross(v, k, k_cross_v);

	float sin_theta = sin(theta);
	float cos_theta = cos(theta);

	coords v_cos_theta;
	coords::mul(v, cos_theta, v_cos_theta);

	coords k_cross_v_sin_theta;
	coords::mul(k_cross_v, sin_theta, k_cross_v_sin_theta);

	float factor = k_dot_v*(1 - cos_theta);
	coords kkvc;
	coords::mul(k, factor, kkvc);

	v_rotated.x = v_cos_theta.x + k_cross_v_sin_theta.x + kkvc.x;
	v_rotated.y = v_cos_theta.y + k_cross_v_sin_theta.y + kkvc.y;
	v_rotated.z = v_cos_theta.z + k_cross_v_sin_theta.z + kkvc.z;
	}

void cmd_test_coords()
	{
	coords xaxis(1, 0, 0);
	coords yaxis(0, 1, 0);
	coords zaxis;
	coords::cross(xaxis, yaxis, zaxis);
	zaxis.logme();
	brk(1);
	}

