#include "myutils.h"
#include "coords.h"

void coords::assert_unit_basis(const coords &X,
						const coords &Y,
						const coords &Z)
	{
	asserta(feq(X.norm(), 1));
	asserta(feq(Y.norm(), 1));
	asserta(feq(Z.norm(), 1));

	asserta(feq(coords::angle_radians(X, Y), PI/2));
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

void calculate_theta_phi(const coords &A, const coords &B,
						 const coords &C, const coords &D,
						 float &theta_radians, float &phi_radians)
	{
	coords ab, bc, cd;
	coords::sub(B, A, ab);
	coords::sub(C, B, bc);
	coords::sub(D, C, cd);

	coords xaxis, yaxis, zaxis;
	bc.normalize(xaxis);
	coords::cross(bc, ab, zaxis);
	zaxis.normalize();
	coords::cross(xaxis, yaxis, zaxis);
	coords::assert_unit_basis(xaxis, yaxis, zaxis);

	theta_radians = coords::angle_radians(bc, cd);
	coords vert_bcd;
	coords::cross(bc, cd, vert_bcd);
	phi_radians = coords::angle_radians(zaxis, vert_bcd);
	float check_angle_radians = coords::angle_radians(zaxis, cd);
	if (check_angle_radians < PI/2)
		phi_radians = -phi_radians;
	}

//def calculate_D(A, B, C, theta_rad, phi_rad):
//	ab = subtract(B, A);
//	bc = subtract(C, B);
//	unit_ab = normalize(ab);
//	unit_bc = normalize(bc);
//	vert_abc = cross_product(ab, bc);
//	unit_vert_abc = normalize(vert_abc);
//	cd1 = rotate_around_vector(unit_bc, unit_vert_abc, theta_rad);
//	unit_cd1 = normalize(cd1);
//	cd = rotate_around_vector(unit_cd1, unit_bc, phi_rad);
//	unit_cd = normalize(cd);
//	D1 = add(C, unit_cd);
//	theta2_rad, phi2_rad = calculate_theta_phi(A, B, C, D1)
//	assert feq(theta2_rad, theta_rad)
//	assert feq(phi2_rad, phi_rad)
//	cd_calpha = multiply(unit_cd, 3.81)
//	D = add(C, cd_calpha)
//	dist_CD = get_dist(C, D)
//	assert feq(dist_CD, 3.81)
//	return D
void calculate_D(const coords &A, const coords &B, const coords &C,
				 float &theta_radians, float &phi_radians,
				 coords &D)
	{
	coords ab, bc;
	coords::sub(B, A, ab);
	coords::sub(C, B, bc);

	coords unit_ab, unit_bc;
	coords::normalize(ab, unit_ab);
	coords::normalize(bc, unit_bc);

	coords vert_abc, unit_vert_abc;
	coords::cross(ab, bc, vert_abc);
	coords::normalize(vert_abc, unit_vert_abc);

	coords cd1, unit_cd1;
	coords::rotate(unit_bc, unit_vert_abc, theta_radians, cd1);
	coords::normalize(cd1, unit_cd1);

	coords cd, unit_cd;
	coords::rotate(unit_cd1, unit_bc, phi_radians, cd);
	coords::normalize(cd, unit_cd);

	coords cd_alpha;
	coords::mul(unit_cd, 3.81f, cd_alpha);

	coords::add(C, cd_alpha, D);
#if DEBUG
	{
	assert(feq(cd_alpha.norm(), 3.81f));
	float d = coords::dist(C, D);
	asserta(feq(d, 3.81f));
	}
#endif
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

