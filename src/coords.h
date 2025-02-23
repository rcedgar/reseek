#pragma once

#include "mypi.h"
static const float PI_plus_epsilon = 3.142f;

static float degrees(float radians)
	{
	return radians*180.0f/PI;
	}

static float radians(float degrees)
	{
	return degrees*PI/180.0f;
	}

class coords
	{
public:
	float x;
	float y;
	float z;

public:
	coords()
		{
		x = FLT_MAX;
		y = FLT_MAX;
		z = FLT_MAX;
		}

	coords(float _x, float _y, float _z)
		{
		x = _x;
		y = _y;
		z = _z;
		}

	coords(const coords &B)
		{
		x = B.x;
		y = B.y;
		z = B.z;
		}

	coords &operator=(const coords &B)
		{
		x = B.x;
		y = B.y;
		z = B.z;
		return *this;
		}

	coords &operator+(const coords &B)
		{
		x += B.x;
		y += B.y;
		z += B.z;
		return *this;
		}

	coords &operator-(const coords &B)
		{
		x -= B.x;
		y -= B.y;
		z -= B.z;
		return *this;
		}

	coords &operator*=(float a)
		{
		x *= a;
		y *= a;
		z *= a;
		return *this;
		}

	void asserteq(const coords &B) const
		{
		assert(feq(x, B.x));
		assert(feq(y, B.y));
		assert(feq(z, B.z));
		}

	float norm() const
		{
		return sqrt(x*x + y*y + z*z);
		}

	coords &normalize()
		{
		float mod = norm();
		x /= mod;
		y /= mod;
		z /= mod;
		return *this;
		}

	const coords &normalize(coords &result) const
		{
		float mod = norm();
		result.x = x/mod;
		result.y = y/mod;
		result.z = z/mod;
		return result;
		}

	void logme(const char *name = "")
		{
		if (*name != 0)
			Log("%s = ", name);
		Log(" x=%.3g y= %.3g z=%.3g\n", x, y, z);
		}

public:
	static coords &cross(const coords &A, const coords &B,
						 coords &result)
		{
		result.x = A.y*B.z - A.z*B.y;
		result.y = A.z*B.x - A.x*B.z;
		result.z = A.x*B.y - A.y*B.x;
#if DEBUG
		{
		float angleAr = angle_radians(A, result);
		float angleBr = angle_radians(B, result);
		asserta(feq(angleAr, PI/2));
		asserta(feq(angleBr, PI/2));
		}
#endif
		return result;
		}

	static float dot(const coords &A, const coords &B)
		{
		float xx = A.x*B.x;
		float yy = A.y*B.y;
		float zz = A.z*B.z;
		return xx + yy + zz;
		}

	static float dist(const coords &A, const coords &B)
		{
		float dx = A.x - B.x;
		float dy = A.y - B.y;
		float dz = A.z - B.z;
		return sqrt(dx*dx + dy*dy + dz*dz);
		}

	static void mul(const coords &v, float a, coords &result)
		{
		result.x = v.x*a;
		result.y = v.y*a;
		result.z = v.z*a;
		}

	static void add(const coords &A, const coords &B, coords &result)
		{
		result.x = A.x + B.x;
		result.y = A.y + B.y;
		result.z = A.z + B.z;
		}

	static void sub(const coords &A, const coords &B, coords &result)
		{
		result.x = A.x - B.x;
		result.y = A.y - B.y;
		result.z = A.z - B.z;
		}

	static void normalize(const coords &A, coords &result)
		{
		result = A;
		result.normalize();
		}

// Angle between vectors in range 0 to pi
	static float angle_radians(const coords &A, const coords &B)
		{
		float modA = A.norm();
		float modB = B.norm();
		assert(modA > 0);
		assert(modB > 0);
		float dotAB = dot(A, B);
		float theta = acos(dotAB/(modA*modB));
		assert(theta >= 0 and theta <= PI_plus_epsilon);
		return theta;
		}

// Angle between vectors in range 0 to 180
	static float angle_degrees(const coords &A, const coords &B)
		{
		return degrees(angle_radians(A, B));
		}

	static void assert_unit_basis(const coords &X,
						   const coords &Y,
						   const coords &Z);

	static void rotate(const coords &v, const coords &axis, float theta,
				coords &v_rotated);
	};
