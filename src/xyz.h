#pragma once

float dot(const float a[3], const float b[3]);
void transform(const float t[3], const float u[3][3],
  const float x[3], float x_transformed[3]);

static const float PI = 3.1415926535f;
static const int X = 0;
static const int Y = 1;
static const int Z = 2;

extern const uint MotifLVec[3];

void LogVec(const string &Msg, const vector<float> &v);
void LogMx(const string &Msg, const vector<vector<float> > &Mx);
float GetMxDeterminant(const vector<vector<float> > &Mx);
void InvertMx(const vector<vector<float> > &Mx,
  vector<vector<float> > &InvMx);
void MulMxVec(const vector<vector<float> > &Mx,
  const vector<float> &Vec,
  vector<float> &Result);
void MulMx(
  const vector<vector<float> > &A,
  const vector<vector<float> > &B,
  vector<vector<float> > &Prod);
void GetBasisR(const vector<vector<float> > &Basis,
  vector<vector<float> > &R);
void LogMx(const string &Msg, const vector<vector<float> > &Mx);
void RotateMx(const vector<vector<float> > &Mx,
  uint Axis, float Theta, vector<vector<float> > &RotatedMx);
void CrossProduct(
  const vector<float> &a,
  const vector<float> &b,
  vector<float> &Prod);
void GetTriForm(
  vector<vector<float> > &MotifCoords,
  vector<float> &t,
  vector<vector<float> > &R);

void GetTriangleCentroid(const vector<vector<float> > &MotifCoords,
  vector<float> &CentroidCoords);
void GetTriangleBasis(const vector<vector<float> > &MotifCoords,
  vector<float> &CentroidCoords, vector<vector<float> > &Basis);

static inline uint GetOtherAxis_i(uint Axis) { return (Axis+1)%3;  }
static inline uint GetOtherAxis_j(uint Axis) { return (Axis+2)%3;  }

static inline void Resize3(vector<float> &v) { v.resize(3); }

static inline void Resize3x3(vector<vector<float> > &Mx)
	{
	Mx.resize(3);
	for (uint i = 0; i < 3; ++i)
		Mx[i].resize(3);
	}

static inline void MulVecScalar(const vector<float> &v,
  float s, vector<float> &r)
	{
	asserta(SIZE(v) == 3);
	r.resize(3);
	r[0] = v[0]*s;
	r[1] = v[1]*s;
	r[2] = v[2]*s;
	}

static inline float GetDist2(
  const vector<float> &Pt1,
  const vector<float> &Pt2)
	{
	assert(Pt1.size() == 3);
	assert(Pt2.size() == 3);
	float dx = Pt1[X] - Pt2[X];
	float dy = Pt1[Y] - Pt2[Y];
	float dz = Pt1[Z] - Pt2[Z];
	float d2 = dx*dx + dy*dy + dz*dz;
	return d2;
	}

static inline float GetDist(
  const vector<float> &Pt1,
  const vector<float> &Pt2)
	{
	assert(Pt1.size() == 3);
	assert(Pt2.size() == 3);
	float dx = Pt1[X] - Pt2[X];
	float dy = Pt1[Y] - Pt2[Y];
	float dz = Pt1[Z] - Pt2[Z];
	float d2 = dx*dx + dy*dy + dz*dz;
	float d = sqrt(d2);
	return d;
	}

static inline float GetDist3D(
  float x1, float y1, float z1,
  float x2, float y2, float z2)
	{
	float dx = x1 - x2;
	float dy = y1 - y2;
	float dz = z1 - z2;
	float d2 = dx*dx + dy*dy + dz*dz;
	float d = sqrt(d2);
	return d;
	}

static inline float GetDist2_Mxij(const vector<vector<float> > &Mx,
  uint i, uint j)
	{
	float xi = Mx[i][X];
	float yi = Mx[i][Y];
	float zi = Mx[i][Z];

	float xj = Mx[j][X];
	float yj = Mx[j][Y];
	float zj = Mx[j][Z];

	float dx = xi - xj;
	float dy = yi - yj;
	float dz = zi - zj;

	float d2 = dx*dx + dy*dy + dz*dz;
	return d2;
	}

static inline float GetDist_Mxij(const vector<vector<float> > &Mx,
  uint i, uint j)
	{
	float xi = Mx[i][X];
	float yi = Mx[i][Y];
	float zi = Mx[i][Z];

	float xj = Mx[j][X];
	float yj = Mx[j][Y];
	float zj = Mx[j][Z];

	float dx = xi - xj;
	float dy = yi - yj;
	float dz = zi - zj;

	float d2 = dx*dx + dy*dy + dz*dz;
	float d = sqrt(d2);
	return d;
	}

static inline float GetMod_xyz(float x, float y, float z)
	{
	float d2 = x*x + y*y + z*z;
	float Mod = sqrt(d2);
	return Mod;
	}

static inline float GetMod_Vec(const vector<float> &v)
	{
	float x = v[X];
	float y = v[Y];
	float z = v[Z];
	float Mod = GetMod_xyz(x, y, z);
	return Mod;
	}

static inline void NormalizeVec(vector<float> &v)
	{
	float Mod = GetMod_Vec(v);
	assert(Mod > 0);
	v[X] /= Mod;
	v[Y] /= Mod;
	v[Z] /= Mod;
	assert(feq(GetMod_Vec(v), 1));
	}

static inline float GetMod_Mxi(const vector<vector<float> > &Mx,
  uint i)
	{
	float x = Mx[i][X];
	float y = Mx[i][Y];
	float z = Mx[i][Z];
	float Mod = GetMod_xyz(x, y, z);
	return Mod;
	}

static inline float GetTheta3D(
  float xi, float yi, float zi,
  float xj, float yj, float zj)
	{
	float DotProd = xi*xj + yi*yj + zi*zj;
	float Modi = GetMod_xyz(xi, yi, zi);
	float Modj = GetMod_xyz(xj, yj, zj);
	if (fabs(Modi*Modj) < 1e-6)
		return 0;
	float cos_theta = DotProd/(Modi*Modj);
	asserta(cos_theta >= -1.02 && cos_theta <= 1.02);
	if (cos_theta < -1)
		cos_theta = -1;
	else if (cos_theta > 1)
		cos_theta = 1;
	float theta = acos(cos_theta);
	return theta;
	}

static inline float GetTheta3D(
  const vector<float> &vi,
  const vector<float> &vj)
	{
	asserta(SIZE(vi) == 3);
	asserta(SIZE(vj) == 3);
	return GetTheta3D(vi[0], vi[1], vi[2], vj[0], vj[1], vj[2]);
	}

static inline void Sub_Vecs(const vector<float> &vi,
  const vector<float> &vj, vector<float> &Diff)
	{
	Resize3(Diff);
	Diff[X] = vi[X] - vj[X];
	Diff[Y] = vi[Y] - vj[Y];
	Diff[Z] = vi[Z] - vj[Z];
	}

static inline void Add_Vecs(const vector<float> &vi,
  const vector<float> &vj, vector<float> &Sum)
	{
	Resize3(Sum);
	Sum[X] = vi[X] + vj[X];
	Sum[Y] = vi[Y] + vj[Y];
	Sum[Z] = vi[Z] + vj[Z];
	}

static inline float GetTheta_Vecs(const vector<float> &vi,
  const vector<float> &vj)
	{
	float xi = vi[X];
	float yi = vi[Y];
	float zi = vi[Z];

	float xj = vj[X];
	float yj = vj[Y];
	float zj = vj[Z];

	float theta = GetTheta3D(xi, yi, zi, xj, yj, zj);
	return theta;
	}

static inline float GetTheta_Mxij(const vector<vector<float> > &Mx,
  uint i, uint j)
	{
	float theta = GetTheta_Vecs(Mx[i], Mx[j]);
	return theta;
	}

static inline float degrees(float Radians) { return Radians*180.0f/PI; }

static inline float degrees_0_to_360(float Radians)
	{
	float Deg = Radians*180.0f/PI;
	Deg = fmodf(Deg, 360);
	if (Deg < 0)
		Deg += 360;
	assert(Deg >= 0 && Deg < 360);
	return Deg;
	}

void GetIdentityMx(vector<vector<float> > &Mx);

void XFormPt(
  const vector<float> &Pt,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<float> &XPt);

void XFormMx(
  const vector<vector<float> > &Mx,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<vector<float> > &XMx);

void XFormPt(
  const vector<float> &Pt,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<float> &XPt);

void XFormMx(
  const vector<vector<float> > &Mx,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<vector<float> > &XMx);

#if DEBUG
static void AssertCanonicalUnitBasis(const vector<vector<float> > &Basis)
	{
	assert(SIZE(Basis) == 3);
	assert(SIZE(Basis[0]) == 3);
	assert(SIZE(Basis[1]) == 3);
	assert(SIZE(Basis[2]) == 3);

	assert(feq(Basis[0][X], 1));
	assert(feq(Basis[0][Y], 0));
	assert(feq(Basis[0][Z], 0));

	assert(feq(Basis[1][X], 0));
	assert(feq(Basis[1][Y], 1));
	assert(feq(Basis[1][Z], 0));

	assert(feq(Basis[2][X], 0));
	assert(feq(Basis[2][Y], 0));
	assert(feq(Basis[2][Z], 1));
	}
#else
#define AssertCanonicalUnitBasis(x)	0
#endif // DEBUG

static void AssertUnitBasisA(const vector<vector<float> > &Basis)
	{
// check length of basis vectors is 1 and angles are PI/2
	for (uint k = 0; k < 3; ++k)
		{
		float Modk = GetMod_Mxi(Basis, k);
		assert(Modk >= 0.98 && Modk <= 1.02);

		uint Axis_i = GetOtherAxis_i(k);
		uint Axis_j = GetOtherAxis_j(k);

		float theta = GetTheta_Mxij(Basis, Axis_i, Axis_j);
		assert(feq(theta, PI/2));
		}
	}

#if DEBUG
static void AssertUnitBasis(const vector<vector<float> > &Basis)
	{
// check length of basis vectors is 1 and angles are PI/2
	for (uint k = 0; k < 3; ++k)
		{
		float Modk = GetMod_Mxi(Basis, k);
		assert(Modk >= 0.98 && Modk <= 1.02);

		uint Axis_i = GetOtherAxis_i(k);
		uint Axis_j = GetOtherAxis_j(k);

		float theta = GetTheta_Mxij(Basis, Axis_i, Axis_j);
		assert(feq(theta, PI/2));
		}
	}
#else
#define AssertUnitBasis(x)	0
#endif // DEBUG

#if DEBUG
static void AssertSameLengths(const vector<vector<float> > &Mx1,
  const vector<vector<float> > &Mx2)
	{
	for (uint k = 0; k < 3; ++k)
		{
		float Modk1 = GetMod_Mxi(Mx1, k);
		float Modk2 = GetMod_Mxi(Mx2, k);
		assert(feq(Modk1, Modk2));
		}
	}
#else
#define AssertSameLengths(x, y)	0
#endif // DEBUG

#if DEBUG
static void AssertSameAngles(const vector<vector<float> > &Mx1,
  const vector<vector<float> > &Mx2)
	{
	for (uint i = 0; i < 3; ++i)
		{
		for (uint j = i + 1; j < 3; ++j)
			{
			float Theta1 = GetTheta_Mxij(Mx1, i, j);
			float Theta2 = GetTheta_Mxij(Mx2, i, j);
			asserta(feq(Theta1, Theta2));
			}
		}
	}
#else
#define AssertSameAngles(x, y)	0
#endif // DEBUG

#if DEBUG
static void AssertMx3D(const vector<vector<float> > &Mx)
	{
	asserta(SIZE(Mx) == 3);
	for (uint i = 0; i < 3; ++i)
		asserta(SIZE(Mx[i]) == 3);
	}
#else
#define AssertMx3D(x)	0
#endif
