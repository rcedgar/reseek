#include "myutils.h"
#include "xyz.h"

void LogVec(const string &Msg, const vector<float> &v)
	{
	Log("%s(%.2f, %.2f, %.2f)\n", Msg.c_str(), v[X], v[Y], v[Z]);
	}

void LogMx(const string &Msg, const vector<vector<float> > &Mx)
	{
	Log("\n");
	Log("%-10.10s    x         y         z\n", Msg.c_str());
	for (int Motif = 0; Motif < 3; ++Motif)
		{
		Log("[%c]  ", "ABC"[Motif]);
		for (int Axis = 0; Axis < 3; ++Axis)
			{
			float v = Mx[Motif][Axis];
			Log("  %8.2f", v);
			}
		Log("\n");
		}

	int A = 0;
	int B = 1;
	int C = 2;
	float AB = GetDist3D(
	  Mx[A][X], Mx[A][Y], Mx[A][Z], 
	  Mx[B][X], Mx[B][Y], Mx[B][Z]);

	float AC = GetDist3D(
	  Mx[A][X], Mx[A][Y], Mx[A][Z], 
	  Mx[C][X], Mx[C][Y], Mx[C][Z]);

	float BC = GetDist3D(
	  Mx[B][X], Mx[B][Y], Mx[B][Z], 
	  Mx[C][X], Mx[C][Y], Mx[C][Z]);
	Log("AB=%.2f, BC=%.2f, AC=%.2f\n", AB, BC, AC);
	}

void GetIdentityMx(vector<vector<float> > &Mx)
	{
	Resize3x3(Mx);

	Mx[0][0] = 1;
	Mx[0][1] = 0;
	Mx[0][2] = 0;

	Mx[1][0] = 0;
	Mx[1][1] = 1;
	Mx[1][2] = 0;

	Mx[2][0] = 0;
	Mx[2][1] = 0;
	Mx[2][2] = 1;
	}

void CrossProduct(
  const vector<float> &a,
  const vector<float> &b,
  vector<float> &Prod)
	{
	Resize3(Prod);
	Prod[X] = a[Y]*b[Z] - a[Z]*b[Y];
	Prod[Y] = a[Z]*b[X] - a[X]*b[Z];
	Prod[Z] = a[X]*b[Y] - a[Y]*b[X];
	}

void XFormMx(
  const vector<vector<float> > &Mx,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<vector<float> > &XMx)
	{
	Resize3x3(XMx);
	for (uint i = 0; i < 3; ++i)
		XFormPt(Mx[i], t, R, XMx[i]);
	}

void XFormPt(
  const vector<float> &Pt,
  const vector<float> &t,
  const vector<vector<float> > &R,
  vector<float> &XPt)
	{
	Resize3(XPt);

	vector<float> Tmp;
	Resize3(Tmp);
	for (uint i = 0; i < 3; ++i)
		Tmp[i] = Pt[i] + t[i];

	for (uint i = 0; i < 3; ++i)
		XPt[i] = R[0][i]*Tmp[0] + R[1][i]*Tmp[1] + R[2][i]*Tmp[2];

	}

void RotateMx(const vector<vector<float> > &Mx,
  uint Axis, float Theta, vector<vector<float> > &RotatedMx)
	{
	RotatedMx.clear();
	RotatedMx.resize(3);

	uint OtherAxis_i = GetOtherAxis_i(Axis);
	uint OtherAxis_j = GetOtherAxis_j(Axis);

	const float cos_Theta = cos(Theta);
	const float sin_Theta = sin(Theta);

	for (uint k = 0; k < 3; ++k)
		{
		RotatedMx[k].resize(3);

		RotatedMx[k][Axis] = Mx[k][Axis];

		float Coord_i = Mx[k][OtherAxis_i];
		float Coord_j = Mx[k][OtherAxis_j];

		RotatedMx[k][OtherAxis_i] = Coord_i*cos_Theta - Coord_j*sin_Theta;
		RotatedMx[k][OtherAxis_j] = Coord_i*sin_Theta + Coord_j*cos_Theta;
		}

	AssertSameLengths(Mx, RotatedMx);
	AssertSameAngles(Mx, RotatedMx);
	}

float dot(const float a[3], const float b[3])
	{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
	}

void transform(const float t[3], const float u[3][3],
  const float x[3], float x_transformed[3])
	{
	x_transformed[0] = t[0] + dot(u[0], x);
	x_transformed[1] = t[1] + dot(u[1], x);
	x_transformed[2] = t[2] + dot(u[2], x);
	}
