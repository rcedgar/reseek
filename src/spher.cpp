#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"
#include "spher.h"

/***
Spherical coordinates
https://en.wikipedia.org/wiki/Spherical_coordinate_system

(r, theta, phi)
	theta is polar angle to z axis 0 .. 180 degrees (latitude)
	phi is azimuthal angle to x axis 0 .. 360 degrees (longitude)

	r		= sqrt(x^2 + y^2 + z^2)
	theta	= arccos(z/r)
	phi		= sgn(y) arccos(x / sqrt(x^2 + y^2))

	x = r sin(theta) cos(phi)
	y = r sin(theta) sin(phi)
	z = r cos(theta)
***/

void CartToSpher(double x, double y, double z,
  double &r, double &theta, double &phi)
	{
	r = sqrt(x*x + y*y + z*z);
	if (r < 0.1)
		{
		theta = 0;
		phi = 0;
		return;
		}

	theta = acos(z/r);
	assert(theta >= 0 && theta < PI);
	double sign_y = (y >= 0 ? 1 : -1);
	if (fabs(x) < 0.1 && fabs(y) < 0.1)
		phi = 0;
	else
		{
		double d = sqrt(x*x + y*y);
		phi = sign_y*acos(x/d);
		}
	}

void CartToSpher(const vector<double> &Pt,
  double &r, double &theta, double &phi)
	{
	CartToSpher(Pt[X], Pt[Y], Pt[Z], r, theta, phi);
	}

void SpherToCart(double r, double theta, double phi,
  double &x, double &y, double &z)
	{
	x = r*sin(theta)*cos(phi);
	y = r*sin(theta)*sin(phi);
	z = r*cos(theta);
	}

// Number of bins on the axis
static const uint g_Lat_N = 8;
static const uint g_Long_N = 16;
static double MaxDist = 20;

static void TestCoords(double x, double y, double z)
	{
	double theta, r, phi;
	CartToSpher(x, y, z, r, theta, phi);

	double theta_deg = degrees_0_to_360(theta);
	double phi_deg = degrees_0_to_360(phi);

	double x2, y2, z2;
	SpherToCart(r, theta, phi, x2, y2, z2);

	Log("\n");
	Log("  x  = %8.3f, y  = %.1f, z  = %.1f\n", x, y, z);
	Log("  x2 = %8.3f, y2 = %.1f, z2 = %.1f\n", x2, y2, z2);
	Log("  r = %8.3f, theta = %.1f, phi = %.1f\n", r, theta_deg, phi_deg);

	asserta(feq(x2, x));
	asserta(feq(y2, y));
	asserta(feq(z2, z));
	}

#if 0
static void Old_FindCavity(const PDBChain &Chain)
	{
	const uint L = Chain.GetSeqLength();

	vector<vector<double> > MinDistMx;
	vector<vector<uint> > MinPosMx;

	MinDistMx.resize(g_Lat_N);
	MinPosMx.resize(g_Lat_N);
	for (uint Lat = 0; Lat < g_Lat_N; ++Lat)
		{
		MinDistMx[Lat].resize(g_Long_N, DBL_MAX);
		MinPosMx[Lat].resize(g_Long_N, UINT_MAX);
		}

	for (uint Pos = 0; Pos < L; ++Pos)
		{
		vector<double> Pt;
		Chain.GetPt(Pos, Pt);
#if 0
		TestCoords(Pt[X], Pt[Y], Pt[Z]);
#endif
		uint Lat, Long;
		double r;
		GetLatLongTick(Pt[X], Pt[Y], Pt[Z], g_Lat_N, g_Long_N, Lat, Long, r);
		if (r > MaxDist)
			continue;

		char c = Chain.m_Seq[Pos];
#if 0
		Log("[%4u]", Pos);
		Log("  %c", c);
		Log("  %8.3f", Pt[X]);
		Log("  %8.3f", Pt[Y]);
		Log("  %8.3f", Pt[Z]);
		Log("  %3u", Lat);
		Log("  %3u", Long);
		Log("  %9.3f", r);
		Log("\n");
#endif
		if (r < MinDistMx[Lat][Long])
			{
			MinDistMx[Lat][Long] = r;
			MinPosMx[Lat][Long] = Pos;
			}
		}

	Log("\n\n");
	Log("Max Lat %u, Lon %u\n", g_MaxLat, g_MaxLong);

	uint Count = 0;
	for (uint Lat = 0; Lat < g_Lat_N; ++Lat)
		{
		for (uint Long = 0; Long < g_Long_N; ++Long)
			{
			uint Pos = MinPosMx[Lat][Long];
			if (Pos == UINT_MAX)
				continue;

			++Count;
			double r = MinDistMx[Lat][Long];
			vector<double> Pt;
			Chain.GetPt(Pos, Pt);

			char c = Chain.m_Seq[Pos];
			Log("[%4u]", Pos);
			Log("  %c", c);
			Log("  %8.3f", Pt[X]);
			Log("  %8.3f", Pt[Y]);
			Log("  %8.3f", Pt[Z]);
			Log("  %3u", Lat);
			Log("  %3u", Long);
			Log("  %9.3f", r);
			Log("\n");
			}
		}
	Log("\n");
	Log("%u points found\n", Count);

	for (uint Lat = 0; Lat < g_Lat_N; ++Lat)
		{
		for (uint Long = 0; Long < g_Long_N; ++Long)
			{
			uint Pos = MinPosMx[Lat][Long];
			if (Pos == UINT_MAX)
				continue;
			Log("select resi %u\n", Pos+1);
			Log("color magenta, sele\n");
			Log("deselect\n");
			}
		}
	}
#endif // 0
