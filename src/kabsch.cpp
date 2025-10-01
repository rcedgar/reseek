#include "myutils.h"
#include "pdbchain.h"

/***
Based on Kabsch() function in TM-align source code v20220412.
https://zhanggroup.org/TM-align/
https://zhanggroup.org/TM-align/TMalign.cpp
Y Zhang, J Skolnick. Nucl Acids Res 33, 2302-9 (2005)
***/

/**************************************************************************
Implemetation of Kabsch algoritm for finding the best rotation matrix
---------------------------------------------------------------------------
x     - x(i,m) are coordinates of atom m in set x            (input)
y     - y(i,m) are coordinates of atom m in set y            (input)
n     - n is number of atom pairs                            (input)
rms   - sum of w*(ux+t-y)**2 over all atom pairs             (output)
u     - u(i,j) is   rotation  matrix for best superposition  (output)
t     - t(i)   is translation vector for best superposition  (output)
**************************************************************************/
float Kabsch(
  const float * const *x, 
  const float * const *y, int n,
  float t[3], float u[3][3])
	{
	int i, j, m, m1, l, k;
	float e0, rms1, d, h, g;
	float cth, sth, sqrth, p, det, sigma;
	float xc[3], yc[3];
	float a[3][3], b[3][3], r[3][3], e[3], rr[6], ss[6];
	float sqrt3 = 1.73205080756888f, tol = 0.01f;
	int ip[] = { 0, 1, 3, 1, 2, 4, 3, 4, 5 };
	int ip2312[] = { 1, 2, 0, 1 };
	asserta(n > 0);

	int a_failed = 0, b_failed = 0;
	float epsilon = 0.00000001f;

	//initializtation
	float rms = 0;
	rms1 = 0;
	e0 = 0;
	float c1[3], c2[3];
	float s1[3], s2[3];
	float sx[3], sy[3], sz[3];
	for (i = 0; i < 3; i++)
		{
		s1[i] = 0.0;
		s2[i] = 0.0;

		sx[i] = 0.0;
		sy[i] = 0.0;
		sz[i] = 0.0;
		}

	for (i = 0; i < 3; i++)
		{
		xc[i] = 0.0;
		yc[i] = 0.0;
		t[i] = 0.0;
		for (j = 0; j < 3; j++)
			{
			u[i][j] = 0.0;
			r[i][j] = 0.0;
			a[i][j] = 0.0;
			if (i == j)
				{
				u[i][j] = 1.0;
				a[i][j] = 1.0;
				}
			}
		}

	//compute centers for vector sets x, y
	for (i = 0; i < n; i++)
		{
		for (j = 0; j < 3; j++)
			{
			c1[j] = x[i][j];
			c2[j] = y[i][j];

			s1[j] += c1[j];
			s2[j] += c2[j];
			}

		for (j = 0; j < 3; j++)
			{
			sx[j] += c1[0] * c2[j];
			sy[j] += c1[1] * c2[j];
			sz[j] += c1[2] * c2[j];
			}
		}
	for (i = 0; i < 3; i++)
		{
		xc[i] = s1[i] / n;
		yc[i] = s2[i] / n;
		}

	for (int mm = 0; mm < n; mm++)
		for (int nn = 0; nn < 3; nn++)
			e0 += (x[mm][nn] - xc[nn]) * (x[mm][nn] - xc[nn]) +
			(y[mm][nn] - yc[nn]) * (y[mm][nn] - yc[nn]);

	for (j = 0; j < 3; j++)
		{
		r[j][0] = sx[j] - s1[0] * s2[j] / n;
		r[j][1] = sy[j] - s1[1] * s2[j] / n;
		r[j][2] = sz[j] - s1[2] * s2[j] / n;
		}

	//compute determinat of matrix r
	det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1])\
		- r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0])\
		+ r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0]);
	sigma = det;

	//compute tras(r)*r
	m = 0;
	for (j = 0; j < 3; j++)
		{
		for (i = 0; i <= j; i++)
			{
			rr[m] = r[0][i] * r[0][j] + r[1][i] * r[1][j] + r[2][i] * r[2][j];
			m++;
			}
		}

	float spur = (rr[0] + rr[2] + rr[5]) / 3.0f;
	float cof = (((((rr[2] * rr[5] - rr[4] * rr[4]) + rr[0] * rr[5])\
		- rr[3] * rr[3]) + rr[0] * rr[2]) - rr[1] * rr[1]) / 3.0f;
	det = det * det;

	for (i = 0; i < 3; i++) e[i] = spur;

	if (spur > 0)
		{
		d = spur * spur;
		h = d - cof;
		g = (spur * cof - det) / 2.0f - spur * h;

		if (h > 0)
			{
			sqrth = sqrt(h);
			d = h * h * h - g * g;
			if (d < 0.0) d = 0.0;
			d = atan2(sqrt(d), -g) / 3.0f;
			cth = sqrth * cos(d);
			sth = sqrth * sqrt3 * sin(d);
			e[0] = (spur + cth) + cth;
			e[1] = (spur - cth) + sth;
			e[2] = (spur - cth) - sth;

				{//compute a                
				for (l = 0; l < 3; l = l + 2)
					{
					d = e[l];
					ss[0] = (d - rr[2]) * (d - rr[5]) - rr[4] * rr[4];
					ss[1] = (d - rr[5]) * rr[1] + rr[3] * rr[4];
					ss[2] = (d - rr[0]) * (d - rr[5]) - rr[3] * rr[3];
					ss[3] = (d - rr[2]) * rr[3] + rr[1] * rr[4];
					ss[4] = (d - rr[0]) * rr[4] + rr[1] * rr[3];
					ss[5] = (d - rr[0]) * (d - rr[2]) - rr[1] * rr[1];

					if (fabs(ss[0]) <= epsilon) ss[0] = 0.0;
					if (fabs(ss[1]) <= epsilon) ss[1] = 0.0;
					if (fabs(ss[2]) <= epsilon) ss[2] = 0.0;
					if (fabs(ss[3]) <= epsilon) ss[3] = 0.0;
					if (fabs(ss[4]) <= epsilon) ss[4] = 0.0;
					if (fabs(ss[5]) <= epsilon) ss[5] = 0.0;

					if (fabs(ss[0]) >= fabs(ss[2]))
						{
						j = 0;
						if (fabs(ss[0]) < fabs(ss[5])) j = 2;
						}
					else if (fabs(ss[2]) >= fabs(ss[5])) j = 1;
					else j = 2;

					d = 0.0;
					j = 3 * j;
					for (i = 0; i < 3; i++)
						{
						k = ip[i + j];
						a[i][l] = ss[k];
						d = d + ss[k] * ss[k];
						}


					//if( d > 0.0 ) d = 1.0 / sqrt(d);
					if (d > epsilon) d = 1.0f / sqrtf(d);
					else d = 0.0;
					for (i = 0; i < 3; i++) a[i][l] = a[i][l] * d;
					}//for l

				d = a[0][0] * a[0][2] + a[1][0] * a[1][2] + a[2][0] * a[2][2];
				if ((e[0] - e[1]) > (e[1] - e[2]))
					{
					m1 = 2;
					m = 0;
					}
				else
					{
					m1 = 0;
					m = 2;
					}
				p = 0;
				for (i = 0; i < 3; i++)
					{
					a[i][m1] = a[i][m1] - d * a[i][m];
					p = p + a[i][m1] * a[i][m1];
					}
				if (p <= tol)
					{
					p = 1.0;
					for (i = 0; i < 3; i++)
						{
						if (p < fabs(a[i][m])) continue;
						p = fabs(a[i][m]);
						j = i;
						}
					k = ip2312[j];
					l = ip2312[j + 1];
					p = sqrt(a[k][m] * a[k][m] + a[l][m] * a[l][m]);
					if (p > tol)
						{
						a[j][m1] = 0.0;
						a[k][m1] = -a[l][m] / p;
						a[l][m1] = a[k][m] / p;
						}
					else a_failed = 1;
					}//if p<=tol
				else
					{
					p = 1.0f / sqrtf(p);
					for (i = 0; i < 3; i++) a[i][m1] = a[i][m1] * p;
					}//else p<=tol  
				if (a_failed != 1)
					{
					a[0][1] = a[1][2] * a[2][0] - a[1][0] * a[2][2];
					a[1][1] = a[2][2] * a[0][0] - a[2][0] * a[0][2];
					a[2][1] = a[0][2] * a[1][0] - a[0][0] * a[1][2];
					}
				}//if(mode!=0)       
			}//h>0

			{
			//compute b
			for (l = 0; l < 2; l++)
				{
				d = 0.0;
				for (i = 0; i < 3; i++)
					{
					b[i][l] = r[i][0] * a[0][l] +
						r[i][1] * a[1][l] + r[i][2] * a[2][l];
					d = d + b[i][l] * b[i][l];
					}
				//if( d > 0 ) d = 1.0 / sqrt(d);
				if (d > epsilon) d = 1.0f / sqrtf(d);
				else d = 0.0;
				for (i = 0; i < 3; i++) b[i][l] = b[i][l] * d;
				}
			d = b[0][0] * b[0][1] + b[1][0] * b[1][1] + b[2][0] * b[2][1];
			p = 0.0;

			for (i = 0; i < 3; i++)
				{
				b[i][1] = b[i][1] - d * b[i][0];
				p += b[i][1] * b[i][1];
				}

			if (p <= tol)
				{
				p = 1.0;
				for (i = 0; i < 3; i++)
					{
					if (p < fabs(b[i][0])) continue;
					p = fabs(b[i][0]);
					j = i;
					}
				k = ip2312[j];
				l = ip2312[j + 1];
				p = sqrt(b[k][0] * b[k][0] + b[l][0] * b[l][0]);
				if (p > tol)
					{
					b[j][1] = 0.0;
					b[k][1] = -b[l][0] / p;
					b[l][1] = b[k][0] / p;
					}
				else b_failed = 1;
				}//if( p <= tol )
			else
				{
				p = 1.0f / sqrtf(p);
				for (i = 0; i < 3; i++) b[i][1] = b[i][1] * p;
				}
			if (b_failed != 1)
				{
				b[0][2] = b[1][0] * b[2][1] - b[1][1] * b[2][0];
				b[1][2] = b[2][0] * b[0][1] - b[2][1] * b[0][0];
				b[2][2] = b[0][0] * b[1][1] - b[0][1] * b[1][0];
				//compute u
				for (i = 0; i < 3; i++)
					for (j = 0; j < 3; j++)
						u[i][j] = b[i][0] * a[j][0] +
						b[i][1] * a[j][1] + b[i][2] * a[j][2];
				}

			//compute t
			for (i = 0; i < 3; i++)
				t[i] = ((yc[i] - u[i][0] * xc[0]) - u[i][1] * xc[1]) -
				u[i][2] * xc[2];
			}//if(mode!=0 && a_failed!=1)
		}//spur>0
	//compute rms
	for (i = 0; i < 3; i++)
		{
		if (e[i] < 0) e[i] = 0;
		e[i] = sqrt(e[i]);
		}
	d = e[2];
	if (sigma < 0.0) d = -d;
	d = (d + e[1]) + e[0];

	rms = (e0 - d) - d;
	if (rms < 0.0) rms = 0.0;

	return rms;
	}

float Kabsch(const PDBChain &ChainA, const PDBChain &ChainB,
  uint LoA, uint LoB, const string &Path,
  float t[3], float u[3][3])
	{
	uint ColCount = SIZE(Path);
	float **x = myalloc(float *, ColCount);
	float **y = myalloc(float *, ColCount);
	uint M = 0;
	uint PosA = LoA;
	uint PosB = LoB;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		switch (Path[Col])
			{
		case 'M':
			{
			x[M] = myalloc(float, 3);
			y[M] = myalloc(float, 3);
			vector<float> PtA;
			vector<float> PtB;
			ChainA.GetPt(PosA, PtA);
			ChainB.GetPt(PosB, PtB);
			x[M][0] = PtA[0];
			y[M][0] = PtB[0];
			x[M][1] = PtA[1];
			y[M][1] = PtB[1];
			x[M][2] = PtA[2];
			y[M][2] = PtB[2];
			++PosA;
			++PosB;
			++M;
			break;
			}

		case 'D':
			++PosA;
			break;

		case 'I':
			++PosB;
			break;

		default:
			asserta(false);
			}
		}
	float RMS = Kabsch(x, y, M, t, u);
	for (uint i = 0; i < M; ++i)
		{
		myfree(x[i]);
		myfree(y[i]);
		}
	myfree(x);
	myfree(y);
	asserta(M > 0);
	return RMS/M;
	}

#if 0
#include "xyz.h"

static void Test(const float t_in[3],
  const float u_in[3][3], uint n)
	{
	Log("t_in = %.1f %.1f %.1f\n", t_in[0], t_in[1], t_in[2]);
	Log("uin0 = %7.3f %7.3f %7.3f\n", u_in[0][0], u_in[0][1], u_in[0][2]);
	Log("uin1 = %7.3f %7.3f %7.3f\n", u_in[1][0], u_in[1][1], u_in[1][2]);
	Log("uin2 = %7.3f %7.3f %7.3f\n", u_in[2][0], u_in[2][1], u_in[2][2]);

	float **x = myalloc(float *, n);
	float **y = myalloc(float *, n);

	for (uint i = 0; i < n; ++i)
		{
		x[i] = myalloc(float, 3);
		y[i] = myalloc(float, 3);

		x[i][0] = float(randu32()%10);
		x[i][1] = float(randu32()%10);
		x[i][2] = float(randu32()%10);

		transform(t_in, u_in, x[i], y[i]);
		}

	float t[3];
	float u[3][3];

	float rms = Kabsch(x, y, n, t, u);

	Log(" rms = %.2f\n", rms);
	Log("   t = %.1f %.1f %.1f\n", t[0], t[1], t[2]);
	Log("  u0 = %7.3f %7.3f %7.3f\n", u[0][0], u[0][1], u[0][2]);
	Log("  u1 = %7.3f %7.3f %7.3f\n", u[1][0], u[1][1], u[1][2]);
	Log("  u2 = %7.3f %7.3f %7.3f\n", u[2][0], u[2][1], u[2][2]);

	Log(" %8.8s,  %8.8s,  %8.8s ", "x1", "y1", "z1");
	Log("   %8.8s,  %8.8s,  %8.8s ", "x2", "y2", "z2");
	Log("   %8.8s,  %8.8s,  %8.8s ", "Kx", "Ky", "Kz");
	Log("\n");
	for (uint i = 0; i < n; ++i)
		{
		float x_transformed[3];
		transform(t, u, x[i], x_transformed);

		Log("(%8.1f,  %8.1f,  %8.1f)", x[i][0], x[i][1], x[i][2]);
		Log("  (%8.1f,  %8.1f,  %8.1f)", y[i][0], y[i][1], y[i][2]);
		Log("  (%8.1f,  %8.1f,  %8.1f)", x_transformed[0], x_transformed[1], x_transformed[2]);
		Log("\n");
		}
	}

void cmd_test()
	{
	float t[3] = { 1, 2, 3 };
	float u[3][3] =
		{
		{ 1, 0, 0 },
		{ 0, 1, 0 },
		{ 0, 0, 1 }
		};

// https://en.wikipedia.org/wiki/Rotation_matrix

	float theta = 1;
	float c = cos(theta);
	float s = sin(theta);

	u[0][0] = 1;		u[0][1] = 0;		u[0][1] = 0;
	u[1][0] = 0;		u[1][1] = c;		u[1][1] = -s;
	u[1][0] = 0;		u[1][1] = s;		u[1][1] = c;

	Test(t, u, 4);
	}
#endif // 0