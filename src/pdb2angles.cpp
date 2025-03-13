#include "myutils.h"
#include "coords.h"
#include "chainreader2.h"

// notebooks/2025_01_30_fake_backbone_chain_angles.docx
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

	coords::cross(xaxis, zaxis, yaxis);
	yaxis.normalize();

	coords::assert_unit_basis(xaxis, yaxis, zaxis);

	theta_radians = coords::angle_radians(bc, cd);
	coords vert_bcd;
	coords::cross(bc, cd, vert_bcd);
	phi_radians = coords::angle_radians(zaxis, vert_bcd);
	float check_angle_radians = coords::angle_radians(zaxis, cd);
	if (check_angle_radians < PI/2)
		phi_radians = -phi_radians;
	}

// notebooks/2025_01_30_fake_backbone_chain_angles.docx
void calculate_D(const coords &A, const coords &B, const coords &C,
				 float theta_radians, float phi_radians,
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

static float err(float a)
	{
	if (fabs(a) < 0.1)
		return 0;
	return a;
	}

void GetAngles(const PDBChain &Chain, uint Pos,
					float &theta_deg, float &phi_deg)
	{
	theta_deg = FLT_MAX;
	phi_deg = FLT_MAX;
	if(Pos < 3)
		return;

	const vector<float> &Xs = Chain.m_Xs;
	const vector<float> &Ys = Chain.m_Ys;
	const vector<float> &Zs = Chain.m_Zs;
	const uint L = Chain.GetSeqLength();
	coords A, B, C, D;
	Chain.GetCoords(Pos-3, A);
	Chain.GetCoords(Pos-2, B);
	Chain.GetCoords(Pos-1, C);
	Chain.GetCoords(Pos-0, D);
	float X = D.x;
	float Y = D.y;
	float Z = D.z;

	float dCD = coords::dist(C, D);

	float theta_rad, phi_rad;
	calculate_theta_phi(A, B, C, D, theta_rad, phi_rad);

	theta_deg = degrees(theta_rad);
	phi_deg = degrees(phi_rad);
	}

static void ChainAngles(FILE *f, const PDBChain &Chain)
	{
	if (f == 0)
		return;
	const vector<float> &Xs = Chain.m_Xs;
	const vector<float> &Ys = Chain.m_Ys;
	const vector<float> &Zs = Chain.m_Zs;
	const uint L = Chain.GetSeqLength();
	const string &Seq = Chain.m_Seq;
	for (uint i = 3; i < L; ++i)
		{
		coords A, B, C, D;
		Chain.GetCoords(i-3, A);
		Chain.GetCoords(i-2, B);
		Chain.GetCoords(i-1, C);
		Chain.GetCoords(i-0, D);
		float X = D.x;
		float Y = D.y;
		float Z = D.z;

		float dCD = coords::dist(C, D);

		float theta_rad, phi_rad;
		calculate_theta_phi(A, B, C, D, theta_rad, phi_rad);

		float theta_deg = degrees(theta_rad);
		float phi_deg = degrees(phi_rad);

		coords D2;
		calculate_D(A, B, C, theta_rad, phi_rad, D2);

		float eX = err(D2.x - X);
		float eY = err(D2.y - Y);
		float eZ = err(D2.z - Z);

		float e = coords::dist(D, D2);

		char c = Seq[i];
		fprintf(f, "%c", c);
		fprintf(f, "\t%.1f", X);
		fprintf(f, "\t%.1f", D2.x);
		fprintf(f, "\t%.1f", eX);

		fprintf(f, "\t%.1f", Y);
		fprintf(f, "\t%.1f", D2.y);
		fprintf(f, "\t%.1f", eY);

		fprintf(f, "\t%.1f", Z);
		fprintf(f, "\t%.1f", D2.z);
		fprintf(f, "\t%.1f", eZ);

		fprintf(f, "\t%.1f", dCD);
		fprintf(f, "\t%.1f", e);

		fprintf(f, "\t%.1f", theta_deg);
		fprintf(f, "\t%.1f", phi_deg);

		fprintf(f, "\n");
		}
	}

void cmd_pdb2angles()
	{
	ChainReader2 CR;
	CR.Open(g_Arg1);
	FILE *fOut = CreateStdioFile(opt_output);

	if (fOut != 0)
		fprintf(fOut, "aa	X	X'	eX	Y	Y'	eY	Z	Z'	eZ	dCD	e	theta	phi\n");
	for (;;)
		{
		PDBChain *Chain = CR.GetNext();
		if (Chain == 0)
			break;
		ChainAngles(fOut, *Chain);
		}
	CloseStdioFile(fOut);
	}
