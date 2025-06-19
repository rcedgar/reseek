#include "myutils.h"
#include "dssaligner.h"
#include "pdbchain.h"

float GetSelfRevScore(DSSAligner &DA, DSS &D, const PDBChain &Chain,
					  const vector<vector<byte> > &Profile,
					  const vector<byte> *ptrMuLetters,
					  const vector<uint> *ptrMuKmers);

static void TraceAln1(const DSSParams &Params,
					  const PDBChain &Q, const PDBChain &T)
	{
	Log("\n______________________________________________\n");
	Log("Q>%s(%u)\n", Q.m_Label.c_str(), Q.GetSeqLength());
	Log("T>%s(%u)\n", T.m_Label.c_str(), T.GetSeqLength());

	DSS D;
	DSSAligner DA;
	D.SetParams(Params);
	DA.SetParams(Params);

	vector<vector<byte> > ProfileQ;
	vector<vector<byte> > ProfileT;

	vector<byte> MuLettersQ;
	vector<byte> MuLettersT;

	vector<uint> MuKmersQ;
	vector<uint> MuKmersT;

	D.Init(Q);
	D.GetProfile(ProfileQ);
	D.GetMuLetters(MuLettersQ);
	D.GetMuKmers(MuLettersQ, MuKmersQ, Params.m_MKFPatternStr);

	D.Init(T);
	D.GetProfile(ProfileT);
	D.GetMuLetters(MuLettersT);
	D.GetMuKmers(MuLettersT, MuKmersT, Params.m_MKFPatternStr);

	float SelfRevScoreQ = GetSelfRevScore(DA, D, Q, ProfileQ, &MuLettersQ, &MuKmersQ);
	float SelfRevScoreT = GetSelfRevScore(DA, D, T, ProfileT, &MuLettersT, &MuKmersT);
	Log("SelfRevScoreQ=%.1f\n", SelfRevScoreQ);
	Log("SelfRevScoreT=%.1f\n", SelfRevScoreT);

	DA.SetQuery(Q, &ProfileQ, &MuLettersQ, &MuKmersQ, SelfRevScoreQ);
	DA.SetTarget(T, &ProfileT, &MuLettersT, &MuKmersT, SelfRevScoreT);
	DA.AlignQueryTarget();
	bool DoMKF = DA.DoMKF();
	Log("Path=(%u)%.10s...\n", SIZE(DA.m_Path), DA.m_Path.c_str());
	float E = DA.m_EvalueA;
	if (E > 1e5)
		Log("EvalueA=%.3g\n", E);
	else
		Log("EvalueA=%.1f\n", E);
	Log("AlnFwdScore=%.3g\n", DA.m_AlnFwdScore);
	Log("DoMKF=%c\n", tof(DoMKF));
	if (DoMKF)
		{
		Log("m_MKF.BestChainScore=%d\n", DA.m_MKF.m_BestChainScore);
		Log("m_XDropScore=%.1f\n", DA.m_XDropScore);
		}
	Log("Omega=%.1f\n", Params.m_Omega);
	Log("DoMuFilter=%c\n", tof(Params.m_Omega>0));
	Log("MuFilterOk=%c\n", tof(DA.MuFilter()));
	}

void cmd_tracealn()
	{
	asserta(optset_db);
	vector<PDBChain *> Qs, Ts;
	ReadChains(g_Arg1, Qs);
	ReadChains(opt(db), Ts);

	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast);

	const uint NQ = SIZE(Qs);
	const uint NT = SIZE(Ts);
	for (uint iq = 0; iq < NQ; ++iq)
		{
		const PDBChain &Q = *Qs[iq];
		for (uint it = 0; it < NT; ++it)
			{
			const PDBChain &T = *Ts[it];
			TraceAln1(Params, Q, T);
			}
		}
	}
