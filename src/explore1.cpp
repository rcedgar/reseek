#include "myutils.h"
#include "sweeper.h"
#include "sort.h"

const uint ITERS1 = 999999;
static uint TRIES = 8;

static float GetDelta(const string &Name)
	{
	if (Name == "GapOpen")
		return 1.1f;
	if (Name == "DALIw")
		return 1.05f;
	return 1.2f;
	}

static float GetZ(const string &Name)
	{
	if (Name == "Bias")
		return 0.02f;
	return 0.001f;
	}

bool Sweeper::Explore1(DSSParams &Params,
  uint Idx, float Delta, float Z, uint MaxTries)
	{
	const string &ParamName = m_ParamNames[Idx];
	bool AnyBetter = false;
	DSSParams TryParams = Params;
	for (uint Try = 0; Try < TRIES; ++Try)
		{
		float OldValue = Params.GetParam(ParamName);
		float NewValue = OldValue*Delta + Z;
		TryParams.SetParam(ParamName, NewValue, false);
		TryParams.NormalizeWeights();
		TryParams.ApplyWeights();

		ProgressLog("  %s %.4g => %.4g\n", ParamName.c_str(), OldValue, NewValue);

		uint SavedBestScore = m_BestScore;
		uint Score = Run(TryParams);
		if (Score <= SavedBestScore)
			return AnyBetter;

		Params = TryParams;
		AnyBetter = true;
		}
	return AnyBetter;
	}

void cmd_explore1()
	{
	const string &CalFN = g_Arg1;
	SCOP40Bench SB;
	asserta(optset_benchlevel);
	DSSParams Params;
	SB.ReadChains(CalFN, "");

	Params.SetFromCmdLine();
	Params.m_DBSize = (float) SB.m_ChainCount;

	Sweeper S;
	if (optset_fev)
		{
		S.m_fFev = CreateStdioFile(opt_fev);
		setbuf(S.m_fFev, 0);
		}

	asserta(optset_params);
	Split(opt_params, S.m_ParamNames, '_');

	const uint N = SIZE(S.m_ParamNames);
	vector<float> Values;
	vector<float> ParamDeltas;
	vector<float> ParamZs;
	for (uint i = 0; i < N; ++i)
		{
		const string &Name = S.m_ParamNames[i];
		float Value = Params.GetParam(Name);
		float Delta = GetDelta(Name);
		float Z = GetZ(Name);
		Values.push_back(Value);
		ParamDeltas.push_back(Delta);
		ParamZs.push_back(Z);
		}
	vector<float> FirstValues = Values;

	S.Run(Params);
	for (uint Loop = 0;; ++Loop)
		{
		uint Improvements = 0;
		for (uint Idx = 0; Idx < N; ++Idx)
			{
			const string &Name = S.m_ParamNames[Idx];
			float Delta = ParamDeltas[Idx];
			float Z = ParamZs[Idx];
			ProgressLog("\n=== [%u] Idx %u/%u %s (delta %.3g, Z %.3g) === %u improves\n",
			  Loop+1, Idx+1, N, Name.c_str(), Delta, Z, Improvements);
			bool BetterUp = S.Explore1(Params, Idx, Delta, Z, TRIES);
			if (BetterUp)
				{
				++Improvements;
				continue;
				}
			bool BetterDn = S.Explore1(Params, Idx, 1.0f/Delta, -Z, TRIES);
			if (BetterDn)
				++Improvements;
			}
		if (Improvements < 2)
			break;

		for (uint i = 0; i < N; ++i)
			{
			float Delta = ParamDeltas[i];
			float Z = ParamZs[i];
			ParamDeltas[i] = sqrtf(Delta);
			ParamZs[i] = Z*0.8f;
			}
		}
	CloseStdioFile(S.m_fFev);
	ProgressLog("\n");
	uint First = S.m_FirstScore;
	uint Best = S.m_BestScore;
	int d = int(Best) - int(First);
	double Pct = GetPct(fabs(d), First);
	ProgressLog("Score %u => %u (%+.2f%%)\n", First, Best, Pct);
	for (uint i = 0; i < N; ++i)
		{
		const string &Name = S.m_ParamNames[i];
		float FirstValue = FirstValues[i];
		float Value = Params.GetParam(Name);
		if (Value == FirstValue)
			ProgressLog("%s :: %.6g (unchanged)\n", Name.c_str(), Value);
		else
			ProgressLog("%s %.6g => %.6g\n", Name.c_str(), FirstValue, Value);
		}
	}
