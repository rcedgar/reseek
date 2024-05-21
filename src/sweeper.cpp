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

uint Sweeper::Run(const DSSParams &Params)
	{
	m_SB.m_Params = &Params;
	m_SB.Run();
	uint Score = m_SB.GetSens1stFP();
	if (m_FirstScore == UINT_MAX)
		m_FirstScore = Score;
	int d = int(Score) - int(m_BestScore);
	bool Improve = (d > 0);
	if (Improve)
		m_BestScore = Score;
	int d2 = int(m_BestScore) - int(m_FirstScore);

	if (m_fFev != 0)
		{
		fprintf(m_fFev, "%4u\t", Score);
		m_SB.m_Params->ToFev(m_fFev);
		}

	Params.ToFev(g_fLog);
	ProgressLog(" Sens %u [%u %+d]", Score, m_BestScore, d);
	if (Improve)
		ProgressLog(" +++ ");
	ProgressLog(" (%+.2f%%)\n", GetPct(d2, m_FirstScore));

	//m_Scores.push_back(Score);
	return Score;
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

	Sweeper S;
	if (optset_fev)
		{
		S.m_fFev = CreateStdioFile(opt_fev);
		setbuf(S.m_fFev, 0);
		}

	DSSParams Params;
	Params.SetFromCmdLine();
	S.m_SB.Setup(Params);
	S.m_SB.ReadChains(CalFN);

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
