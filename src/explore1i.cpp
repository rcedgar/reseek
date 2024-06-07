#include "myutils.h"
#include "sweeper.h"
#include "sort.h"

const uint ITERS1 = 999999;
static uint TRIES = 8;

bool Sweeper::Explore1i(DSSParams &Params, uint Idx, bool Plus, uint MaxTries)
	{
	const string &ParamName = m_ParamNames[Idx];
	bool AnyBetter = false;
	DSSParams TryParams = Params;
	for (uint Try = 0; Try < TRIES; ++Try)
		{
		int OldValue = Params.GetIntParam(ParamName);
		int NewValue = OldValue + (Plus ? 1 : -1);
		TryParams.SetIntParam(ParamName, NewValue);

		ProgressLog("  %s %d => %d\n", ParamName.c_str(), OldValue, NewValue);

		uint SavedBestScore = m_BestScore;
		uint Score = Run(TryParams);
		if (Score <= SavedBestScore)
			return AnyBetter;

		Params = TryParams;
		AnyBetter = true;
		}
	return AnyBetter;
	}

void cmd_explore1i()
	{
	const string &CalFN = g_Arg1;
	Sweeper S;
	if (optset_fev)
		{
		S.m_fFev = CreateStdioFile(opt_fev);
		setbuf(S.m_fFev, 0);
		}
	SCOP40Bench &SB = S.m_SB;
	asserta(optset_benchlevel);
	DSSParams Params;
	SB.ReadChains(CalFN, "");

	opt_sensitive = true;
	optset_sensitive = true;

	Params.SetFromCmdLine();
	Params.m_DBSize = (float) SB.m_ChainCount;
	Params.m_ParaComboGapOpen = 3;
	Params.m_ParaComboGapExt = 1;
	SB.Setup(Params);

	Params.m_ComboScoreOnly = true;
	SB.m_ScoresAreEvalues = false;

	//asserta(optset_params);
	//Split(opt_params, S.m_ParamNames, '_');
	S.m_ParamNames.push_back("ParaComboGapOpen");
	S.m_ParamNames.push_back("ParaComboGapExt");

	vector<int> Values;
	const uint N = SIZE(S.m_ParamNames);
	for (uint i = 0; i < N; ++i)
		{
		const string &Name = S.m_ParamNames[i];
		int Value = Params.GetIntParam(Name);
		Values.push_back(Value);
		}
	vector<int> FirstValues = Values;

	DSSParams DefaultParams;
	DefaultParams.SetFromCmdLine(true);
	DefaultParams.m_ComboScoreOnly = true;
	uint SFFP = S.Run(DefaultParams);

	ProgressLog("%u\topen=%d\text=%d\t@SFFP@\n",
	  SFFP, DefaultParams.m_ParaComboGapOpen, DefaultParams.m_ParaComboGapExt);

	for (int iOpen = 2; iOpen <= 7; ++iOpen)
		{
		for (int iExt = 0; iExt <= 3; ++iExt)
			{
			Params.m_ParaComboGapOpen = iOpen;
			Params.m_ParaComboGapExt = iExt;
			uint SFFP = S.Run(Params);
			ProgressLog("%u\topen=%d\text=%d\t@SFFP@\n", SFFP, iOpen, iExt);
			}
		}
	return;

	S.Run(Params);
	for (uint Loop = 0; ; ++Loop)
		{
		uint Improvements = 0;
		for (uint Idx = 0; Idx < N; ++Idx)
			{
			const string &Name = S.m_ParamNames[Idx];
			ProgressLog("\n=== [%u] Idx %u/%u %s === %u improves\n",
			  Loop+1, Idx+1, N, Name.c_str(), Improvements);
			bool BetterUp = S.Explore1i(Params, Idx, true, TRIES);
			if (BetterUp)
				{
				++Improvements;
				continue;
				}
			bool BetterDn = S.Explore1i(Params, Idx, false, TRIES);
			if (BetterDn)
				++Improvements;
			}
		if (Improvements < 1)
			break;
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
		int FirstValue = FirstValues[i];
		int Value = Params.GetIntParam(Name);
		if (Value == FirstValue)
			ProgressLog("%s :: %d (unchanged)\n", Name.c_str(), Value);
		else
			ProgressLog("%s %d => %d\n", Name.c_str(), FirstValue, Value);
		}
	}
