#include "myutils.h"
#include "scop40bench.h"

//  m_Evalue_Gumbel_mu = 2.5;
//  m_Evalue_Gumbel_beta = 0.613f;
static float AdjustTS(uint Method, float TS,
  float mu1, float beta1, float mu2, float beta2)
	{
	switch (Method)
		{
	case 0:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 1:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS - dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 2:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + 0.5f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 3:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS - 0.5f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 4:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS - 0.25f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 5:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + 0.25f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 6:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu2;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 7:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu2;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS - dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 8:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu2;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + 0.5f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 9:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu2;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS - 0.5f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 10:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu2;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS - 0.25f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 11:
		{
		if (TS <= 0.001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu2;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + 0.25f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 12:
		{
		if (TS <= 0.0001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu2;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + 0.1f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 13:
		{
		if (TS <= 0.0001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS + 0.1f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	case 14:
		{
		if (TS <= 0.0001)
			return TS;
		const float gmu = 2.5f;
		const float gbeta = 0.613f;
		float dmu = gmu - mu1;
		float logTS = logf(TS);
		float Adjusted_logTS = logTS - 0.1f*dmu;
		float Adjusted_TS = expf(Adjusted_logTS);
		return Adjusted_TS;
		}

	default:
		asserta(false);
		}
	return FLT_MAX;
	}

void cmd_adjust()
	{
	asserta(optset_input);
	asserta(optset_input2);

	SCOP40Bench SB;
	SB.ReadBit(g_Arg1);
	SB.LoadDB(opt_input);
	SB.BuildDomSFIndexesFromDBChainLabels();
	const uint DomCount = SB.GetDBChainCount();

	vector<float> mus(DomCount, FLT_MAX);
	vector<float> betas(DomCount, FLT_MAX);
	{
	string Line;
	vector<string> Fields;
	FILE *f = OpenStdioFile(opt_input2);
	while (ReadLineStdioFile(f, Line))
		{
		Split(Line, Fields, '\t');
		const string &Label = Fields[0];
		string Dom;
		SCOP40Bench::GetDomFromLabel(Label, Dom);
		map<string, uint>::iterator iter = SB.m_DomToIdx.find(Dom);
		asserta(iter != SB.m_DomToIdx.end());
		uint DomIdx = iter->second;
		asserta(DomIdx < DomCount);
		float mu = float(StrToFloat(Fields[1]));
		float beta = float(StrToFloat(Fields[2]));
		mus[DomIdx] = mu;
		betas[DomIdx] = beta;
		}
	}

	DSSParams Params;
	Params.SetFromCmdLine(DomCount);

	ProgressLog("Unchanged TSs\n");
	SB.m_ScoresAreEvalues = false;
	SB.m_Level = "sf";
	SB.SetStats(0.1f);
	SB.WriteSummary();

	vector<float> TSs = SB.m_Scores;
	vector<float> Es;
	const uint HitCount = SIZE(SB.m_DomIdx1s);
	for (uint i = 0; i < HitCount; ++i)
		{
		float TS = TSs[i];
		float E = Params.GetEvalue(TS);
		Es.push_back(E);
		}

	ProgressLog("\nDefault E-values\n");
	SB.m_Scores = Es;
	SB.m_ScoresAreEvalues = true;
	SB.SetStats(0.1f);
	SB.WriteSummary();

	for (uint Method = 13; Method <= 14; ++Method)
		{
		vector<float> AdjustedTSs;
		for (uint i = 0; i < HitCount; ++i)
			{
			uint Dom1 = SB.m_DomIdx1s[i];
			uint Dom2 = SB.m_DomIdx2s[i];
			float TS = TSs[i];

			float mu1 = mus[Dom1];
			float mu2 = mus[Dom2];

			float beta1 = betas[Dom1];
			float beta2 = betas[Dom2];

			float AdjustedTS = AdjustTS(Method, TS, mu1, beta1, mu2, beta2);
			AdjustedTSs.push_back(AdjustedTS);
			}

		ProgressLog("\nAdjusted TSs %u\n", Method);
		SB.m_Scores = AdjustedTSs;
		SB.m_ScoresAreEvalues = false;
		SB.SetStats(0.1f);
		SB.WriteSummary();
		}
	}
