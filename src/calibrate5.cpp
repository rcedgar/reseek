#include "myutils.h"
#include "chainreader2.h"
#include "dbsearcher.h"

/***
* Calibrate DB following similar method to MASM in src/scop40s.
* Specify reversed query on cmdline.
* This is first step to report TS distribution for each query.
***/

class C5_Searcher : public DBSearcher
	{
public:
	vector<vector<float> > m_TSVec;
	vector<vector<float> > m_DPVec;
	map<string, uint> m_DBLabelToChainIdx;

public:
	virtual void OnSetup()
		{
		const uint N = GetDBChainCount();
		m_DPVec.resize(N);
		m_TSVec.resize(N);
		for (uint i = 0; i < N; ++i)
			{
			const string &DBLabel = m_DBChains[i]->m_Label;
			m_DBLabelToChainIdx[DBLabel] = i;
			}
		}

	virtual void OnAln(DSSAligner &DA, bool Up)
		{
		const string &LabelB = DA.m_ChainB->m_Label;
		map<string, uint>::const_iterator iter =
		  m_DBLabelToChainIdx.find(LabelB);
		asserta(iter != m_DBLabelToChainIdx.end());
		uint ChainIndexA = iter->second;
		m_DPVec[ChainIndexA].push_back(DA.m_AlnFwdScore);
		m_TSVec[ChainIndexA].push_back(DA.m_TestStatisticB);
		}

	void WriteVecs(const string &FN1, const string &FN2) const
		{
		FILE *f1 = CreateStdioFile(FN1);
		FILE *f2 = CreateStdioFile(FN2);
		const uint N = GetDBChainCount();
		asserta(SIZE(m_DPVec) == N);
		asserta(SIZE(m_TSVec) == N);
		float MaxScore = 0;
		string MaxLabel;
		for (uint i = 0; i < N; ++i)
			{
			const string &DBLabel = m_DBChains[i]->m_Label;
			fprintf(f1, "%s", DBLabel.c_str());
			const vector<float> &v1 = m_TSVec[i];
			for (uint j = 0; j < SIZE(v1); ++j)
				{
				float Score = v1[j];
				fprintf(f1, "\t%.3g", Score);
				if (Score > MaxScore)
					{
					MaxScore = Score;
					MaxLabel = DBLabel;
					}
				}
			fprintf(f1, "\n");

			fprintf(f2, "%s", DBLabel.c_str());
			const vector<float> &v2 = m_DPVec[i];
			for (uint j = 0; j < SIZE(v2); ++j)
				{
				float Score = v2[j];
				fprintf(f2, "\t%.3g", Score);
				if (Score > MaxScore)
					{
					MaxScore = Score;
					MaxLabel = DBLabel;
					}
				}
			fprintf(f2, "\n");
			}
		CloseStdioFile(f1);
		CloseStdioFile(f2);
		ProgressLog("Max score %.3g >%s\n", 
		  MaxScore, MaxLabel.c_str());
		}
	};

void cmd_calibrate5()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt_db;

	opt_verysensitive = true;
	optset_verysensitive = true;

	C5_Searcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast, SCOP40_DBSIZE);
	DBS.m_Params = &Params;

	DBS.LoadDB(DBFN);

	Params.m_DBSize = (float) DBS.GetDBSize();
	DBS.Setup();

	ChainReader2 QCR;
	QCR.Open(QFN);
	DBS.RunQuery(QCR);
	DBS.WriteVecs(opt_calib_output5a, opt_calib_output5b);
	}
