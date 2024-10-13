#include "myutils.h"
#include "chainreader2.h"
#include "dbsearcher.h"
#include "quarts.h"

// Measure binned TS distribution of database vs. reversed database.

class C5_Searcher : public DBSearcher
	{
public:
	vector<float> m_TSVec;

public:
	virtual void OnSetup()
		{
		const uint N = GetDBChainCount();
		m_TSVec.resize(N);
		}

	virtual void OnAln(DSSAligner &DA, bool Up)
		{
		asserta(Up);
		m_TSVec.push_back(DA.m_TestStatisticA);
		}

	void Write(const string &FN) const
		{
		FILE *f = CreateStdioFile(FN);
		QuartsFloat QF;
		GetQuartsFloat(m_TSVec, QF);
		QF.LogMe();
		QF.WriteMe(f);
		QF.WriteMe(stderr);
		CloseStdioFile(f);
		}
	};

void cmd_calibrate7()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt_db;

	opt_verysensitive = true;
	optset_verysensitive = true;

	C5_Searcher DBS;
	DSSParams Params;
	Params.SetFromCmdLine(10000);
	DBS.m_Params = &Params;

	DBS.LoadDB(DBFN);

	Params.m_DBSize = (float) DBS.GetDBSize();
	DBS.Setup();

	ChainReader2 QCR;
	QCR.Open(QFN);
	DBS.RunQuery(QCR);
	DBS.Write(opt_output);
	}
