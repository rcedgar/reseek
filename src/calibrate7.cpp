#include "myutils.h"
#include "chainreader2.h"
#include "dbsearcher.h"
#include "quarts.h"
#include "binner.h"

// Min=-0.127, LoQ=-0.0482, Med=-0.0254, HiQ=-0.00619, Max=0.15, Avg=-0.0145
// Measure binned TS distribution of database vs. reversed database.
class C7_Searcher : public DBSearcher
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
		if (DA.m_NewTestStatisticA > 0)
			m_TSVec.push_back(DA.m_NewTestStatisticA);
		}
	};

void cmd_calibrate7()
	{
	const string &QFN = g_Arg1;
	const string &DBFN = opt(db);

	opt(verysensitive) = true;
	optset_verysensitive = true;

	C7_Searcher DBS;
	DSSParams Params;
	Params.SetDSSParams(DM_DefaultFast);
	DBS.m_Params = &Params;

	DBS.LoadDB(DBFN);

	DBS.Setup();

	ChainReader2 QCR;
	QCR.Open(QFN);
	DBS.RunQuery(QCR);

	FILE *f = CreateStdioFile(opt(output));
	QuartsFloat QF;
	GetQuartsFloat(DBS.m_TSVec, QF);
	QF.LogMe();
	QF.WriteMe(f);
	QF.WriteMe(stderr);

	Binner B(DBS.m_TSVec, 21, 0.0f, 0.1f);
	B.ToTsv(f);
	B.AccumToTsv(f);
	B.AccumToTsvReverse(f);
	CloseStdioFile(f);
	}
