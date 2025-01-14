#include "myutils.h"
#include "dss.h"
#include "seqdb.h"
#include "museqsource.h"

void MuFilter(const DSSParams &Params,
			  SeqDB &QueryDB,
			  SeqSource &DBSS,
			  const string &OutputFileName);

void PostMuFilter(const DSSParams &Params,
				  const string &MuFilterTsvFN,
				  const string &QueryCAFN,
				  const string &DBBCAFN,
				  float MaxEvalue,
				  const string &HitsFN);

void cmd_bigsearch()
	{
	asserta(optset_db);
	const string &QueryFN = g_Arg1;
	const string &DBFN = string(opt_db);

	if (!EndsWith(DBFN, ".bca"))
		Die(".bca format required for -db");

	DSSParams Params;
	Params.SetFromCmdLine(10000);
	const string &PatternStr = Params.m_PatternStr;
	asserta(PatternStr == "111");

	string MuFilterTsvFN;
	GetTmpFileName(MuFilterTsvFN);
	Log("MuFilterTsvFN=%s\n", MuFilterTsvFN.c_str());

	MuSeqSource QSS;
	MuSeqSource DBSS;
	QSS.Open(QueryFN, Params);
	if (optset_input2)
		{
		Progress("open %s\n", opt_input2);
		DBSS.OpenFasta(opt_input2);
		}
	else
		DBSS.Open(DBFN, Params);

	SeqDB MuQueryDB;
	MuQueryDB.FromSS(QSS);

	float MaxEvalue = 10;
	if (optset_evalue)
		MaxEvalue = (float) opt_evalue;

	MuFilter(Params, MuQueryDB, DBSS, MuFilterTsvFN);
	PostMuFilter(Params, MuFilterTsvFN, QueryFN, DBFN, MaxEvalue, opt_output);

	if (!opt_keeptmp)
		DeleteStdioFile(MuFilterTsvFN);
	}
