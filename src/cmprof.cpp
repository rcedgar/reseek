#include "myutils.h"
#include "cmprof.h"
#include "alpha.h"

double GetNormal(double Mu, double Sigma, double x)
	{
	static double TWOPI = (2.0*3.1415926535);
	static double FACTOR = 1.0/sqrt(TWOPI);

	double a = (x - Mu)/Sigma;
	double y = (FACTOR/Sigma)*exp(-0.5*a*a);
	return y;
	}

void CMProf::GetDistMx(const PDBChain &Chain, const vector<uint> &PosVec,
  vector<vector<double> > &DistMx)
	{
	DistMx.clear();
	const uint N = SIZE(PosVec);
	DistMx.resize(N);
	for (uint i = 0; i < N; ++i)
		DistMx[i].resize(N, DBL_MAX);

	for (uint i = 0; i < N; ++i)
		{
		uint Posi = PosVec[i];
		DistMx[i][i] = 0;
		if (Posi == UINT_MAX)
			continue;
		for (uint j = 0; j < N; ++j)
			{
			uint Posj = PosVec[j];
			if (Posj == UINT_MAX)
				continue;
			double d = Chain.GetDist(Posi, Posj);
			DistMx[i][j] = d;
			DistMx[j][i] = d;
			}
		}
	}

void CMProf::MxToFile(FILE *f, const string &Name,
  const vector<vector<double> > &Mx) const
	{
	if (f == 0)
		return;

	const uint CoreColCount = GetCoreColCount();
	asserta(SIZE(Mx) == CoreColCount);
	for (uint i = 0; i < CoreColCount; ++i)
		{
		fprintf(f, "%s\t%u", Name.c_str(), i);
		for (uint j = 0; j <= i; ++j)
			{
			double x = Mx[i][j];
			if (i == j)
				asserta(x == 0);
			asserta(Mx[j][i] == Mx[i][j]);
			fprintf(f, "\t%.3g", x);
			}
		fprintf(f, "\n");
		}
	}

void CMProf::ToFile(const string &FileName) const
	{
	if (FileName == "")
		return;
	FILE *f = CreateStdioFile(FileName);
	asserta(f != 0);
	const uint ColCount = GetColCount();
	const uint CoreColCount = GetCoreColCount();
	fprintf(f, "CMProf\t%u\n", CoreColCount);
	for (uint i = 0; i < ColCount; ++i)
		fprintf(f, "%c", m_ColIsCore[i] ? '1' : '0');
	fprintf(f, "\n");
	MxToFile(f, "mean", m_MeanDistMx);
	MxToFile(f, "stddev", m_StdDevs);
	CloseStdioFile(f);
	}

void CMProf::MxFromFile(FILE *f, string &Name, uint CoreColCount,
  vector<vector<double> > &Mx)
	{
	Mx.resize(CoreColCount);
	for (uint i = 0; i < CoreColCount; ++i)
		Mx[i].resize(CoreColCount, DBL_MAX);

	string Line;
	vector<string> Fields;
	for (uint i = 1; i < CoreColCount; ++i)
		{
		bool Ok = ReadLineStdioFile(f, Line);
		if (!Ok)
			Die("Premature EOF in CMProf file");
		Split(Line, Fields, '\t');
		asserta(SIZE(Fields) == i+3);
		if (i == 1)
			Name = Fields[0];
		else
			asserta(Fields[0] == Name);
		asserta(StrToUint(Fields[1]) == i);

		for (uint j = 0; j <= i; ++j)
			{
			double Value = StrToFloat(Fields[j+2]);
			if (i == j)
				asserta(Value == 0);
			Mx[i][j] = Value;
			Mx[j][i] = Value;
			}
		}
	}

void CMProf::FromFile(FILE *f)
	{
	Clear();
	string Line;
	vector<string> Fields;
	bool Ok = ReadLineStdioFile(f, Line);
	if (!Ok)
		Die("Premature EOF in CMProf file");
	Split(Line, Fields, '\t');
	if (SIZE(Fields) != 2 || Fields[0] != "CMProf")
		Die("Invalid CMProf file (hdr)");
	const uint CoreColCount = StrToUint(Fields[1]);

	string Name;
	MxFromFile(f, Name, CoreColCount, m_MeanDistMx);
	asserta(Name == "mean");

	MxFromFile(f, Name, CoreColCount, m_StdDevs);
	asserta(Name == "stddev");
	}

void CMProf::FromFile(const string &FileName)
	{
	asserta(FileName != "");
	FILE *f = OpenStdioFile(FileName);
	FromFile(f);
	CloseStdioFile(f);
	}

bool CMProf::TrainChain(const PDBChain &Q)
	{
	const string &Seq = Q.m_Seq;
	const string Label = Q.m_Label;
	const uint L = SIZE(Seq);
	map<string, uint>::const_iterator p = m_UngappedSeqToIdx.find(Seq);
	if (p == m_UngappedSeqToIdx.end())
		return false;

	uint SeqIdx = p->second;
	const string &Row = m_MSA->GetSeq(SeqIdx);
	const uint ColCount = GetColCount();
	vector<uint> PosVec;
	uint Pos = 0;
	for (uint Col = 0; Col < ColCount; ++Col)
		{
		if (m_ColIsCore[Col])
			{
			char c = Row[Col];
			if (isgap(c))
				PosVec.push_back(UINT_MAX);
			else
				PosVec.push_back(Pos++);
			}
		}
	vector<vector<double> > DistMx;
	GetDistMx(Q, PosVec, DistMx);

	m_DistMxVec.push_back(DistMx);
	return true;
	}

void CMProf::FinalizeTrain()
	{
	const uint CoreColCount = GetCoreColCount();

	m_MeanDistMx.clear();
	m_StdDevs.clear();

	m_MeanDistMx.resize(CoreColCount);
	m_StdDevs.resize(CoreColCount);
	for (uint i = 0; i < CoreColCount; ++i)
		{
		m_MeanDistMx[i].resize(CoreColCount, DBL_MAX);
		m_StdDevs[i].resize(CoreColCount, DBL_MAX);
		}

	for (uint i = 0; i < CoreColCount; ++i)
		{
		m_MeanDistMx[i][i] = 0;
		m_StdDevs[i][i] = 0;
		for (uint j = 0; j < i; ++j)
			{
			double Mean, StdDev;
			GetMeanStdDev(i, j, Mean, StdDev);

			m_MeanDistMx[i][j] = Mean;
			m_MeanDistMx[j][i] = Mean;

			m_StdDevs[i][j] = StdDev;
			m_StdDevs[j][i] = StdDev;
			}
		}
	}

void CMProf::GetMeanStdDev(uint i, uint j,
  double &Mean, double &StdDev) const
	{
	Mean = DBL_MAX;
	StdDev = DBL_MAX;
	const uint N = SIZE(m_DistMxVec);
	asserta(N > 0);
	double Sum = 0;
	uint n = 0;
	for (uint k = 0; k < N; ++k)
		{
		asserta(i < SIZE(m_DistMxVec[k]));
		asserta(j < SIZE(m_DistMxVec[k][i]));
		double d = m_DistMxVec[k][i][j];
		if (d == DBL_MAX)
			continue;
		++n;
		Sum += d;
		}
	if (n == 0)
		return;

	Mean = (n == 0 ? 0 : double(Sum)/n);

	double Sumd2 = 0;
	for (uint k = 0; k < N; ++k)
		{
		double d = m_DistMxVec[k][i][j];
		if (d == DBL_MAX)
			continue;
		double d2 = (d - Mean)*(d - Mean);
		Sumd2 += d2;
		}
	StdDev = sqrt(Sumd2/n);
	}

void CMProf::SetMSA(const SeqDB &MSA)
	{
	m_ColIsCore.clear();
	m_CoreCols.clear();
	m_UngappedSeqToIdx.clear();
	m_MSA = &MSA;
	const uint SeqCount = MSA.GetSeqCount();
	if (SeqCount <= 2)
		Die("MSA must have > 2 sequences");
	const uint ColCount = MSA.GetColCount();

	double MaxGapPct = (optset_maxgappct ? opt_maxgappct/100.0 : 50);
	double MaxGapFract = MaxGapPct/100.0;
	uint MinLetters = uint((1 - MaxGapFract)*SeqCount + 1);
	if (MinLetters < 2)
		MinLetters = 2;
	ProgressLog("Max gap pct %.1f, min %u letters/col\n",
	  MaxGapPct, MinLetters);

	for (uint ColIndex = 0; ColIndex < ColCount; ++ColIndex)
		{
		uint LetterCount = MSA.GetLetterCount(ColIndex);
		bool IsCore = (LetterCount >= MinLetters);
		if (IsCore)
			m_CoreCols.push_back(ColIndex);
		m_ColIsCore.push_back(IsCore);
		}

	for (uint SeqIdx = 0; SeqIdx < SeqCount; ++SeqIdx)
		{
		string UngappedSeq;
		MSA.GetSeq_StripGaps(SeqIdx, UngappedSeq, true);
		const string &Label = MSA.GetLabel(SeqIdx);
		m_UngappedSeqToIdx[UngappedSeq] = SeqIdx;
		}
	}
