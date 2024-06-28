#include "myutils.h"
#include "pdbchain.h"
#include "abcxyz.h"

void PDBChain::FromCal(const string &FileName)
	{
	vector<string> Lines;
	ReadLinesFromFile(FileName, Lines);
	FromCalLines(Lines);
	}

void PDBChain::FromCalLines(const vector<string> &Lines)
	{
	Clear();

	if (Lines.empty())
		return;

	const string &FirstLine = Lines[0];
	m_Label = FirstLine.substr(1, string::npos);

/***
>102l
M       43.619  -1.924  8.869
N       40.445  -0.876  10.670
I       38.254  2.240   11.220
F       40.340  3.621   14.036
***/
	const uint N = SIZE(Lines);
	vector<string> Fields;
	for (uint LineNr = 1; LineNr < N; ++LineNr)
		{
		const string &Line = Lines[LineNr];
		Split(Line, Fields, '\t');
		if (Fields.size() != 4 || Fields[0].size() != 1)
			Die("Invalid .cal record '%s'", Line.c_str());

		char aa = Fields[0][0];
		float X = StrToFloatf(Fields[1]);
		float Y = StrToFloatf(Fields[2]);
		float Z = StrToFloatf(Fields[3]);

		m_Seq.push_back(aa);
		m_Xs.push_back(X);
		m_Ys.push_back(Y);
		m_Zs.push_back(Z);
		}
	}

