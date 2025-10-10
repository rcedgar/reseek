#include "myutils.h"
#include "scop40bench.h"
#include <set>

/***
==> 3dblastaln <==
d12asa_ d12asa_ 1981
d12asa_ d1nnha_ 428
d12asa_ d1eova2 406

==> 3dblastswaln <==
                    <<< blank line!!  
d12asa_ d12asa_ 1981
d12asa_ d1nnha_ 428

==> cealn <==
d12asa_ d12asa_ 8.03 0.00 1.00
d12asa_ d1nnha_ 6.81 2.40 0.81
d12asa_ d1b8aa2 6.70 2.76 0.75

==> cleswaln <==
d12asa_ d12asa_ 15284
d12asa_ d1b8aa2 5739
d12asa_ d1nnha_ 5620

==> dalialn <==
d12asa_ d12asa_ 58.0 0.0
d12asa_ d1nnha_ 27.0 2.7
d12asa_ d1b8aa2 25.1 2.7

==> foldseekaln <==
      0       1     2     3   4 5   6     7 8     9        10     11
d1a1xa_	d1a1xa_	0.000	106	105	0	1	106	1	106	2.549E-163	1045
d1a1xa_	d1jsga_	0.000	108	105	0	1	106	4	111	9.137E-93	622
d1a1xa_	d3saoa_	0.000	77	76	0	26	102	27	103	2.448E-21	188

==> mmseqsaln <==
d12asa_ d12asa_ 674
d12asa_ d2gz4a1 39
d12asa_ d1b8aa2 30

==> tmaln <==
d12asa_ d12asa_ 1.0
d12asa_ d1b8aa2 0.75576
d12asa_ d1nnha_ 0.751037

==> tmfastaln <==
d12asa_ d12asa_ 1.0000 1.0000
d12asa_ d1b8aa2 0.7558 0.7394
d12asa_ d1nnha_ 0.7510 0.8309
***/

void SCOP40Bench::ReadHits(const string &FN)
	{
	uint QueryFieldNr = 0;
	uint TargetFieldNr = 1;
	uint ScoreFieldNr = 2;

	string Algo;
	GetStemName(FN, Algo);

	FILE *f = OpenStdioFile(FN);
	uint64 FileSize = GetStdioFileSize64(f);
	Progress("Reading hits %s\n", Algo.c_str());
	uint HitCount = 0;
	string Line;
	vector<string> Fields;
	uint BadLineCount = 0;
	set<string> NotFound;
	time_t t0 = time(0);
	while (ReadLineStdioFile(f, Line))
		{
		if (HitCount > 0 && HitCount%1000 == 0)
			{
			time_t t1 = time(0);
			if (t1 > t0)
				{
				uint64 Pos = GetStdioFilePos64(f);
				Progress("Hits %.2f%%  %s\r",
						 GetPct(double(Pos), double(FileSize)),
						 IntToStr(HitCount));
				t0 = t1;
				}
			}
		for (uint i = 0; i < SIZE(Line); ++i)
			if (Line[i] == ' ')
				Line[i] = '\t';
		Split(Line, Fields, '\t');
	// 3dblastswaln has blank line
		uint FieldCount = SIZE(Fields);
		if (FieldCount <= ScoreFieldNr)
			{
			++BadLineCount;
			continue;
			}
		string Label1 = Fields[QueryFieldNr];
		string Label2 = Fields[TargetFieldNr];
		uint DomIdx1 = GetDomIdx(Label1, false);
		uint DomIdx2 = GetDomIdx(Label2, false);
		float Score = (float) StrToFloat(Fields[ScoreFieldNr]);
		m_DomIdx1s.push_back(DomIdx1);
		m_DomIdx2s.push_back(DomIdx2);
		m_Scores.push_back(Score);
		++HitCount;
		}
	uint NotFoundCount = SIZE(NotFound);
	ProgressLog("%u hits, %u bad lines %s, %u unknown domains\n",
	  HitCount, BadLineCount, Algo.c_str(), NotFoundCount);
	for (set<string>::const_iterator iter = NotFound.begin();
	  iter != NotFound.end(); ++iter)
		Log("NOTFOUND %s\n", iter->c_str());
	}
