#pragma once

#include "alpha.h"

void ByteSeqToHexIntStr(
	uint AlphaSize,
	const vector<byte> &Letters,
	string &Str);


template<class t> void GetHexIntSeq(
	uint AlphaSize, const string &Seq, vector<t> &IntSeq)
	{
	IntSeq.clear();
	const uint SeqL = SIZE(Seq);
	const uint L = (AlphaSize < 35 ? SeqL : SeqL/2);
	IntSeq.reserve(L);
	if (AlphaSize < 36)
		{
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			char c = Seq[Pos];
			byte Letter = g_CharToLetterMu[c];
			asserta(Letter < AlphaSize);
			IntSeq.push_back(Letter);
			}
		}
	else
		{
		char Tmp[3];
		Tmp[2] = 0;
		for (uint Pos = 0; Pos < L; ++Pos)
			{
			Tmp[0] = Seq[2*Pos];
			Tmp[1] = Seq[2*Pos+1];
			char *endptr;
			uint Letter32 = strtoul(Tmp, &endptr, 16);
			byte Letter = byte(Letter32);
			asserta(uint(Letter) == Letter32);
			IntSeq.push_back(Letter);
			}
		}
	}

template<class t> void ReadHexIntSeqs(
	uint AlphaSize,
	const string &FN,
	vector<vector<t> > &IntSeqs,
	vector<string> &Labels,
	map<string, uint> &LabelToSeqIdx)
	{
	IntSeqs.clear();
	Labels.clear();
	LabelToSeqIdx.clear();

	vector<string> Lines;
	ReadLinesFromFile(FN, Lines);
	const uint LineCount = SIZE(Lines);
	string Seq;
	string Label;
	vector<t> IntSeq;
	vector<string> Fields;
	for (uint LineIdx = 0; LineIdx < LineCount; ++LineIdx)
		{
		const string &Line = Lines[LineIdx];
		if (Line.empty())
			continue;
		if (Line[0] == '>')
			{
			if (!Seq.empty())
				{
				uint SeqIdx = SIZE(IntSeqs);
				GetHexIntSeq(AlphaSize, Seq, IntSeq);
				IntSeqs.push_back(IntSeq);
				asserta(!Label.empty());
				Labels.push_back(Label);
				LabelToSeqIdx[Label] = SeqIdx;
				Seq.clear();
				}
			Label = Line.substr(1, string::npos);
			if (!opt(keepscopid))
				{
				Split(Label, Fields, '/');
				Label = Fields[0];
				}
			}
		else
			Seq += Line;
		}
	if (!Seq.empty())
		{
		uint SeqIdx = SIZE(IntSeqs);
		GetHexIntSeq(AlphaSize, Seq, IntSeq);
		IntSeqs.push_back(IntSeq);
		asserta(!Label.empty());
		LabelToSeqIdx[Label] = SeqIdx;
		Labels.push_back(Label);
		}
	}
