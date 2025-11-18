#include "myutils.h"
#include "fastaseqsource.h"
#include "seqinfo.h"
#include "alpha.h"

bool FastaFileIsNucleo(FILE *f);

FASTASeqSource::FASTASeqSource()
	{
	}

FASTASeqSource::~FASTASeqSource()
	{
	}

bool FASTASeqSource::GetIsNucleo()
	{
	FILE *f = m_LR.m_f;
	asserta(f != 0);
	return FastaFileIsNucleo(f);
	}

// Caller must own memory because SeqSource may be shared
// between threads, so SeqInfo should be thread-private.
bool FASTASeqSource::GetNextLo(SeqInfo *SI)
	{
	if (m_LR.m_EOF)
		return false;

	bool TruncLabels = opt(trunclabels);
// Outer for loop just to allow skipping of empty sequences
	for (;;)
		{
	// Special case at start of file
		if (m_LR.m_LineNr == 0)
			{
			bool Ok = ReadLine();
			if (!Ok)
				return false;
			}

		unsigned SeqIndex = m_SeqCount;
		SI->Init(SeqIndex);

		const char *Line = m_LineBuff.Data;
		unsigned n = m_LineBuff.Size;
		if (n == 0)
			{
			bool Ok = ReadLine();
			if (!Ok)
				return false;
			asserta(n > 0);
			}
		if (Line[0] != '>')
			Die("Bad FASTA file %s, expected '>' in line %u",
			  GetFileNameC(), m_LR.m_LineNr);
		SI->AllocLabel(n);
		char *Label = SI->m_LabelBuffer;
		for (unsigned i = 1; i < n; ++i)
			{
			byte c = Line[i];
			if (TruncLabels && isspace(c))
				{
				Label[i-1] = 0;
				break;
				}
			Label[i-1] = c;
			}
		Label[n-1] = 0;

		unsigned SeqLength = 0;
		for (;;)
			{
			bool Ok = ReadLine();
			if (!Ok)
				break;
			const char *Line = m_LineBuff.Data;
			unsigned n = m_LineBuff.Size;
			if (n > 0 && Line[0] == '>')
				break;
			SI->m_L = SeqLength;
			SI->AllocL(SeqLength + n);
			byte *Seq = SI->m_SeqBuffer;
			for (unsigned i = 0; i < n; ++i)
				{
				byte c = (byte) Line[i];
				if (isspace(c))
					continue;
				if (c == '-' || c == '.')
					{
					if (m_StripGaps)
						continue;
					}
				else if (!isalpha(c) && c != '*')
					{
					if (m_AllowDigits && isdigit(c))
						Seq[SeqLength++] = c;
					else
						BadByte(c);
					continue;
					}
				Seq[SeqLength++] = c;
				}
			}

		SI->m_L = SeqLength;
		if (SeqLength > 0)
			return true;
		else
			{
			Warning("Empty sequence at line %u in FASTA file %s, label >%s",
			  GetLineNr(), GetFileNameC(), SI->m_Label);
			continue;
			}
		}
	}
