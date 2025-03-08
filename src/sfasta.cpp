#include "sfasta.h"
#include "alpha.h"
#include "timing.h"

void SeqToFasta(FILE *f, const char *Label, const char *Seq, unsigned L,
  uint ROWLEN)
	{
	if (f == 0)
		return;
	if (L == 0)
		return;

	if (Label != 0)
		fprintf(f, ">%s\n", Label);
	unsigned BlockCount = (L + ROWLEN - 1)/ROWLEN;
	for (unsigned BlockIndex = 0; BlockIndex < BlockCount; ++BlockIndex)
		{
		unsigned From = BlockIndex*ROWLEN;
		unsigned To = From + ROWLEN;
		if (To >= L)
			To = L;
		for (unsigned Pos = From; Pos < To; ++Pos)
			fputc(Seq[Pos], f);
		fputc('\n', f);
		}
	}

void SeqToFasta(FILE *f, const string &Label, const string &Seq,
  uint BLOCKLEN)
	{
	SeqToFasta(f, Label.c_str(), Seq.c_str(), SIZE(Seq), BLOCKLEN);
	}

bool FastaFileIsNucleo(FILE *f)
	{
	unsigned SampleSize = 1024;
	uint64 CurrPos = GetStdioFilePos64(f);
	uint64 FileSize = GetStdioFileSize64(f);

	SetStdioFilePos64(f, 0);
	char lastc = '\n';
	bool InLabel = false;
	unsigned NucleoCount = 0;
	unsigned LetterCount = 0;
	for (uint64 Pos = 0; Pos < FileSize; ++Pos)
		{
		char c;
		ReadStdioFile(f, &c, 1);
		if (c == '>' && (lastc == '\n' || lastc == '\r'))
			InLabel = true;
		else if (InLabel && (c == '\n' || c == '\r'))
			InLabel = false;
		else if (!InLabel && isalpha(c))
			{
			++LetterCount;
			if (g_IsNucleoChar[c])
				++NucleoCount;
			if (LetterCount >= SampleSize)
				break;
			}
		lastc = c;
		}

	bool IsNucleo = (LetterCount > 0 && double(NucleoCount)/double(LetterCount) > 0.9);
	SetStdioFilePos64(f, CurrPos);
	return IsNucleo;
	}

static unsigned GetMaxPoly(const char *Seq, unsigned L)
	{
	char CurrChar = Seq[0];
	unsigned Start = 0;
	unsigned MaxLen = 1;
	for (unsigned i = 1; i < L; ++i)
		{
		char c = Seq[i];
		if (c != CurrChar || i+1 == L)
			{
			unsigned Len = i - Start;
			if (Len > MaxLen)
				MaxLen = Len;
			CurrChar = c;
			Start = i;
			}
		}
	return MaxLen;
	}

SFasta::SFasta()
	{
	m_FileName = "";
	m_File = 0;
	m_Buffer = 0;
	m_BufferSize = 0;
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_FilePos = 0;
	m_FileSize = 0;
	m_Label = 0;
	m_SeqLength = 0;
	m_TooLongCount = 0;
	m_LongestLength = 0;
	m_IsNucleo = false;
	m_IsNucleoSet = false;
	}

SFasta::~SFasta()
	{
	Clear();
	}

void SFasta::Clear()
	{
	myfree(m_Buffer);
	if (m_File != 0)
		CloseStdioFile(m_File);

	m_FileName = "";
	m_File = 0;
	m_Buffer = 0;
	m_BufferSize = 0;
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_FilePos = 0;
	m_FileSize = 0;
	m_Label = 0;
	m_SeqLength = 0;
	m_SeqIndex = UINT_MAX;
	m_AllowGaps = false;
	m_IsNucleo = false;
	m_IsNucleoSet = false;
	m_TooLongCount = 0;
	m_LongestLength = 0;
	m_BadCharToCount.clear();
	}

void SFasta::LogMe() const
	{
	Log("\n");
	Log("SFasta::LogMe()\n");
	Log("FileName=%s\n", m_FileName.c_str());
	Log("FileSize=%u\n", (unsigned) m_FileSize);
	Log("FilePos=%u\n", (unsigned) m_FilePos);
	Log("BufferSize=%u\n", m_BufferSize);
	Log("BufferPos=%u\n", m_BufferOffset);
	Log("BufferBytes=%u\n", m_BufferBytes);
	if (m_Label == 0)
		Log("Label=NULL\n");
	else
		Log("Label=%s\n", m_Label);
	Log("SeqLength=%u\n", m_SeqLength);
	}

const char *SFasta::GetNextSeq()
	{
	for (;;)
		{
		const char *Seq = GetNextSeqLo();
		if (Seq == 0)
			{
			if (m_TooLongCount > 0)
				Warning("%u long sequences (-maxseqlength %u, longest %u) discarded from %s",
				  m_TooLongCount, opt(maxseqlength), m_LongestLength, m_FileName.c_str());

			if (!m_BadCharToCount.empty())
				{
				string BadLetterString;
				unsigned Total = 0;
				Log("\n");
				Log("Bad char  Count\n");
				Log("--------  -----\n");
				for (map<char, unsigned>::const_iterator p = m_BadCharToCount.begin();
				  p != m_BadCharToCount.end(); ++p)
					{
					char c = p->first;
					if (isprint(c))
						BadLetterString += c;
					else
						Psa(BadLetterString, "[%02x]", c);
					unsigned Count = p->second;
					Total += Count;
					if (isprint(c))
						Log("%8c  %5u\n", c, Count);
					else
						Log("    0x%02x  %5u\n", c, Count);
					}
				Warning("%u total invalid letters in FASTA {%s}",
				  Total, BadLetterString.c_str());
				}
			return 0;
			}
		if (m_SeqLength > opt(maxseqlength) && opt(maxseqlength) != 0)
			{
			if (m_LongestLength == 0 || m_SeqLength > m_LongestLength)
				m_LongestLength = m_SeqLength;
			++m_TooLongCount;
			continue;
			}
		return Seq;
		}
	}

const char *SFasta::GetNextSeqLo()
	{
// End of cache?
	if (m_BufferOffset == m_BufferBytes)
		{
	// End of file?
		if (m_FilePos == m_FileSize)
			return 0;
		FillCache();
		}

	StartTimer(SF_GetNextSeq);
	asserta(m_Buffer[m_BufferOffset] == '>');
	m_Label = (char *) (m_Buffer + m_BufferOffset + 1);
	
	char *ptr = 0;
	for (unsigned i = m_BufferOffset; i < m_BufferSize; ++i)
		{
		char c = m_Buffer[i];
		if (c == '\n' || c == '\r')
			{
			ptr = m_Buffer + i;
			break;
			}
		}
	asserta(ptr != 0);

	if (opt(trunclabels))
		{
		for (char *p = m_Label; *p; ++p)
			if (isspace(*p))
				{
				*p = 0;
				break;
				}
		}
	else
		{
		for (char *p = m_Label; *p; ++p)
			{
			if (*p == '\t')
				*p = ' ';
			else if (*p == '\r' || *p == '\n')
				{
				*p = 0;
				char NextChar = *(p+1);
				if (NextChar == '\r' || NextChar == '\n')
					++p;
				break;
				}
			}
		}

// ptr points to end-of-line.
// Move to start of sequence data.
	char *Seq = ++ptr;

// Delete white space in-place
	char *To = ptr;
	m_BufferOffset = (unsigned) (ptr - m_Buffer);
	while (m_BufferOffset < m_BufferBytes)
		{
		char c = m_Buffer[m_BufferOffset];
		if (c == '>')
			{
			char prevc = '\n';
			if (m_BufferOffset > 0)
				prevc = m_Buffer[m_BufferOffset-1];
			if (prevc == '\n' || prevc == '\r' || prevc == 0)
				break;
			}
		++m_BufferOffset;
		if (isgap(c))
			{
			if (m_AllowGaps)
				*To++ = c;
			continue;
			}
		if (g_IsSeqChar[c])
			*To++ = c;
		else if (isspace(c))
			continue;
		else
			{
			if (m_BadCharToCount.find(c) == m_BadCharToCount.end())
				m_BadCharToCount[c] = 1;
			else
				m_BadCharToCount[c] += 1;
			if (m_IsNucleo)
				*To++ = 'N';
			else
				*To++ = 'X';
			continue;
			}
		}
	m_SeqLength = unsigned(To - Seq);

	if (m_SeqIndex == UINT_MAX)
		m_SeqIndex = 0;
	else
		++m_SeqIndex;

	EndTimer(SF_GetNextSeq);
	return Seq;
	}

void SFasta::Open(const string &FileName)
	{
	Clear();
	m_FileName = FileName;
	m_File = OpenStdioFile(FileName);
	m_BufferSize = opt(sfasta_buff_bytes);
	//m_Buffer = myalloc<char>(m_BufferSize);
	m_Buffer = myalloc(char, m_BufferSize);
	m_FileSize = GetStdioFileSize64(m_File);
	}

void SFasta::Rewind()
	{
	m_BufferOffset = 0;
	m_BufferBytes = 0;
	m_FilePos = 0;
	m_SeqIndex = UINT_MAX;
	m_TooLongCount = 0;
	}

bool SFasta::GetIsNucleo()
	{
	if (m_IsNucleoSet)
		return m_IsNucleo;

	if (m_File == 0)
		Die("FASTASeqSource::GetBuffer, not open");

	bool FastaFileIsNucleo(FILE *f);
	m_IsNucleo = FastaFileIsNucleo(m_File);
	m_IsNucleoSet = true;

	return m_IsNucleo;
	}

void SFasta::FillCache()
	{
	StartTimer(SF_FillCache);
	asserta(m_FilePos < m_FileSize);

	uint64 BytesToRead64 = m_FileSize - m_FilePos;

	bool FinalBuffer = true;
	if (BytesToRead64 > m_BufferSize)
		{
		FinalBuffer = false;
		BytesToRead64 = m_BufferSize;
		}

	asserta(BytesToRead64 <= UINT32_MAX);
	uint32 BytesToRead = (uint32) BytesToRead64;
	SetStdioFilePos64(m_File, m_FilePos);
	ReadStdioFile(m_File, m_Buffer, BytesToRead);
	if (m_Buffer[0] != '>')
		{
		if (m_FilePos == 0)
			Die("Input is not FASTA file");
		else
			Die("SFasta::FillCache() failed, expected '>'");
		}

	m_BufferOffset = 0;

// If last buffer in file, done
	if (FinalBuffer)
		{
		m_BufferBytes = BytesToRead;
		m_FilePos += BytesToRead;
		EndTimer(SF_FillCache);
		return;
		}

// If not last buffer, truncate any partial sequence
// at end of buffer. Search backwards to find last '>'.
	char *ptr = m_Buffer + BytesToRead - 1;
	while (ptr > m_Buffer)
		{
		if (ptr[0] == '>' && (ptr[-1] == '\n' || ptr[-1] == '\r'))
			break;
		--ptr;
		}

	if (ptr == m_Buffer)
		{
		LogMe();
		if (*ptr != '>')
			{
	// No '>' found.
	// This might techincally be legal FASTA if the entire
	// buffer is white space, but strange if not the last buffer
	// in the file, so quit anyway.
			Die("Failed to find '>' (pos=%u, bytes=%u)",
			  (unsigned) m_FilePos, BytesToRead);
			}
		else
			{
	// Entire buffer is one sequence which may be truncated.
			Die("Sequence too long (pos=%u, bytes=%u, -sfasta_buff_bytes %u)",
			  (unsigned) m_FilePos, BytesToRead, opt(sfasta_buff_bytes));
			}
		}

	asserta(*ptr == '>');

	m_BufferBytes = unsigned(ptr - m_Buffer);
	m_FilePos += m_BufferBytes;

	EndTimer(SF_FillCache);
	}

unsigned SFasta::GetPctDoneX10() const
	{
	if (m_FilePos == 0 || m_FileSize == 0)
		return 0;

	assert(m_FilePos >= m_BufferBytes);
	uint64 BufferStart = m_FilePos - m_BufferBytes;
	uint64 BufferPos = BufferStart + m_BufferOffset;

	unsigned iPctX10 = unsigned(10.0*double(BufferPos)*100.0/double(m_FileSize));
	if (iPctX10 == 0)
		return 1;
	if (iPctX10 >= 999)
		return 998;
	return iPctX10;
	}

#if	TEST
void TestSFasta()
	{
	SFasta SF;
	SF.Open(opt(input));

	if (opt(verbose))
		{
		Log("  Index   Length  Label\n");
		Log("-------  -------  -----\n");
		}

	unsigned Index = 0;
	unsigned SeqCount = 0;
	double LetterCount = 0.0;
	ProgressStep(0, 1000, "Reading");
	for (;;)
		{
		const char *Seq = SF.GetNextSeq();
		if (Seq == 0)
			break;
		ProgressStep(SF.GetPctDoneX10(), 1000, "Reading");
		const char *Label = SF.GetLabel();
		unsigned L = SF.GetSeqLength();
		++SeqCount;
		LetterCount += L;

		if (opt(verbose))
			{
			Log(">%7u  %7u  '%s'\n", Index, L, Label);
			Log("+%7.7s  %7.7s  \"%*.*s\"\n", "", "", L, L, Seq);
			}

		++Index;
		}
	ProgressStep(999, 1000, "Reading");

	Progress("%u seqs, %s letters\n", SeqCount, FloatToStr(LetterCount));
	Log("%u seqs, %s letters\n", SeqCount, FloatToStr(LetterCount));
	}
#endif // TEST
