#ifndef alpha_h
#define alpha_h

#include <limits.h>
#include <string>

using namespace std;

const unsigned char INVALID_LETTER = 0xff;
const unsigned char INVALID_CHAR = '?';
const unsigned BAD_WORD = UINT_MAX;

extern unsigned char g_AminoAcidChars[];
extern unsigned char g_CharToLetterAmino[];
extern unsigned char g_CharToLetterAminoStop[];
extern unsigned char g_LetterToCharAmino[];
extern unsigned char g_CharToLetterNucleo[];
extern unsigned char g_LetterToCharNucleo[];
extern unsigned char g_CodonWordToAminoLetter[];
extern unsigned char g_CodonWordToAminoChar[];
extern unsigned char g_CharToCompChar[];
extern unsigned char g_CharToCompLetter[];
extern unsigned char g_CharToLetterMu[];
extern unsigned char g_LetterToCharMu[];

extern bool g_IsAminoChar[];
extern bool g_IsNucleoChar[];
extern bool g_IsACGTU[];
extern bool g_IsSeqChar[];

extern float g_AminoFreqs[];

extern unsigned g_CharToLetterRed[];
extern unsigned char g_LetterToCharRed[];
extern unsigned g_RedAlphaSize;

void LogRedAlphaRed();
void ReadRedAlphaFromFile(const string &FileName);
unsigned char GetAminoCharFrom3NucChars(unsigned char c1, unsigned char c2,
  unsigned char c3);

extern unsigned char g_Nucleo_CharToBit[256];
extern unsigned char g_IUPAC_CharToBits[256];

static inline bool IUPAC_Eq(unsigned char Char, unsigned char CharOrWildcard)
	{
	unsigned char Bit = g_Nucleo_CharToBit[Char]; 
	unsigned char Bits = g_IUPAC_CharToBits[CharOrWildcard]; 
	if ((Bit & Bits) == 0)
		return false;
	return true;
	}

static inline bool AminoLetterIsStartCodon(unsigned char Letter)
	{
	return Letter == 10;
	}

static inline bool AminoLetterIsStopCodon(unsigned char Letter)
	{
	return Letter == 20;
	}

const char *WordToStr(unsigned Word, unsigned WordLength, bool Nucleo);
const char *WordToStrNucleo(unsigned Word, unsigned WordLength);
const char *WordToStrAmino(unsigned Word, unsigned WordLength);
const char *WordToStrAmino2(unsigned Word, unsigned WordLength, char *Str);
unsigned StrToWordAmino(const char *Str, unsigned WordLength);

static inline bool isgap(unsigned char c)
	{
	return c == '-' || c == '.';
	}

static inline bool isunaligned(unsigned char c)
	{
	return islower(c) || c == '.';
	}

static inline bool isaligned(unsigned char c)
	{
	return isupper(c) || c == '-';
	}

#endif // alpha_h
