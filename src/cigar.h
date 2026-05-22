#pragma once

// NOTE -- sometimes parasail_cigar_decode returns
// cig_rev->beg_query=0, cig_rev->beg_ref=0 and
// a CIGAR string which begins with Ds or Is
// This is a quirk more than a bug, it's a non-
// standard way to represent local alignment

void parasail_cigar_to_path(
	const string &para_cigar,
	uint para_loi, uint para_loj,
	uint &loi, uint &loj,
	string &path);

void LocalPathToCIGAR(const char *Path, uint LoQ, uint LoR, string &CIGAR,
  bool FlipDI);
void PathToCIGAR(const char *Path, string &CIGAR, bool FlipDI = false);
void CIGARGetOps(const string &CIGAR, string &Ops, vector<uint> &Lengths);
const char *LocalCIGARToPath(const string &CIGAR, string &Path,
	uint &LoQ, uint &LoR, bool FlipDI);
const char *CIGARToPath(const string &CIGAR, string &Path, bool FlipDI = false);
void CIGARToLs(const string &CIGAR, uint &QL, uint &TL);
void PathToLs(const string &Path, uint &QL, uint &TL);
void ExpandParaCigar(const string &s, string &Path);
void ExpandParaCigar_reverseDI(const string &s, string &Path);
void GetPathCounts(const string &Path, uint &M, uint &D, uint &I);

uint find_closest_point(
	const string &cigar,
	uint loQ, uint loT,
	uint LQ, uint LT,
	uint posQ, uint posT,
	uint &closest_posQ, uint &closest_posT);
