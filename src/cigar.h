#pragma once

void LocalPathToCIGAR(const char *Path, uint LoQ, uint LoR, string &CIGAR,
  bool FlipDI);
void PathToCIGAR(const char *Path, string &CIGAR, bool FlipDI = false);
void CIGARGetOps(const string &CIGAR, string &Ops, vector<uint> &Lengths);
const char *LocalCIGARToPath(const string &CIGAR, string &Path,
  uint &LoQ, uint &LoR, bool FlipDI);
const char *CIGARToPath(const string &CIGAR, string &Path);
void CIGARToLs(const string &CIGAR, uint &QL, uint &TL);
void PathToLs(const string &Path, uint &QL, uint &TL);
void ExpandParaCigar(const string &s, string &Path);
void ExpandParaCigar_reverseDI(const string &s, string &Path);
