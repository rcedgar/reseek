#pragma once

extern FILE *g_fTsv;
extern FILE *g_fAln;
extern FILE *g_fFasta2;

void OpenOutputFiles();
void CloseOutputFiles();
