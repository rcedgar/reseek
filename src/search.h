#pragma once

uint MuFilter(const DSSParams &Params,
			  SeqDB &QueryDB,
			  SeqSource &DBSS,
			  const string &OutputFileName);

void PostMuFilter(const DSSParams &Params,
				  const string &MuFilterTsvFN,
				  const string &QueryCAFN,
				  const string &DBBCAFN,
				  float MaxEvalue,
				  const string &HitsFN);
