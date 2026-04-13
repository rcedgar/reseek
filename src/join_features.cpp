#include "myutils.h"
#include "lookup.h"
#include "fastbench.h"

static void read_all(
	const string &fn,
	const lookup &look,
	uint fldidxq,
	uint fldidxt,
	uint fldidxvalue,
	float *values)
	{
	ProgressLog("%s\n", fn.c_str());
	uint maxfldidx = max(max(fldidxq, fldidxt), fldidxvalue);
	FILE *f = OpenStdioFile(fn);
	string line;
	vector<string> flds;
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() > maxfldidx);
		const string &q = flds[fldidxq];
		const string &t = flds[fldidxt];
		const string &v = flds[fldidxvalue];
		uint domidxq = look.get_domidx(q);
		uint domidxt = look.get_domidx(t);
		float value = StrToFloatf(v);
		uint k = look.
			get_pair_idx_upper_triangle_with_diagonal(domidxq, domidxt);
		values[k] = value;
		}
	CloseStdioFile(f);

	const uint npair =
		look.get_pair_count_upper_triangle_with_diagonal();
	uint nmiss = 0;
	for (uint i = 0; i < npair; ++i)
		if (values[i] == FLT_MAX) ++nmiss;
	if (nmiss > 0)
		ProgressLog("%s: %u missing\n", fn.c_str(), nmiss);
	}

static void read_self(
	const string &fn,
	const lookup &look,
	float *values)
	{
	ProgressLog("%s\n", fn.c_str());
	FILE *f = OpenStdioFile(fn);
	string line;
	vector<string> flds;
	while (ReadLineStdioFile(f, line))
		{
		Split(line, flds, '\t');
		asserta(flds.size() == 2);
		const string &q = flds[0];
		const string &v = flds[1];
		uint domidxq = look.get_domidx(q);
		float value = StrToFloatf(v);
		values[domidxq] = value;
		}
	CloseStdioFile(f);

	const uint ndom = look.get_ndom();
	uint nmiss = 0;
	for (uint i = 0; i < ndom; ++i)
		if (values[i] == FLT_MAX) ++nmiss;
	if (nmiss > 0)
		ProgressLog("%s: %u missing\n", fn.c_str(), nmiss);
	}

void cmd_join_features()
	{
	asserta(optset_lookup);
	asserta(optset_output);
	asserta(optset_output2);

	FILE *fout = CreateStdioFile(opt(output));
	FILE *fout2 = CreateStdioFile(opt(output2));

	vector<string> lines;
	ReadLinesFromFile(g_Arg1, lines);
	const uint nfile = uint(lines.size());

	const string &dir = opt(filesdir);

	lookup look;
	look.from_tsv(opt(lookup));

	const uint ndom = look.get_ndom();
	const uint npair =
		look.get_pair_count_upper_triangle_with_diagonal();

	vector<string> namevec;
	vector<float *> valuevec;
	vector<bool> selfvec;

	vector<string> flds;
	for (uint fileidx = 0; fileidx < nfile; ++fileidx)
		{
		//  self	nu		selfrev_nu.tsv
		//	all		entropy	entropy.tsv		2	3	1
		const string &line = lines[fileidx];
		if (line == "" || line[0] == '#')
			continue;
		Split(line, flds, '\t');
		asserta(flds.size() >= 3);
		const string &cat = flds[0];
		const string &name = flds[1];
		const string &fn = dir + "/" + flds[2];
		if (cat == "all")
			{
			float *values = myalloc(float, npair);
			for (uint i = 0; i < npair; ++i) values[i] = FLT_MAX;
			asserta(flds.size() == 6);
			uint fldidxq = StrToUint(flds[3]) - 1;
			uint fldidxt = StrToUint(flds[4]) - 1;
			uint fldidxvalue = StrToUint(flds[5]) - 1;
			read_all(fn, look, fldidxq, fldidxt, fldidxvalue, values);
			namevec.push_back(name);
			valuevec.push_back(values);
			selfvec.push_back(false);
			}
		else if (cat == "self")
			{
			float *values = myalloc(float, ndom);
			for (uint i = 0; i < ndom; ++i) values[i] = FLT_MAX;
			asserta(flds.size() == 3);
			read_self(fn, look, values);
			namevec.push_back(name);
			valuevec.push_back(values);
			selfvec.push_back(true);
			}
		else
			Die("cat=%s", cat.c_str());
		}

	const uint nname = uint(namevec.size());
	fprintf(fout, "query+target");
	uint32_t nf = 0;
	string hdr;
	for (uint nameidx = 0; nameidx < nname; ++nameidx)
		{
		bool self = selfvec[nameidx];
		if (self)
			{
			string sq = "q_" + namevec[nameidx];
			string st = "t_" + namevec[nameidx];

			fprintf(fout, "\t%s", sq.c_str());
			fprintf(fout, "\t%s", st.c_str());
			nf += 2;

			hdr += sq + ";";
			hdr += st + ";";
			}
		else
			{
			const string &s = namevec[nameidx];
			fprintf(fout, "\t%s", s.c_str());
			nf += 1;

			hdr += s + ";";
			}
		}
	fprintf(fout, "\tTP");
	fprintf(fout, "\n");

	asserta(hdr.size() < 100);
	hdr.resize(100);

	WriteStdioFile(fout2, &nf, sizeof(nf));
	WriteStdioFile(fout2, hdr.c_str(), 100);

	uint counter = 0;
	for (uint idxq = 0; idxq < ndom; ++idxq)
		{
		const string &q = look.get_dom(idxq);
		for (uint idxt = idxq; idxt < ndom; ++idxt)
			{
			vector<float> vv;
			uint k = look.
				get_pair_idx_upper_triangle_with_diagonal(idxq, idxt);
			ProgressStep(counter++, npair, "Writing output");

			const string &t = look.get_dom(idxt);

			fprintf(fout, "%s", q.c_str());
			fprintf(fout, "+%s", t.c_str());

			for (uint nameidx = 0; nameidx < nname; ++nameidx)
				{
				bool self = selfvec[nameidx];
				if (self)
					{
					float vq = valuevec[nameidx][idxq];
					float vt = valuevec[nameidx][idxt];

					fprintf(fout, "\t%.4g", vq);
					fprintf(fout, "\t%.4g", vt);

					vv.push_back(vq);
					vv.push_back(vt);
					}
				else
					{
					float v = valuevec[nameidx][k];
					fprintf(fout, "\t%.4g", v);

					vv.push_back(v);
					}
				}

			bool is_tp = look.is_tp_k(k);
			fprintf(fout, "\t%d", int(is_tp));
			fprintf(fout, "\n");

			asserta(vv.size() == nf);
			WriteStdioFile(fout2, vv.data(), nf*sizeof(float));
			}
		}

	CloseStdioFile(fout);
	CloseStdioFile(fout2);
	}

float *read_join_data(
	const string &fn,
	const lookup &look,
	vector<string> &names)
	{
	const uint npair =
		look.get_pair_count_upper_triangle_with_diagonal();
	FILE *f = OpenStdioFile(fn);
	uint32_t nf;
	ReadStdioFile(f, &nf, sizeof(nf));
	ProgressLog("nf=%u\n", nf);
	char hdr[101];
	ReadStdioFile(f, hdr, 100);
	hdr[100] = 0;
	Split(hdr, names, ';');
	const uint nname = uint(names.size());
	ProgressLog("%u names\n", nname);
	unordered_map<string, uint> name2idx;
	for (size_t i = 0; i < nname; ++i)
		{
		const string &name = names[i];
		name2idx[name] = uint(i);
		ProgressLog("[%2u]  %s\n" , i, name.c_str());
		}

	uint nfloat = npair*nname;
	float *data = myalloc(float, nfloat);
	uint bytes = nfloat*sizeof(float);
	ProgressLog("Reading %s bytes... ", IntToStr(bytes));
	ReadStdioFile(f, data, bytes);
	ProgressLog("ok\n");
	CloseStdioFile(f);
	return data;
	}

void cmd_join_stats()
	{
	asserta(optset_lookup);
	lookup look;
	look.from_tsv(opt(lookup));
	vector<string> names;
	float *data = read_join_data(g_Arg1, look, names);

	FastBench FB;
	FB.m_scores_are_evalues = opt(scores_are_evalues);
	FB.ReadLookup(opt(lookup));
	FB.Alloc();

	uint nf = uint(names.size());
	const uint npair =
		look.get_pair_count_upper_triangle_with_diagonal();
	for (uint i = 0; i < npair; ++i)
		{
		float mega = data[nf*i + 4];
		FB.m_Scores[i] = mega;
		}
	ProgressLog("Sorting...");
	FB.SetScoreOrder();
	ProgressLog("\n");
	FB.Bench();
	}
