#include "myutils.h"
#include "desc_db.h"

#include <regex>
#include <unordered_set>

uint FindPrime(uint Min, uint Max);

static_assert(sizeof(desc_db_slot) == 16, "desc_db_slot");

uint32 desc_db::HashLabel(const char *s)
	{
	uint32 h = 2166136261u;
	for (; *s; ++s)
		{
		h ^= (byte) *s;
		h *= 16777619u;
		}
	if (h == 0)
		h = 1;
	return h;
	}

void desc_db::Clear()
	{
	if (m_f != 0)
		{
		CloseStdioFile(m_f);
		m_f = 0;
		}
	myfree(m_Slots);
	m_Slots = 0;
	m_n = 0;
	m_TableSize = 0;
	m_RecordsOff = 0;
	m_RecordsBytes = 0;
	m_Labels.clear();
	m_Descs.clear();
	}

void desc_db::AllocSlots(uint32 TableSize)
	{
	asserta(TableSize > 0);
	m_TableSize = TableSize;
	m_Slots = myalloc(desc_db_slot, TableSize);
	memset(m_Slots, 0, sizeof(desc_db_slot)*TableSize);
	}

void desc_db::InsertSlot(uint32 hash, uint32 rec_bytes, uint64 rec_off)
	{
	asserta(m_Slots != 0);
	asserta(m_TableSize > 0);
	asserta(hash != 0);
	uint32 i = hash % m_TableSize;
	for (uint32 n = 0; n < m_TableSize; ++n)
		{
		desc_db_slot &s = m_Slots[i];
		if (s.hash == 0)
			{
			s.hash = hash;
			s.rec_bytes = rec_bytes;
			s.rec_off = rec_off;
			return;
			}
		i = i + 1;
		if (i == m_TableSize)
			i = 0;
		}
	Die("desc_db hash table full");
	}

void desc_db::FromTsv(const string &fn)
	{
	Clear();

	unordered_set<string> seen;
	string line;
	FILE *f = OpenStdioFile(fn);
	ProgressFileInit(f, "Reading %s", fn.c_str());
	while (ReadLineStdioFile(f, line))
		{
		ProgressFileStep();
		if (line.empty() || StartsWith(line, "#"))
			continue;
		size_t tab = line.find('\t');
		if (tab == string::npos)
			Die("desc_db: expected label<TAB>description: %s", line.c_str());
		string label = line.substr(0, tab);
		string desc = line.substr(tab + 1);
		if (label.empty())
			Die("desc_db: empty label");
		if (seen.find(label) != seen.end())
			Die("desc_db: duplicate label '%s'", label.c_str());
		seen.insert(label);
		m_Labels.push_back(label);
		m_Descs.push_back(desc);
		}
	ProgressFileDone();
	CloseStdioFile(f);

	if (SIZE(m_Labels) > UINT32_MAX/4)
		Die("desc_db: too many entries");
	m_n = uint32(SIZE(m_Labels));
	asserta(SIZE(m_Descs) == m_n);

	uint32 TableSize = 1;
	if (m_n > 0)
		{
		uint Lo = m_n*2;
		uint Hi = Lo + m_n/4;
		TableSize = FindPrime(Lo, Hi);
		}
	AllocSlots(TableSize);

	m_RecordsOff = uint64(DESC_DB_HEADER_BYTES) +
	  uint64(m_TableSize)*sizeof(desc_db_slot);
	m_RecordsBytes = 0;
	for (uint32 i = 0; i < m_n; ++i)
		{
		uint64 rec = uint64(m_Labels[i].size()) + 1 +
		  uint64(m_Descs[i].size()) + 1;
		if (rec > UINT32_MAX)
			Die("desc_db: record too large for label '%s'",
			  m_Labels[i].c_str());
		uint32 rec_bytes = uint32(rec);
		uint64 rec_off = m_RecordsOff + m_RecordsBytes;
		InsertSlot(HashLabel(m_Labels[i].c_str()), rec_bytes, rec_off);
		m_RecordsBytes += rec;
		}
	}

void desc_db::ToFile(const string &fn) const
	{
	asserta(fn != "");
	asserta(m_Slots != 0);
	asserta(SIZE(m_Labels) == m_n);
	asserta(SIZE(m_Descs) == m_n);

	FILE *f = CreateStdioFile(fn);
	uint32 Magic = DESC_DB_MAGIC;
	uint32 Version = DESC_DB_VERSION;
	uint32 Reserved[2] = { 0, 0 };

	WriteStdioFile(f, &Magic, sizeof(Magic));
	WriteStdioFile(f, &Version, sizeof(Version));
	WriteStdioFile(f, &m_n, sizeof(m_n));
	WriteStdioFile(f, &m_TableSize, sizeof(m_TableSize));
	WriteStdioFile(f, &m_RecordsOff, sizeof(m_RecordsOff));
	WriteStdioFile(f, &m_RecordsBytes, sizeof(m_RecordsBytes));
	WriteStdioFile(f, Reserved, sizeof(Reserved));

	uint64 SlotBytes = uint64(m_TableSize)*sizeof(desc_db_slot);
	asserta(SlotBytes <= UINT32_MAX);
	WriteStdioFile(f, m_Slots, uint32(SlotBytes));

	uint64 rec_off = GetStdioFilePos64(f);
	asserta(rec_off == m_RecordsOff);
	for (uint32 i = 0; i < m_n; ++i)
		{
		const string &label = m_Labels[i];
		const string &desc = m_Descs[i];
		WriteStdioFile(f, label.c_str(), uint32(label.size() + 1));
		WriteStdioFile(f, desc.c_str(), uint32(desc.size() + 1));
		}
	uint64 rec_end = GetStdioFilePos64(f);
	asserta(rec_end == m_RecordsOff + m_RecordsBytes);

	WriteStdioFile(f, &Magic, sizeof(Magic));
	CloseStdioFile(f);

	const uint64 SlotRAM = uint64(m_TableSize)*sizeof(desc_db_slot);
	ProgressLog("Wrote desc_db %s  n=%u  table=%u  slots=%s  records=%llu\n",
	  fn.c_str(), m_n, m_TableSize, MemBytesToStr(double(SlotRAM)),
	  (unsigned long long) m_RecordsBytes);
	}

void desc_db::FromFile(const string &fn)
	{
	Clear();
	asserta(fn != "");

	m_f = OpenStdioFile(fn);
	const uint64 FileSize = GetStdioFileSize64(m_f);
	if (FileSize < DESC_DB_HEADER_BYTES + sizeof(uint32))
		Die("desc_db: file too small %s", fn.c_str());

	uint32 Magic = 0;
	uint32 Version = 0;
	uint32 Reserved[2];
	ReadStdioFile(m_f, &Magic, sizeof(Magic));
	asserta(Magic == DESC_DB_MAGIC);
	ReadStdioFile(m_f, &Version, sizeof(Version));
	asserta(Version == DESC_DB_VERSION);
	ReadStdioFile(m_f, &m_n, sizeof(m_n));
	ReadStdioFile(m_f, &m_TableSize, sizeof(m_TableSize));
	ReadStdioFile(m_f, &m_RecordsOff, sizeof(m_RecordsOff));
	ReadStdioFile(m_f, &m_RecordsBytes, sizeof(m_RecordsBytes));
	ReadStdioFile(m_f, Reserved, sizeof(Reserved));

	asserta(m_TableSize > 0);
	asserta(m_n <= m_TableSize);
	uint64 SlotBytes = uint64(m_TableSize)*sizeof(desc_db_slot);
	asserta(m_RecordsOff == uint64(DESC_DB_HEADER_BYTES) + SlotBytes);
	asserta(FileSize == m_RecordsOff + m_RecordsBytes + sizeof(uint32));

	AllocSlots(m_TableSize);
	asserta(SlotBytes <= UINT32_MAX);
	ReadStdioFile(m_f, m_Slots, uint32(SlotBytes));

	uint32 Magic2 = 0;
	ReadStdioFile64(m_f, m_RecordsOff + m_RecordsBytes, &Magic2, sizeof(Magic2));
	asserta(Magic2 == DESC_DB_MAGIC);
	}

bool desc_db::Get(const char *label, string &desc) const
	{
	desc.clear();
	asserta(label != 0);
	asserta(m_f != 0);
	asserta(m_Slots != 0);
	asserta(m_TableSize > 0);

	const uint32 h = HashLabel(label);
	uint32 i = h % m_TableSize;
	for (uint32 n = 0; n < m_TableSize; ++n)
		{
		const desc_db_slot &s = m_Slots[i];
		if (s.hash == 0)
			return false;
		if (s.hash == h)
			{
			if (s.rec_bytes < 2)
				Die("desc_db: bad rec_bytes");
			string rec;
			rec.resize(s.rec_bytes);
			ReadStdioFile64(m_f, s.rec_off, &rec[0], s.rec_bytes);
			if (rec[s.rec_bytes - 1] != 0)
				Die("desc_db: record not NUL-terminated");
			if (strcmp(rec.c_str(), label) == 0)
				{
				const char *d = rec.c_str() + strlen(rec.c_str()) + 1;
				asserta(d < rec.c_str() + s.rec_bytes);
				desc = d;
				return true;
				}
			}
		i = i + 1;
		if (i == m_TableSize)
			i = 0;
		}
	return false;
	}

void cmd_createdesc()
	{
	asserta(optset_output);
	desc_db DB;
	DB.FromTsv(g_Arg1);
	DB.ToFile(opt(output));
	}

static void DumpRecords(const desc_db &DB, FILE *fout)
	{
	vector<uint32> idxs;
	idxs.reserve(DB.m_n);
	for (uint32 i = 0; i < DB.m_TableSize; ++i)
		{
		if (DB.m_Slots[i].hash != 0)
			idxs.push_back(i);
		}
	asserta(SIZE(idxs) == DB.m_n);

	sort(idxs.begin(), idxs.end(),
	  [&DB](uint32 a, uint32 b)
		{
		return DB.m_Slots[a].rec_off < DB.m_Slots[b].rec_off;
		});

	string rec;
	for (uint32 k = 0; k < SIZE(idxs); ++k)
		{
		const desc_db_slot &s = DB.m_Slots[idxs[k]];
		rec.resize(s.rec_bytes);
		ReadStdioFile64(DB.m_f, s.rec_off, &rec[0], s.rec_bytes);
		const char *label = rec.c_str();
		const char *desc = label + strlen(label) + 1;
		fputs(label, fout);
		fputc('\t', fout);
		fputs(desc, fout);
		fputc('\n', fout);
		}
	}

void cmd_append_desc()
	{
	asserta(optset_output);
	asserta(optset_desc);
	default_opt(tfield, 2);
	if (opt(tfield) == 0)
		Die("-tfield must be >= 1");
	const uint tidx = opt(tfield) - 1;

	regex labeledit_re;
	const bool use_labeledit = optset_labeledit;
	if (use_labeledit)
		{
		try
			{
			labeledit_re.assign(opt(labeledit));
			}
		catch (const regex_error &)
			{
			Die("Invalid -labeledit regex");
			}
		if (labeledit_re.mark_count() != 1)
			Die("-labeledit must have exactly one capturing group");
		}

	desc_db DB;
	DB.FromFile(opt(desc));

	FILE *fin = OpenStdioFile(g_Arg1);
	FILE *fout = CreateStdioFile(opt(output));
	uint bad = 0;
	string line;
	vector<string> flds;
	string desc;
	smatch m;
	ProgressFileInit(fin, "Appending descriptions");
	while (ReadLineStdioFile(fin, line))
		{
		ProgressFileStep();
		if (line.empty() || StartsWith(line, "#"))
			{
			fputs(line.c_str(), fout);
			fputc('\n', fout);
			continue;
			}
		Split(line, flds, '\t');
		if (SIZE(flds) <= tidx)
			{
			++bad;
			fputs(line.c_str(), fout);
			fputc('\n', fout);
			continue;
			}
		string key = flds[tidx];
		if (use_labeledit)
			{
			if (!regex_search(key, m, labeledit_re) || m[1].length() == 0)
				{
				++bad;
				fputs(line.c_str(), fout);
				fputc('\n', fout);
				continue;
				}
			key = m[1].str();
			}
		if (!DB.Get(key.c_str(), desc))
			desc = "(missing description)";
		flds[tidx] = key + " " + desc;
		for (uint i = 0; i < SIZE(flds); ++i)
			{
			if (i > 0)
				fputc('\t', fout);
			fputs(flds[i].c_str(), fout);
			}
		fputc('\n', fout);
		}
	ProgressFileDone();
	CloseStdioFile(fin);
	CloseStdioFile(fout);
	if (bad > 0)
		Warning("%u bad line%s (too few fields or -labeledit no match)",
		  bad, bad == 1 ? "" : "s");
	}

void cmd_desc_lookup()
	{
	if (!optset_label && !optset_output)
		Die("Must set -label or -output");

	desc_db DB;
	DB.FromFile(g_Arg1);

	if (optset_label)
		{
		string desc;
		if (!DB.Get(opt(label), desc))
			Die("label not found '%s'", opt(label));
		fputs(desc.c_str(), stdout);
		fputc('\n', stdout);
		}

	if (optset_output)
		{
		FILE *fout = CreateStdioFile(opt(output));
		DumpRecords(DB, fout);
		CloseStdioFile(fout);
		}
	}
