#pragma once

class tabbedlines
	{
public:
	vector<string> m_lines;
	vector<string> m_flds;
	size_t m_linenr = 0;

public:
	tabbedlines() {}							// for writing
	tabbedlines(const vector<string> &lines)	// for reading
		{
		m_lines = lines;
		m_linenr = 0;
		}

	void get_eof() const
		{
		asserta(m_linenr == SIZE(m_lines));
		}

	const string &get()
		{
		for (;;)
			{
			asserta(m_linenr < m_lines.size());
			const string& line = m_lines[m_linenr++];
			if (StartsWith(line, "#"))
				continue;
			Split(line, m_flds, '\t');
			return line;
			}
		}

	const string &get_str(const string &fld0)
		{
		get();
		asserta(SIZE(m_flds) == 2);
		asserta(m_flds[0] == fld0);
		return m_flds[1];
		}

	uint get_int(const string &fld0)
		{
		get();
		asserta(SIZE(m_flds) == 2);
		asserta(m_flds[0] == fld0);
		return StrToUint(m_flds[1]);
		}

	float get_float(const string &fld0)
		{
		get();
		asserta(SIZE(m_flds) == 2);
		asserta(m_flds[0] == fld0);
		return StrToFloatf(m_flds[1]);
		}

	void get_int_vec(const string &fld0, uint n, vector<uint> &v)
		{
		v.clear();
		get();
		asserta(SIZE(m_flds) == n+1);
		asserta(m_flds[0] == fld0);
		for (uint i = 0; i < n; ++i)
			v.push_back(StrToUint(m_flds[i+2]));
		}

	uint8_t *get_int8_vec(const string &fld0, uint &n)
		{
		get();
		uint nfld = SIZE(m_flds);
		asserta(nfld > 2);
		asserta(m_flds[0] == fld0);
		n = nfld - 2;
		uint8_t *v = myalloc(uint8_t, n);
		for (uint i = 0; i < n; ++i)
			{
			uint value = StrToUint(m_flds[i+2]);
			asserta(value <= UINT8_MAX);
			v[i] = value;
			}
		return v;
		}

	uint16_t *get_int16_vec(const string &fld0, uint &n)
		{
		get();
		uint nfld = SIZE(m_flds);
		asserta(nfld > 2);
		asserta(m_flds[0] == fld0);
		n = nfld - 2;
		uint16_t *v = myalloc(uint16_t, n);
		for (uint i = 0; i < n; ++i)
			{
			uint value = StrToUint(m_flds[i+2]);
			asserta(value <= UINT16_MAX);
			v[i] = value;
			}
		return v;
		}

	uint *get_int_vec(const string &fld0, uint &n)
		{
		get();
		uint nfld = SIZE(m_flds);
		asserta(nfld > 2);
		n = nfld - 2;
		asserta(m_flds[0] == fld0);
		uint *v = myalloc(uint, n);
		for (uint i = 0; i < n; ++i)
			v[i] = StrToUint(m_flds[i+2]);
		return v;
		}

	void get_signed_int_vec(const string &fld0, uint n,
		vector<int> &v)
		{
		v.clear();
		get();
		asserta(SIZE(m_flds) == n+1);
		asserta(m_flds[0] == fld0);
		for (uint i = 0; i < n; ++i)
			v.push_back(StrToInt(m_flds[i+1]));
		}

	int* get_signed_int_flat_vec(const string &fld0, uint &n)
		{
		get();
		asserta(m_flds[0] == fld0);
		asserta(SIZE(m_flds) > 2);
		n = StrToUint(m_flds[1]);
		int *v = myalloc(int, n);
		asserta(SIZE(m_flds) == n+2);
		for (uint i = 0; i < n; ++i)
			v[i] = StrToInt(m_flds[i+2]);
		return v;
		}

	float *get_float_vec(const string &name, uint &n)
		{
		get();
		asserta(m_flds[0] == name);
		uint nfld = uint(m_flds.size());
		asserta(nfld > 2);
		n = nfld - 2;
		float *v = myalloc(float, n);
		for (uint i = 0; i < n; ++i)
			v[i] = StrToFloatf(m_flds[2+i]);
		return v;
		}

	void get_str_vec(const string &name, vector<string> &v)
		{
		get();
		asserta(m_flds[0] == name);
		uint nfld = uint(m_flds.size());
		asserta(nfld > 2);
		uint n = nfld - 2;
		v.clear();
		v.reserve(n);
		for (uint i = 0; i < n; ++i)
			{
			const string &s = m_flds[2+i];
			v.push_back(s);
			}
		}

	void put_float_flat_square_mx(
		const string &name, uint n, const float *v)
		{
		put_int(name, n);
		put_float_flat_square_mx(n, v);
		}

	void put_float_flat_square_mx(uint n, const float *v)
		{
		for (uint i = 0; i < n; ++i)
			{
			string line;
			Ps(line, "%u", i);
			for (uint j = 0; j < n; ++j)
				{
				float x = v[i*n + j];
				if (x == FLT_MAX)
					Psa(line, "\t*");
				else
					Psa(line, "\t%.4g", x);
				}
			m_lines.push_back(line);
			}
		}

	float *get_float_flat_square_mx(const string &name, uint &n)
		{
		n = get_int(name);
		float *v = myalloc(float, n*n);
		get_float_flat_square_mx(n, v);
		return v;
		}

	void get_float_flat_square_mx(uint n, float *v)
		{
		for (uint i = 0; i < n; ++i)
			{
			get();
			asserta(SIZE(m_flds) == n+1);
			uint i2 = StrToUint(m_flds[0]);
			asserta(i2 == i);
			for (uint j = 0; j < n; ++j)
				v[i*n + j] = (float) StrToFloat(m_flds[j+1]);
			}
		}

	void get_double_flat_square_mx(uint n, double *v)
		{
		for (uint i = 0; i < n; ++i)
			{
			get();
			asserta(SIZE(m_flds) == n+1);
			uint i2 = StrToUint(m_flds[0]);
			asserta(i2 == i);
			for (uint j = 0; j < n; ++j)
				v[i*n + j] = StrToFloat(m_flds[j+1]);
			}
		}

	uint16_t* get_int16_flat_vec(const string &fld0, uint n)
		{
		get();
		asserta(m_flds[0] == fld0);
		asserta(SIZE(m_flds) > 2);
		uint n2 = StrToUint(m_flds[1]);
		asserta(n2 == n);
		uint16_t *v = myalloc(uint16_t, n);
		asserta(SIZE(m_flds) == n+2);
		for (uint i = 0; i < n; ++i)
			v[i] = StrToInt(m_flds[i+2]);
		return v;
		}

	const string &fld(uint i) const
		{
		asserta(i < SIZE(m_flds));
		return m_flds[i];
		}

	void put_float(const string &name, float x)
		{
		string line;
		if (x == FLT_MAX)
			Ps(line, "%s\tFLT_MAX", name.c_str(), x);
		else
			Ps(line, "%s\t%.4g", name.c_str(), x);
		m_lines.push_back(line);
		}

	void put_int(const string &name, uint i)
		{
		string line;
		if (i == UINT_MAX)
			Ps(line, "%s\tUINT_MAX", name.c_str(), i);
		else
			Ps(line, "%s\t%u", name.c_str(), i);
		m_lines.push_back(line);
		}

	void put_int8_vec(const string &name, const uint8_t *v, uint n)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), n);
		for (uint i = 0; i < n; ++i)
			Psa(line, "\t%u", v[i]);
		m_lines.push_back(line);
		}

	void put_int16_vec(const string &name, const uint16_t *v, uint n)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), n);
		for (uint i = 0; i < n; ++i)
			Psa(line, "\t%u", v[i]);
		m_lines.push_back(line);
		}

	void put_int_vec(const string &name, const vector<uint> &v)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), SIZE(v));
		for (uint i = 0; i < SIZE(v); ++i)
			Psa(line, "\t%u", v[i]);
		m_lines.push_back(line);
		}

	void put_int_vec(const string &name, uint32_t *v, size_t n)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), n);
		for (uint i = 0; i < n; ++i)
			Psa(line, "\t%u", v[i]);
		m_lines.push_back(line);
		}

	void put_float_vec(const string &name, const float *v, size_t n)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), n);
		for (uint i = 0; i < n; ++i)
			Psa(line, "\t%.4g", v[i]);
		m_lines.push_back(line);
		}

	void put_str_vec(const string &name, const vector<string> &v)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), SIZE(v));
		for (uint i = 0; i < SIZE(v); ++i)
			Psa(line, "\t%s", v[i].c_str());
		m_lines.push_back(line);
		}

	void put_signed_int_vec(const string &name, const vector<int> &v)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), SIZE(v));
		for (uint i = 0; i < SIZE(v); ++i)
			Psa(line, "\t%d", v[i]);
		m_lines.push_back(line);
		}

	void put_int16_flat_vec(const string &name, const uint16_t *v, uint n)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), n);
		for (uint i = 0; i < n; ++i)
			Psa(line, "\t%u", v[i]);
		m_lines.push_back(line);
		}

	void put_int_flat_vec(const string &name, const uint *v, uint n)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), n);
		for (uint i = 0; i < n; ++i)
			Psa(line, "\t%u", v[i]);
		m_lines.push_back(line);
		}

	void put_signed_int_flat_vec(const string &name, const int *v, uint n)
		{
		string line;
		Ps(line, "%s\t%u", name.c_str(), n);
		for (uint i = 0; i < n; ++i)
			Psa(line, "\t%d", v[i]);
		m_lines.push_back(line);
		}

	void to_tsv(FILE *f) const
		{
		if (f == 0)
			return;
		for (uint i = 0; i < SIZE(m_lines); ++i)
			fprintf(f, "%s\n", m_lines[i].c_str());
		}

	void to_tsv(const string &fn) const
		{
		FILE *f = CreateStdioFile(fn);
		to_tsv(f);
		CloseStdioFile(f);
		}

	static void to_tsv(const string &fn, const vector<string> &lines)
		{
		FILE *f = CreateStdioFile(fn);
		to_tsv(f, lines);
		CloseStdioFile(f);
		}
	
	static void to_tsv(FILE *f, const vector<string> &lines)
		{
		if (f == 0)
			return;
		for (uint i = 0; i < SIZE(lines); ++i)
			fprintf(f, "%s\n", lines[i].c_str());
		}
	};
