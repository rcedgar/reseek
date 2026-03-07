#pragma once

#include "flat_base.h"

static const size_t RESERVE_CHAIN_LENGTH = 400;
static const uint M = 64;

class flat_chain
	{
public:
	string m_label;
	chainxyz_t *m_xyz;
	chainaa_t *m_aa;
	vector<string> m_lines;

	flat_chain()
		{
		m_label = "(undef)";
		m_xyz = nullptr;
		m_aa = nullptr;
		}

	flat_chain(
		const string &label,
		const vector<char> aas,
		const vector<float> &Xs,
		const vector<float> &Ys,
		const vector<float> &Zs);

	~flat_chain()
		{
		clear();
		}

	void clear()
		{
		m_lines.clear();
		down0(m_xyz);
		down0(m_aa);
		}

	void set_xyz(const vector<float> &Xs,
		const vector<float> &Ys, const vector<float> &Zs);
	void set_aa(const vector<char> &aas);
	bool from_pdb_lines(const string &label,
		const vector<string> &lines, bool save_lines);
	void to_fasta(const string &fn) const;
	void to_fasta(FILE *f) const;
	void to_cal(const string &fn) const;
	void to_cal(FILE *f) const;

	uint32_t get_length() const
		{
		assert(m_aa); return m_aa->m_size;
		}

	void get_ics(uint i, uint16_t &icx, uint16_t &icy, uint16_t &icz) const
		{
		const uint16_t *data = m_xyz->m_data;
		uint k = 3*i;
		uint ic_x = data[k];
		uint ic_y = data[k+1];
		uint ic_z = data[k+2];
		}

	void get_coords(uint i, float &x, float &y, float &z) const
		{
		assert(m_xyz != 0);
		const uint16_t *data = m_xyz->m_data;
		uint k = 3*i;
		uint ic_x = data[k];
		uint ic_y = data[k+1];
		uint ic_z = data[k+2];
		x = ic2coord(ic_x);
		y = ic2coord(ic_y);
		z = ic2coord(ic_z);
		}

	uint16_t slow_sd(uint i, uint j) const
		{
		assert(i < m_aa->m_size);
		assert(j < m_aa->m_size);
		const uint16_t *data = m_xyz->m_data;

		int32_t dx = int32_t(data[3*i]) - int32_t(data[3*j]);
		int32_t dy = int32_t(data[3*i+1]) - int32_t(data[3*j+1]);
		int32_t dz = int32_t(data[3*i+2]) - int32_t(data[3*j+2]);
		return dx*dx + dy*dy + dz*dz;
		}

	float slow_float_dist(uint i, uint j) const
		{
		assert(i < m_aa->m_size);
		assert(j < m_aa->m_size);
		const uint16_t *data = m_xyz->m_data;

		int32_t dx = int32_t(data[3*i]) - int32_t(data[3*j]);
		int32_t dy = int32_t(data[3*i+1]) - int32_t(data[3*j+1]);
		int32_t dz = int32_t(data[3*i+2]) - int32_t(data[3*j+2]);
		
		float d2 = float(dx*dx + dy*dy + dz*dz);
		float d = sqrtf(d2)/10;
		return d;
		}

	char get_aa(uint i) const
		{
		assert(m_xyz != 0);
		assert(i < m_aa->m_size);
		return m_aa->m_data[i];
		}

public:
	static uint16_t coord2ic(float X) { return uint16_t((X + 1000)*10 + 0.5); }
	static float ic2coord(uint16_t IC) { return float(IC/10.0f) - 1000; }
	};

void read_flat_chains(const string &fn, vector<flat_chain *> &chains);
