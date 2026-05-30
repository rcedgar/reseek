#pragma once

class PDBChain;
#include "flat_chain.h"
#include "query_data.h"
#include "chaq.h"

// Binary C-alpha
class BCAData
	{
public:
	vector<string> m_Labels;
	vector<uint64_t> m_Offsets; // start of IC vector in file
	vector<uint32_t> m_SeqLengths;
	bool m_HasNuSequences = false;
	string m_FN;
	FILE *m_f = 0;
	bool m_Writing = false;
	bool m_Reading = false;
	uint64 m_SeqLengthsPos64 = UINT64_MAX;
	uint64 m_LabelDataSize64 = UINT64_MAX;
	sid_t *m_distmx = 0;
	chaq_vecs2 m_cv;
	uint8_t *m_scratch_buffer = 0;
	uint m_scratch_buffer_bytes = 0;
	uint8_t *m_codeseq_nu = 0;
	uint m_maxL = 4000;
	mutable mutex m_ReadLock;

public:
	void Clear();
	void Create(const string &FN, bool WithNu = false);
	void Open(const string &FN);
	void WriteChain(const PDBChain &Chain);

	const string &GetLabel(uint idx) const
		{
		asserta(idx < m_Labels.size());
		return m_Labels[idx];
		}

	uint GetSeqLength(uint idx) const
		{
		asserta(idx < m_SeqLengths.size());
		return m_SeqLengths[idx];
		}

	void ReadChain(uint64 ChainIdx, PDBChain &Chain) const;
	flat_chain_t* read_flat_chain(uint64 ChainIdx) const;
	void write_flat_chain(const flat_chain_t *chain);
	void Close();
	uint GetChainCount() const { return SIZE(m_Labels); }
	uint64 GetSeqOffset(uint64 ChainIdx) const;
	uint GetSeqLength(uint64 ChainIdx) const;
	void append_codeseq_nu(const flat_chain_t *chain);
	uint read_codeseq_nu(
		uint8_t *codeseq_nu, uint idx, uint buffer_length) const;
	uint64 get_offset_ICs(uint idx) const;
	uint64 get_offset_aaseq(uint idx) const;
	uint64 get_offset_nuseq(uint idx) const;
	void make_kappa_codeseqs(
		uint8_t ***ptr_kappa_codeseqs,
		uint **ptr_lengths) const;
	void make_nu_and_kappa_codeseqs(
		uint8_t ***ptr_nu_codeseqs,
		uint8_t ***ptr_kappa_codeseqs,
		uint **ptr_lengths) const;
	query_data *get_query_data(const flat_params &params,
		uint idx);
	query_data **get_query_data_vec(const flat_params &params);

private:
	void CloseWriter();
	void CloseReader();
	};

const uint32_t BCA_MAGIC = 0xBCABCA;
const uint32_t BCB_MAGIC = 0xBCBBCB;
