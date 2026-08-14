#pragma once

#define DESC_DB_MAGIC		0x44455343	// 'DESC'
#define DESC_DB_VERSION		1
#define DESC_DB_HEADER_BYTES	40

/***
On-disk layout (native endian):
	uint32 magic ('DESC')
	uint32 version (1)
	uint32 n
	uint32 table_size
	uint64 records_off
	uint64 records_bytes
	uint32 reserved[2]
	Slot[table_size]		// loaded into RAM
	uint8  records[...]		// NOT loaded; label\\0 desc\\0 per entry
	uint32 magic
***/

struct desc_db_slot
	{
	uint32 hash;		// 0 = empty; real hash 0 stored as 1
	uint32 rec_bytes;	// bytes at rec_off (label\0 desc\0)
	uint64 rec_off;
	};

class desc_db
	{
public:
	FILE *m_f = 0;
	uint32 m_n = 0;
	uint32 m_TableSize = 0;
	uint64 m_RecordsOff = 0;
	uint64 m_RecordsBytes = 0;
	desc_db_slot *m_Slots = 0;

// Convert-time only
	vector<string> m_Labels;
	vector<string> m_Descs;

public:
	~desc_db()
		{
		Clear();
		}

	void Clear();
	void FromTsv(const string &fn);
	void ToFile(const string &fn) const;
	void FromFile(const string &fn);
	bool Get(const char *label, string &desc) const;

private:
	static uint32 HashLabel(const char *s);
	void AllocSlots(uint32 TableSize);
	void InsertSlot(uint32 hash, uint32 rec_bytes, uint64 rec_off);
	};
