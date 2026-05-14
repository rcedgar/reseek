#pragma once

class scratch_mem
	{
public:
	uint8_t *mem = 0;
	size_t size = 0;
	uint8_t *ptr = 0;
	bool owner = false;

public:
	scratch_mem(size_t n)
		{
		asserta(n > 0);
		mem = myalloc64(uint8_t, n);
		size = n;
		ptr = mem;
		owner = true;
		}

	scratch_mem(uint8_t *buffer, size_t n)
		{
		asserta(n > 0);
		mem = buffer;
		size = n;
		ptr = buffer;
		owner = false;
		}

	~scratch_mem()
		{
		if (owner)
			myfree(mem);
		}

	template<class t> t *get(uint n)
		{
		uint8_t *tmp_ptr = ptr;
		ptr += n*sizeof(t);
		asserta(size_t(ptr - mem) <= size);
		return (t *) tmp_ptr;
		}
	};
