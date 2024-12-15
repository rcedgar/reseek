#pragma once
#if 0
// Modified from original here:
// https://learn.microsoft.com/en-us/cpp/standard-library/allocators?view=msvc-170

#include <stdlib.h> //size_t, malloc, free
#include <new> // bad_alloc, bad_array_new_length
#include <memory>

template <class T>
struct Myallocator
	{
	typedef T value_type;
	Myallocator() noexcept {} //default ctor not required by C++ Standard Library

	// A converting copy constructor:
	template<class U> Myallocator(const Myallocator<U> &) noexcept {}
	template<class U> bool operator==(const Myallocator<U> &) const noexcept
		{
		return true;
		}
	template<class U> bool operator!=(const Myallocator<U> &) const noexcept
		{
		return false;
		}
	T *allocate(const size_t n) const;
	void deallocate(T *const p, size_t) const noexcept;
	};

template <class T>
T *Myallocator<T>::allocate(const size_t n) const
	{
	if(n == 0)
		{
		return nullptr;
		}
	if(n > static_cast<size_t>(-1) / sizeof(T))
		{
		throw std::bad_array_new_length();
		}
	void *const pv = malloc(n * sizeof(T));
	if(!pv) { throw std::bad_alloc(); }
	return static_cast<T *>(pv);
	}

template<class T>
void Myallocator<T>::deallocate(T *const p, size_t) const noexcept
	{
	free(p);
	}
#endif // 0