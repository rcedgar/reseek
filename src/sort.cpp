#if 0 // @@DELETE
#include "myutils.h"
#include "sort.h"
#include <algorithm>

template<class T> class CmpVecAsc
	{
public:
	const vector<T> *Vec;

	bool operator() (unsigned i, unsigned j) const
		{
		return (*Vec)[i] < (*Vec)[j];
		}
	};

template<class T> class CmpPtrAsc
	{
public:
	const T *Vec;

	bool operator() (unsigned i, unsigned j) const
		{
		return Vec[i] < Vec[j];
		}
	};

template<class T> class CmpVecDesc
	{
public:
	const vector<T> *Vec;

	bool operator() (unsigned i, unsigned j) const
		{
		return (*Vec)[i] > (*Vec)[j];
		}
	};

template<class T> class CmpPtrDesc
	{
public:
	const T *Vec;

	bool operator() (unsigned i, unsigned j) const
		{
		return Vec[i] > Vec[j];
		}
	};

void SortDescending(const vector<unsigned> &Values, unsigned N, unsigned *Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpVecDesc<unsigned> CV;
	CV.Vec = &Values;
	sort(Order, Order + N, CV);
	EndTimer(Sort);
	}

void SortAscending(const vector<unsigned> &Values, unsigned N, unsigned *Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpVecAsc<unsigned> CV;
	CV.Vec = &Values;
	sort(Order, Order + N, CV);
	EndTimer(Sort);
	}

void SortDescending(const uint16 *Values, unsigned N, unsigned *Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpPtrDesc<uint16> CV;
	CV.Vec = Values;
	sort(Order, Order + N, CV);
	EndTimer(Sort);
	}

void SortDescending(const unsigned *Values, unsigned N, unsigned *Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpPtrDesc<unsigned> CV;
	CV.Vec = Values;
	sort(Order, Order + N, CV);
	EndTimer(Sort);
	}

void SortAscending(const unsigned *Values, unsigned N, unsigned *Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpPtrAsc<unsigned> CV;
	CV.Vec = Values;
	sort(Order, Order + N, CV);
	EndTimer(Sort);
	}

void SortDescending(const unsigned *Values, unsigned N, vector<unsigned> &Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpPtrDesc<unsigned> CV;
	CV.Vec = Values;
	sort(Order.begin(), Order.end(), CV);
	EndTimer(Sort);
	}

void SortDescending(const vector<unsigned> &Values, vector<unsigned> &Order)
	{
	StartTimer(Sort);
	const unsigned N = SIZE(Values);
	Range(Order, N);
	CmpVecDesc<unsigned> CV;
	CV.Vec = &Values;
	sort(Order.begin(), Order.end(), CV);
	EndTimer(Sort);
	}

void SortAscending(const vector<unsigned> &Values, vector<unsigned> &Order)
	{
	StartTimer(Sort);
	const unsigned N = SIZE(Values);
	Range(Order, N);
	CmpVecAsc<unsigned> CV;
	CV.Vec = &Values;
	sort(Order.begin(), Order.end(), CV);
	EndTimer(Sort);
	}
///
void SortDescending(const float *Values, unsigned N, unsigned *Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpPtrDesc<float> CV;
	CV.Vec = Values;
	sort(Order, Order + N, CV);
	EndTimer(Sort);
	}

void SortAscending(const float *Values, unsigned N, unsigned *Order)
	{
	StartTimer(Sort);
	Range(Order, N);
	CmpPtrAsc<float> CV;
	CV.Vec = Values;
	sort(Order, Order + N, CV);
	EndTimer(Sort);
	}

void SortDescending(const vector<float> &Values, vector<unsigned> &Order)
	{
	StartTimer(Sort);
	const unsigned N = SIZE(Values);
	Range(Order, N);
	CmpVecDesc<float> CV;
	CV.Vec = &Values;
	sort(Order.begin(), Order.end(), CV);
	EndTimer(Sort);
	}

void SortAscending(const vector<float> &Values, vector<unsigned> &Order)
	{
	StartTimer(Sort);
	const unsigned N = SIZE(Values);
	Range(Order, N);
	CmpVecAsc<float> CV;
	CV.Vec = &Values;
	sort(Order.begin(), Order.end(), CV);
	EndTimer(Sort);
	}
#endif