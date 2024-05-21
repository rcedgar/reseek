#ifndef sort_h
#define sort_h

#include "timing.h"

void SortDescending(const unsigned *Values, unsigned N, unsigned *Order);
void SortAsccending(const unsigned *Values, unsigned N, unsigned *Order);
void SortDescending(const unsigned *Values, unsigned N, vector<unsigned> &Order);
void SortAscending(const unsigned *Values, unsigned N, vector<unsigned> &Order);
void SortAscending(const unsigned *Values, unsigned N, unsigned *Order);
void SortDescending(const vector<unsigned> &Values, vector<unsigned> &Order);
void SortAscending(const vector<unsigned> &Values, vector<unsigned> &Order);
void SortDescending(const float *Values, unsigned N, unsigned *Order);
void SortAscending(const float *Values, unsigned N, unsigned *Order);
void SortDescending(const vector<float> &Values, vector<unsigned> &Order);
void SortAscending(const vector<float> &Values, vector<unsigned> &Order);

inline void Range(vector<unsigned> &v, unsigned N)
	{
	v.clear();
	v.reserve(N);
	for (unsigned i = 0; i < N; ++i)
		v.push_back(i);
	}

inline void Range(unsigned *v, unsigned N)
	{
	for (unsigned i = 0; i < N; ++i)
		v[i] = i;
	}

template<class T, bool Desc> void QuickSortInPlaceRecurse(T *Values, int left, int right)
	{
	int i = left;
	int j = right;
	int Mid = (left + right)/2;
	T pivot = Values[Mid];

	while (i <= j)
		{
		if (Desc)
			{
			while (Values[i] > pivot)
				i++;
			while (Values[j] < pivot)
				j--;
			}
		else
			{
			while (Values[i] < pivot)
				i++;
			while (Values[j] > pivot)
				j--;
			}

		if (i <= j)
			{
			swap(Values[i], Values[j]);
			i++;
			j--;
			}
		}

	if (left < j)
		QuickSortInPlaceRecurse<T, Desc>(Values, left, j);

	if (i < right)
		QuickSortInPlaceRecurse<T, Desc>(Values, i, right);
	}

template<class T, bool Desc> void QuickSortOrderRecurse(const T *Values, int left, int right, unsigned *Order)
	{
	int i = left;
	int j = right;
	int Mid = (left + right)/2;
	T pivot = Values[Order[Mid]];

	while (i <= j)
		{
		if (Desc)
			{
			while (Values[Order[i]] > pivot)
				i++;
			while (Values[Order[j]] < pivot)
				j--;
			}
		else
			{
			while (Values[Order[i]] < pivot)
				i++;
			while (Values[Order[j]] > pivot)
				j--;
			}

		if (i <= j)
			{
			swap(Order[i], Order[j]);
			i++;
			j--;
			}
		}

	if (left < j)
		QuickSortOrderRecurse<T, Desc>(Values, left, j, Order);

	if (i < right)
		QuickSortOrderRecurse<T, Desc>(Values, i, right, Order);
	}

template<class T> void QuickSortInPlace(T *Values, unsigned N)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortInPlace);
	QuickSortInPlaceRecurse<T, false>(Values, 0, int(N-1));
	EndTimer(QuickSortInPlace);
	}

template<class T> void QuickSortInPlaceDesc(T *Values, unsigned N)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortInPlaceDesc);
	QuickSortInPlaceRecurse<T, true>(Values, 0, int(N-1));
	EndTimer(QuickSortInPlaceDesc);
	}

template<class T> void QuickSortOrder(const T *Values, unsigned N, unsigned *Order)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortOrder);
	Range(Order, N);
	QuickSortOrderRecurse<T, false>(Values, 0, int(N-1), Order);
	EndTimer(QuickSortOrder);
	}

template<class T> void QuickSortOrderDesc(const T *Values, unsigned N, unsigned *Order)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortOrderDesc);
	Range(Order, N);
	QuickSortOrderRecurse<T, true>(Values, 0, int(N-1), Order);
	EndTimer(QuickSortOrderDesc);
	}

template<class T> void QuickSortSubset(const T *Values, unsigned N, unsigned *Subset)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortSubset);
	QuickSortOrderRecurse<T, false>(Values, 0, int(N-1), Subset);
	EndTimer(QuickSortSubset);
	}

template<class T> void QuickSortSubsetDesc(const T *Values, unsigned N, unsigned *Subset)
	{
	if (N == 0)
		return;
	asserta(N < INT_MAX);

	StartTimer(QuickSortSubsetDesc);
	QuickSortOrderRecurse<T, true>(Values, 0, int(N-1), Subset);
	EndTimer(QuickSortSubsetDesc);
	}

#endif // sort_h
