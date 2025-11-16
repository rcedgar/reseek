#pragma once

template <class T>
class Binner
	{
private:
	const vector<T> *m_ptrValues;
	T m_MinValue;
	T m_MaxValue;
	uint32_t m_BinCount;
	vector<uint32_t> m_Bins;

public:
	~Binner() {}

// Private methods
private:
	void Count()
		{
		const vector<T> &Values = *m_ptrValues;
		for (size_t i = 0; i < Values.size(); ++i)
			{
			T Value = Values[i];
			uint32_t Bin = ValueToBin(Value);
			asserta(Bin < m_Bins.size());
			m_Bins[Bin] += 1;
			}
		}

// Constructors
public:
	Binner()
		{
		m_ptrValues = 0;
		m_MinValue = 0;
		m_MaxValue = 0;
		m_BinCount = 0;
		}

	Binner(const vector<T> &Values, uint32_t BinCount, T MinValue, T MaxValue) : Binner()
		{
		m_MinValue = MinValue;
		m_MaxValue = MaxValue;
		SetValues(Values);
		SetBinCount(BinCount);
		Count();
		}

	Binner(const vector<T> &Values, uint32_t BinCount, T MinValue) : Binner()
		{
		m_MinValue = MinValue;
		SetValues(Values);
		SetMaxFromValues();
		SetBinCount(BinCount);
		m_MinValue = MinValue;
		Count();
		}

	Binner(const vector<T> &Values, uint32_t BinCount) : Binner()
		{
		SetValues(Values);
		SetMinFromValues();
		SetMaxFromValues();
		SetBinCount(BinCount);
		Count();
		}

// Public methods
public:
	T GetMax() const { return m_MaxValue; }
	T GetMin() const { return m_MinValue; }

	void SetBinCount(uint32_t n)
		{
		m_Bins.clear();
		m_BinCount = n;
		m_Bins.resize(n);
		}

	void SetValues(const vector<T> &Values)
		{
		m_ptrValues = &Values;
		}

	void SetMin(T MinValue)
		{
		m_MinValue = MinValue;
		}

	void SetMax(T MaxValue)
		{
		m_MaxValue = MaxValue;
		}

	void SetMinFromValues()
		{
		const vector<T> &Values = *m_ptrValues;
		size_t ValueCount = Values.size();
		for (size_t i = 0; i < ValueCount; ++i)
			{
			T Value = Values[i];
			if (i == 0)
				m_MinValue = Value;
			else
				m_MinValue = min(Value, m_MinValue);
			}
		}

	void SetMaxFromValues()
		{
		const vector<T> &Values = *m_ptrValues;
		size_t ValueCount = Values.size();
		for (size_t i = 0; i < ValueCount; ++i)
			{
			T Value = Values[i];
			if (i == 0)
				m_MaxValue = Value;
			else
				m_MaxValue = max(Value, m_MaxValue);
			}
		}

	uint32_t ValueToBin(T Value) const
		{
		if (Value > m_MaxValue)
			Value = m_MaxValue;
		if (Value < m_MinValue)
			Value = m_MinValue;
		T Range = m_MaxValue - m_MinValue;
		asserta(Range > 0);
		T r = (Value - m_MinValue)/Range;
		asserta(r >= 0 && r <= 1);
		uint32_t Bin = uint32_t(r*(m_BinCount-1));
		asserta(Bin < m_BinCount);
		return Bin;
		}

	uint32_t GetCount(uint32_t Bin) const
		{
		asserta(Bin < SIZE(m_Bins));
		uint32_t n = m_Bins[Bin];
		return n;
		}

	T GetBinSize() const
		{
		asserta(m_MaxValue > m_MinValue);
		asserta(m_BinCount > 0);
		T BinSize = (m_MaxValue - m_MinValue)/m_BinCount;
		return BinSize;
		}

	T GetBinLo(uint32_t Bin) const
		{
		T BinSize = GetBinSize();
		T BinLo = m_MinValue + Bin*BinSize;
		return BinLo;
		}

	T GetBinHi(uint32_t Bin) const
		{
		T BinSize = GetBinSize();
		T BinHi = m_MinValue + (Bin+1)*BinSize;
		return BinHi;
		}

	T GetBinMid(uint32_t Bin) const
		{
		T BinSize = GetBinSize();
		T BinMid = m_MinValue + Bin*BinSize + BinSize/2;
		return BinMid;
		}

	uint32_t GetBinCount() const { return m_BinCount; }

	const vector<uint32_t> &GetBins() { return m_Bins; }

	void ValueToStr(T x, string &s) const
		{
		double d = double(x);
		Ps(s, "%.4g", d);
		}

	void ToTsv(const string &FileName) const
		{
		if (FileName == "")
			return;
		FILE *f = CreateStdioFile(FileName);
		ToTsv(f);
		CloseStdioFile(f);
		}

	void AccumToTsv(FILE *f, const string &Msg = "") const
		{
		vector<uint> AccumBins;
		GetAccumBins(AccumBins);
		asserta(SIZE(AccumBins) == m_BinCount);
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			T Mid = GetBinMid(Bin);
			string s;
			ValueToStr(Mid, s);
			uint n = AccumBins[Bin];
			if (n == 0)
				fprintf(f, "%u\t%s\t", Bin, s.c_str());
			else
				fprintf(f, "%u\t%s\t%u", Bin, s.c_str(), n);
			if (Msg != "")
				fprintf(f, "\t%s", Msg.c_str());
			fprintf(f, "\n");
			}
		}

	void AccumToTsvReverse(FILE *f, const string &Msg = "") const
		{
		vector<uint> AccumBins;
		GetAccumBinsReverse(AccumBins);
		asserta(SIZE(AccumBins) == m_BinCount);
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			T Mid = GetBinMid(Bin);
			string s;
			ValueToStr(Mid, s);
			uint n = AccumBins[Bin];
			if (n == 0)
				fprintf(f, "%u\t%s\t", Bin, s.c_str());
			else
				fprintf(f, "%u\t%s\t%u", Bin, s.c_str(), n);
			if (Msg != "")
				fprintf(f, "\t%s", Msg.c_str());
			fprintf(f, "\n");
			}
		}

	void ToTsv(FILE *f, const string &Msg = "") const
		{
		if (f == 0)
			return;
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			T Mid = GetBinMid(Bin);
			string s;
			ValueToStr(Mid, s);
			uint n = m_Bins[Bin];
			fprintf(f, "%u\t%s\t%u", Bin, s.c_str(), n);
			if (Msg != "")
				fprintf(f, "\t%s", Msg.c_str());
			fprintf(f, "\n");
			}
		}

	void ToHist(FILE *f) const
		{
		if (f == 0)
			return;
		const uint m = GetMaxCount();
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			T Mid = GetBinMid(Bin);
			string s;
			ValueToStr(Mid, s);
			uint n = m_Bins[Bin];
			string h;
			uint w = (n*80)/m;
			for (uint i = 0; i < w; ++i)
				h += '*';
			fprintf(f, "%10.10s  %10u  %s", s.c_str(), n, h.c_str());
			fprintf(f, "\n");
			}
		}

	void GetAccumBins(vector<uint> &AccumBins) const
		{
		AccumBins.clear();
		uint Sum = 0;
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			uint n = m_Bins[Bin];
			Sum += n;
			AccumBins.push_back(Sum);
			}
		}

	void GetAccumBinsReverse(vector<uint> &AccumBins) const
		{
		AccumBins.clear();
		AccumBins.resize(m_BinCount);
		uint Sum = 0;
		for (uint i = 0; i < m_BinCount; ++i)
			{
			uint Bin = m_BinCount - i - 1;
			uint n = m_Bins[Bin];
			Sum += n;
			AccumBins[Bin] = Sum;
			}
		}

	void DeleteAbove(T Value)
		{
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			T Mid = GetBinMid(Bin);
			if (Mid > Value)
				m_Bins[Bin] = 0;
			}
		}

	uint GetTotalCount() const
		{
		uint Sum = 0;
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			uint n = m_Bins[Bin];
			Sum += n;
			}
		return Sum;
		}

	uint GetMaxCount() const
		{
		uint Max = 0;
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			uint n = m_Bins[Bin];
			Max = max(Max, n);
			}
		return Max;
		}

// Cutoff which eliminates the largest TopN values.
	T GetCutoff_TopN(uint TopN) const
		{
		uint Total = GetTotalCount();
		uint Sum = 0;
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			uint n = m_Bins[Bin];
			Sum += n;
			if (Sum + TopN >= Total)
				return GetBinMid(Bin);
			}
		return m_MaxValue;
		}

// Cutoff C such that Value < C contains Fract*Total
	T GetCutoff_Fract(double Fract) const
		{
		asserta(Fract >= 0 && Fract <= 1.0);
		uint Total = GetTotalCount();
		uint Sum = 0;
		for (uint Bin = 0; Bin < m_BinCount; ++Bin)
			{
			uint n = m_Bins[Bin];
			Sum += n;
			if (Sum >= uint(Total*Fract))
				return GetBinMid(Bin);
			}
		return m_MaxValue;
		}
	};
