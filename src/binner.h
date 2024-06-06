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
	Binner()
		{
		m_ptrValues = 0;
		m_MinValue = 0;
		m_MaxValue = 0;
		m_BinCount = 0;
		}
	~Binner() {}

public:
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

	uint32_t ValueToBin(float Value) const
		{
		if (Value > m_MaxValue)
			Value = m_MaxValue;
		if (Value < m_MinValue)
			Value = m_MinValue;
		float Range = m_MaxValue - m_MinValue;
		asserta(Range > 0);
		float r = (Value - m_MinValue)/Range;
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

	T GetBinMid(uint32_t Bin) const
		{
		T BinSize = GetBinSize();
		T BinMid = m_MinValue + Bin*BinSize + BinSize/2;
		return BinMid;
		}

	uint32_t GetBinCount() const { return m_BinCount; }

	const vector<uint32_t> &GetBins() { return m_Bins; }

	Binner(const vector<T> &Values, uint32_t BinCount, T MinValue) : Binner()
		{
		m_MinValue = MinValue;
		SetValues(Values);
		SetMaxFromValues();
		SetBinCount(BinCount);
		m_MinValue = MinValue;
		for (size_t i = 0; i < Values.size(); ++i)
			{
			T Value = Values[i];
			uint32_t Bin = ValueToBin(Value);
			asserta(Bin < m_Bins.size());
			m_Bins[Bin] += 1;
			}
		}
	};
