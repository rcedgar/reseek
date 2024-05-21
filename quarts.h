#ifndef quarts_h
#define quarts_h

struct Quarts
	{
	unsigned Min;
	unsigned LoQ;
	unsigned Med;
	unsigned HiQ;
	unsigned Max;
	unsigned Total;
	float Avg;

	Quarts()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		}

	void Clear()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		}

	void LogMe() const
		{
		Log("Min=%u", Min);
		Log(", LoQ=%u", LoQ);
		Log(", Med=%u", Med);
		Log(", HiQ=%u", HiQ);
		Log(", Max=%u", Max);
		Log(", Avg=%.1f", Avg);
		Log("\n");
		}
	};

struct QuartsFloat
	{
	float Min;
	float LoQ;
	float Med;
	float HiQ;
	float Max;
	float Total;
	float Avg;
	float StdDev;

	QuartsFloat()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		StdDev = 0;
		}

	void Clear()
		{
		Min = 0;
		LoQ = 0;
		Med = 0;
		HiQ = 0;
		Max = 0;
		Total = 0;
		Avg = 0;
		StdDev = 0;
		}

	void LogMe() const
		{
		Log("Min=%.3g", Min);
		Log(", LoQ=%.3g", LoQ);
		Log(", Med=%.3g", Med);
		Log(", HiQ=%.3g", HiQ);
		Log(", Max=%.3g", Max);
		Log(", Avg=%.3g", Avg);
		Log(", StdDev=%.3g", StdDev);
		Log("\n");
		}
	};

void GetQuarts(const vector<unsigned> &v, Quarts &Q);
void GetQuartsFloat(const vector<float> &v, QuartsFloat &Q);

#endif // quarts_h
