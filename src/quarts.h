#ifndef quarts_h
#define quarts_h

struct Quarts
	{
	unsigned N;
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

	void WriteMe(FILE *f) const
		{
		if (f == 0)
			return;
		fprintf(f, "Min=%u", Min);
		fprintf(f, ", LoQ=%u", LoQ);
		fprintf(f, ", Med=%u", Med);
		fprintf(f, ", HiQ=%u", HiQ);
		fprintf(f, ", Max=%u", Max);
		fprintf(f, ", Avg=%.1f", Avg);
		fprintf(f, "\n");
		}

	void LogMe() const
		{
		Log("N=%u", N);
		Log(", Min=%u", Min);
		Log(", LoQ=%u", LoQ);
		Log(", Med=%u", Med);
		Log(", HiQ=%u", HiQ);
		Log(", Max=%u", Max);
		Log(", Avg=%3g", Avg);
		Log("\n");
		}
	};

struct QuartsFloat
	{
	uint N;
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

	void WriteMe(FILE *f) const
		{
		if (f == 0)
			return;
		fprintf(f, "Min=%.3g", Min);
		fprintf(f, ", LoQ=%.3g", LoQ);
		fprintf(f, ", Med=%.3g", Med);
		fprintf(f, ", HiQ=%.3g", HiQ);
		fprintf(f, ", Max=%.3g", Max);
		fprintf(f, ", Avg=%.3g", Avg);
		fprintf(f, "\n");
		}

	void ProgressLogMe() const
		{
		ProgressLog("N=%u", N);
		ProgressLog(", Min=%.3g", Min);
		ProgressLog(", LoQ=%.3g", LoQ);
		ProgressLog(", Med=%.3g", Med);
		ProgressLog(", HiQ=%.3g", HiQ);
		ProgressLog(", Max=%.3g", Max);
		ProgressLog(", Avg=%.3g", Avg);
		ProgressLog(", StdDev=%.3g", StdDev);
		ProgressLog("\n");
		}

	void LogMe() const
		{
		Log("N=u", N);
		Log(", Min=%.3g", Min);
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
void GetQuarts(const unsigned *v, uint N, Quarts &Q);
void GetQuartsFloat(const vector<float> &v, QuartsFloat &Q);

#endif // quarts_h
