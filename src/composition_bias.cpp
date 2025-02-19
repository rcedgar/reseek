#include "myutils.h"

#if 0
void CalcLocalBiasCorrection_3Di(const byte *Seq, uint L, int W, float Scale,
									vector<float> &BiasVec,
									vector<int8_t> &BiasVec8)
	{
	extern int8_t threedi_substmx[20][20];
	extern float threedi_letter_freqs[20];

	BiasVec.clear();
	BiasVec8.clear();
	BiasVec.reserve(L);
	BiasVec8.reserve(L);

	//const int W = BIAS_WINDOW;
	const int N = int(L);
	for (int i = 0; i < N; ++i)
		{
        const int minPos = max(0, (i - W/2));
        const int maxPos = min(N, (i + W/2));
        const int w = maxPos - minPos;
		int sumSubScores = 0;
		byte Letter_i = Seq[i]%20;
		const int8_t *row_i = threedi_substmx[Letter_i];
		for (int j = minPos; j < maxPos; ++j)
			{
			if (i == j)
				continue;
			byte Letter_j = Seq[j]%20;
			sumSubScores += row_i[Letter_j];
			}
		float deltaS_i = -float(sumSubScores)/w;
		for (int Letter = 0; Letter < 20; ++Letter)
			deltaS_i += threedi_letter_freqs[Letter]*row_i[Letter];
		float Bias = deltaS_i*Scale;
		int8_t Bias8 = int8_t(Bias < 0 ? Bias/4 - 0.5 : Bias/4 + 0.5);
		BiasVec.push_back(Bias);
		BiasVec8.push_back(Bias8);
		}
	}
#endif

void CalcLocalBiasCorrection_Mu(const byte *Seq, uint L, int W, float Scale,
									vector<float> &BiasVec,
									vector<int8_t> &BiasVec8)
	{
	BiasVec.clear();
	BiasVec.reserve(L);
	BiasVec8.reserve(L);
	BiasVec.resize(L);
	BiasVec8.resize(L);
	}