#include "myutils.h"

// Utility functions derived from foldseek source code
// https://github.com/steineggerlab/foldseek
// GPL3 license
// van Kempen M, Kim S, Tumescheit C, Mirdita M, Lee J, Gilchrist CLM, Soding J,
// and Steinegger M. Fast and accurate protein structure search with Foldseek.
// Nature Biotechnology (2023)
// doi:10.1038/s41587-023-01773-0

void LogCoords16(const char *mem, uint chainLength)
	{
	uint bufferSize = 3*chainLength;
	const char *data = mem;
	int32_t start;
	memcpy(&start, data, sizeof(int32_t));
	data += sizeof(int32_t);
	Log("X[0] = %.1f\n", start/1000.0);

	int32_t diffSum = 0;
	int16_t intDiff = 0;
	for(size_t i = 1; i < chainLength; ++i)
		{
		memcpy(&intDiff, data, sizeof(int16_t));
		data += sizeof(int16_t);
		diffSum += intDiff;
		//buffer[i] = (start + diffSum) / 1000.0f;
		Log("X[%d] = %.1f\n", i, (start + diffSum) / 1000.0);
		}

	diffSum = 0;
	memcpy(&start, data, sizeof(int32_t));
	data += sizeof(int32_t);
	//buffer[chainLength] = start / 1000.0f;
	Log("Y[0] = %.1f\n", start/1000.0);

	for(size_t i = chainLength + 1; i < 2 * chainLength; ++i)
		{
		memcpy(&intDiff, data, sizeof(int16_t));
		data += sizeof(int16_t);
		diffSum += intDiff;
		//buffer[i] = (start + diffSum) / 1000.0f;
		Log("Y[%d] = %.1f\n", i-chainLength, (start + diffSum) / 1000.0);
		}

	diffSum = 0;
	memcpy(&start, data, sizeof(int32_t));
	data += sizeof(int32_t);
	//buffer[2 * chainLength] = start / 1000.0f;
	Log("Z[0] = %.1f\n", start/1000.0);

	for(size_t i = 2 * chainLength + 1; i < 3 * chainLength; ++i)
		{
		memcpy(&intDiff, data, sizeof(int16_t));
		data += sizeof(int16_t);
		diffSum += intDiff;
		//buffer[i] = (start + diffSum) / 1000.0f;
		Log("Z[%d] = %.1f\n", i-2*chainLength, (start + diffSum) / 1000.0);
		}
	}

// Convert database buffer mem to float coordinates x,y,z
// float* Coordinate16::read(const char* mem, size_t chainLength, size_t entryLength)
// src/commons/Coordinate16.h
// Caller must free the returned float array.
float *GetCoordsFromMem(const char *mem, uint chainLength, uint entryLength)
	{
	if(entryLength >= (chainLength * 3) * sizeof(float))
		return (float *)mem;

	uint bufferSize = 3*chainLength;
	float *buffer = myalloc(float, bufferSize);
	const char *data = mem;
	int32_t start;
	memcpy(&start, data, sizeof(int32_t));
	data += sizeof(int32_t);

	int32_t diffSum = 0;
	buffer[0] = start / 1000.0f;
	int16_t intDiff = 0;
	for(size_t i = 1; i < chainLength; ++i)
		{
		memcpy(&intDiff, data, sizeof(int16_t));
		data += sizeof(int16_t);
		diffSum += intDiff;
		buffer[i] = (start + diffSum) / 1000.0f;
		}

	diffSum = 0;
	memcpy(&start, data, sizeof(int32_t));
	data += sizeof(int32_t);
	buffer[chainLength] = start / 1000.0f;
	for(size_t i = chainLength + 1; i < 2 * chainLength; ++i)
		{
		memcpy(&intDiff, data, sizeof(int16_t));
		data += sizeof(int16_t);
		diffSum += intDiff;
		buffer[i] = (start + diffSum) / 1000.0f;
		}

	diffSum = 0;
	memcpy(&start, data, sizeof(int32_t));
	data += sizeof(int32_t);
	buffer[2 * chainLength] = start / 1000.0f;
	for(size_t i = 2 * chainLength + 1; i < 3 * chainLength; ++i)
		{
		memcpy(&intDiff, data, sizeof(int16_t));
		data += sizeof(int16_t);
		diffSum += intDiff;
		buffer[i] = (start + diffSum) / 1000.0f;
		}
	return buffer;
	}

char *CoordsToMem(const float *coords, int L, uint &membytes)
	{
	membytes = 3*sizeof(int32_t) + 3*(L-1)*sizeof(int16_t);
	char *buffer = myalloc(char, membytes);
	int16_t *ptr = (int16_t *) buffer;

	for (int Axis = 0; Axis < 3; ++Axis)
		{
		int16_t DiffSum = 0;
		int32_t first = (int32_t)(coords[L*Axis] * 1000);
		memcpy(ptr, &first, sizeof(int32_t));
		ptr += 2;

		int32_t last = first;
		for (size_t i = 1; i < L; ++i)
			{
			int32_t curr = (int32_t)(coords[L*Axis+i] * 1000);
			int32_t diff32 = int32_t(curr) - int32_t(last);
			int16_t diff16 = int16_t(diff32);
			if (int32_t(diff16) == diff32)
				*ptr++ = diff16;
			else
				{
			// overflow
				free(buffer);
				return 0;
				}
			last = curr;
			}
		}
	size_t bytes2 = (char *) ptr - buffer;
    return buffer;
	}
