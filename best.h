#pragma once

static inline void Best5(float v, float w, float x, float y, float z,
  char cv, char cw, char cx, char cy, char cz,
  float &BestScore, char &BestChar)
	{
	if (v >= w && v >= x && v >= y && v >= z)
		{
		BestScore = v;
		BestChar = cv;
		}
	else if (w >= v && w >= x && w >= y && w >= z)
		{
		BestScore = w;
		BestChar = cw;
		}
	else if (x >= v && x >= w && x >= y && x >= z)
		{
		BestScore = x;
		BestChar = cx;
		}
	else if (y >= v && y >= w && y >= x && y >= z)
		{
		BestScore = y;
		BestChar = cy;
		}
	else
		{
		BestScore = z;
		BestChar = cz;
		}
	}

static inline void Best4(float w, float x, float y, float z,
  char cw, char cx, char cy, char cz,
  float &BestScore, char &BestChar)
	{
	if (w >= x && w >= y && w >= z)
		{
		BestScore = w;
		BestChar = cw;
		}
	else if (x >= w && x >= y && x >= z)
		{
		BestScore = x;
		BestChar = cx;
		}
	else if (y >= w && y >= x && y >= z)
		{
		BestScore = y;
		BestChar = cy;
		}
	else
		{
		BestScore = z;
		BestChar = cz;
		}
	}

static inline void Best3(float x, float y, float z, char cx, char cy, char cz,
  float &BestScore, char &BestChar)
	{
	if (x >= y && x >= z)
		{
		BestScore = x;
		BestChar = cx;
		}
	else if (y >= x && y >= z)
		{
		BestScore = y;
		BestChar = cy;
		}
	else
		{
		BestScore = z;
		BestChar = cz;
		}
	}

static inline void Best2(float x, float y, char cx, char cy,
  float &BestScore, char &BestChar)
	{
	if (x >= y)
		{
		BestScore = x;
		BestChar = cx;
		}
	else
		{
		BestScore = y;
		BestChar = cy;
		}
	}
