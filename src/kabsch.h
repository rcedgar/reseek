#pragma once

/** Least-squares rigid superposition (TM-align Kabsch): maps x -> u*x + t ~ y. */
double Kabsch(
	const double * const *x,
	const double * const *y, int n,
	double t[3], double u[3][3]);
