#pragma once

void CartToSpher(double x, double y, double z,
  double &r, double &theta, double &phi);

void CartToSpher(const vector<double> &Pt,
  double &r, double &theta, double &phi);

void SpherToCart(double r, double theta, double phi,
  double &x, double &y, double &z);
