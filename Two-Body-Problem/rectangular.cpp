#include<iostream>
#include<cmath>

const double R = 5.20336301;
const double a = 5.20336301;
const double e = 0.04839266;
const double theta0 = 0.0;
const double m0 = 1.0;
const double m1 = 0.9;
const double G = 1.0;

const double b = a - e;
const double mu = (m0 * m1) / (m0 + m1);
const double L = R * m0;
const double C = G * m0 * m1;
const double E = (e * e - 1) * mu * C * C / (2 * L * L);
