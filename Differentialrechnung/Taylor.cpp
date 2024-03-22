#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>

const double InitialX = 0.0;
const double InitialY = 1.0;
const double StepSize = 0.001;

double StepSize2 = StepSize * StepSize;
double StepSize3 = StepSize2 * StepSize;

double analitcySolution(double x) {
	return -1.0 / ( (x * x) * 0.25 + x * 0.5 - 1);
}

double firstDerivative(double x, double y) {
	return 0.5 * (1 + x) * y * y;
}

double secondDerivative(double x, double y) {
	return 0.5 * y * y + (1 + x) * y * firstDerivative(x, y);
}

double thirdDerivative(double x, double y) {
	return 2 * y * firstDerivative(x, y) + (1 + x) * firstDerivative(x, y) * firstDerivative(x, y) + (1 + x) * y * secondDerivative(x, y);
}

double taylorSolution(double x, double y, double stepSize) {
	return x + stepSize * firstDerivative(x, y) + 0.5 * StepSize2 * secondDerivative(x, y) + (1.0 / 6) * StepSize3 * thirdDerivative(x, y);
}

double sumTaylorSolution(double initalx, double initialy, double stepSize, double upperLimit) {
	double x = initalx;
	double y = initialy;
	for (int i = 1; i < upperLimit; ++i) {
		x += stepSize;
		y += taylorSolution(x, y, stepSize);

	}
	return y;
}

int main() {
	
}
