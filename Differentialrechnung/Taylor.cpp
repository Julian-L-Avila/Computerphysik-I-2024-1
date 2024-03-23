#include<iostream>
#include<fstream>
#include<cstdlib>
#include<cmath>

const double InitialX = 0.0; // Initial Condition x
const double InitialY = 1.0; // Inital Condition y
const double StepSize = 0.001; // Stepsize 1 (h) for Taylor Method
const double StepSize2 = 0.01; // Stepsize 2 (h) for Taylor Method

// Analytic solution of differenital equation
// y' = (1 + x)y² / 2
double analyticSolution(double x) {
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
	double stepSizeSquare = stepSize * stepSize;
	double stepSizeCube = stepSizeSquare * stepSize;
	return x + stepSize * firstDerivative(x, y) + 0.5 * stepSizeSquare * secondDerivative(x, y) + (1.0 / 6) * stepSizeCube * thirdDerivative(x, y);
}

void  saveDataToFile(double initalx, double initialy, double stepSize, double upperLimit, double lowerLimit, const std::string& filename) {
	double x = initalx;
	double yTaylorh = initialy;
	double y = analyticSolution(x);
	double N = (upperLimit - lowerLimit) / stepSize;

	std::ofstream datafile(filename);

	for (int i = 1; i < N; ++i) {
		datafile << x << '\t' << y << '\t' << yTaylorh << std::endl;
		yTaylorh += taylorSolution(x, yTaylorh, stepSize);
		y = analyticSolution(x);
		x += stepSize;
	}
	datafile.close();
}

int main() {
	saveDataToFile(InitialX, InitialY, StepSize, 1.0, 0.0, "results-h0.001.dat");	

	return 0;
}
