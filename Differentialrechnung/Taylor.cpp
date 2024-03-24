#include<iostream>
#include<fstream>
#include<iomanip>
#include<cstdlib>

const double InitialX = 0.0; // Initial Condition x
const double InitialY = 1.0; // Inital Condition y
const double UpperLimit = 1.0; // Upper Limit Interval (b)
const double h1 = 0.001; // Stepsize 1 (h) for Taylor Method
const double h2 = 0.01; // Stepsize 2 (h) for Taylor Method
const double h3 = 0.1; // Steptize 3 (h) for Taylor Method
int desiredPrecision = 10;

// Differential Equation
// y' = (1 + x)y² / 2

double firstDerivative(double x, double y) {
	return 0.5 * (1 + x) * y * y;
}

double secondDerivative(double x, double y) {
	return 0.5 * y * y + (1 + x) * y * firstDerivative(x, y);
}

double thirdDerivative(double x, double y) {
	return 2 * y * firstDerivative(x, y) + (1 + x) * firstDerivative(x, y) * firstDerivative(x, y) + (1 + x) * y * secondDerivative(x, y);
}

double taylorSolution(double x, double y, double stepSize, int order) {
	double stepSizeSquare = stepSize * stepSize;
	double stepSizeCube = stepSizeSquare * stepSize;
	switch (order) {
		case 1:
			stepSizeSquare = 0.0;
			stepSizeCube = 0.0;
			break;
		case 2:
			stepSizeCube = 0.0;
			break;
		defult:
			break;
	};
		
	return y + stepSize * firstDerivative(x, y) + 0.5 * stepSizeSquare * secondDerivative(x, y) + (1.0 / 6) * stepSizeCube * thirdDerivative(x, y);
}

void  saveDataToFile(double initalx, double initialy, double stepSize, double upperLimit, int order, const std::string& filename) {
	double x = initalx;
	double yTaylorh = initialy;
	double N = (upperLimit - initalx) / stepSize;

	std::ofstream datafile(filename);

	for (int i = 1; i <= N + 1; ++i) {
		datafile << std::setprecision(desiredPrecision) << x << '\t' << yTaylorh << std::endl;
		yTaylorh = taylorSolution(x, yTaylorh, stepSize, order);
		x += stepSize;
	}
	datafile.close();
}

int main() {
	saveDataToFile(InitialX, InitialY, h1, UpperLimit, 3, "results-h0.001.dat");	
	saveDataToFile(InitialX, InitialY, h2, UpperLimit, 3, "results-h0.01.dat");
	saveDataToFile(InitialX, InitialY, h3, UpperLimit, 3, "results-h0.1.dat");
	saveDataToFile(InitialX, InitialY, h1, UpperLimit, 1, "results-o1h0.001.dat");	
	saveDataToFile(InitialX, InitialY, h2, UpperLimit, 1, "results-o1h0.01.dat");
	saveDataToFile(InitialX, InitialY, h3, UpperLimit, 1, "results-o1h0.1.dat");

	std::ofstream scriptFile("taylormethod.gnu");

	scriptFile << "set term pdfcairo\n"
		<< "set output 'Taylor-Method.pdf'\n"
		<< "set multiplot layout 1,2 tit 'Taylor Method'\n"
		<< "set grid \n"
		<< "set key l t\n"
		<< "set xrange [0:1]\n"
		<< "set auto y\n"
		<< "set tit 'Taylor Method with Order 3'\n"
		<< "f(x) = - 1 / ( x ** 2 /4 + x /2 - 1) \n"
		<< "set xlabel 'x'\n"
		<< "set ylabel 'y(x)'\n"
		<< "p f(x) tit 'Exact', 'results-h0.001.dat' u 1:2 w l dt 2 tit 'h = 0.001', 'results-h0.01.dat' u 1:2 w l dt 3 tit 'h = 0.01', 'results-h0.1.dat' u 1:2 w l dt 4 tit 'h = 0.1'\n"
		<< "set grid\n"
		<< "set key l t\n"
		<< "set xrange [0:1]\n"
		<< "set auto y\n"
		<< "set tit 'Taylor Method with Order 1'\n"
		<< "set xlabel 'x'\n"
		<< "set ylabel 'y(x)'\n"
		<< "p f(x) tit 'Exact', 'results-o1h0.001.dat' u 1:2 w l dt 2 tit 'h = 0.001', 'results-o1h0.01.dat' u 1:2 w l dt 3 tit 'h = 0.01', 'results-o1h0.1.dat' u 1:2 w l dt 4 tit 'h = 0.1'\n"
		<< "unset multiplot"
		<< std::endl;
	scriptFile.close();

	system("gnuplot -p 'taylormethod.gnu'");
	return 0;
}
