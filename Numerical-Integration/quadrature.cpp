#include<iostream>
#include<cmath>
#include<fstream>

const double LowerLimit = 0.0;
const double UpperLimit = 1.0;
const int MaxN = 10000;
const int NPlot = 3;
const double ExactIntegral = (199 * std::sin(20) + 20 * std::cos(20)) / 40;

// Function to integrate
double function(double x) {
	return 100 * std::pow(x, 2) * std::cos(20 * x);
}

// Approximate integral using rectangular quadrature
double rectangularIntegral(double lowerLimit, double upperLimit, int N) {
	double stepSize = (upperLimit - lowerLimit) / N;
	double sum = 0.0;

	for (int i = 0; i < N; ++i) {
		double center = lowerLimit + i * stepSize;
		sum += function(center);
	}
	return stepSize * sum;
}

double trapezoidIntegral(double lowerLimit, double upperLimit, int N) {
	double stepSize = (upperLimit - lowerLimit) / N;
	double sum = (function(lowerLimit) + function(upperLimit)) / 2;
	
	for (int i = 1; i < N; ++i) {
		double center = lowerLimit + i * stepSize;
		sum += function(center);
	}
	return stepSize * sum;
}

double simpsonIntegral(double lowerLimit, double upperLimit, int N) {
	double stepSize = (upperLimit - lowerLimit) / N;
	double h2 = stepSize / 2;
	double sum = function(lowerLimit) + function(upperLimit);

	for (int i =0; i < N; ++i) {
		double center = lowerLimit + i * stepSize;
		sum += 4 * function(center + h2) + 2 * function(center);
	}
	return h2 * sum / 3;
}

int main() {

	std::ofstream datafile("results.dat");

	for (int N = 1; N <= MaxN; ++N){
		double integralrect = rectangularIntegral(LowerLimit, UpperLimit, N);
		double integraltrap = trapezoidIntegral(LowerLimit, UpperLimit, N);
		double integralsimpson = simpsonIntegral(LowerLimit, UpperLimit, N);

		datafile << N << '\t' << ExactIntegral << '\t' << integralrect << '\t' << integraltrap << '\t' << integralsimpson << std::endl;
	}
	datafile.close();

	std::ofstream scriptFile1("numerical-integration.gnu");
	scriptFile1 << "set term pdfcairo\n"
		<< "set output 'Approximate-Integral.pdf'\n"
		<< "set grid\n"
		<< "set xlabel 'N'\n"
		<< "set ylabel 'Numerical Approx'\n"
		<< "set logscale x\n"
		<< "p 'results.dat' u 1:2 w l tit 'Exact Integral', 'results.dat' u 1:3 w l tit 'Rectangular Approx', 'results.dat' u 1:4 w l tit 'Trapezoid Approx.', 'results.dat' u 1:5 w l tit 'Simpson Approx.'\n";
	scriptFile1.close();

	std::ofstream datafile2("dataRectangular.dat");

	double hPlot = (UpperLimit - LowerLimit) / NPlot;

	for (double x = LowerLimit; x <= UpperLimit; x += hPlot) {
		datafile2 << x << '\t' << function(x) << std::endl;
	}

	datafile2.close();

	std::ofstream scriptFile2("rectangles.gnu");
	scriptFile2 << "set term pdfcairo\n"
		<< "set output 'Rectangles.pdf'\n"
		<< "set grid\n"
		<< "set xlabel 'x'\n"
		<< "set ylabel 'f(x)'\n"
		<< "f(x) = 100 * (x ** 2) * cos(20* x)\n"
		<< "p f(x) tit 'Function', 'dataRectangular.dat' u 1:2 w boxes tit 'Rectangles'\n";
	scriptFile2.close();

	system("gnuplot -p 'numerical-integration.gnu'");
	system("gnuplot -p 'rectangles.gnu'");
	std::cout << "Data saved and Ploted" << std::endl;


	return 0;
}
