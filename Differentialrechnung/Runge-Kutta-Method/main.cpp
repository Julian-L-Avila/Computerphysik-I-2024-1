#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>
#include <string>

const double kInitialt = 1.0;
const double kInitialy = 1.0;
const double kFinalt   = 5.0;

const double kStepSize         = 1e-1;
const double kDesiredPrecision = 20;

double FirstDerivative(double t, double y) {
	return t * t / y;
}

double AnalyticSolution(double t) {
	return std::sqrt(2 * (t * t * t / 3.0 + 1.0 / 6.0));
}

double RungeKuttaMethod(double previous_t, double previous_y, double step_size, double(*derivative)(double, double)) {
	double k1, k2, k3, k4;

	k1 = derivative(previous_t, previous_y);
	k2 = derivative(previous_t + 0.5 * step_size, previous_y + 0.5 * step_size * k1);
	k3 = derivative(previous_t + 0.5 * step_size, previous_y + 0.5 * step_size * k2);
	k4 = derivative(previous_t + step_size, previous_y + step_size * k3);

	return previous_y + step_size * (k1 + 2 * k2 + 2 * k3 + k4) / 6.0;
}

double Error(double real_value, double value) {
	return std::abs((real_value - value) / real_value) * 100;
}

void IterationLoop(const std::string& method_name, double initial_t, double final_t, double initial_y, double step_size, double (*method)(double, double, double, double (*)(double, double))) {
	double previous_y = initial_y;
	double analytic_y = initial_y;
	double y_error    = 0.0;

	std::ofstream datafile("dat-runge-kutta.dat");

	datafile << "# Time (s)" << '\t' << "Analytic" << '\t' << "Approx"  << '\t' << "Relative Percentage Error" << '\n'
		<< std::setprecision(kDesiredPrecision)
		<< initial_t << '\t' << analytic_y << '\t' << previous_y << '\t' << y_error << '\n';

	for (double t = initial_t + step_size; t <= final_t; t += step_size) {

		analytic_y = AnalyticSolution(t);
		previous_y = method(t, previous_y, step_size, FirstDerivative);
		y_error    = Error(analytic_y, previous_y);

		datafile << std::setprecision(kDesiredPrecision)
			<< t << '\t' << analytic_y << '\t' << previous_y << '\t' << y_error << '\n';
	}

	datafile.close();
}

int main() {
	system("clear");

	IterationLoop("runge-kutta", kInitialt, kFinalt, kInitialy, kStepSize, RungeKuttaMethod);

	std::string plotname = "plot-runge-kutta.gnu";

	std::ofstream plotfile(plotname);

	plotfile << "set term pdfcairo" << '\n'
		<< "set output 'runge-kutta.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set ylabel 'v [ms^{-1}]'" << '\n'
		<< "set xlabel 't [s]'" << '\n'
		<< "set auto xy" << '\n'
		<< "p 'dat-runge-kutta.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Approx.'" << '\n'
		<< "set output 'error-runge-kutta.pdf'" << '\n'
		<< "set ylabel 'Error %'" << '\n'
		<< "set logscale y" << '\n'
		<< "set auto xy" << '\n'
		<< "p 'dat-runge-kutta.dat' u 1:4 w l tit 'Error'" << std::endl;

	plotfile.close();

	plotname = "gnuplot -p " + plotname;
	const char* plotname_c = plotname.c_str();
	system(plotname_c);

	return 0;
}
