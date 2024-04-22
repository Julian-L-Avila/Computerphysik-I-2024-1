#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>
#include <map>
#include <string>

const double kInitialTime      = 0.0;
const double kInitialVelocity  = 0.0;
const double kInitialPosition  = 1e3;
const double kMass             = 70.0;
const double kDeltaCoefficient = 0.8;
const double kAreaCrossSection = 0.6;

const double kGravityAcceleration = 9.81;
const double kAirDensity  = 1.293;

const double kDragCoefficient  = kAirDensity * kAreaCrossSection * kDeltaCoefficient / 2;
const double kTerminalVelocity = std::sqrt(kMass * kGravityAcceleration / kDragCoefficient);

const double kStepSize      = 1e-2;
const int kDesiredPrecision = 20;

double FirstDerivative(double t, double v) {
	return (kDragCoefficient * v * v / kMass) - kGravityAcceleration;
}

double SecondDerivative(double t, double v) {
	return (2 * kDragCoefficient * v * FirstDerivative(t, v));
}

double ThirdDerivative(double t, double v) {
	return 2 * kDragCoefficient * ((FirstDerivative(t, v) * FirstDerivative(t, v)) + v * SecondDerivative(t, v));
}

double AnalyticSolutionVelocity(double t) {
	double term1 = (kInitialVelocity - kTerminalVelocity) * exp(kGravityAcceleration * (t - kInitialTime) / kTerminalVelocity);
	double term2 = (kInitialVelocity + kTerminalVelocity) * exp(-kGravityAcceleration * (t - kInitialTime) / kTerminalVelocity);

	return -kTerminalVelocity * ((term1 + term2) / (term1 - term2));
}

double AnalyticSolutionPosition(double t) {
	double term1 = (kTerminalVelocity - kInitialTime) * exp((2 * kGravityAcceleration * t) / kTerminalVelocity);
	double term2 = (kTerminalVelocity + kInitialTime) * exp((2 * kGravityAcceleration * kInitialTime) / kTerminalVelocity);
	double term3 = (kTerminalVelocity + kInitialVelocity) * exp((2 * kGravityAcceleration * kInitialTime) / kTerminalVelocity);
	double term4 = (kTerminalVelocity - kInitialVelocity);

	return -(kTerminalVelocity * (kTerminalVelocity * log(std::abs(term1 + term2)) - kGravityAcceleration * t - kTerminalVelocity * log(std::abs(term3 + term4)))) / kGravityAcceleration;
}

double TaylorMethod(double previous_t, double previous_y, double step_size, double (*derivative)(double, double)) {
	double step_size_square = step_size * step_size;
	double step_size_cube = step_size_square * step_size;

	return previous_y + step_size * FirstDerivative(previous_t, previous_y) + 0.5 * step_size_square * SecondDerivative(previous_t, previous_y)
		+ (1.0 / 6.0) * step_size_cube * ThirdDerivative(previous_t, previous_y);
}

double EulerMethod(double previous_t, double previous_y, double step_size, double (*derivative)(double, double)) {
	return previous_y + step_size * derivative(previous_t, previous_y);
}

double HeunMethod(double previous_t, double previous_y, double step_size, double (*derivative)(double, double)) {
	double euler_y = EulerMethod(previous_t, previous_y, step_size, derivative);

	return previous_y + 0.5 * step_size * (derivative(previous_t, previous_y) + derivative(previous_t + step_size, euler_y));
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

double VoidFunction(double previous_t, double previous_y, double step_size, double (*derivative)(double, double)) {
	return AnalyticSolutionVelocity(previous_t);
}

void IterationLoop(const std::string& method_name, double initial_time, double initial_velocity, double initial_position, double step_size, double (*method)(double, double, double, double (*)(double, double))) {
	double previous_velocity = initial_velocity;
	double analytic_velocity = initial_velocity;
	double previous_position = initial_position;
	double velocity_error    = 0.0;

	std::string velocity_datafile_name = "dat-velocity-" + method_name + ".dat";
	std::ofstream velocity_datafile(velocity_datafile_name);

	std::cout << '\n' << "Chosen method has started ..." << std::endl;

	velocity_datafile << "# Data for " << method_name << '\n'
		<< "# Time (s)" << '\t' << "Analytic Velocity (ms^{-1})" << '\t' << "Approx. Velocity (ms^{-1})"  << '\t' << "Relative Percentage Error" << '\n'
		<< std::setprecision(kDesiredPrecision)
		<< initial_time << '\t' << analytic_velocity << '\t' << previous_velocity << '\t' << velocity_error << '\n';

	for (double time = initial_time + step_size; previous_position >= 0; time += step_size) {
		previous_position = kInitialPosition + AnalyticSolutionPosition(time);

		if (previous_position <= 0) {
			std::cout << "Trooper has landed at a time = " << time << '\n'
				<< '\n' << "Data has been saved." << std::endl;
			break;
		}

		analytic_velocity = AnalyticSolutionVelocity(time);
		previous_velocity = method(time, previous_velocity, step_size, FirstDerivative);
		velocity_error    = Error(analytic_velocity, previous_velocity);

		velocity_datafile << std::setprecision(kDesiredPrecision)
			<< time << '\t' << analytic_velocity << '\t' << previous_velocity << '\t' << velocity_error << '\n';
	}

	velocity_datafile.close();
}

int MenuOption () {
	std::cout << "This script solves the parachute problem's differential equation with different methods, please select one: " << '\n'
		<< "0. Analytic Solution" << '\n'
		<< "1. Taylor's Method" << '\n'
		<< "2. Euler's Method" << '\n'
		<< "3. Heun's Method" << '\n'
		<< "4. Runge-Kutta's Method" << '\n';

	int option;
	std::cout << "Press key: ";
	std::cin >> option;

	return option;
}

int main() {
	std::map<int, std::string> name_option = {
		{0, "Analytic"},
		{1, "Taylor"},
		{2, "Euler"},
		{3, "Heun"},
		{4, "Runge-Kutta"}
	};

	std::map<std::string, double (*)(double, double, double, double (*)(double, double))> method_lookup = {
			{"Analytic", VoidFunction},
			{"Taylor", TaylorMethod},
			{"Euler", EulerMethod},
			{"Heun", HeunMethod},
			{"Runge-Kutta", RungeKuttaMethod}
	};

	system("clear");

	int option = MenuOption();
	if (option < 0 || option > name_option.size()) {
		std::cout << "No method was chosen." << '\n' << "Exit." << std::endl;
		return 1;
	}

	std::string& method_name = name_option[option];

	IterationLoop(method_name, kInitialTime, kInitialVelocity, kInitialPosition, kStepSize, method_lookup[method_name]);

	std::string plotname = "plot-" + method_name + ".gnu";

	std::ofstream plotfile(plotname);

	if (option == 0) {
		plotfile << "set term pdfcairo" << '\n'
			<< "set output 'velocity-" << method_name << ".pdf'" << '\n'
			<< "set grid" << '\n'
			<< "set ylabel 'v [ms^{-1}]'" << '\n'
			<< "set xlabel 't [s]'" << '\n'
			<< "set auto xy" << '\n'
			<< "set tit 'Analytic Velocity'" << '\n'
			<< "p 'dat-velocity-" << method_name << ".dat' u 1:2 w l tit 'Analytic'" << '\n'
			<< "unset term" << '\n'
			<< "rep" << '\n'
			<< "pause -1" << '\n';

		plotfile.close();
	}

	plotfile << "set term pdfcairo" << '\n'
		<< "set output 'velocity-" << method_name << ".pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set ylabel 'v [ms^{-1}]'" << '\n'
		<< "set xlabel 't [s]'" << '\n'
		<< "set auto xy" << '\n'
		<< "set tit '" << method_name << " Method'" << '\n'
		<< "p 'dat-velocity-" << method_name << ".dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Approx.'" << '\n'
		<< "unset term" << '\n'
		<< "rep" << '\n'
		<< "pause -1" << '\n'
		<< "set term pdfcairo" << '\n'
		<< "set output 'error-velocity-" << method_name << ".pdf'" << '\n'
		<< "set ylabel 'Error %'" << '\n'
		<< "set logscale y" << '\n'
		<< "set auto xy" << '\n'
		<< "set tit 'Error with " << method_name << " Method'" << '\n'
		<< "p 'dat-velocity-" << method_name << ".dat' u 1:4 w l tit 'Error'" << '\n'
		<< "unset term" << '\n'
		<< "rep" << '\n'
		<< "pause -1" << '\n';

	plotfile.close();

	plotname = "gnuplot -p " + plotname;
	const char* plotname_c = plotname.c_str();
	system(plotname_c);

	return 0;
}
