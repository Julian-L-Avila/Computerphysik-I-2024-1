/* Solution to harmonic differential equation for "Fisica Computacional I"
 *
 * Made by:
 * Julian L. Avila - 20212107030
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <filesystem>

const double kInitialTime     = 0.0;
const double kFinalTime       = 10.0;
const double kInitialPosition = 5.0;
const double kInitialVelocity = 0.0;

double spring_constant          = 5.0;
double mass                     = 10.0;
double natural_frequency_square = spring_constant / mass;
double natural_frequency        = sqrt(natural_frequency_square);

const double kStepSize      = 1e-2;
const int kDesiredPrecision = 20;

double FirstVelocityDerivative(double time, double velocity, double position) {
	return -natural_frequency_square * position;
}

double FirstPositionDerivative(double time, double position, double velocity) {
	return velocity;
}

double AnalyticSolutionVelocity(double t) {
	return 0.0;
}

double AnalyticSolutionPosition(double t) {
	return 0.0;
}

double EulerMethod(double previous_t, double previous_x, double previous_y, double step_size, double (*Derivative)(double, double, double)) {
	return previous_x + step_size * Derivative(previous_t, previous_x,  previous_y);
}

void EulerLoop(const std::string& mass, double initial_time, double final_time, double initial_position, double initial_velocity, double step_size, double (*Derivativex)(double, double, double), double (*Derivativey)(double, double, double)) {
	double previous_x, previous_v;
	previous_x = initial_position;
	previous_v = initial_velocity;

	std::string path_file_name = "./Approx-Data/euler-" + mass + ".dat";

	std::ofstream datafile(path_file_name);

	datafile << "# Euler Data" << '\n'
		<< "#Time (s) \t Position (m) \t Velocity (ms^-1)" << '\n'
		<< std::setprecision(kDesiredPrecision)
		<< initial_time << '\t' << initial_position << '\t' << initial_velocity << '\n';

	for (double t = initial_time; t <= final_time; t += step_size) {
	double x = EulerMethod(t, previous_x, previous_v, step_size, Derivativex);
	double v = EulerMethod(t, previous_v, previous_x, step_size, Derivativey);

	datafile << t << '\t' << x << '\t' << v << '\n';
	previous_x = x;
	previous_v = v;
	}

	datafile.close();
}

/* double HeunMethod(double previous_t, double previous_y, double step_size, double (*derivative)(double, double)) { */
/* 	double euler_y = EulerMethod(previous_t, previous_y, step_size, derivative); */

/* 	return previous_y + 0.5 * step_size * (derivative(previous_t, previous_y) + derivative(previous_t + step_size, euler_y)); */
/* } */

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

int main() {
	EulerLoop("200g", kInitialTime, kFinalTime, kInitialPosition, kInitialVelocity, kStepSize, FirstPositionDerivative, FirstVelocityDerivative);
	return 0;
}
