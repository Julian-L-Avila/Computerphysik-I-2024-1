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

double EulerMethod(double previous_t, double previous_x, double previous_y, double (*Derivative)(double, double, double)) {
	return previous_x + kStepSize * Derivative(previous_t, previous_x,  previous_y);
}

void EulerLoop(const std::string& mass, double initial_time, double final_time, double initial_position, double initial_velocity, double (*Derivativex)(double, double, double), double (*Derivativey)(double, double, double)) {
	double previous_x, previous_v;
	previous_x = initial_position;
	previous_v = initial_velocity;

	std::string path_file_name = "./Approx-Data/euler-" + mass + ".dat";

	std::ofstream datafile(path_file_name);

	datafile << "# Euler Data" << '\n'
		<< "#Time (s) \t Position (m) \t Velocity (ms^-1)" << '\n'
		<< std::setprecision(kDesiredPrecision) << std::fixed
		<< initial_time << '\t' << initial_position << '\t' << initial_velocity << '\n';

	for (double t = initial_time + kStepSize; t <= final_time; t += kStepSize) {
	double x = EulerMethod(t, previous_x, previous_v, Derivativex);
	double v = EulerMethod(t, previous_v, previous_x, Derivativey);

	datafile << t << '\t' << x << '\t' << v << '\n';
	previous_x = x;
	previous_v = v;
	}

	datafile.close();
}

void HeunLoop(const std::string& mass, double initial_time, double final_time, double initial_position, double initial_velocity, double (*Derivativex)(double, double, double), double (*Derivativey)(double, double, double)) {

	double previous_x, previous_v;
	previous_x = initial_position;
	previous_v = initial_velocity;

	std::string path_file_name = "./Approx-Data/heun-" + mass + ".dat";

	std::ofstream datafile(path_file_name);

	datafile << "# Heun Data" << '\n'
		<< "#Time (s) \t Position (m) \t Velocity (ms^-1)" << '\n'
		<< std::setprecision(kDesiredPrecision) << std::fixed
		<< initial_time << '\t' << initial_position << '\t' << initial_velocity << '\n';

	for (double t = initial_time + kStepSize; t <= final_time; t += kStepSize) {
		double euler_x = EulerMethod(t, previous_x, previous_v, Derivativex);
		double euler_y = EulerMethod(t, previous_v, previous_x, Derivativey);

		double x = previous_x + 0.5 * kStepSize * (Derivativex(t, previous_x, previous_v) + Derivativex(t + kStepSize, euler_x, euler_y));
		double v = previous_v + 0.5 * kStepSize * (Derivativey(t, previous_v, previous_x) + Derivativey(t, euler_y, euler_x));

		datafile << t << '\t' << x << '\t' << v << '\n';
		previous_x = x;
		previous_v = v;
	}

	datafile.close();
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

int main() {
	EulerLoop("200g", kInitialTime, kFinalTime, kInitialPosition, kInitialVelocity, FirstPositionDerivative, FirstVelocityDerivative);
	HeunLoop("200g", kInitialTime, kFinalTime, kInitialPosition, kInitialVelocity, FirstPositionDerivative, FirstVelocityDerivative);
	return 0;
}
