/* Solution to harmonic differential equation for "Fisica Computacional I"
 *
 * Made by:
 * Julian L. Avila - 20212107030
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

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

double FirstVelocityDerivative(double velocity, double position) {
	return -natural_frequency_square * position;
}

double FirstPositionDerivative(double position, double velocity) {
	return velocity;
}

double AnalyticSolutionVelocity(double t) {
	return 0.0;
}

double AnalyticSolutionPosition(double t) {
	return 0.0;
}

double EulerMethod(double previous_x, double previous_y, double (*Derivative)(double, double)) {
	return previous_x + kStepSize * Derivative(previous_x,  previous_y);
}

void EulerLoop(const std::string& mass, double initial_time, double final_time, double initial_position, double initial_velocity) {
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
	double x = EulerMethod(previous_x, previous_v, FirstPositionDerivative);
	double v = EulerMethod(previous_v, previous_x, FirstVelocityDerivative);

	datafile << t << '\t' << x << '\t' << v << '\n';
	previous_x = x;
	previous_v = v;
	}

	datafile.close();
}

void HeunLoop(const std::string& mass, double initial_time, double final_time, double initial_position, double initial_velocity) {

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
		double euler_x = EulerMethod(previous_x, previous_v, FirstPositionDerivative);
		double euler_y = EulerMethod(previous_v, previous_x, FirstPositionDerivative);

		double x = previous_x + 0.5 * kStepSize * (FirstPositionDerivative(previous_x, previous_v) + FirstPositionDerivative(euler_x, euler_y));
		double v = previous_v + 0.5 * kStepSize * (FirstVelocityDerivative(previous_v, previous_x) + FirstVelocityDerivative(euler_y, euler_x));

		datafile << t << '\t' << x << '\t' << v << '\n';
		previous_x = x;
		previous_v = v;
	}

	datafile.close();
}

void RungeKuttaLoop(const std::string& mass, double initial_time, double final_time, double initial_position, double initial_velocity) {

	double previous_x, previous_v;
	previous_x = initial_position;
	previous_v = initial_velocity;

	std::string path_file_name = "./Approx-Data/rungekutta-" + mass + ".dat";

	std::ofstream datafile(path_file_name);

	datafile << "# RungeKutta Data" << '\n'
		<< "#Time (s) \t Position (m) \t Velocity (ms^-1)" << '\n'
		<< std::setprecision(kDesiredPrecision) << std::fixed
		<< initial_time << '\t' << initial_position << '\t' << initial_velocity << '\n';

	for (double t = initial_time + kStepSize; t <= final_time; t += kStepSize) {
		long double k1 = FirstPositionDerivative(previous_x, previous_v);
		long double m1 = FirstVelocityDerivative(previous_v, previous_x);
		long double k2 = FirstPositionDerivative(previous_x + 0.5 * k1 * kStepSize, previous_v + 0.5 * m1 * kStepSize);
		long double m2 = FirstVelocityDerivative(previous_v + 0.5 * m1 * kStepSize, previous_x + 0.5 * k1 * kStepSize);
		long double k3 = FirstPositionDerivative(previous_x + 0.5 * k2 * kStepSize, previous_v + 0.5 * m2 * kStepSize);
		long double m3 = FirstVelocityDerivative(previous_v + 0.5 * m2 * kStepSize, previous_x + 0.5 * k2 * kStepSize);
		long double k4 = FirstPositionDerivative(previous_x + k3 * kStepSize, previous_v + m3 * kStepSize);
		long double m4 = FirstVelocityDerivative(previous_v + m3, previous_x + k4);

		long double x = previous_x + kStepSize * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
		long double v = previous_v + kStepSize * (m1 + 2.0 * (m2 + m3) + m4) / 6.0;

		datafile << t << '\t' << x << '\t' << v << '\n';
		previous_x = x;
		previous_v = v;
	}

	datafile.close();
}

double Error(double real_value, double value) {
	return std::abs((real_value - value) / real_value) * 100;
}

int main() {
	EulerLoop("200g", kInitialTime, kFinalTime, kInitialPosition, kInitialVelocity);
	HeunLoop("200g", kInitialTime, kFinalTime, kInitialPosition, kInitialVelocity);
	RungeKuttaLoop("200g", kInitialTime, kFinalTime, kInitialPosition, kInitialVelocity);
	return 0;
}
