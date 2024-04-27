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

const double kInitialTime      = 0.0;

const double kStepSize      = 1e-2;
const int kDesiredPrecision = 20;

double FirstDerivative(double t, double v) {
	return 0.0;
}

double AnalyticSolutionVelocity(double t) {
	return 0.0;
}

double AnalyticSolutionPosition(double t) {
	return 0.0;
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

void IterationLoop(const std::string& method_name, double initial_time, double initial_velocity, double initial_position, double step_size, double (*method)(double, double, double, double (*)(double, double))) {
}

int main() {
	return 0;
}
