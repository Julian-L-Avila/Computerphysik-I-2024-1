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
const double kFinalTime       = 20.0;
const double kInitialPosition = 5.0;
const double kInitialVelocity = 0.0;
const double kSpringConstant  = 0.1;

// Must go after user choice
double mass                     = 0.1;
double natural_frequency_square = kSpringConstant / mass;
double natural_frequency        = sqrt(natural_frequency_square);
//

const double kStepSize      = 1e-2;
const int kDesiredPrecision = 20;

long double FirstVelocityDerivative(long double velocity, long double position) {
	return -natural_frequency_square * position;
}

long double FirstPositionDerivative(long double position, long double velocity) {
	return velocity;
}

long double AnalyticSolutionVelocity(double t) {
	return -natural_frequency * kInitialPosition * sin(natural_frequency * t);
}

long double AnalyticSolutionPosition(double t) {
	return kInitialPosition * cos(natural_frequency * t);
}

long double Energy(double mass, long double x, long double v) {
	return 0.5 * mass * v * v + 0.5 * kSpringConstant * x * x;
}

long double EulerMethod(long double previous_x, long double previous_y, long double (*Derivative)(long double, long double)) {
	return previous_x + kStepSize * Derivative(previous_x,  previous_y);
}

void InitialLoop(const std::string& mass_string, const std::string& method_name, double initial_time, double final_time, double initial_position, double initial_velocity, void(*MethodLoop) (double&, long double&, long double&, long double&, long double&, long double&)) {
	long double previous_x, previous_v, energy, x, v;
	previous_x = initial_position;
	previous_v = initial_velocity;
	energy = Energy(mass, previous_x, previous_v);

	std::string path_file_name = "./Approx-Data/" + method_name + "-" + mass_string + ".dat";

	std::ofstream datafile(path_file_name);

	datafile << "# " + method_name + "Data" << '\n'
		<< "#Time (s) \t Position (m) \t Velocity (ms^-1) \t Energy (J)" << '\n'
		<< std::setprecision(kDesiredPrecision) << std::fixed
		<< initial_time << '\t' << initial_position << '\t' << initial_velocity << '\t' << energy << '\n';

	for (double t = initial_time + kStepSize; t <= final_time; t += kStepSize) {
		MethodLoop(t, previous_x, previous_v, energy, x, v);

		datafile << t << '\t' << x << '\t' << v << '\t' << energy << '\n';
		previous_x = x;
		previous_v = v;
	}

	datafile.close();
}

void EulerLoop(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	x = EulerMethod(previous_x, previous_v, FirstPositionDerivative);
	v = EulerMethod(previous_v, previous_x, FirstVelocityDerivative);
	energy = Energy(mass, previous_x, previous_v);
}

void HeunLoop(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	long double euler_x = EulerMethod(previous_x, previous_v, FirstPositionDerivative);
	long double euler_v = EulerMethod(previous_v, previous_x, FirstVelocityDerivative);

	x = previous_x + 0.5 * kStepSize * (FirstPositionDerivative(previous_x, previous_v) + FirstPositionDerivative(euler_x, euler_v));
	v = previous_v + 0.5 * kStepSize * (FirstVelocityDerivative(previous_v, previous_x) + FirstVelocityDerivative(euler_v, euler_x));
	energy = Energy(mass, previous_x, previous_v);
}

void RungeKuttaLoop(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	long double k1 = FirstPositionDerivative(previous_x, previous_v);
	long double m1 = FirstVelocityDerivative(previous_v, previous_x);
	long double k2 = FirstPositionDerivative(previous_x + 0.5 * k1 * kStepSize, previous_v + 0.5 * m1 * kStepSize);
	long double m2 = FirstVelocityDerivative(previous_v + 0.5 * m1 * kStepSize, previous_x + 0.5 * k1 * kStepSize);
	long double k3 = FirstPositionDerivative(previous_x + 0.5 * k2 * kStepSize, previous_v + 0.5 * m2 * kStepSize);
	long double m3 = FirstVelocityDerivative(previous_v + 0.5 * m2 * kStepSize, previous_x + 0.5 * k2 * kStepSize);
	long double k4 = FirstPositionDerivative(previous_x + k3 * kStepSize, previous_v + m3 * kStepSize);
	long double m4 = FirstVelocityDerivative(previous_v + kStepSize * m3, previous_x + kStepSize * k4);

	x = previous_x + kStepSize * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
	v = previous_v + kStepSize * (m1 + 2.0 * (m2 + m3) + m4) / 6.0;
	energy = Energy(mass, x, v);
}

void AnalyticLoop(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	x = AnalyticSolutionPosition(t);
	v = AnalyticSolutionVelocity(t);
}

double Error(long double real_value, long double value) {
	return std::abs((real_value - value) / real_value) * 100;
}

int main() {
	return 0;
}
