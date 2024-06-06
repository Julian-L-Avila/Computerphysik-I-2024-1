/* Pendulum and Double Pendulum Simulation for Computational Physics 2nd Exam
 * UD FCJ - FCMN - Physics
 * Made by:
 * - Julian Avila - 20212107030
 */

#include <iostream>
#include <cmath>
#include <vector>

long double gravity_acceleration, natural_frequency, natural_frequency_square,
		Amplitude, Shift;

std::vector<long double> StateVariables;

struct VariablesAt0 {
	double angle_1;
	double angle_velocity_1;
	double length_1;
	double mass_1;
	double angle_2;
	double angle_velocity_2;
	double length_2;
	double mass_2;
} InitialConditions;

void GetConstantsLinear(VariablesAt0& InitialConditions);
long double AngleDerivative(long double& angle_velocity);
long double VelocityDerivativeLinear(long double& angle,
		double& natural_frequency_square);
long double AnalyticAngleLinear(double t);
long double AnalyticVelocityLinear(double t);
long double EnergyLinear(long double& angle, long double& angle_velocity,
		double& mass);
long double VelocityDerivative(long double& angle,
		double& natural_frequency_square);
long double EllipticIntegralFirstKind(int& degree, double& angle);

int main() {
	natural_frequency = std::sqrt(9.81 / 1.0);
	double angle_1 = 20.0 * M_PI / 180.0;
	int b = 10000;
	long double a = EllipticIntegralFirstKind(b, angle_1);
	std::cout << a;
}

// Pendulum (Linear)

void GetConstantsLinear(VariablesAt0& InitialConditions) {
	natural_frequency_square = gravity_acceleration / InitialConditions.length_1;
	natural_frequency = std::sqrt(natural_frequency_square);
	Amplitude = InitialConditions.angle_1;
	Shift = InitialConditions.angle_velocity_1 / natural_frequency;
}

long double AngleDerivative(long double& angle_velocity) {
	return angle_velocity;
}

long double VelocityDerivativeLinear(long double& angle,
		double& natural_frequency_square) {
	return - natural_frequency_square * angle;
}

long double AnalyticAngleLinear(double t) {
	return Amplitude * std::cos(natural_frequency * t) +
		Shift * std::sin(natural_frequency * t);
}

long double AnalyticVelocityLinear(double t) {
	return Shift * natural_frequency * std::cos(natural_frequency * t) -
		Amplitude * natural_frequency * std::sin(natural_frequency * t);
}

long double EnergyLinear(long double& angle, long double& angle_velocity,
		double& mass) {
	long double y_position = - InitialConditions.length_1 * std::cos(angle);

	long double K = InitialConditions.length_1 * InitialConditions.length_1 *
		angle_velocity * angle_velocity / 4.0;
	long double V = mass * gravity_acceleration * y_position;

	return K + V;
}

// Pendulum (Non-Linear)

long double VelocityDerivative(long double& angle,
		double& natural_frequency_square) {
	return - natural_frequency_square * std::sin(angle);
}

long double EllipticIntegralFirstKind(int& degree, double& angle) {
	long double k = std::sin(angle / 2.0);
	long double sum = 1.0;
	long double product = 1.0;
	long double product_k = 1.0;

	for (int i = 0; i <= degree + 1; i+=1) {
		product *= (2.0 * i + 1.0) / (2.0 * i + 2.0);
		product_k *= k * k;
		sum += product * product * product_k;
	}

	return 4.0 * sum * M_PI / (2.0 * natural_frequency);
}
