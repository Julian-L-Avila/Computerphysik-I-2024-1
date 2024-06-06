/* Pendulum and Double Pendulum Simulation for Computational Physics 2nd Exam
 * UD FCJ - FCMN - Physics
 * Made by:
 * - Julian Avila - 20212107030
 */

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

long double gravity_acceleration, natural_frequency, natural_frequency_square,
		Amplitude, Shift;
int DesiredPrecision;

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
long double EllipticIntegralFirstKind(int& degree, long double& angle,
		long double& modulus);
long double PeriodNotLinear(long double& angle, int& degree);
long double SinusAmplitudis(int& degree, long double& angle,
		long double& modulus, double& t);
long double AnalyticAngleNonLinear(int& degree, double& t,
		long double& initial_angle);

int main() {
}

// Pendulum (Linear)

void GetConstants(VariablesAt0& InitialConditions) {
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

long double EllipticIntegralFirstKind(int& degree, long double& angle,
		long double& modulus) {
	long double sum = 1.0;
	long double product = 1.0;
	long double product_k = 1.0;

	for (int i = 0; i <= degree + 1; i+=1) {
		product *= (2.0 * i + 1.0) / (2.0 * i + 2.0);
		product_k *= modulus * modulus;
		sum += product * product * product_k;
	}

	return sum * M_PI / 2.0;
}

long double PeriodNotLinear(int& degree, long double& angle,
		long double& modulus) {
	return 4 * EllipticIntegralFirstKind(degree, angle, modulus) / natural_frequency;
}

long double SinusAmplitudis(int& degree, long double& angle,
		long double& modulus, double& t) {
	long double K = EllipticIntegralFirstKind(degree, angle, modulus);
	K -= natural_frequency * t;
	return boost::math::jacobi_sn(modulus, K);
}

long double AnalyticAngleNonLinear(int& degree, double& t,
		long double& initial_angle) {
	long double modulus = std::sin(initial_angle / 2.0);
	long double theta = modulus * SinusAmplitudis(degree, initial_angle,
			modulus, t);
	return 2.0 * std::asin(theta);
}
