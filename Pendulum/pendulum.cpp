/* Pendulum and Double Pendulum Simulation for Computational Physics 2nd Exam
 * UD FCJ - FCMN - Physics
 * Made by:
 * - Julian Avila - 20212107030
 */

#include <boost/math/special_functions/math_fwd.hpp>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

long double gravity_acceleration, natural_frequency, natural_frequency_square,
		Amplitude, Shift, modulus, modulus_square, elliptic_integral;

int approx_degree;

const int kDesiredPrecision = 20;

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

void LinearConstans(VariablesAt0& InitialConditions);
void NonLinearConstants(VariablesAt0 InitialConditions);
void GetConstants(VariablesAt0& InitialConditions, std::vector<int> systems);
long double AngleDerivative(long double& angle_velocity);
long double VelocityDerivativeLinear(long double& angle,
		double& natural_frequency_square);
long double AnalyticAngleLinear(double t);
long double AnalyticVelocityLinear(double t);
long double EnergySimple(long double& angle, long double& angle_velocity,
		double& mass);
long double VelocityDerivative(long double& angle,
		double& natural_frequency_square);
long double EllipticIntegralFirstKind(long double& angle);
long double PeriodNotLinear(long double& angle);
long double SinusAmplitudis(long double& angle, double& t);
long double CosinusAmplitudis(long double& angle, double& t);
long double DeltaAmplitudis(long double& angle, double& t);
long double AnalyticAngleNonLinear(int& degree, double& t,
		long double& initial_angle);
long double AnalyticVelocityNonLinear(double& t,
		long double& initial_angle);

int main() {
}

void LinearConstans(VariablesAt0& InitialConditions) {
	natural_frequency_square = gravity_acceleration / InitialConditions.length_1;
	natural_frequency = std::sqrt(natural_frequency_square);
	Amplitude = InitialConditions.angle_1;
	Shift = InitialConditions.angle_velocity_1 / natural_frequency;
}

void NonLinearConstants(VariablesAt0 InitialConditions) {
	modulus = std::sin(InitialConditions.angle_1 / 2.0);
	modulus_square = modulus * modulus;
	elliptic_integral = EllipticIntegralFirstKind(Amplitude);
}

void GetConstants(VariablesAt0& InitialConditions, std::vector<int> systems) {
	if (systems[0] == 1) {
		LinearConstans(InitialConditions);
	};

	if (systems[1] == 1) {
		NonLinearConstants(InitialConditions);
	}

	if (systems[2] == 1) {
	}
}

// Pendulum (Linear)

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

long double EnergySimple(long double& angle, long double& angle_velocity,
		double& mass, double& length) {
	long double y_position = - length * std::cos(angle);

	long double K = length * length *
		angle_velocity * angle_velocity / 4.0;
	long double V = mass * gravity_acceleration * y_position;

	return K + V;
}

// Pendulum (Non-Linear)

long double VelocityDerivative(long double& angle,
		double& natural_frequency_square) {
	return - natural_frequency_square * std::sin(angle);
}

long double EllipticIntegralFirstKind(long double& angle) {
	long double sum = 1.0;
	long double product = 1.0;
	long double product_k = 1.0;

	for (int i = 0; i <= approx_degree + 1; i+=1) {
		product *= (2.0 * i + 1.0) / (2.0 * i + 2.0);
		product_k *= modulus_square;
		sum += product * product * product_k;
	}

	return sum * M_PI / 2.0;
}

long double PeriodNotLinear(long double& angle) {
	return 4 * EllipticIntegralFirstKind(angle) / natural_frequency;
}

long double SinusAmplitudis(long double& angle, double& t) {
	elliptic_integral -= natural_frequency * t;
	return boost::math::jacobi_sn(modulus_square, elliptic_integral);
}

long double CosinusAmplitudis(long double& angle, double& t) {
	elliptic_integral -= natural_frequency * t;
	return boost::math::jacobi_cn(modulus_square, elliptic_integral);
}

long double DeltaAmplitudis(long double& angle, double& t) {
	elliptic_integral -= natural_frequency * t;
	return boost::math::jacobi_dn(modulus_square, elliptic_integral);
}

long double AnalyticAngleNonLinear(double& t,
		long double& initial_angle) {
	long double theta = modulus * SinusAmplitudis(initial_angle, t);
	return 2.0 * std::asin(theta);
}

long double AnalyticVelocityNonLinear(double& t,
		long double& initial_angle) {
	long double cn = CosinusAmplitudis(initial_angle, t);
	long double dn = DeltaAmplitudis(initial_angle, t);
	return - 2.0 * natural_frequency * modulus * cn * dn;
}

// Double Pendulum

long double AngleDerivativeDouble(double& t, long double& angle,
		long double& angle_velocity) {
	return angle_velocity;
}

long double VelocityDerivativeDouble1(double& t, long double& angle_1,
		long double& angle_2, long double& angle_velocity_1,
		long double& angle_velocity_2, double& mass_1, double& mass_2,
		double& length_1, double& length_2) {
	long double angle_diff = angle_1 - angle_2;
	long double sin_of_diff = std::sin(angle_diff);
	long double angle_velocity_1_square = angle_velocity_1 * angle_velocity_1;
	long double angle_velocity_2_square = angle_velocity_2 * angle_velocity_2;
	long double numerator = 18.0 * mass_2 * std::cos(angle_diff) *
		(gravity_acceleration * std::sin(angle_2) - length_1 *
		angle_velocity_1_square * sin_of_diff) - 12.0 *
		(mass_2 * length_2 * angle_velocity_2_square * sin_of_diff +
		(mass_1 + 2.0 * mass_2) * gravity_acceleration * std::sin(angle_1));
	long double denominator = length_1 * (15.0 * mass_2 + 8.0 * mass_1 -
		9.0 * mass_2 * std::cos(2.0 * angle_diff));
	return numerator / denominator;
}

long double VelocityDerivativeDouble2(double& t, long double& angle_1,
		long double& angle_2, long double& angle_velocity_1,
		long double& angle_velocity_2, double& mass_1, double& mass_2,
		double& length_1, double& length_2) {
	long double angle_diff = angle_1 - angle_2;
	long double sin_of_diff = std::sin(angle_diff);
	long double angle_velocity_1_square = angle_velocity_1 * angle_velocity_1;
	long double angle_velocity_2_square = angle_velocity_2 * angle_velocity_2;
	long double numerator = 36.0 * std::cos(angle_diff) *
		(mass_2 * length_2 * angle_velocity_2_square * sin_of_diff +
		(mass_1 + 2.0 * mass_2) * gravity_acceleration * std::sin(angle_1)) -
		4.0 * (mass_1 + 3.0 * mass_2) * (gravity_acceleration * std::sin(angle_2) -
		length_1 * angle_velocity_1_square + sin_of_diff);
	long double denominator = length_2 * (15.0 * mass_2 + 8.0 * mass_1 - 9.0 *
			mass_2 * std::cos(2.0 * angle_diff));
	return numerator / denominator;
}

long double EnergyDouble(long double& angle_1, long double& angle_2,
		long double& angle_velocity_1, long double& angle_velocity_2,
		double& mass_1, double& mass_2, double& length_1, double& length_2) {
	long double k_1 = mass_1 * length_1 *
		(length_1 * angle_velocity_1 * angle_velocity_1 / 3.0 -
			gravity_acceleration * std::cos(angle_1));
	long double k_2 = mass_2 * length_2 * (length_2 * (angle_velocity_1 *
			angle_velocity_1 + angle_velocity_2 * angle_velocity_2 / 3.0 +
			angle_velocity_1 * angle_velocity_2 * std::cos(angle_1 - angle_2)) -
			gravity_acceleration * (2.0 * std::cos(angle_1) - std::cos(angle_2)));
	return k_1 + k_2;
}
