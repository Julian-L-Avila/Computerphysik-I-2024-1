/* Pendulum and Double Pendulum Simulation for Computational Physics 2nd Exam
 * UD FCJ - FCMN - Physics
 * Made by:
 * - Julian Avila - 20212107030
 */

#include <iostream>
#include <cmath>

double initial_angle, initial_velocity, natural_frequency,
			natural_frequency_square;

int main() {
}

// Pendulum (Linear)
// TODO

long double AngleDerivativeLinear(long double& velocity) {
	return velocity;
}

long double VelocityDerivativeLinear(long double& angle,
		double& natural_frequency_square) {

	return - natural_frequency_square * angle;
}

long double AnalyticAngleLinear(double t) {
	double A = initial_angle;
	long double B = initial_velocity / natural_frequency;

	return A * std::cos(natural_frequency * t) + B * std::sin(natural_frequency * t);
}

long double AnalyticVelocityLinear(double t) {
	double A = initial_angle;
	long double B = initial_velocity / natural_frequency;

	return B * natural_frequency * std::cos(natural_frequency * t) -
		A * natural_frequency * std::sin(natural_frequency * t);
}
