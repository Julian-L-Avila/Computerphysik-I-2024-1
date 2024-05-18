/* Pendulum and Double Pendulum Simulation for Computational Physics 2nd Exam
 * UD FCJ - FCMN - Physics
 * Made by:
 * - Julian Avila - 20212107030
 */

#include <iostream>
#include <cmath>

double initial_angle, initial_angle_velocity, natural_frequency,
			natural_frequency_square, length1;

int main() {
}

// Pendulum (Linear)
// TODO Add A and B constants in code

long double AngleDerivativeLinear(long double& angle_velocity) {
	return angle_velocity;
}

long double VelocityDerivativeLinear(long double& angle,
		double& natural_frequency_square) {

	return - natural_frequency_square * angle;
}

long double AnalyticAngleLinear(double t) {
	double A = initial_angle;
	long double B = initial_angle_velocity / natural_frequency;

	return A * std::cos(natural_frequency * t) +
		B * std::sin(natural_frequency * t);
}

long double AnalyticVelocityLinear(double t) {
	double A = initial_angle;
	long double B = initial_angle_velocity / natural_frequency;

	return B * natural_frequency * std::cos(natural_frequency * t) -
		A * natural_frequency * std::sin(natural_frequency * t);
}

long double EnergyLinear(long double& angle, long double& angle_velocity) {
	long double x_position = length1 * std::sin(angle);
	long double y_position = - length1 * std::cos(angle);
	long double x_velocity = - angle_velocity * y_position;
	long double y_velocity = angle_velocity * x_position;
}
