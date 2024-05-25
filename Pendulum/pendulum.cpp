/* Pendulum and Double Pendulum Simulation for Computational Physics 2nd Exam
 * UD FCJ - FCMN - Physics
 * Made by:
 * - Julian Avila - 20212107030
 */

#include <cmath>

long double gravity_acceleration, natural_frequency, natural_frequency_square,
		Amplitude, Shift;

struct StateVariables {
	double angle_1;
	double angle_velocity_1;
	double length_1;
	double mass_1;
	double angle_2;
	double angle_velocity_2;
	double length_2;
	double mass_2;
} InitialConditions;

int main() {
}

// Pendulum (Linear)

void GetConstantsLinear(StateVariables& InitialConditions) {
	natural_frequency_square = gravity_acceleration / InitialConditions.length_1;
	natural_frequency = std::sqrt(natural_frequency_square);
	Amplitude = InitialConditions.angle_1;
	Shift = InitialConditions.angle_velocity_1 / natural_frequency;
}

long double AngleDerivativeLinear(long double& angle_velocity) {
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
