#include <iostream>
#include <cmath>

const double Initial_Time       = 0.0;
const double Initial_Velocity   = 0.0;
const double Initial_Position   = 1e3;
const double Final_Time         = 1e10;
const double Mass               = 70.0;
const double Delta_Coefficient  = 0.8;
const double Area_Cross_Section = 0.6;

const double Acceleration = 9.81;
const double Air_Density  = 1.293;

const double Drag_Coefficient  = Air_Density * Area_Cross_Section * Delta_Coefficient / 2;
const double Terminal_Velocity = std::sqrt(Mass * Acceleration / Drag_Coefficient);

const double Step_Size     = 1e-2;
const int desiredPrecision = 10;

double first_derivative(double t, double v) {
	return (Drag_Coefficient * v * v / Mass) - Acceleration;
}

double second_derivative(double t, double v) {
	return (2 * Drag_Coefficient * v * first_derivative(t, v));
}

double third_derivative(double t, double v) {
	return 0;
}

double analytic_solution_velocity(double t) {
	return -Terminal_Velocity * (((Initial_Velocity - Terminal_Velocity) * exp(Acceleration * (t - Initial_Time) / Terminal_Velocity)) + ((Initial_Velocity + Terminal_Velocity) * exp(-Acceleration * (t - Initial_Time) / Terminal_Velocity))) / (((Initial_Velocity - Terminal_Velocity) * exp(Acceleration * (t - Initial_Time) / Terminal_Velocity)) - ((Initial_Velocity + Terminal_Velocity) * exp(-Acceleration * (t - Initial_Time) / Terminal_Velocity)));
}

double taylor_method(double previous_t, double previous_y, double step_size, double (*derivative)(double t, double y)) {
	double step_size_square = step_size * step_size;
	double step_size_cube = step_size_square * step_size;

	return previous_y + step_size * first_derivative(previous_t, previous_y) + 0.5 * step_size_square * second_derivative(previous_t, previous_y) + (1.0 / 6) * step_size_cube * third_derivative(previous_t, previous_y);
}

double euler_method(double previous_t, double previous_y, double step_size, double (*derivative)(double t, double y)) {
	return previous_y + step_size * derivative(previous_t, previous_y);
}

double heun_method(double previous_t, double previous_y, double step_size, double (*derivative)(double t, double y)) {
	double euler_y = euler_method(previous_t, previous_y, step_size, derivative);
	return previous_y + 0.5 * step_size * (derivative(previous_t, previous_y) + derivative(previous_t + step_size, euler_y));
}

void iteration_loop(double initial_t, double initial_y, double step_size, double (*method)(double previous_t, double previous_y, double step_size, double (*derivative)(double t, double y))) {
}

int menu_option () {
	std::cout << "This script solves the parachute problem with different methods, please select one: " << '\n'
		<< "1. Taylor's Method" << '\n'
		<< "2. Euler's Method" << '\n'
		<< "3. Heun's Method" << std::endl;
	
	int option;
	std::cin >> option;

	return option;
}

void menu() {
	switch (menu_option()) {
		case 1:
	}
}
