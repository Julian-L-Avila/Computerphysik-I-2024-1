#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>
#include <map>
#include <string>

const double Initial_Time       = 0.0;
const double Initial_Velocity   = 0.0;
const double Initial_Position   = 1e3;
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
	double term1 = (Initial_Velocity - Terminal_Velocity) * exp(Acceleration * (t - Initial_Time) / Terminal_Velocity);
	double term2 = (Initial_Velocity + Terminal_Velocity) * exp(-Acceleration * (t - Initial_Time) / Terminal_Velocity);

	return -Terminal_Velocity * ((term1 + term2) / (term1 - term2));
}

double analytic_solution_position(double t) {
	double term1 = (Terminal_Velocity - Initial_Time) * exp((2 * Acceleration * t) / Terminal_Velocity);
	double term2 = (Terminal_Velocity + Initial_Time) * exp((2 * Acceleration * Initial_Time) / Terminal_Velocity);
	double term3 = (Terminal_Velocity + Initial_Velocity) * exp((2 * Acceleration * Initial_Time) / Terminal_Velocity);
	double term4 = (Terminal_Velocity - Initial_Velocity);

	return -(Terminal_Velocity * (Terminal_Velocity * log(std::abs(term1 + term2)) - Acceleration * t - Terminal_Velocity * log(std::abs(term3 + term4)))) / Acceleration;
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

// TODO Implement Runge Kutta Method
void runge_kutta_method() {
}

void iteration_loop(const std::string& method_name, double initial_time, double initial_velocity, double initial_position, double step_size, double (*method)(double, double, double, double (*)(double, double))) {
	std::cout << '\n' << "Chosen method has started ..." << std::endl;
	double previous_velocity = initial_velocity;
	double analytic_velocity = initial_velocity;
	double previous_position = initial_position;

	const std::string velocity_datafile_name = "dat-velocity-" + method_name + ".dat";

	std::ofstream velocity_datafile(velocity_datafile_name);

	velocity_datafile << "# Data for " << method_name << '\n'
		<< "# Time (s)" << '\t' << "Analytic Velocity (ms^{-1})" << "Approx. Velocity (ms^{-1})" << '\n'
		<< std::setprecision(desiredPrecision)
		<< initial_time << '\t' << analytic_velocity << '\t' << previous_velocity << '\n';

	for (double time = initial_time + step_size; previous_position >= 0; time += step_size) {
		previous_position = Initial_Position + analytic_solution_position(time);


		if (previous_position <= 0) {
			std::cout << "Trooper has landed at a time = " << time << '\n'
				<< '\n' << "Data has been saved." << std::endl;
			break;
		}

		analytic_velocity = analytic_solution_velocity(time);
		previous_velocity = method(time, previous_velocity, step_size, first_derivative);

		velocity_datafile << std::setprecision(desiredPrecision)
			<< time << '\t' << analytic_velocity << '\t' << previous_velocity << '\n';
	}

	velocity_datafile.close();
}

int menu_option () {
	std::cout << "This script solves the parachute problem's differential equation with different methods, please select one: " << '\n'
		<< "1. Taylor's Method" << '\n'
		<< "2. Euler's Method" << '\n'
		<< "3. Heun's Method" << '\n'
		<< "4. Runge-Kutta's Method" << '\n';

	int option;
	std::cout << "Press key: ";
	std::cin >> option;

	return option;
}

int main() {
	std::map<int, std::string> name_option = {
		{1, "taylor"},
		{2, "euler"},
		{3, "heun"},
	};

	std::map<std::string, double (*)(double, double, double, double (*)(double, double))> method_lookup = {
			{"taylor", taylor_method},
			{"euler", euler_method},
			{"heun", heun_method},
	};

	system("clear");

	int option = menu_option();
	if (option < 1 || option > name_option.size()) {
		std::cout << "No method was chosen." << '\n' << "Exit." << std::endl;
		return 1;
	}

	const std::string& method_name = name_option[option];

	iteration_loop(method_name, Initial_Time, Initial_Velocity, Initial_Position, Step_Size, method_lookup[method_name]);

	std::string plotname = "plot-" + method_name + ".gnu";
	std::string outputname = "velocity-" + method_name + ".pdf";

	std::ofstream plotfile(plotname);

	plotfile << "set term pdfcairo" << '\n'
		<< "set output '" << outputname << "'" << '\n'
		<< "set grid" << '\n'
		<< "set ylabel 'v [ms^{-1}]'" << '\n'
		<< "set xlabel 't [s]'" << '\n'
		<< "set auto xy" << '\n'
		<< "p 'dat-velocity-" << method_name << ".dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Approx.'" << std::endl;

	plotname = "gnuplot -p " + plotname;
	const char* plotname_c = plotname.c_str();
	system(plotname_c);

	return 0;
}
