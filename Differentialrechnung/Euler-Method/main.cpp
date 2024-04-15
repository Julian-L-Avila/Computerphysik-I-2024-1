#include <iomanip>
#include <iostream>
#include <cmath>
#include <ostream>
#include <fstream>

const double Initial_Time       = 0.0;
const double Initial_Velocity   = 0.0;
const double Final_Time         = 10.0;
const double Mass               = 70.0;
const double Delta_Coefficient  = 0.8;
const double Area_Cross_Section = 0.6;

const double Acceleration = 9.81;
const double Air_Density  = 1.293;

const double Drag_Coefficient  = Air_Density * Area_Cross_Section * Delta_Coefficient / 2;
const double Terminal_Velocity = std::sqrt(Mass * Acceleration / Drag_Coefficient);

const double Step_Size = 0.5;
const int desiredPrecision = 10;

double analytic_solution_velocity(double t) {
	return Terminal_Velocity * (((Initial_Velocity - Terminal_Velocity) * exp(Acceleration * (t - Initial_Time) / Terminal_Velocity)) + ((Initial_Velocity + Terminal_Velocity) * exp(-Acceleration * (t - Initial_Time) / Terminal_Velocity))) / (((Initial_Velocity - Terminal_Velocity) * exp(Acceleration * (t - Initial_Time) / Terminal_Velocity)) - ((Initial_Velocity + Terminal_Velocity) * exp(-Acceleration * (t - Initial_Time) / Terminal_Velocity)));
}

double first_derivative_velocity(double t, double v) {
	return -(Drag_Coefficient * v * v / Mass) + Acceleration;
}

double euler_method(double t, double previous_y, double stepsize, double (*derivative)(double t, double previous_y)) {
	return previous_y + stepsize * derivative(t, previous_y);
}

int main() {
	double analytic_velocity = Initial_Velocity;
	double approx_velocity = Initial_Velocity;

	system("clear");

	std::cout << "Starting script" << '\n';

	std::ofstream datafile("dat-velocity-euler.dat");
	datafile << "Time (s)" << '\t' << "Velocity Analytic (m/s)" << '\t' << "Velocity Euler (m/s)" << '\n'
		<< std::setprecision(desiredPrecision) << Initial_Time << '\t'
		<< analytic_velocity << '\t'
		<< approx_velocity
		<< std::endl;

	for (double t = Initial_Time + Step_Size; t <= Final_Time; t += Step_Size) {
		analytic_velocity = analytic_solution_velocity(t);
		approx_velocity = euler_method(t, approx_velocity, Step_Size, first_derivative_velocity);

		datafile << std::setprecision(desiredPrecision) << t << '\t'
			<< analytic_velocity << '\t'
			<< approx_velocity << '\t'
			<< '\n';
	}

	datafile.close();

	std::cout << "Data saved" << '\n'
		<< "Plotting velocity comparison" << '\n';

	std::ofstream velocity_plot("plot-velocity.gnu");

	velocity_plot << "set term pdfcairo" << '\n'
		<< "set output 'velocity.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set auto xy" << '\n'
		<< "set tit 'Velocity vs Time'" << '\n'
		<< "set ylabel 'v(t) [ms^(-1)]'" << '\n'
		<< "set xlabel 't [s]' " << '\n'
		<< "p 'dat-velocity-euler.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Euler Method'" << '\n';

	velocity_plot.close();

	std::cout << "Done" << std::endl;

	return 0;
}
