#include <iomanip>
#include <iostream>
#include <cmath>
#include <ostream>
#include <fstream>

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

double analytic_solution_velocity(double t) {
	return -Terminal_Velocity * (((Initial_Velocity - Terminal_Velocity) * exp(Acceleration * (t - Initial_Time) / Terminal_Velocity)) + ((Initial_Velocity + Terminal_Velocity) * exp(-Acceleration * (t - Initial_Time) / Terminal_Velocity))) / (((Initial_Velocity - Terminal_Velocity) * exp(Acceleration * (t - Initial_Time) / Terminal_Velocity)) - ((Initial_Velocity + Terminal_Velocity) * exp(-Acceleration * (t - Initial_Time) / Terminal_Velocity)));
}

double analytic_solution_possition(double t) {
	return -(Terminal_Velocity * (Terminal_Velocity * log(std::abs((Terminal_Velocity - Initial_Time) * exp((2 * Acceleration * t)/ Terminal_Velocity) + (Terminal_Velocity + Initial_Time) * exp((2* Acceleration * Initial_Time) / Terminal_Velocity))) - Acceleration * t - Terminal_Velocity * log(std::abs((Terminal_Velocity + Initial_Velocity) * exp((2 * Acceleration * Initial_Time) / Terminal_Velocity ) + Terminal_Velocity - Initial_Velocity))))/ Acceleration;
}

double first_derivative_velocity(double t, double v) {
	return (Drag_Coefficient * v * v / Mass) - Acceleration;
}

double euler_method(double t, double previous_y, double stepsize, double (*derivative)(double t, double previous_y)) {
	return previous_y + stepsize * derivative(t, previous_y);
}

double simpson_integral(double lower_limit, double upper_limit, double step_size, double rectangle_number, double (*function)(double t)) {
	step_size = (upper_limit - lower_limit) / rectangle_number;
	double halft_step_size = 0.5 * step_size;
	double sum = function(lower_limit) + function(upper_limit);

	for (double t = lower_limit; t <= upper_limit - step_size; t += step_size) {
		sum += 4.0 * function(t + halft_step_size) + 2.0 * function(t);
	}
	return halft_step_size * sum / 3.0;
}

double error(double real_value, double approx_value) {
	return std::abs((real_value - approx_value) / real_value) * 100;
}

int main() {
	double analytic_velocity = Initial_Velocity;
	double approx_velocity   = Initial_Velocity;
	double analytic_position = Initial_Position;
	double approx_position   = Initial_Position;

	system("clear");

	std::cout << "Starting script" << '\n';

	std::ofstream datafile("dat-velocity.dat");
	datafile << "Time (s)" << '\t' << "Velocity Analytic (m/s)" << '\t' << "Velocity Euler (m/s)" << '\n'
		<< std::setprecision(desiredPrecision) << Initial_Time << '\t'
		<< analytic_velocity << '\t'
		<< approx_velocity
		<< std::endl;

	std::ofstream position_datafile("dat-position.dat");
	position_datafile << "Time (s)" << '\t' << "Position Analytic (m)" << '\t' << "Position Euler (m)" << '\n'
		<< std::setprecision(desiredPrecision) << Initial_Time << '\t'
		<< analytic_position << '\t'
		<< approx_position
		<< std::endl;

	std::ofstream errorfile("dat-error.dat");
	errorfile << "Time (s)" << '\t' << "Error Velocity (ms^(-1))" << '\t' << "Error Position (m)" << '\n';

	for (double t = Initial_Time + Step_Size; t <= Final_Time; t += Step_Size) {
		analytic_velocity = analytic_solution_velocity(t);
		approx_velocity = euler_method(t, approx_velocity, Step_Size, first_derivative_velocity);
		analytic_position = Initial_Position + analytic_solution_possition(t);
		approx_position += simpson_integral(t, t + Step_Size, Step_Size, 100, analytic_solution_velocity);

		if (analytic_position <= 0.0) {
			std::cout << "Landed" << std::endl;
			break;
		}

		datafile << std::setprecision(desiredPrecision) << t << '\t'
			<< analytic_velocity << '\t'
			<< approx_velocity << '\t'
			<< '\n';

		position_datafile << std::setprecision(desiredPrecision) << t << '\t'
			<< analytic_position << '\t'
			<< approx_position << '\t'
			<< '\n';

		errorfile << std::setprecision(desiredPrecision) << t << '\t'
			<< error(analytic_velocity, approx_velocity) << '\t'
			<< error(analytic_position, approx_position) << '\t'
			<< '\n';
	}

	datafile.close();
	position_datafile.close();
	errorfile.close();

	std::cout << "Velocity Data saved" << '\n'
		<< "Plotting velocity comparison" << '\n';

	std::ofstream velocity_plot("plot-velocity.gnu");

	velocity_plot << "set term pdfcairo" << '\n'
		<< "set output 'velocity.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set auto xy" << '\n'
		<< "set tit 'Velocity vs Time'" << '\n'
		<< "set ylabel 'v(t) [ms^{-1}]'" << '\n'
		<< "set xlabel 't [s]' " << '\n'
		<< "p 'dat-velocity.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Euler Method'" << '\n';

	velocity_plot.close();

	system("gnuplot -p 'plot-velocity.gnu'");

	std::cout << "Velocity plotted" << '\n';

	std::ofstream position_plot("plot-position.gnu");

	position_plot << "set term pdfcairo" << '\n'
		<< "set output 'position.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set auto xy" << '\n'
		<< "set tit 'Position vs Time'" << '\n'
		<< "set xlabel 'x(t) [m]'" << '\n'
		<< "set ylabel 't [s]'" << '\n'
		<< "p 'dat-position.dat' u 1:2 w l tit 'Analytic', '' u 1:3 w l tit 'Simpson Method'" << '\n';

	position_plot.close();

	system("gnuplot -p 'plot-position.gnu'");

	std::cout << "Position plotted" << '\n';

	std::ofstream error_plot("plot-error.gnu");

	error_plot << "set term pdfcairo" << '\n'
		<< "set output 'error.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set logscale y" << '\n'
		<< "set tit 'Velocity and Position Error'" << '\n'
		<< "set xlabel 't [s]'" << '\n'
		<< "set ylabel 'Error %'" << '\n'
		<< "p 'dat-error.dat' u 1:2 w l tit 'Velocity Error', '' u 1:3 w l tit 'Position Error'" << '\n';

	error_plot.close();

	system("gnuplot -p 'plot-error.gnu'");

	std::cout << "Done" << std::endl;

	return 0;
}
