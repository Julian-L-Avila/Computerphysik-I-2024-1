/* Solution to harmonic differential equation for "Fisica Computacional I"
 *
 * Made by:
 * Bryan Martinez   - 20212107008
 * Laura Y. Herrera - 20212107011
 * Julian L. Avila  - 20212107030
 * Juan S. Acuña    - 20212107034
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <iomanip>

const double kInitialTime      = 0.0;
const double kInitialVelocity  = 0.0;
const double kSpringConstant   = 19.49;
const double kStepSize         = 0.01;
const int    kDesiredPrecision = 20;

double mass, natural_frequency, natural_frequency_square, final_time, initial_position;
std::string mass_as_string, method_name;

long double FirstVelocityDerivative(long double velocity, long double position);
long double FirstPositionDerivative(long double position, long double velocity);
long double AnalyticSolutionVelocity(double t);
long double AnalyticSolutionPosition(double t);
long double Energy(long double x, long double v);
void AnalyticImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v);
long double EulerMethod(long double previous_x, long double previous_y,
													long double (*Derivative)(long double, long double));
void EulerImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v);
void HeunImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v);
void RungeKuttaImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v);
long double AbsolutePercentageError(long double real_value, long double value);
void DataLoop(std::string& mass_as_string, std::string& method_name,
									double initial_time, double final_time, double initial_position,
									double initial_velocity, void(*MethodLoop) (double&, long double&, long double&, long double&, long double&, long double&));
std::string MassChoice();
void MethodData(std::string& mass_as_string);
void PlotData(std::string& mass_as_string, std::string& method_name);

int main() {
	int stop = 0;
	system("clear");

	std::cout << ":::::::::::::::::::::::::::::::WELCOME TO A MASS-SPRING SIMULATION:::::::::::::::::::::::::::::::" << '\n';

	while (stop == 0) {
	mass_as_string = MassChoice();

	method_name = "Analytic";
	DataLoop(mass_as_string, method_name, kInitialTime, final_time, initial_position, kInitialVelocity, AnalyticImplementation);
	MethodData(mass_as_string);

	PlotData(mass_as_string, method_name);
	std::string plotfile_name = "./Plot-Data/" + method_name + "-" + mass_as_string + ".gnu";
	system(("gnuplot -p " + plotfile_name).c_str());
	std::string plot_pdf = "./Plot-Data/plot-" + method_name + "-" + mass_as_string + ".pdf";
	system(("atril " + plot_pdf + "&").c_str());

	std::cout << "Press any key to continue or 'q' to quit: ";
	char choice;
	std::cin >> choice;
	if (choice == 'q' || choice == 'Q') {
		std::cout << "End of script" << std::endl;
		system("clear");
		break;
	}
	system("clear");
	}

	return 0;
}

long double FirstVelocityDerivative(long double velocity, long double position) {
	return -natural_frequency_square * position;
}

long double FirstPositionDerivative(long double position, long double velocity) {
	return velocity;
}

long double AnalyticSolutionVelocity(double t) {
	return -natural_frequency * initial_position * sin(natural_frequency * t);
}

long double AnalyticSolutionPosition(double t) {
	return initial_position * cos(natural_frequency * t);
}

long double Energy(long double x, long double v) {
	return 0.5 * mass * v * v + 0.5 * kSpringConstant * x * x;
}

long double EulerMethod(long double previous_x, long double previous_y, long double (*Derivative)(long double, long double)) {
	return previous_x + kStepSize * Derivative(previous_x,  previous_y);
}

void AnalyticImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	x = AnalyticSolutionPosition(t);
	v = AnalyticSolutionVelocity(t);
}

void EulerImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	x = EulerMethod(previous_x, previous_v, FirstPositionDerivative);
	v = EulerMethod(previous_v, previous_x, FirstVelocityDerivative);
}

void HeunImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	long double euler_x = EulerMethod(previous_x, previous_v, FirstPositionDerivative);
	long double euler_v = EulerMethod(previous_v, previous_x, FirstVelocityDerivative);

	x = previous_x + 0.5 * kStepSize * (FirstPositionDerivative(previous_x, previous_v) + FirstPositionDerivative(euler_x, euler_v));
	v = previous_v + 0.5 * kStepSize * (FirstVelocityDerivative(previous_v, previous_x) + FirstVelocityDerivative(euler_v, euler_x));
}

void RungeKuttaImplementation(double& t, long double& previous_x, long double& previous_v, long double& energy, long double& x, long double& v) {
	long double k1 = FirstPositionDerivative(previous_x, previous_v);
	long double m1 = FirstVelocityDerivative(previous_v, previous_x);
	long double k2 = FirstPositionDerivative(previous_x + 0.5 * k1 * kStepSize, previous_v + 0.5 * m1 * kStepSize);
	long double m2 = FirstVelocityDerivative(previous_v + 0.5 * m1 * kStepSize, previous_x + 0.5 * k1 * kStepSize);
	long double k3 = FirstPositionDerivative(previous_x + 0.5 * k2 * kStepSize, previous_v + 0.5 * m2 * kStepSize);
	long double m3 = FirstVelocityDerivative(previous_v + 0.5 * m2 * kStepSize, previous_x + 0.5 * k2 * kStepSize);
	long double k4 = FirstPositionDerivative(previous_x + k3 * kStepSize, previous_v + m3 * kStepSize);
	long double m4 = FirstVelocityDerivative(previous_v + kStepSize * m3, previous_x + kStepSize * k4);

	x = previous_x + kStepSize * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
	v = previous_v + kStepSize * (m1 + 2.0 * (m2 + m3) + m4) / 6.0;
}

long double AbsolutePercentageError(long double actual, long double estimated) {
	return std::abs((actual - estimated) / actual) * 100.0;
}

void DataLoop(std::string& mass_as_string, std::string& method_name, double initial_time, double final_time, double initial_position, double initial_velocity, void(*MethodLoop) (double&, long double&, long double&, long double&, long double&, long double&)) {
	long double previous_x, previous_v, energy, x, v;
	previous_x = initial_position;
	previous_v = initial_velocity;
	energy = Energy(previous_x, previous_v);

	natural_frequency_square   = kSpringConstant / mass;
	natural_frequency          = sqrt(natural_frequency_square);

	std::string path_file_name = "./Approx-Data/" + method_name + "-" + mass_as_string + ".dat";

	std::ofstream datafile(path_file_name);

	datafile << "# " + method_name + "Data" << '\n'
		<< "#Time (s) \t Position (m) \t Velocity (ms^-1) \t Energy (J)" << '\n'
		<< std::setprecision(kDesiredPrecision) << std::fixed
		<< initial_time << '\t' << initial_position << '\t' << initial_velocity << '\t'
		<< energy << '\n';

	for (double t = initial_time + kStepSize; t <= final_time; t += kStepSize) {
		MethodLoop(t, previous_x, previous_v, energy, x, v);
		energy = Energy(previous_x, previous_v);

		datafile << t << '\t' << x << '\t' << v << '\t' << energy << '\n';
		previous_x = x;
		previous_v = v;
	}

	datafile.close();
}

std::string MassChoice() {
	char mass_option;

	std::cout << "Select the mass value of the object to obtain the experimental, analytical and numerical data of:" << '\n'
		<< '\t' << "1. 100 g" << '\n'
		<< '\t' << "2. 200 g" << '\n'
		<< '\t' << "3. 250 g" << '\n'
		<< '\t' << "4. 270 g" << '\n'
		<< '\t' << "5. 280 g" << '\n';

	std::cin >> mass_option;

	switch (mass_option) {
		case '1':
			mass = 0.1;
			initial_position = -0.118;
			final_time = 11.180;
			return "100";
		case '2':
			mass = 0.2;
			initial_position = -7.53e-2;
			final_time = 12.110;
			return "200";
		case '3':
			mass = 0.25;
			initial_position = -1.187e-2;
			final_time = 10.610;
			return "250";
		case '4':
			mass = 0.27;
			initial_position = -5.644e-2;
			final_time = 10.160;
			return "270";
		case '5':
			mass = 0.28;
			initial_position = -5.300e-2;
			final_time = 8.740;
			return "280";
		default:
			std::cout << "Invalid input." << '\n';
			return MassChoice();
			break;
	}
}

void MethodData(std::string& mass_as_string) {
	char method_choice;

	std::cout << "Select numerical method to compare with experimental and analytical data:" << '\n'
		<< '\t' << "1. Euler" << '\n'
		<< '\t' << "2. Heun" << '\n'
		<< '\t' << "3. Runge-Kutta (Order 4)" << '\n';

	std::cin >> method_choice;

	switch (method_choice) {
		case '1':
			method_name = "Euler";
			DataLoop(mass_as_string, method_name, kInitialTime, final_time, initial_position, kInitialVelocity, EulerImplementation);
			break;
		case '2':
			method_name = "Heun";
			DataLoop(mass_as_string, method_name, kInitialTime, final_time, initial_position, kInitialVelocity, HeunImplementation);
			break;
		case '3':
			method_name = "Runge-Kutta";
			DataLoop(mass_as_string, method_name, kInitialTime, final_time, initial_position, kInitialVelocity, RungeKuttaImplementation);
			break;
		default:
			std::cout << "Invalid input." << '\n';
			return MethodData(mass_as_string);
			break;
	}
}

void PlotData(std::string& mass_as_string, std::string& method_name) {
	std::string plot_file_name = "./Plot-Data/" + method_name + "-" + mass_as_string + ".gnu";

	std::ofstream plotfile(plot_file_name);

	plotfile << "set grid" << '\n'
		<< "set xlabel 't [s]'" << '\n'
		<< "set auto xy" << '\n'
		<< "set term pdf" << '\n'
		<< "set output './Plot-Data/plot-" << method_name << "-" << mass_as_string << ".pdf'" << '\n';

	plotfile << "set tit 'Position'" << '\n'
		<< "set ylabel 'x [m]'" << '\n'
		<< "p './Approx-Data/" << method_name << "-" << mass_as_string << ".dat' u 1:2 w l tit '" << method_name << "',"
		<< "'./Approx-Data/Analytic-" << mass_as_string << ".dat' u 1:2 w l tit 'Analytic',"
		<< "'./Experimental-Data/experimental-" << mass_as_string << ".dat' u 1:2 w l tit 'Experimental'" << '\n';

	plotfile << "set tit 'Velocity'" << '\n'
		<< "set ylabel 'v [ms^{-1}]'" << '\n'
		<< "p './Approx-Data/" << method_name << "-" << mass_as_string << ".dat' u 1:3 w l tit '" << method_name << "',"
		<< "'./Approx-Data/Analytic-" << mass_as_string << ".dat' u 1:3 w l tit 'Analytic',"
		<< "'./Experimental-Data/experimental-" << mass_as_string << ".dat' u 1:3 w l tit 'Experimental'" << '\n';

	long double E = Energy(initial_position, kInitialVelocity);

	plotfile << "set tit 'Energy'" << '\n'
		<< "set ylabel 'E [J]'" << '\n'
		<< "E = " << E << '\n'
		<< "p './Approx-Data/" << method_name << "-" << mass_as_string << ".dat' u 1:4 w l tit '" << method_name << "',"
		<< "'./Approx-Data/Analytic-" << mass_as_string << ".dat' u 1:4 w l tit 'Analytic',"
		<< "E tit 'Experimental'" << '\n'
		<< "exit";
}
