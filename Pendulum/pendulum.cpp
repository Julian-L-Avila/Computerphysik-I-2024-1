/* Pendulum and Double Pendulum Simulation for Computational Physics 2nd Exam
 * UD FCJ - FCMN - Physics
 * Made by:
 * - Julian Avila - 20212107030
 * - Laura Herrera - 20212107011
 * - Bryan Martinez - 2021210708
 * - Juan Acuna - 20212107034
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/jacobi_elliptic.hpp>

long double step_size, gravity_acceleration, natural_frequency,
		natural_frequency_square, Amplitude, Shift, modulus, modulus_square,
		elliptic_integral, final_time;

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
void GetConstants(VariablesAt0& InitialConditions);
long double AngleDerivative(long double& angle_velocity);
long double VelocityDerivativeLinear(long double& angle);
long double AnalyticAngleLinear(double t);
long double AnalyticVelocityLinear(double t);
long double EnergySimple(long double& angle, long double& angle_velocity,
		double& mass);
long double VelocityDerivative(long double& angle);
long double EllipticIntegralFirstKind(long double& angle);
long double PeriodNotLinear(long double& angle);
long double SinusAmplitudis(long double& angle, double& t);
long double CosinusAmplitudis(long double& angle, double& t);
long double DeltaAmplitudis(long double& angle, double& t);
long double AnalyticAngleNonLinear(int& degree, double& t,
		long double& initial_angle);
long double AnalyticVelocityNonLinear(double& t,
		long double& initial_angle);
long double AngleDerivativeDouble(double& t, long double& angle,
		long double& angle_velocity);
long double VelocityDerivativeDouble1(double& t, long double& angle_1,
		long double& angle_2, long double& angle_velocity_1,
		long double& angle_velocity_2);
long double VelocityDerivativeDouble2(double& t, long double& angle_1,
		long double& angle_2, long double& angle_velocity_1,
		long double& angle_velocity_2);
long double EnergyDouble(long double& angle_1, long double& angle_2,
		long double& angle_velocity_1, long double& angle_velocity_2);

int main() {
	GetConstants(InitialConditions);
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

void GetConstants(VariablesAt0& InitialConditions) {
	InitialConditions.angle_1 = 5 * M_PI/180;
	InitialConditions.angle_velocity_1 = 0.0;
	InitialConditions.length_1 = 50.0e-2;
	InitialConditions.mass_1 = 1.0;
	LinearConstans(InitialConditions);
	NonLinearConstants(InitialConditions);
}

// Pendulum (Linear)

long double AngleDerivative(double& t, long double& angle,
		long double& angle_velocity) {
	return angle_velocity;
}

long double VelocityDerivativeLinear(double& t, double& angle_velocity,
		long double& angle) {
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

long double EnergySimple(long double& angle, long double& angle_velocity) {
	double length = InitialConditions.length_1;
	double mass   = InitialConditions.mass_1;

	long double y_position = - length * std::cos(angle);

	long double K = length * length *
		angle_velocity * angle_velocity / 4.0;
	long double V = mass * gravity_acceleration * y_position;

	return K + V;
}

// Pendulum (Non-Linear)

long double VelocityDerivative(double& t, long double& angle,
		long double& angle_velocity) {
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
		long double& angle_velocity_2) {
	double mass_1 = InitialConditions.mass_1;
	double mass_2 = InitialConditions.mass_2;
	double length_1 = InitialConditions.length_1;
	double length_2 = InitialConditions.length_2;

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
		long double& angle_velocity_2) {
	double mass_1 = InitialConditions.mass_1;
	double mass_2 = InitialConditions.mass_2;
	double length_1 = InitialConditions.length_1;
	double length_2 = InitialConditions.length_2;

	long double angle_diff = angle_1 - angle_2;
	long double sin_of_diff = std::sin(angle_diff);
	long double angle_velocity_1_square = angle_velocity_1 * angle_velocity_1;
	long double angle_velocity_2_square = angle_velocity_2 * angle_velocity_2;

	long double numerator = 18.0 * std::cos(angle_diff) *
		(mass_2 * length_2 * angle_velocity_2_square * sin_of_diff +
		(mass_1 + 2.0 * mass_2) * gravity_acceleration * std::sin(angle_1)) -
		12.0 * (mass_1 + 3.0 * mass_2) * (gravity_acceleration * std::sin(angle_2) -
		length_1 * angle_velocity_1_square + sin_of_diff);

	long double denominator = length_2 * (15.0 * mass_2 + 8.0 * mass_1 - 9.0 *
			mass_2 * std::cos(2.0 * angle_diff));

	return numerator / denominator;
}

long double EnergyDouble(long double& angle_1, long double& angle_2,
		long double& angle_velocity_1, long double& angle_velocity_2) {
	double mass_1 = InitialConditions.mass_1;
	double mass_2 = InitialConditions.mass_2;
	double length_1 = InitialConditions.length_1;
	double length_2 = InitialConditions.length_2;

	long double k_1 = mass_1 * length_1 *
		(length_1 * angle_velocity_1 * angle_velocity_1 / 3.0 -
			gravity_acceleration * std::cos(angle_1));

	long double k_2 = mass_2 * length_2 * (length_2 * (angle_velocity_1 *
			angle_velocity_1 + angle_velocity_2 * angle_velocity_2 / 3.0 +
			angle_velocity_1 * angle_velocity_2 * std::cos(angle_1 - angle_2)) -
			gravity_acceleration * (2.0 * std::cos(angle_1) - std::cos(angle_2)));

	return k_1 + k_2;
}

long double EulerMethod(double& t, long double& previous_x,
		long double& previous_y, long double (*Derivative)(double&, long double&,
		long double&)) {
	return t + step_size * Derivative(t, previous_x, previous_y);
	}

void EulerImplementation(double& t, long double& previous_x,
		long double& previous_y, long double& x, long double& y,
		long double (*Derivative)(double&, long double&, long double&),
		long double (*Derivative1)(double&, long double&, long double&)) {
	x = EulerMethod(t, previous_x, previous_y, Derivative);
	y = EulerMethod(t, previous_y, previous_x, Derivative1);
}

void HeunImplementation(double& t, long double& previous_x,
		long double& previous_y, long double& x, long double& y,
		long double (*Derivative)(double&, long double&, long double&),
		long double (*Derivative1)(double&, long double&, long double&)) {
	long double euler_x = EulerMethod(t, previous_x, previous_y, Derivative);
	long double euler_y = EulerMethod(t, previous_y, previous_x, Derivative1);

	x = previous_x + 0.5 * step_size * (Derivative(t, previous_x, previous_y) + Derivative(t, euler_x, euler_y));
	y = previous_y + 0.5 * step_size * (Derivative1(t, previous_y, previous_x) + Derivative1(t, euler_y, euler_x));
}

void RungeKuttaImplementation(double& t, long double& previous_x,
		long double& previous_y, long double& x, long double& y,
		long double (*Derivative)(double&, long double, long double),
		long double (*Derivative1)(double&, long double, long double)) {
	long double k1 = Derivative(t, previous_x, previous_y);
	long double m1 = Derivative1(t, previous_y, previous_x);
	long double k2 = Derivative(t, previous_x + 0.5 * k1 * step_size, previous_y + 0.5 * m1 * step_size);
	long double m2 = Derivative1(t, previous_y + 0.5 * m1 * step_size, previous_x + 0.5 * k1 * step_size);
	long double k3 = Derivative(t, previous_x + 0.5 * k2 * step_size, previous_y + 0.5 * m2 * step_size);
	long double m3 = Derivative1(t, previous_y + 0.5 * m2 * step_size, previous_x + 0.5 * k2 * step_size);
	long double k4 = Derivative(t, previous_x + k3 * step_size, previous_y + m3 * step_size);
	long double m4 = Derivative1(t, previous_y + step_size * m3, previous_x + step_size * k4);

	x = previous_x + step_size * (k1 + 2.0 * (k2 + k3) + k4) / 6.0;
	y = previous_y + step_size * (m1 + 2.0 * (m2 + m3) + m4) / 6.0;
}

long double AbsoluteError(long double actual, long double estimated) {
	return std::abs(actual - estimated);
}

void DataLoop(std::string& mass_as_string, std::string& method_name, double initial_time, double final_time, double initial_position, double initial_velocity, void(*MethodLoop) (double&, long double&, long double&, long double&, long double&, long double&)) {
	long double previous_x, previous_v, energy, x, v;
	previous_x = initial_position;
	previous_v = initial_velocity;
	energy = EnergySimple(previous_x, previous_v);

	std::string path_file_name = "./Approx-Data/" + method_name + '-' + mass_as_string + ".dat";

	std::ofstream datafile(path_file_name);

	datafile << "# " + method_name + "Data" << '\n'
		<< "#Time (s) \t Position (m) \t Velocity (ms^-1) \t Energy (J)" << '\n'
		<< std::setprecision(kDesiredPrecision) << std::fixed
		<< initial_time << '\t' << initial_position << '\t' << initial_velocity << '\t'
		<< energy << '\n';

	for (double t = initial_time + step_size; t <= final_time; t += step_size) {
		MethodLoop(t, previous_x, previous_v, energy, x, v);
		energy = EnergySimple(previous_x, previous_v);

		datafile << t << '\t' << x << '\t' << v << '\t' << energy << '\n';
		previous_x = x;
		previous_v = v;
	}

	datafile.close();
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

void ErrorData(std::string& mass_as_string, std::string& method_name) {
	std::string experimental_file_name = "./Experimental-Data/experimental-" + mass_as_string + ".dat";
	std::string analytical_file_name   = "./Approx-Data/Analytic-" + mass_as_string + ".dat";
	std::string numerical_file_name    = "./Approx-Data/" + method_name + '-' + mass_as_string + ".dat";
	std::string error_file_name        = "./Error-Data/Error-" + method_name + '-' + mass_as_string + ".dat";

	std::ifstream experimental_data(experimental_file_name);
	std::ifstream analytical_data(analytical_file_name);
	std::ifstream numerical_data(numerical_file_name);
	std::ofstream error_data(error_file_name);

	error_data << "# Error Data " << method_name << ' ' << mass_as_string << '\n'
		<< "# Time" << '\t' <<  "Position (Ex-An)" << '\t' << "Position (Ex-Nu)" << '\t' << "Position (An-Nu)" << '\t'
		<< "Velocity (Ex-An)" << '\t' << "Velocity (Ex-Nu)" << '\t' << "Velocity (An-Nu)" << '\t'
		<< "Energy (An-Nu)" << '\n';


	long double x_ex, v_ex, x_an, v_an, e_an, x_nu, v_nu, e_nu, t;

	std::string line_ex, line_an, line_nu;
	while (std::getline(experimental_data, line_ex)) {
		std::getline(analytical_data, line_an);
		std::getline(numerical_data, line_nu);

		std::stringstream ex(line_ex);
		std::stringstream an(line_an);
		std::stringstream nu(line_nu);

		ex >> t >> x_ex >> v_ex;
		an >> t >> x_an >> v_an >> e_an;
		nu >> t >> x_nu >> v_nu >> e_nu;

		error_data << std::setprecision(kDesiredPrecision) << t << '\t'
			<< AbsoluteError(x_ex, x_an) << '\t' << AbsoluteError(x_ex, x_nu) << '\t' << AbsoluteError(x_an, x_nu) << '\t'
			<< AbsoluteError(v_ex, v_an) << '\t' << AbsoluteError(v_ex, v_nu) << '\t' << AbsoluteError(v_an, v_nu) << '\t'
			<< AbsoluteError(e_an, e_nu) << '\n';
	}

	experimental_data.close();
	analytical_data.close();
	numerical_data.close();
	error_data.close();
}

void PlotData(std::string& mass_as_string, std::string& method_name) {
	std::string plot_file_name = "./Plot-Data/" + method_name + '-' + mass_as_string + ".gnu";

	std::ofstream plotfile(plot_file_name);

	plotfile << "set grid" << '\n'
		<< "set xlabel 't [s]'" << '\n'
		<< "set auto xy" << '\n'
		<< "set term pdf" << '\n'
		<< "set output './Plot-Data/plot-" << method_name << '-' << mass_as_string << ".pdf'" << '\n';

	plotfile << "set tit 'Position'" << '\n'
		<< "set ylabel 'x [m]'" << '\n'
		<< "set auto xy" << '\n'
		<< "p './Approx-Data/" << method_name << '-' << mass_as_string << ".dat' u 1:2 w l tit '" << method_name << "',"
		<< "'./Approx-Data/Analytic-" << mass_as_string << ".dat' u 1:2 w l tit 'Analytic',"
		<< "'./Experimental-Data/experimental-" << mass_as_string << ".dat' u 1:2 w l tit 'Experimental'" << '\n';

	plotfile << "set tit 'Velocity'" << '\n'
		<< "set ylabel 'v [ms^{-1}]'" << '\n'
		<< "set auto xy" << '\n'
		<< "p './Approx-Data/" << method_name << '-' << mass_as_string << ".dat' u 1:3 w l tit '" << method_name << "',"
		<< "'./Approx-Data/Analytic-" << mass_as_string << ".dat' u 1:3 w l tit 'Analytic',"
		<< "'./Experimental-Data/experimental-" << mass_as_string << ".dat' u 1:3 w l tit 'Experimental'" << '\n';

	long double E = Energy(initial_position, kInitialVelocity);

	plotfile << "set tit 'Energy'" << '\n'
		<< "set ylabel 'E [J]'" << '\n'
		<< "set auto xy" << '\n'
		<< "E = " << E << '\n'
		<< "p './Approx-Data/" << method_name << '-' << mass_as_string << ".dat' u 1:4 w l tit '" << method_name << "',"
		<< "'./Approx-Data/Analytic-" << mass_as_string << ".dat' u 1:4 w l tit 'Analytic',"
		<< "E tit 'Experimental'" << '\n';

	plotfile << "set tit 'Position Error'" << '\n'
		<< "set ylabel 'Absolute Error'" << '\n'
		<< "set key left" << '\n';

	if (method_name == "Runge-Kutta") {
		plotfile << "set key bottom right" << '\n';
	}

	plotfile << "set log y" << '\n'
		<< "set auto xy" << '\n'
		<< "p './Error-Data/Error-" << method_name << '-' << mass_as_string << ".dat' u 1:2 w l tit 'Exp. vs An.', "
		<< "'' u 1:3 w l tit 'Exp. vs Nu.', "
		<< "'' u 1:4 w l tit 'An vs Nu.'" << '\n';

	plotfile << "set tit 'Velocity Error'" << '\n'
		<< "set auto xy" << '\n'
		<< "p './Error-Data/Error-" << method_name << '-' << mass_as_string << ".dat' u 1:5 w l tit 'Exp. vs An.', "
		<< "'' u 1:6 w l tit 'Exp. vs Nu.', "
		<< "'' u 1:7 w l tit 'An vs Nu.'" << '\n';

	plotfile << "set tit 'Energy Error'" << '\n'
		<< "set auto xy" << '\n'
		<< "p './Error-Data/Error-" << method_name << '-' << mass_as_string << ".dat' u 1:8 w l tit 'An. vs Nu.'" << '\n'
		<< "exit";
}
