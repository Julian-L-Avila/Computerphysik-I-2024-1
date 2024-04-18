#include<iostream>

const double Initial_t = 0.0;
const double Initial_y = 1.0;
const double Final_t   = 3.0;

const double Step_Size = 0.1;

double first_derivative(double t, double y) {
	return t*t + y;
}

double euler_method(double previous_t, double previous_y, double step_size, double (*derivative)(double t, double y)) {
	return previous_y + step_size * derivative(previous_t, previous_y);
}

double heun_method(double previous_t, double previous_y, double step_size, double (*derivative)(double t, double y)) {
	double euler_y = euler_method(previous_t, previous_y, step_size, derivative);
	return previous_y + 0.5 * step_size * (derivative(previous_t, previous_y) + derivative(previous_t + step_size, euler_y));
}

double iteration_loop(double initial_t, double final_t, double initial_y, double step_size, double (*derivative)(double t, double y)) {
	double previous_y = initial_y;
	for (double time = initial_t; time <= final_t; time += step_size) {
		previous_y = heun_method(time, previous_y, step_size, derivative);
	}
	return previous_y;
}

int main() {
	std::cout << iteration_loop(Initial_t, Final_t, Initial_y, Step_Size, first_derivative) << std::endl;
	return 0;
}
