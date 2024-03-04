#include<iostream>
#include<cstdlib>
#include<cmath>

//Made by Julian L. Avila

double x0;
bool printapprox = true;
double tol = 1e-20;
int maxiter = 1000;

// f(x) = e^(2x) + ln(2x)
// f'(x) = 2e^(2x) + 1/x

double f(double x) {
	return std::exp(2 * x) + log(2 * x); 
}

double df(double x) {
	return 2 * std::exp(2 * x) + 1.0 / x;
}

double Newton_Raphson(double f(double), double df(double), double x0, double tol, int maxiter, bool printapprox) {
	if (tol <= 0 || maxiter <= 0) {
		throw std::invalid_argument("Tolerance and maximum iterations must be positive.");
	}

	double x = x0;
	for (int i = 1; i <= maxiter; ++i) {
		double fx = f(x);
		double dfx = df(x);

		if (std::abs(dfx) < tol) {
			throw std::runtime_error("Derivative is zero at iteration " + std::to_string(i));
		}

		x -= fx / dfx;

		if (std::abs(f(x)) < tol || std::abs(x - x0) < tol) {
			return x;
		}

		if (printapprox) {
			std::cout << "Approximation at iteration " << i << ": " << x << std::endl;
		}

		x0 = x;
	}

	throw std::runtime_error("Newton Raphson method failed to converge after " + std::to_string(maxiter) + " iterations.");
}

int main() {
	try {
		std::cout << "Enter the initial guess (x0): ";
		std::cin >> x0;

		double root = Newton_Raphson(f, df, x0, tol, maxiter, printapprox);
		std::cout << "The root is: " << root << std::endl;
	} catch (const std::exception& e) {
		std::cerr << "Error: " << e.what() << std::endl;
		return 1;
	}

	return 0;
}
