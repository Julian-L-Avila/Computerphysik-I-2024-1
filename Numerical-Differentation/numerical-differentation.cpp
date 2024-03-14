// Julian Leonardo Avila Martinez - 20212107030
// Laura Yeraldin Herrera Martinez - 20212107011
#include<iostream>
#include<fstream>
#include<cmath>

// Function
double function(double x) {
	return x * x * x;
}

// Annalytic Derivative
double analyticDerivative(double x) {
	return 3 * x * x;
}

// Forward Differences (2 points)
double forwardDifference(double x, double h) {
	return ((function(x + h) - function(x)) / h);
}

// Backward Differences (2 points)
double backwardDifference(double x, double h) {
	return ((function(x) - function(x - h)) / h);
}

// Central Differences (3 points)
double centralDifference(double x, double h) {
	return ((function(x + h) - function(x - h)) / (2 * h));
}

// Forward Differences (3 points)
double forwardDifference3(double x, double h) {
	return ((4 * function(x + h) - function(x + 2 * h) - 3 * function(x)) / (2 * h));
}

// Backward Differences (3 points)
double backwardDifference3(double x, double h) {
	return ((3 * function(x) - 4 * function(x - h) + function(x - 2 * h)) / (2 * h));
}

// Central Differences (4 points)
double centralDifference4(double x, double h) {
	return ((-function(x + 2 * h) + 8 * function(x + h) - 8 * function(x - h) + function(x - 2 * h)) / (12 * h));
}

int main() {
	const int points = 100;
	const double step = 0.1;

	// Open File
	std::ofstream datafile("derivative.dat");

	for(int i = 0; i < points; i++) {
		double x = static_cast<double>(i) / (points - 1);

		double annalyticResult = analyticDerivative(x);
		double forwardResult = forwardDifference(x, step);
		double backwardResult = backwardDifference(x, step);
		double centralResult = centralDifference(x, step);
		double forwardResult2 = forwardDifference3(x, step);
		double backwardResult2 = backwardDifference3(x, step);
		double centralResult2 = centralDifference4(x, step);

		// Save on File
		datafile << x << "\t" << annalyticResult << "\t" << forwardResult << "\t" << backwardResult << "\t" << centralResult << "\t" << forwardResult2 << "\t" << backwardResult2 << "\t" << centralResult2 << '\n';
	}

	// Closed File
	datafile.close();

	// GNUplot
	std::ofstream gnuplotScript("plotScript");
	gnuplotScript << "set grid" << '\n';
	gnuplotScript << "p 'derivative.dat' u 1:2 w l tit 'Annalytic',"
		<< "'' u 1:3 w l tit 'Forward',"
		<< "'' u 1:4 w l tit 'Backward',"
		<< "'' u 1:5 w l tit 'Central',"
		<< "'' u 1:6 w l tit 'Forward (3 points)',"
		<< "'' u 1:7 w l tit 'Backward (3 points)',"
		<< "'' u 1:8 w l tit 'Central (4 points)',";

	gnuplotScript.close();

	system("gnuplot -p plotScript");

	return 0;
}
