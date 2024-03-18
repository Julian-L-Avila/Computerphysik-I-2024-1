#include<iostream>
#include<cmath>
#include<fstream>

// Function to integrate
double function(double x) {
	return 100 * pow(x, 2) * cos(20 * x);
}

// Approximate integral using rectangular quadrature
double rectangular(double a, double b, int N) {
	double h = (b - a) / N;
	double sum = 0.0;
	for (int i = 0; i < N; ++i) {
		double xi = a + i * h;
		sum += function(xi);
	}
	return h * sum;
}

int main() {

	std::ofstream datafile("resultsRectangular.dat");

	for (int N = 1; N <= 10000; ++N){
		double integral = rectangular(0.0, 1.0, N);
		double exact = (199 * std::sin(20) + 20 * std::cos(20)) / 40;
		
		datafile << N << '\t' << integral << '\t' << exact << std::endl;
	}

	datafile.close();

	std::ofstream scriptFile1("numerical-integration");
	scriptFile1 << "set term pdfcairo\n";
	scriptFile1 << "set output 'Approximate-Integral.pdf'\n";
	scriptFile1 << "set grid\n";
	scriptFile1 << "set xlabel 'N'\n";
	scriptFile1 << "set ylabel 'Numerical Approx'\n";
	scriptFile1 << "set logscale x\n";
	scriptFile1 << "p 'results.dat' u 1:2 w l tit 'Numerical Integral', 'results.dat' u 1:3 w l tit 'Exact Integral'\n";
	scriptFile1.close();

	std::ofstream datafile2("dataRectangular.dat");

	int N_plot = 10;
	double h_plot = 1.0 / N_plot;

	for (double x = 0.0; x <= 1.0; x += h_plot) {
		datafile2 << x << '\t' << function(x) << std::endl;
		datafile2 << x << '\t' << 0.0 << std::endl;
		datafile2 << std::endl;
	}

	datafile2.close();

	std::ofstream scriptFile2("rectangles");
	scriptFile2 << "set term pdfcairo\n";
	scriptFile2 << "set output 'Rectangles.pdf'\n";
	scriptFile2 << "set grid\n";
	scriptFile2 << "set xlabel 'x'\n";
	scriptFile2 << "set ylabel 'f(x)'\n";
	scriptFile2 << "p 'dataRectangular.dat' u 1:2 w l tit 'Function', 'dataRectangular.dat' u 1:2 w boxes tit 'Rectangles'\n";
	scriptFile2.close();

	system("gnuplot -p 'numerical-integration'");
	system("gnuplot -p 'rectangles'");

	return 0;
}
