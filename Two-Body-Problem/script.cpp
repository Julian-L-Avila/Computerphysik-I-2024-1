#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>

const bool   ActiveLog       = true;

const double SemimajorAxis   = 5.0;
const double Mass1           = 1.0; // Solar Masses
const double Mass2           = 0.9; // Solar Masses
const double Eccentricity    = 0.5;

const double RectangleNumber  = 1000;
const int    desiredPrecision = 200;

const double GravitationalConstant = 4 * M_PI * M_PI;
const double ReducedMass           = (Mass1 * Mass2) / (Mass1 + Mass2);
const double PotentialConstant     = GravitationalConstant * Mass1 * Mass2;
const double Energy                = PotentialConstant / (-2 * SemimajorAxis);
const double AngularMomentum       = std::sqrt((Eccentricity * Eccentricity - 1) * ReducedMass * PotentialConstant * -SemimajorAxis);
const double SquareAngularMomentum = AngularMomentum * AngularMomentum;
const double InitialDistance       = SemimajorAxis * (1 - Eccentricity * Eccentricity) / (1 + Eccentricity);
const double FinalDistance         = SemimajorAxis * (1 - Eccentricity * Eccentricity) / (1 - Eccentricity);

const double StepSize = Eccentricity / RectangleNumber;

long double analyticFunction(double position) {
	return - std::asin((1.0 / position - ReducedMass * GravitationalConstant / SquareAngularMomentum) * std::sqrt(SquareAngularMomentum * SquareAngularMomentum / (2 * ReducedMass * Energy * SquareAngularMomentum + ReducedMass * ReducedMass * PotentialConstant * PotentialConstant)));
}

long double positionIntegralFunction(double position) {
	return AngularMomentum / position * std::sqrt(2 * ReducedMass * position * (Energy * position + PotentialConstant) - SquareAngularMomentum);
}

long double rectangularMethod(double lowerLimit, double upperLimit, double stepSize, long double (*function)(double x)) {
	stepSize = (upperLimit - lowerLimit) / RectangleNumber;
	long double sum = 0.0;

	for (double r = lowerLimit; r <= upperLimit; r += stepSize) {
		sum += function(r);
	}
	return stepSize * sum;
}

long double trapezoidMethod(double lowerLimit, double upperLimit, double stepSize, long double (*function)(double x)) {
	stepSize = (upperLimit - lowerLimit) / RectangleNumber;
	long double sum = 0.5 * (function(lowerLimit) + function(upperLimit));

	for (int i = 1; i <= RectangleNumber; i++) {
		double r = lowerLimit + i * stepSize;
		sum += function(r);
	}
	return stepSize * sum;
}

long double simpsonMethod(double lowerLimit, double upperLimit, double stepSize, long double (*function)(double x)) {
	stepSize = (upperLimit - lowerLimit) / RectangleNumber;
	double halfStepSize = 0.5 * stepSize;
	long double sum = function(lowerLimit) + function(upperLimit);

	for (double r = lowerLimit; r <= upperLimit - stepSize; r += stepSize) {
		sum += 4 * function(r + halfStepSize) + 2 * function(r);
	}
	return halfStepSize * sum / 3.0;
}

long double gaussMethod(double lowerLimit, double upperLimit, double stepSize, long double (*function)(double x)) {
	stepSize = (upperLimit - lowerLimit) / RectangleNumber;
	double innerCoefficient1 = (upperLimit - lowerLimit) * 0.5;
	double innerCoefficient2 = (upperLimit + lowerLimit) * 0.5;
	double r1 = 1.0 / std::sqrt(3);
	double r2 = -r1;

	return innerCoefficient1 * (function(innerCoefficient1 * r2 + innerCoefficient2) + function(innerCoefficient1 * r1 + innerCoefficient2));
}

void log() {
	std::time_t currentTime = std::time(nullptr);

  std::string timeString = std::ctime(&currentTime);

	std::ofstream datafile("script.log");
	
	datafile << std::setprecision(desiredPrecision)
		<< "Time and Date of Simulation              = " << timeString << '\n'
		<< "-------------------------- Constants --------------------------" << '\n'
		<< "Initial distances (AU)                   = " << InitialDistance << '\n'
		<< "Eccentricity                             = " << Eccentricity << '\n'
		<< "Mass Object 1 (Solar Mass)               = " << Mass1 << '\n'
		<< "Mass Object 2 (Solar Mass)               = " << Mass2 << '\n'
		<< "Orbital Energy (AU² * Solar Mass / yr²)  = " << Energy << '\n'
		<< "Angular Momentum (AU² * Solar Mass / yr) = " << AngularMomentum << '\n'
		<< '\n'
		<< "Number of rectangles                     = " << RectangleNumber << '\n'
		<< "Desired Precision                        = " << desiredPrecision << '\n'
		<< '\n'
		<< "Gravitational Constant                   = " << GravitationalConstant << '\n'
		<< "Reduced Mass                             = " << ReducedMass << '\n'
		<< "Potential Gravitational Constant         = " << PotentialConstant << '\n'
		<< '\n'
		<< "Step size of numerical method            = " << StepSize / RectangleNumber << '\n'
		<< '\n'
		<< "Problem Part            = " << 2 * ReducedMass * InitialDistance * (Energy * InitialDistance + PotentialConstant) << '\n'
		<< "Angular Momentum Square = " << SquareAngularMomentum
		<< std::endl;

	datafile.close();
}

int main() {
	double finalfowardr, finalthetaRectanglar, finalthetaTrapezoid, finalthetaSimpson, finalthetaGauss;

	std::ofstream datafileR("datrectangular.dat");
	std::ofstream errorR("daterrorR.dat");
	std::ofstream datafileT("dattrapezoid.dat");
	std::ofstream errorT("daterrorT.dat");
	std::ofstream datafileS("datsimpson.dat");
	std::ofstream errorS("daterrorS.dat");
	std::ofstream datafileG("datgauss.dat");
	std::ofstream errorG("daterrorG.dat");

	std::ofstream datafileA("datanalytic.dat");

	std::cout << "Started Script" << std::endl;

	if (ActiveLog) {
		log();
		std::cout << ".log file created" << std::endl;
	}

	std::cout << "Initial r = " << InitialDistance << std::endl;

	for (double r = InitialDistance; r <= FinalDistance; r += StepSize) {
		long double deltathetaRectangular, deltathetaTrapezoid, deltathetaSimpson, deltathetaGauss, deltatheta;

		deltathetaRectangular += rectangularMethod(r, r + StepSize, StepSize, positionIntegralFunction);
		deltathetaTrapezoid   += trapezoidMethod(r, r + StepSize, StepSize, positionIntegralFunction);
		deltathetaSimpson     += simpsonMethod(r, r + StepSize, StepSize, positionIntegralFunction);
		deltathetaGauss       += gaussMethod(r, r + StepSize, StepSize, positionIntegralFunction);

		deltatheta = analyticFunction(r + StepSize);

		if (std::isnan(deltathetaRectangular)) {
			std::cout << "Error at r = " << r << '\n'
				<< "Problem Part = " << 2 * ReducedMass * r * (Energy * r + PotentialConstant) << std::endl;
			break;
		}

		datafileR << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaRectangular << '\t' 
			<< analyticFunction(deltathetaRectangular) << '\t'
			<< (r + StepSize) * std::cos(deltathetaRectangular) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaRectangular) 
			<< std::endl;

		errorR << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t'
			<< deltathetaRectangular - analyticFunction(r + StepSize)
			<< std::endl;

		datafileT << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaTrapezoid << '\t' 
			<< analyticFunction(deltathetaTrapezoid) << '\t'
			<< (r + StepSize) * std::cos(deltathetaTrapezoid) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaTrapezoid) 
			<< std::endl;

		errorT << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t'
			<< deltathetaTrapezoid - analyticFunction(r + StepSize)
			<< std::endl;

		datafileS << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaSimpson << '\t' 
			<< analyticFunction(deltathetaSimpson) << '\t'
			<< (r + StepSize) * std::cos(deltathetaSimpson) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaSimpson) 
			<< std::endl;

		errorS << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t'
			<< deltathetaSimpson - analyticFunction(r + StepSize)
			<< std::endl;

		datafileG << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaGauss << '\t' 
			<< analyticFunction(deltathetaGauss) << '\t'
			<< (r + StepSize) * std::cos(deltathetaGauss) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaGauss) 
			<< std::endl;

		errorG << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t'
			<< deltathetaGauss - analyticFunction(r + StepSize)
			<< std::endl;

		datafileA << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t'
			<< deltatheta << '\t'
			<< (r + StepSize) * std::cos(deltatheta) << '\t'
			<< (r + StepSize) * std::sin(deltatheta)
			<< std::endl;

	}

	datafileR.close();
	datafileT.close();
	datafileS.close();
	datafileG.close();
	datafileA.close();

	errorR.close();
	errorT.close();
	errorS.close();
	errorG.close();

	std::ofstream analyticplot("analytic.gnu");

	analyticplot << "set term pdfcairo" << '\n'
		<< "set output 'analytic.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set xlabel 'x (AU)'" << '\n'
		<< "set ylabel 'y (AU)'" << '\n'
		<< "set tit 'Analytic Solution'" << '\n'
		<< "p 'datanalytic.dat' u 3:4 w l tit 'Analytic solution'" 
		<< std::endl;

	analyticplot.close();

	std::ofstream rectangularplot("rectangular.gnu");

	rectangularplot << "set term pdfcairo" << '\n'
		<< "set output 'rectangular.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set xlabel 'x (AU)'" << '\n'
		<< "set ylabel 'y (AU)'" << '\n'
		<< "set tit 'Rectangular Rule'" << '\n'
		<< "p 'datrectangular.dat' u 4:5 w l tit 'Rectangular Approx.'" 
		<< std::endl;

	rectangularplot.close();

	std::ofstream trapezoidplot("trapezoid.gnu");

	trapezoidplot << "set term pdfcairo" << '\n'
		<< "set output 'trapezoid.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set xlabel 'x (AU)'" << '\n'
		<< "set ylabel 'y (AU)'" << '\n'
		<< "set tit 'Trapezoidal Rule'" << '\n'
		<< "p 'dattrapezoid.dat' u 4:5 w l tit 'Trapezoid Approx.'" 
		<< std::endl;

	trapezoidplot.close();

	std::ofstream simpsonplot("simpson.gnu");

	simpsonplot << "set term pdfcairo" << '\n'
		<< "set output 'simpson.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set xlabel 'x (AU)'" << '\n'
		<< "set ylabel 'y (AU)'" << '\n'
		<< "set tit 'Simpson Rule'" << '\n'
		<< "p 'datsimpson.dat' u 4:5 w l tit 'Simpson Approx.'" 
		<< std::endl;

	simpsonplot.close();

	std::ofstream gaussplot("gauss.gnu");

	gaussplot << "set term pdfcairo" << '\n'
		<< "set output 'gauss.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set xlabel 'x (AU)'" << '\n'
		<< "set ylabel 'y (AU)'" << '\n'
		<< "set tit 'Gaussian Quadrature'" << '\n'
		<< "p 'datgauss.dat' u 4:5 w l tit 'Gaussian Approx.'" 
		<< std::endl;

	gaussplot.close();

	system("gnuplot -p 'analytic.gnu'");
	system("gnuplot -p 'rectangular.gnu'");
	system("gnuplot -p 'trapezoid.gnu'");
	system("gnuplot -p 'simpson.gnu'");
	system("gnuplot -p 'gauss.gnu'");

	return 0;
}
