#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>

const bool   ActiveLog       = true;
const double SemimajorAxis   = 5.0;
const double Mass1           = 1.0; // Solar Masses
const double Mass2           = 0.9; // Solar Masses
const double Eccentricity    = 0.5;

const double RectangleNumber  = 10000;
const int    desiredPrecision = 10;

const double GravitationalConstant = 4 * M_PI * M_PI;
const double ReducedMass           = (Mass1 * Mass2) / (Mass1 + Mass2);
const double PotentialConstant     = GravitationalConstant * Mass1 * Mass2;
const double Energy                = PotentialConstant / (-2 * SemimajorAxis);
const double AngularMomentum       = std::sqrt((Eccentricity * Eccentricity - 1) * ReducedMass * PotentialConstant * PotentialConstant / (2 * Energy));
const double SquareAngularMomentum = AngularMomentum * AngularMomentum;
const double InitialDistance       = SemimajorAxis / (1 + Eccentricity);

const double StepSize = Eccentricity / RectangleNumber;

double positionIntegralFunction(double position) {
	return 1 / position * std::sqrt(2 * ReducedMass * position * (Energy * position + PotentialConstant) - SquareAngularMomentum);
}

double rectangularMethod(double lowerLimit, double upperLimit, double stepSize, double (*function)(double x)) {
	stepSize   = stepSize / RectangleNumber;
	double sum = 0.0;

	for (double r = lowerLimit; r <= upperLimit; r += stepSize) {
		sum += function(r);
	}
	return stepSize * sum;
}

double trapezoidMethod(double lowerLimit, double upperLimit, double stepSize, double (*function)(double x)) {
	stepSize = stepSize / RectangleNumber;
	double sum = 0.5 * (function(lowerLimit) + function(upperLimit));

	for (double r = lowerLimit + stepSize; r <= upperLimit - stepSize; r += stepSize) {
		sum += function(r);
	}
	return stepSize * sum;
}

double simpsonMethod(double lowerLimit, double upperLimit, double stepSize, double (*function)(double x)) {
	stepSize = stepSize / RectangleNumber;
	double halfStepSize = 0.5 * stepSize;
	double sum = function(lowerLimit) + function(upperLimit);

	for (double r = lowerLimit; r <= upperLimit - stepSize; r += stepSize) {
		sum += 4 * function(r + halfStepSize) + 2 * function(r);
	}
	return halfStepSize * sum / 3.0;
}

double gaussMethod(double lowerLimit, double upperLimit, double stepSize, double (*function)(double x)) {
	stepSize = stepSize / RectangleNumber;
	double innerCoefficient1 = (upperLimit - lowerLimit) * 0.5;
	double innerCoefficient2 = (upperLimit + lowerLimit) * 0.5;
	double r1 = 1.0 / std::sqrt(3);
	double r2 = -r1;

	return innerCoefficient1 * (function(innerCoefficient1 * r2 + innerCoefficient2) + function(innerCoefficient1 * r1 + innerCoefficient2));
}

void log() {
	std::time_t currentTime = std::time(nullptr);

  std::string timeString = std::ctime(&currentTime);

	std::ofstream datafile("rectangular.log");
	
	datafile << std::setprecision(desiredPrecision)
		<< "Time and Date of Simulation             = " << timeString << '\n'
		<< "-------------------------- Constants --------------------------" << '\n'
		<< "Inital distances (AU)                   = " << InitialDistance << '\n'
    << "Eccentricity                            = " << Eccentricity << '\n'
    << "Mass Object 1 (Solar Mass)              = " << Mass1 << '\n'
    << "Mass Object 2 (Solar Mass)              = " << Mass2 << '\n'
		<< "Obtital Energy (AU² * Solar Mass / yr²) = " << Energy << '\n'
		<< "AngularMomentum (AU² * Solar Mass / yr) = " << AngularMomentum << '\n'
		<< '\n'
		<< "Number of rectangles                    = " << RectangleNumber << '\n'
		<< "Desired Precision                       = " << desiredPrecision << '\n'
		<< '\n'
		<< "Gravitacional Constant                  = " << GravitationalConstant << '\n'
		<< "Reduced Mass                            = " << ReducedMass << '\n'
		<< "Potential Gravitational Constat         = " << PotentialConstant << '\n'
		<< '\n'
		<< "Step size of numerical method           = " << StepSize / RectangleNumber << '\n'
		<< '\n'
		<< "Problem Part = " << 2 * ReducedMass * InitialDistance * (Energy * InitialDistance + PotentialConstant) << '\n'
		<< "Square Mommentum = " << SquareAngularMomentum
		<< std::endl;

	datafile.close();
}

int main() {
	double finalfowardr, finalthetaRectanglar, finalthetaTrapezoid, finalthetaSimpson, finalthetaGauss;

	std::ofstream datafileR("datrectangular.dat");
	std::ofstream datafileT("dattrapezoid.dat");
	std::ofstream datafileS("datsimpson.dat");
	std::ofstream datafileG("datgauss.dat");

	std::cout << "Started Script" << std::endl;

	if (ActiveLog) {
		log();
		std::cout << ".log file created" << std::endl;
	}

	std::cout << "Initial r = " << InitialDistance << std::endl;

	for (double r = InitialDistance; r <= SemimajorAxis / (1 - Eccentricity); r += StepSize) {
		double deltathetaRectangular, deltathetaTrapezoid, deltathetaSimpson, deltathetaGauss;

		deltathetaRectangular += AngularMomentum * rectangularMethod(r, r + StepSize, StepSize, positionIntegralFunction);
		deltathetaTrapezoid   += AngularMomentum * trapezoidMethod(r, r + StepSize, StepSize, positionIntegralFunction);
		deltathetaSimpson     += AngularMomentum * simpsonMethod(r, r + StepSize, StepSize, positionIntegralFunction);
		deltathetaGauss       += AngularMomentum * gaussMethod(r, r + StepSize, StepSize, positionIntegralFunction);

		if (std::isnan(deltathetaRectangular)) {
			std::cout << "Error at r = " << r << '\n'
				<< "Problem Part = " << 2 * ReducedMass * r * (Energy * r + PotentialConstant) << std::endl;
			break;
		}

		datafileR << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaRectangular << '\t' 
			<< (r + StepSize) * std::cos(deltathetaRectangular) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaRectangular) 
			<< std::endl;

		datafileT << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaTrapezoid << '\t' 
			<< (r + StepSize) * std::cos(deltathetaTrapezoid) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaTrapezoid) 
			<< std::endl;

		datafileS << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaSimpson << '\t' 
			<< (r + StepSize) * std::cos(deltathetaSimpson) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaSimpson) 
			<< std::endl;

		datafileG << std::setprecision(desiredPrecision)
			<< r + StepSize << '\t' 
			<< deltathetaGauss << '\t' 
			<< (r + StepSize) * std::cos(deltathetaGauss) << '\t' 
			<< (r + StepSize) * std::sin(deltathetaGauss) 
			<< std::endl;

		finalfowardr         = r;
		finalthetaRectanglar = deltathetaRectangular;
		finalthetaTrapezoid  = deltathetaTrapezoid;
		finalthetaSimpson    = deltathetaSimpson;
		finalthetaGauss      = deltathetaGauss;
	}

	std::cout << "Initial r = " << finalfowardr << std::endl;

	for (double r = finalfowardr;  SemimajorAxis / (1 + Eccentricity) <= r; r -=StepSize) {
		double deltathetaRectangular, deltathetaTrapezoid, deltathetaSimpson, deltathetaGauss;

		deltathetaRectangular -= AngularMomentum * rectangularMethod(r - StepSize, r, StepSize, positionIntegralFunction);

		if (std::isnan(deltathetaRectangular)) {
			std::cout << "Error at r = " << r << '\n'
				<< "Problem Part = " << 2 * ReducedMass * r * (Energy * r + PotentialConstant) << std::endl;
			finalfowardr = r - StepSize;
			break;
		}

		datafileR << std::setprecision(desiredPrecision)
			<< r - StepSize << '\t' 
			<< deltathetaRectangular << '\t' 
			<< (r - StepSize) * std::cos(deltathetaRectangular) << '\t' 
			<< (r - StepSize) * std::sin(deltathetaRectangular) 
			<< std::endl;
	}

	datafileR.close();
	datafileT.close();
	datafileS.close();
	datafileG.close();

	std::ofstream analiticplot("analytic.gnu");

	analiticplot << "set term pdfcairo" << '\n'
		<< "set output 'analytic.pdf'" << '\n'
		<< "set grid" << '\n'
		<< "set parametric" << '\n'
		<< "set xlabel 'x (AU)'" << '\n'
		<< "set ylabel 'y (AU)'" << '\n'
		<< "set tit 'Analytic Solution'" << '\n'
		<< "a = " << SemimajorAxis << '\n'
		<< "e = " << Eccentricity << '\n'
		<< "r(t) = a / (1 + e*cos(t))" << '\n'
		<< "p (r(t)*cos(t)),(r(t)*sin(t)) w l tit 'Analytic solution'" 
		<< std::endl;

	system("gnuplot -p 'analytic.gnu'");
	return 0;
}
