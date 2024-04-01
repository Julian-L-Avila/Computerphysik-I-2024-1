#include <iomanip>
#include <iostream>
#include <cmath>
#include <fstream>

const bool   ActiveLog       = true;
const double InitialDistance = 4.0; // Astronomical Units
const double InitialTheta    = 0.0; // Rad
const double Mass1           = 1.0; // Solar Masses
const double Mass2           = 0.9; // Solar Masses
const double Eccentricity    = 0.9;
const double SemimajorAxis   = 5.0;

const double RectangleNumber  = 10000;
const int    desiredPrecision = 10;

const double GravitationalConstant = 4 * M_PI * M_PI;
const double ReducedMass           = (Mass1 * Mass2) / (Mass1 + Mass2);
const double PotentialConstant     = GravitationalConstant * Mass1 * Mass2;
const double Energy                = PotentialConstant / (-2 * SemimajorAxis);
const double AngularMomentum       = std::sqrt((Eccentricity * Eccentricity - 1) * ReducedMass * PotentialConstant * PotentialConstant / (2 * Energy));
const double SquareAngularMomentum = AngularMomentum * AngularMomentum;

const double StepSize = 5e-3;

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

void log() {
	std::time_t currentTime = std::time(nullptr);

  std::string timeString = std::ctime(&currentTime);

	std::ofstream datafile("rectangular.log");
	
	datafile << std::setprecision(desiredPrecision)
		<< "Time and Date of Simulation             = " << timeString << '\n'
		<< "-------------------------- Constants --------------------------" << '\n'
		<< "Inital distances (AU)                   = " << InitialDistance << '\n'
    << "Initial Angle (rad)                     = " << InitialTheta << '\n'
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
	double finalfowardr;

	std::ofstream datafile("datrectangualr.dat");
	std::cout << "Started Script" << std::endl;

	if (ActiveLog) {
		log();
		std::cout << ".log file created" << std::endl;
	}

	std::cout << "Initial r = " << InitialDistance << std::endl;

	for (double r = InitialDistance + StepSize; r <= 100; r += StepSize) {
		double deltatheta;

		deltatheta += AngularMomentum * rectangularMethod(r - StepSize, r, StepSize, positionIntegralFunction);

		if (std::isnan(deltatheta)) {
			std::cout << "Error at r = " << r << '\n'
				<< "Problem Part = " << 2 * ReducedMass * r * (Energy * r + PotentialConstant) << std::endl;
			finalfowardr = r - StepSize;
			break;
		}

		datafile << std::setprecision(desiredPrecision)
			<< r << '\t' 
			<< deltatheta << '\t' 
			<< r * std::cos(deltatheta) << '\t' 
			<< r * std::sin(deltatheta) 
			<< std::endl;

	}

	for (double r = InitialDistance; 0 <= r; r -=StepSize) {
		double deltatheta;

		deltatheta -= AngularMomentum * rectangularMethod(r - StepSize, r, StepSize, positionIntegralFunction);

		if (std::isnan(deltatheta)) {
			std::cout << "Error at r = " << r << '\n'
				<< "Problem Part = " << 2 * ReducedMass * r * (Energy * r + PotentialConstant) << std::endl;
			finalfowardr = r - StepSize;
			break;
		}

		datafile << std::setprecision(desiredPrecision)
			<< r << '\t' 
			<< deltatheta << '\t' 
			<< r * std::cos(deltatheta) << '\t' 
			<< r * std::sin(deltatheta) 
			<< std::endl;
	}

	datafile.close();
	return 0;
}
