#include <iomanip>
#include<iostream>
#include<cmath>
#include <fstream>

const double InitialDistance       = 5.20336301;
const double SemiMajorAxis         = 5.20336301;
const double Eccentricity          = 0.5;
const double InitialTheta          = 0.0;
const double Mass1                 = 1.0;
const double Mass2                 = 0.9;
const double GravitationalConstant = 4 * M_PI * M_PI;

const double RectangleNumber = 10000;
int desiredPrecision         = 10;

const double SemiMinorAxis     = SemiMajorAxis - Eccentricity;
const double ReducedMass       = (Mass1 * Mass2) / (Mass1 + Mass2);
const double PotentialConstant = GravitationalConstant * Mass1 * Mass2;
const double Energy            = - PotentialConstant / (2*SemiMajorAxis);
const double AngularMomentum   = std::sqrt((std::pow(Eccentricity, 2) - 1) * ReducedMass * std::pow(PotentialConstant, 2) / (2 * Energy));

const double StepSize = Eccentricity / RectangleNumber;

double positionIntegralFunction(double position) {
	return 1 / position * std::sqrt(2 * ReducedMass * position * (Energy * position + PotentialConstant) - std::pow(AngularMomentum, 2));
}

double timeIntegralFunction(double position) {
	return position / std::sqrt(2 * ReducedMass * position * (Energy * position + PotentialConstant) - std::pow(AngularMomentum, 2));
}

double rectangularMethod(double lowerLimit, double upperLimit, double stepSize, double (*function)(double x)) {
	stepSize   = stepSize / RectangleNumber;
	double sum = 0.0;

	for (int i = 0; i < RectangleNumber; ++i) {
		double center  = lowerLimit + i * stepSize;
		sum           += function(center);
	}
	return stepSize * sum;
}

int main() {
	std::ofstream datafile("datrectangualr.dat");

	for (double r = InitialDistance; SemiMinorAxis <= r; r -= StepSize) {
		double deltatheta, deltat, distance;

		deltatheta = AngularMomentum * rectangularMethod(r, SemiMajorAxis, StepSize, positionIntegralFunction);
		deltat     = ReducedMass * rectangularMethod(r, SemiMajorAxis, StepSize, timeIntegralFunction);
		distance   = std::sqrt(AngularMomentum * deltat / (ReducedMass * deltatheta));

		datafile << std::setprecision(desiredPrecision) << r << '\t' << deltatheta << '\t' << deltat << '\t' << distance << '\t' << r * std::cos(deltatheta) << '\t' << r * std::sin(deltatheta) << std::endl;

	}
	datafile.close();
	return 0;
}
