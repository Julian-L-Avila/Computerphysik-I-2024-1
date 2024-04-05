#include<iomanip>
#include<iostream>
#include<fstream>

const double desiredPrecision = 10;

int main() {
	std::string file1Name = "dattrapezoid1s.dat";
	std::string file2Name = "dattrapezoid1l.dat";

	std::ifstream file1(file1Name);
	std::ifstream file2(file2Name);
	std::ofstream file3("error.dat");

	double value1, value2, error;

	while (file1 >> value1 && file2 >> value2) {
		std::string dummy;
		file3 << std::setprecision(desiredPrecision) << value1 << '\t';
		file1 >> dummy;
		file2 >> dummy;

		error = std::abs(value1 - value2) / (value1);

		file3 << std::setprecision(desiredPrecision) << error << '\n';
		file1 >> dummy >> dummy;
		file2 >> dummy >> dummy;
	}

	file1.close();
	file2.close();
	file3.close();
}
