#include <iostream>
#include <cstdlib>

int option, num1, num2;

bool isEven(int number) {
    return number % 2 == 0;
}

int add(int x, int y) {
    return x + y;
}

int subtract(int x, int y) {
    return x - y;
}

int multiply(int x, int y) {
    return x * y;
}

double divide(int x, int y) {
    if (y == 0) {
	    throw std::runtime_error("Division by zero");
    }
    return static_cast<double>(x) / y;
}

int modulo(int x, int y) {
    if (y == 0) {
	    throw std::runtime_error("Modulo by 0");
    }
    return x % y;
}

int main() {
	system("clear");
	std::cout << "Welcome to a basic Calculator in C++" << '\n';

	std::cout << "Enter a whole number x: " << std::endl;
	std::cin >> num1;
	if (isEven(num1)) {
		std::cout << num1 << " is even." << std::endl;
	} else {
		std::cout << num1 << " is odd." << std::endl;
	}

	std::cout << "\nEnter a whole number y: " << std::endl;
	std::cin >> num2;

	system("clear");
	std::cout << "\n x = " << num1 << " and y = " << num2 << '\n';
	std::cout << "::::::::::::::::::::MENU::::::::::::::::::::" << '\n';
	std::cout << "\t 1. Addition" << '\n';
	std::cout << "\t 2. Subtraction" << '\n';
	std::cout << "\t 3. Multiplication" << '\n';
	std::cout << "\t 4. Division" << '\n';
	std::cout << "\t 5. Modulo" << '\n';
	std::cout << "\t 6. Exit" << '\n';
	
	std::cout << "Choose an option: " << std::endl;
	std::cin >> option;

	switch (option) {
		case 1:
			std::cout << num1 << " + " << num2 << " = " << add(num1, num2) << std::endl;
			break;
		case 2:
			std::cout << num1 << " - (" << num2 << ") = " << subtract(num1, num2) << std::endl;
			break;
		case 3:
			std::cout << num1 << "(" << num2 << ") = " << multiply(num1, num2) << std::endl;
			break;
		case 4:
			std::cout << num1 << "(" << num2 << ")^-1 = " << divide(num1, num2) << std::endl;
			break;
		case 5:
			std::cout << num1 << " % " << num2 << " = " << modulo(num1, num2) << std::endl;
			break;
		default:
			std::cout << "Bye Bye" << '\n';
			std::cout << "Made by Julian L. Avila" << std::endl;
			exit(2);
    }

    return 0;
}
