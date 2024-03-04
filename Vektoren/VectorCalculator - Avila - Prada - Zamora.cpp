// Made by Julian Avila - Adriano Patada - Jose Zamora

#include<iostream>
#include<cstdlib>
#include<cmath>
#include<functional>

struct Vec3 {
	double x, y, z;
};

Vec3 inputVec3() {
    Vec3 v;
    std::cout << "Introduzca las componentes del vector (x, y, z): \n";
    std::cin >> v.x >> v.y >> v.z;
    return v;
}
Vec3 componentWise(Vec3 a, Vec3 b, std::function<double(double, double)> op) {
	return {op(a.x, b.x) , op(a.y, b.y), op(a.z, b.z)};
}

Vec3 hadamardProduct(Vec3 a, Vec3 b) {
	return componentWise(a, b, [](double x, double y) { return x * y;});
}

Vec3 scalarProduct(Vec3 a, double s) {
	return hadamardProduct(a, {s, s, s});
}

Vec3 addition(Vec3 a, Vec3 b) {
	return componentWise(a, b, std::plus<double>());
}

Vec3 substraction(Vec3 a, Vec3 b) {
	return componentWise(a, b, std::minus<double>());
}

Vec3 crossProduct(Vec3 a, Vec3 b) {
	return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

double innerProduct(Vec3 a, Vec3 b) {
	Vec3 axb = hadamardProduct(a, b);
	return (axb.x + axb.y + axb.z);
}

double norm(Vec3 a) {
	return std::sqrt(innerProduct(a, a));
}

std::string printVec3(Vec3 a) {
	return ("(" + std::to_string(a.x) + ", " + std::to_string(a.y) + ", " + std::to_string(a.z) + ")");
}

std::string printBiVec3(Vec3 a) {
	return ("(" + std::to_string(a.z) + " e12, " + std::to_string(-1 * a.y) + " e13, " + std::to_string(a.x) + " e23)");
}

std::string geometricProduct(Vec3 a, Vec3 b) {
	return (std::to_string(innerProduct(a, b)) + " + " + printBiVec3(crossProduct(a, b)));
}

int menu() {
	int option;
	std::cout << "::::::::::::::MENU::::::::::::::" << '\n';
	std::cout << "Seleccione una opcion\n\t1. Norma \n\t2. Suma de dos vectores \n\t3. Resta de dos vectores \n\t4. Producto por un escalar \n\t5. Producto Interior \n\t6. Producto Extererior \n\t7. Producto Cruz \n\t8. Producto Geometrico \n\t9. Producto de Hadamard \n\t10. Salir / Creditos" << std::endl;
	std::cin >> option;
	return option;
}

int main() {
	Vec3 v1, v2;
	int option;
	system("clear");
	std::cout << "Bienvenido a la calculadora de operaciones de vectores" << '\n';
	option = menu();
	if (option == 0 || option > 9) {
		std::cout << "Gracias y Adios - Made by J.L. Avila , M. A. Parada, J. L. Zamora" << std::endl;
		return 1;
	} else if (option == 1 || option == 4) {
		v1 = inputVec3();
		system("clear");
		std::cout << "v1 = " << printVec3(v1) << '\n' << std::endl;
	} else {
		v1 = inputVec3();
		v2 = inputVec3();
		system("clear");
		std::cout << "v1 = " << printVec3(v1) << " y v2 = " << printVec3(v2) << '\n'  << std::endl;
	}
	switch (option){
		case 1:
			std::cout << "La norma de " << printVec3(v1) << " es: " << norm(v1) << std::endl;
			break;
		case 2:
			std::cout << "La suma de los vectores es: " << printVec3(addition(v1, v2)) << std::endl;
			break;
		case 3:
			std::cout << "La resta entre los dos vectores es: " << printVec3(substraction(v1, v2)) << std::endl;
			break;
		case 4:
			double i;
			std::cout << "Ingrese el escalar" << std::endl;
			std::cin >> i;
			system("clear");
			std::cout << "El producto escalar de " << i << " el vector es: " << printVec3(scalarProduct(v1, i)) <<std::endl;
			break;
		case 5:
			std::cout << "El producto interior entre los vectores es: " << innerProduct(v1, v2) << std::endl;
			break;
		case 6:
			std::cout << "El producto exterior entre los vectores es: " << printBiVec3(crossProduct(v1, v2)) << std::endl;
			break;
		case 7:
			std::cout << "El producto cruz entre los vectores es:" << printVec3(crossProduct(v1, v2)) << std::endl;
			break;
		case 8:
			std::cout << "El producto geometrico entre los vectores es: " << geometricProduct(v1, v2) << std::endl;
			break;
		case 9:
			std::cout << "El producto de Hadamard entre los vectores es: " << printVec3(hadamardProduct(v1, v2)) << std::endl;
			break;
	}	
	return 0;
}
