#include<iostream>
#include<cstdlib>
#include<cmath>
#include<functional>

struct Vec3 {
	double x, y, z;
};

Vec3 componentWise(Vec3 a, Vec3 b, auto op) {
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

Vec3 exteriorProduct(Vec3 a, Vec3 b) {
	return {a.x * b.y - a.y * b.x, a.x * b.z - a.z * b.x, a.y * b.z - a.z * b.y};
}

Vec3 crossProduct(Vec3 a, Vec3 b) {
	Vec3 avb = hadamardProduct(exteriorProduct(a, b), {1, -1, 1});
	return {avb.z, avb.y, avb.x};
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
	return ("(" + std::to_string(a.x) + " e12, " + std::to_string(a.y) + " e13, " + std::to_string(a.z) + " e23)");
}

std::string geometricProduct(Vec3 a, Vec3 b) {
	return (std::to_string(innerProduct(a, b)) + " + " + printBiVec3(exteriorProduct(a, b)));
}

int main() {
	Vec3 v1 = {1, 2, 3};
	Vec3 v2 = {2, 4, 3};

	std::cout << norm(v1) << std::endl;
	std::cout << printVec3(addition(v1, v2)) << std::endl;
	std::cout << printVec3(substraction(v1, v2)) << std::endl;
	std::cout << printVec3(scalarProduct(v1, 2)) << std::endl;
	std::cout << printVec3(hadamardProduct(v1, v2)) << std::endl;
	std::cout << innerProduct(v1, v2) << std::endl;
	std::cout << printBiVec3(exteriorProduct(v1, v2)) << std::endl;
	std::cout << printVec3(crossProduct(v1, v2)) << std::endl;
	std::cout << geometricProduct(v1, v2) << std::endl;

	return 0;
}
