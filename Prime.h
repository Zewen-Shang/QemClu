#pragma once

#include <Eigen/dense>
#include <vector>

using namespace Eigen;
using namespace std;

typedef Matrix<double, 4, 1> Vector4d;
typedef Matrix<double, 4, 4> Matrix4d;

class Vertex {
public:
	Vector4d pos;
	int cellId;
	Vertex(Vector4d pos);
	Vertex() {};
};

class Tetra {
public:
	int vertexs[4];
	Tetra(int vertexs[4]);
	Tetra() {};
}; 

class Coff
{
public:
	Coff();
	Coff(Matrix4d A, Vector4d b, double c);
	Matrix4d A;
	Vector4d b;
	double c;
	Coff operator+(const Coff& other) {
		return Coff(this->A + other.A, this->b + other.b, this->c + other.c);
	}

	Coff operator*(const double x) {
		return Coff(this->A * x, this->b * x, this->c * x);
	}
	double getError(Vector4d x);
};


vector<Vector4d> smtOrth(vector<Vector4d>A);