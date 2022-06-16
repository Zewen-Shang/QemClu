#include "Prime.h"

Vertex::Vertex(Vector4d pos):pos(pos)
{
}

Tetra::Tetra(int vertexs[4])
{
	for (int i = 0; i < 4; i++) {
		this->vertexs[i] = vertexs[i];
	}
}



Coff::Coff()
{
	A = Matrix4d::Zero();
	b = Vector4d::Zero();
	c = 0;
}

Coff::Coff(Matrix4d A, Vector4d b, double c) :A(A), b(b), c(c)
{
}

double Coff::getError(Vector4d x)
{
	return (x.transpose() * A * x)(0) + 2 * (b.transpose() * x)(0) + c;
}

vector<Vector4d> smtOrth(vector<Vector4d>A) {
	vector<Vector4d>B;
	for (int i = 0; i < A.size(); i++) {
		Vector4d t = A[i];
		for (int j = 0; j < i; j++) {
			t -= A[i].dot(B[j]) / B[j].dot(B[j]) * B[j];
		}
		t.normalized();
		B.push_back(t);
	}
	return B;
}