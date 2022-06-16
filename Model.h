#pragma once

#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkTetra.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkCell.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkDelaunay3D.h>
#include <vtkPolyData.h>

#include <string.h>
#include <iostream>
#include <set>
#include <map>

#include "Prime.h"

using namespace std;

const double eps = 1e-300;

class BoundBox {
public:
	double upBound[3], downBound[3];
	vector<BoundBox> getSons();
	BoundBox(double x0, double x1, double y0, double y1, double z0, double z1);
	BoundBox(){}
	bool include(Vector3d pos);
};

class OctCell {
public:
	int id;
	int sons[8];
	BoundBox box;
	Coff c;
	Vector4d bestPos;
	vector<int>vertexs;
	OctCell(int id, BoundBox box);
	OctCell() {};
};

class Model
{
public:
	double upBound[4], downBound[4];
	Vertex* vertexBuffer;
	int vCnt;
	Tetra* tetraBuffer;
	int tCnt;
	vector<OctCell>cellBuffer;
	map<tuple<int, int, int>, bool>borderMap;
	Model(string fileName);
	void buildOct(int threshold);
	void splitCell(int threshold, int cellId);
	int getCellId(Vector4d pos);
	void addQ();
	void selectBorder(double borderWeight);
	void getBestPos();
	void outputVtk(string fileName);
};
