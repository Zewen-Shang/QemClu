#include "Model.h"

using namespace std;

Model::Model(string fileName)
{
	for (int i = 0; i < 4; i++) {
		upBound[i] = DBL_MIN;
		downBound[i] = DBL_MAX;
	}

	int temp;

	vtkNew<vtkUnstructuredGridReader> reader;
	vtkUnstructuredGrid* unGrid;
	string path = "./input/" + fileName + ".vtk";
	reader->SetFileName(path.c_str());
	reader->Update();
	unGrid = reader->GetOutput();

	vCnt = unGrid->GetPoints()->GetNumberOfPoints();
	vertexBuffer = new Vertex[vCnt + 10];

	string attrName = "attr";

	for (int i = 0; i < vCnt; i++) {
		double* pos = unGrid->GetPoint(i);
		vtkPointData* pointData = unGrid->GetPointData();
		double* attr;
		attr = pointData->GetArray(attrName.c_str())->GetTuple(i);
		for (int j = 0; j < 4; j++) {
			upBound[j] = max(upBound[j], pos[j]);
			downBound[j] = min(downBound[j], pos[j]);
		}
		Vector4d Pos(pos[0], pos[1], pos[2],attr[0]);
		vertexBuffer[i] = Vertex(Pos);
	}
	for (int i = 0; i < 4; i++) {
		upBound[i] += 0.005 * (upBound[i]-downBound[i]);
		downBound[i] -= 0.005 * (upBound[i] - downBound[i]);
	}

	tCnt = unGrid->GetCells()->GetNumberOfCells();
	tetraBuffer = new Tetra[tCnt + 10];

	vtkIdList* ids;
	for (int i = 0; i < tCnt; i++) {
		int id[4];
		ids = unGrid->GetCell(i)->GetPointIds();
		for (int j = 0; j < 4; j++) {
			id[j] = ids->GetId(j);
		}
		tetraBuffer[i] = Tetra(id);

		sort(id, id + 4);
		tuple<int, int, int> tups[4] = {
			tuple<int,int,int>(id[0],id[1],id[2]),
			tuple<int,int,int>(id[0],id[1],id[3]),
			tuple<int,int,int>(id[0],id[2],id[3]),
			tuple<int,int,int>(id[1],id[2],id[3]),
		};

		for (int i = 0; i < 4; i++) {
			if (borderMap.count(tups[i]) == 1) {
				borderMap.erase(tups[i]);
			}
			else {
				borderMap[tups[i]] = 1;
			}
		}
	}
}

void Model::buildOct(int threshold)
{
	BoundBox box(downBound[0], upBound[0], downBound[1], upBound[1], downBound[2], upBound[2]);
	OctCell cell(cellBuffer.size(),box);
	for (int i = 0; i < vCnt; i++) {
		cell.vertexs.push_back(i);
		vertexBuffer[i].cellId = 0;
	}
	cellBuffer.push_back(cell);
	splitCell(threshold,0);
}

void Model::splitCell(int threshold, int cellId)
{
	if (cellBuffer[cellId].vertexs.size() > threshold) {
		OctCell sons[8];
		vector<BoundBox>sonBoxs = cellBuffer[cellId].box.getSons();
		for (int i = 0; i < 8; i++) {
			sons[i].id = cellBuffer.size() + i;
			sons[i].box = sonBoxs[i];
			cellBuffer[cellId].sons[i] = sons[i].id;
		}
		for (int vId : cellBuffer[cellId].vertexs) {
			Vertex& v = vertexBuffer[vId];
			for (int i = 0; i < 8; i++) {
				if (sons[i].box.include(v.pos.block(0,0,3,1))) {
					v.cellId = sons[i].id;
					sons[i].vertexs.push_back(vId);
				}
			}
		}
		for (int i = 0; i < 8; i++) {
			cellBuffer.push_back(sons[i]);
		}

		cellBuffer[cellId].vertexs.clear();

		for (int i = 0; i < 8; i++) {
			splitCell(threshold, sons[i].id);
		}
	}
}



void Model::addQ()
{
	for (int i = 0; i < tCnt; i++) {
		Tetra t = tetraBuffer[i];
		int* ids = t.vertexs;
		vector<Vector4d> es;
		for (int j = 0; j < 3; j++) {
			es.push_back((vertexBuffer[ids[j]].pos - vertexBuffer[ids[3]].pos).normalized());
		}
		es = smtOrth(es);
		Matrix4d A = Matrix4d::Identity();
		for (int j = 0; j < 3; j++) {
			A -= es[j] * es[j].transpose();
		}


		Matrix3d M;
		for (int j = 0; j < 3; j++) {
			M.col(j) = es[j].block(0, 0, 3, 1);
		}
		double deter = abs(M.determinant());

		for (int j = 0; j < 4; j++) {
			Vector4d p = vertexBuffer[ids[j]].pos;
			Vector4d b = -A * vertexBuffer[ids[j]].pos;
			double c = p.transpose() * A * p;
			Coff C(A, b, c);

			cellBuffer[vertexBuffer[ids[j]].cellId].c = cellBuffer[vertexBuffer[ids[j]].cellId].c + C * deter;
		}
	}
}

void Model::selectBorder(double borderWeight)
{
	for (auto i : borderMap) {
		int vPos[3];
		vPos[0] = std::get<0>(i.first);
		vPos[1] = std::get<1>(i.first);
		vPos[2] = std::get<2>(i.first);


		Vector4d e0 = vertexBuffer[vPos[1]].pos - vertexBuffer[vPos[0]].pos, e1 = vertexBuffer[vPos[2]].pos - vertexBuffer[vPos[0]].pos;
		Vector3d a0 = e0.block(0, 0, 3, 1), a1 = e1.block(0, 0, 3, 1);
		double deter = a0.cross(a1).norm();
		e0.normalize();
		e1 = e1 - e0 * (e0.dot(e1));
		e1.normalize();
		Matrix4d A = Matrix4d::Identity();
		A -= (e0 * e0.transpose() + e1 * e1.transpose());
		A *= abs(deter) / 6.0 * borderWeight;

		for (int j = 0; j < 3; j++) {
			Vector4d b = -A * vertexBuffer[vPos[j]].pos;
			double c = vertexBuffer[vPos[j]].pos.transpose() * A * vertexBuffer[vPos[j]].pos;
			//Vector4d b = Vector4d::Zero();
			//double c = 0;
			Coff C(A, b, c);
			cellBuffer[vertexBuffer[vPos[j]].cellId].c = cellBuffer[vertexBuffer[vPos[j]].cellId].c + C;
		}
	}
}

void Model::getBestPos()
{
	int inverse = 0, uninverse = 0;
	for (int i = 0; i < cellBuffer.size(); i++) {
		OctCell& cell = cellBuffer[i];
		if (cell.vertexs.size() == 0) {
			cell.bestPos << 0, 0, 0, 0;
		}
		else {
			if (abs(cell.c.A.determinant()) >= eps) {
				inverse++;
				Matrix4d AI = cell.c.A.inverse();
				cell.bestPos = -AI * cell.c.b;
			}
			else {
				uninverse++;
				cell.bestPos = Vector4d::Zero();
				for (int j : cell.vertexs) {
					cell.bestPos += vertexBuffer[j].pos;
				}
				cell.bestPos /= cell.vertexs.size();
			}
		}

	}
	cout << inverse << "  " << uninverse << endl;
}

void Model::outputVtk(string fileName)
{
	vtkSmartPointer<vtkUnstructuredGrid> unGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
	vtkSmartPointer<vtkCellArray> cellArray = vtkSmartPointer<vtkCellArray>::New();
	vtkNew<vtkPointData> pointData;
	vtkNew<vtkFloatArray> attr;
	string attrName = "attr";

	attr->SetName(attrName.c_str());

	map<tuple<int, int, int, int>, bool>cellMap;
	for (auto cell : cellBuffer) {
		if (cell.vertexs.size() == 0) {
			points->InsertNextPoint(0, 0, 0);
			attr->InsertNextTuple1(0);
		}
		else {
			points->InsertNextPoint(cell.bestPos(0), cell.bestPos(1), cell.bestPos(2));
			attr->InsertNextTuple1(cell.bestPos(3));
		}

	}

	for (int i = 0; i < 3; i++) {
		unGrid->GetPointData()->AddArray(attr);
	}

	vtkNew<vtkPolyData>polyData;
	polyData->SetPoints(points);
	vtkNew<vtkDelaunay3D> del;
	del->SetInputData(polyData);
	del->SetAlpha(0.0105);
	del->Update();
	unGrid = del->GetOutput();

	//for (int i = 0; i < tCnt; i++) {
	//	Tetra* t = &tetraBuffer[i];
	//	set<int>vertexs;
	//	for (int j = 0; j < 4; j++) {
	//		vertexs.insert(vertexBuffer[t->vertexs[j]].cellId);
	//	}
	//	if (vertexs.size() < 4)continue;
	//	tuple<int, int, int, int>four = std::tie(vertexBuffer[t->vertexs[0]].cellId, vertexBuffer[t->vertexs[1]].cellId, vertexBuffer[t->vertexs[2]].cellId, vertexBuffer[t->vertexs[3]].cellId) ;
	//	if (cellMap[four]) {

	//	}
	//	else {
	//		cellMap[four] = 1;
	//		vtkNew<vtkTetra> tetra;
	//		int cnt = 0;
	//		for (auto vId : vertexs) {
	//			tetra->GetPointIds()->SetId(cnt++, vId);
	//		}
	//		cellArray->InsertNextCell(tetra);
	//	}
	//}

	//unGrid->SetPoints(points);
	//unGrid->SetCells(10, cellArray);

	vtkSmartPointer<vtkUnstructuredGridWriter> writer = vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	writer->SetFileName(("./outputVtk/" + fileName + ".vtk").c_str());
	writer->SetInputData(unGrid);
	writer->Write();

}

OctCell::OctCell(int id, BoundBox box):id(id),box(box)
{
}

vector<BoundBox> BoundBox::getSons()
{
	vector<BoundBox> ans;
	double midBound[3];
	for (int i = 0; i < 3; i++) {
		midBound[i] = (upBound[i] + downBound[i]) / 2;
	}
	BoundBox son(downBound[0], midBound[0], downBound[1], midBound[1], downBound[2], midBound[2]);
	ans.push_back(son);
	son = BoundBox(downBound[0], midBound[0], downBound[1], midBound[1], midBound[2], upBound[2]);
	ans.push_back(son);
	son = BoundBox(downBound[0], midBound[0], midBound[1], upBound[1], downBound[2], midBound[2]);
	ans.push_back(son);
	son = BoundBox(downBound[0], midBound[0], midBound[1], upBound[1], midBound[2], upBound[2]);
	ans.push_back(son);

	son = BoundBox(midBound[0], upBound[0], downBound[1], midBound[1], downBound[2], midBound[2]);
	ans.push_back(son);
	son = BoundBox(midBound[0], upBound[0], downBound[1], midBound[1], midBound[2], upBound[2]);
	ans.push_back(son);
	son = BoundBox(midBound[0], upBound[0], midBound[1], upBound[1], downBound[2], midBound[2]);
	ans.push_back(son);
	son = BoundBox(midBound[0], upBound[0], midBound[1], upBound[1], midBound[2], upBound[2]);
	ans.push_back(son);

	return ans;
}

BoundBox::BoundBox(double x0, double x1, double y0, double y1, double z0, double z1)
{
	downBound[0] = x0, upBound[0] = x1;
	downBound[1] = y0, upBound[1] = y1;
	downBound[2] = z0, upBound[2] = z1;
}

bool BoundBox::include(Vector3d pos)
{
	for (int i = 0; i < 3; i++) {
		if (downBound[i] <= pos(i) && upBound[i] >= pos(i)) {
			continue;
		}
		else {
			return false;
		}
	}
	return true;
}
