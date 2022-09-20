#include "FEM_BASE.h"
#include "FEM_RESULTS.h"

using namespace FEM;

//===========================================================//
void NODE2D::SetInternalIndex() {

	//индексация разряженной матрицы
	Index.resize(Matrix.size());

	for (int i(0); i != Cells.size(); ++i)
		for (int j(0); j != 3; ++j) {
			//[6] указывает на this (индексация указывает на этот узел)
			
			Index[i][j] = 6;
			for (int k(0); k != Nodes.size(); ++k)
				if (Nodes[k] == Cells[i]->Nodes[j])
					Index[i][j] = k;
			if (Index[i][j] == 6) 
				IndexSelf[i] = j;
		}
}	
//======================================================================================================================//
void SIDE2D::CreateSide(FemNode2D n1, FemNode2D n2, FemCell2D Cell) {

	//создает элемент границы

	CellLink = Cell;
	CellLink->Xc = Cell->GetBarycenter();
	Xc = (n2->X + n1->X) / 2.0;

	//определяем порядок узлов (контур внешней границы против часовой)
	if (cross(Xc - Cell->Xc, n2->X - n1->X) > 0) {                                               
		NodeBegin = n1;
		NodeEnd = n2;
	}
	else {
		NodeBegin = n2;
		NodeEnd = n1;
	}

	//определение пормали границы
	Norm = NodeEnd->X - NodeBegin->X;
	Norm = normalize(Norm);
	Norm = normal(Norm);
	if (dot(Norm, Xc - Cell->Xc) < 0) Norm = -Norm;

	//определение длины
	Len = distance(NodeBegin->X, NodeEnd->X);
}
//======================================================================================================================//
void SIDE2D::FormCalculate(bool Axis) {

	//расчитывает значения, зависящие от формы

	CellLink->Xc = CellLink->GetBarycenter();
	Xc = (NodeEnd->X + NodeBegin->X) / 2.0;
	Norm = NodeEnd->X - NodeBegin->X;
	Norm = normalize(Norm);
	Norm = normal(Norm);
	Len = distance(NodeBegin->X, NodeEnd->X);

	if (!Axis) {
	
		for (int i = 0; i < 3; i++)
			ShapeMoment1[i] = IntegrateSide([&](dVec2 x) {
			return (CellLink->Form[0][i] + CellLink->Form[1][i] * x[0] + CellLink->Form[2][i] * x[1]);
				}, 50);
	
		for (int i(0); i != 3; ++i)
			for (int j(0); j != 3; ++j)
				ShapeMoment2[i][j] = IntegrateSide([&](dVec2 x) {
				return (CellLink->Form[0][i] + CellLink->Form[1][i] * x[0] + CellLink->Form[2][i] * x[1]) *
					   (CellLink->Form[0][j] + CellLink->Form[1][j] * x[0] + CellLink->Form[2][j] * x[1]);
					   }, 50);
	}	
	else {
	
		for (int i(0); i != 3; ++i)
			ShapeMoment1[i] = IntegrateSide([&](dVec2 x) {
			return (CellLink->Form[0][i] + CellLink->Form[1][i] * x[0] + CellLink->Form[2][i] * x[1]) * x[1];
				}, 50);
	
		for (int i(0); i != 3; ++i)
			for (int j(0); j != 3; ++j)
				ShapeMoment2[i][j] = IntegrateSide([&](dVec2 x) {
				return (CellLink->Form[0][i] + CellLink->Form[1][i] * x[0] + CellLink->Form[2][i] * x[1]) *
					   (CellLink->Form[0][j] + CellLink->Form[1][j] * x[0] + CellLink->Form[2][j] * x[1]) * x[1];
					   }, 50);
	}
}
//======================================================================================================================//
double SIDE2D::IntegrateSide(const std::function<double(dVec2)>& Func, int Num) {

	//интегрирование по элементу границы (численно)
	//Num - число разбиений

	double Sum(0.0);
	//определение шага
	dVec2 dx = (NodeEnd->X - NodeBegin->X) / (double)(Num);
	dVec2 x = this->NodeBegin->X;

	for (int i(0); i != Num; ++i) {
		Sum += Func(x);
		x += dx;
	}

	return (Sum * this->Len / (double)(Num));
}
//=================================================================//
double CELL2D::GetGradient(int i, int k) {

	double Sum(0.0);
	for (int j(0); j != 3; ++j)
		Sum += (this->Nodes[j]->Value[i]) * Form[k][j];
	return Sum;
}
//=====================================================//
double CELL2D::IntegrateElem(const std::function<double(dVec2)>& Func, int Num) {
	
	double Sum(0.0);
	int n(0);
	dVec2 x(this->Xmin);
	dVec2 dx((Xmax - Xmin) / (double)(Num));

	while (x[0] < Xmax[0]) {		
		x[1] = Xmin[1];
		while (x[1] < Xmax[1]) {			
			if (this->PointInside(x)) {
				Sum += Func(x);
				n++;
			}
			x[1] += dx[1];
		}
		x[0] += dx[0];
	}
	return (abs(this->Area) / 2.0 * Sum / (double)(n));
}
//======================================================================================================================//

bool CELL2D::PointInside(dVec2 px) {
	double A[3];
	A[0] = (Nodes[0]->X[0] - px[0]) * (Nodes[1]->X[1] - Nodes[0]->X[1]);
	A[0] -= (Nodes[1]->X[0] - Nodes[0]->X[0]) * (Nodes[0]->X[1] - px[1]);
	A[1] = (Nodes[1]->X[0] - px[0]) * (Nodes[2]->X[1] - Nodes[1]->X[1]);
	A[1] -= (Nodes[2]->X[0] - Nodes[1]->X[0]) * (Nodes[1]->X[1] - px[1]);
	A[2] = (Nodes[2]->X[0] - px[0]) * (Nodes[0]->X[1] - Nodes[2]->X[1]);
	A[2] -= (Nodes[0]->X[0] - Nodes[2]->X[0]) * (Nodes[2]->X[1] - px[1]);

	return (((A[0] >= 0.0) && (A[1] >= 0.0) 
		&& (A[2] >= 0.0)) || ((A[0] <= 0.0) 
			&& (A[1] <= 0.0) && (A[2] <= 0.0)));
}
//======================================================================================================================//
dVec2 CELL2D::GetBarycenter() {
	
	auto x = dVec2_zero;
	for (auto& node : Nodes)
		x += node->X;
	x /= 3.0;
	return x;
}
//======================================================================================================================//
void CELL2D::FormCalculate(bool Axis) {

	//расчет значений зависящих от формы

	Form[0][0] = Nodes[1]->X[0] * Nodes[2]->X[1] - Nodes[1]->X[1] * Nodes[2]->X[0];
	Form[0][1] = Nodes[0]->X[1] * Nodes[2]->X[0] - Nodes[0]->X[0] * Nodes[2]->X[1];
	Form[0][2] = Nodes[0]->X[0] * Nodes[1]->X[1] - Nodes[1]->X[0] * Nodes[0]->X[1];
	Form[1][0] = Nodes[1]->X[1] - Nodes[2]->X[1];
	Form[1][1] = Nodes[2]->X[1] - Nodes[0]->X[1];
	Form[1][2] = Nodes[0]->X[1] - Nodes[1]->X[1];
	Form[2][0] = Nodes[2]->X[0] - Nodes[1]->X[0];
	Form[2][1] = Nodes[0]->X[0] - Nodes[2]->X[0];
	Form[2][2] = Nodes[1]->X[0] - Nodes[0]->X[0];

	//опредение площади элемента (не приведенная)
	Area = det3(1.0, Nodes[0]->X[0], Nodes[0]->X[1],
                1.0, Nodes[1]->X[0], Nodes[1]->X[1],
		        1.0, Nodes[2]->X[0], Nodes[2]->X[1]);


	for (int i(0); i != 3; ++i)
		for (int j(0); j != 3; ++j)
			Form[i][j] /= Area;

	Xc = GetBarycenter();
	Xmin = Nodes[0]->X;
	Xmax = Nodes[0]->X;

	for (auto& n : Nodes) {
		if (n->X[0] < Xmin[0]) Xmin[0] = n->X[0];
		if (n->X[0] > Xmax[0]) Xmax[0] = n->X[0];
		if (n->X[1] < Xmin[1]) Xmin[1] = n->X[1];
		if (n->X[1] > Xmax[1]) Xmax[1] = n->X[1];
	}

	if (!Axis) {

		ShapeMoment0 = abs(Area) / 2.0;

		for (int i(0); i != 3; ++i)
			ShapeMoment1[i] = IntegrateElem([&](dVec2 x) {
			return (Form[0][i] + Form[1][i] * x[0] + Form[2][i] * x[1]); 
				}, 50);

		for (int i(0); i != 3; ++i)
			for (int j(0); j != 3; ++j)
				ShapeMoment2[i][j] = IntegrateElem([&](dVec2 x) {
					return (Form[0][i] + Form[1][i] * x[0] + Form[2][i] * x[1]) *
						   (Form[0][j] + Form[1][j] * x[0] + Form[2][j] * x[1]);
						   }, 50);
	}
	else {

		ShapeMoment0 = IntegrateElem([&](dVec2 x) {return (x[1]); }, 50);

		for (int i(0); i != 3; ++i)
			ShapeMoment1[i] = IntegrateElem([&](dVec2 x) {
			return (Form[0][i] + Form[1][i] * x[0] + Form[2][i] * x[1])*x[1]; 
				    }, 50);

		for (int i(0); i != 3; ++i)
			for (int j(0); j != 3; ++j)
				ShapeMoment2[i][j] = IntegrateElem([&](dVec2 x) {
				return (Form[0][i] + Form[1][i] * x[0] + Form[2][i] * x[1]) *
					   (Form[0][j] + Form[1][j] * x[0] + Form[2][j] * x[1]) *x[1];
			           }, 50);
	}
}