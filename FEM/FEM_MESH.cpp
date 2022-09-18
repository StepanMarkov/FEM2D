#include "FEM_MESH.h"
#include "FEM_RESULTS.h"
#include <omp.h>

using namespace MESH;

void MESH2D::AddTriangle(dVec2 Xhead) {

	//создание треугольного элемента

	FemSide2D Side1(new SIDE2D());
	FemSide2D Side2(new SIDE2D());
	FemSide2D Side3(new SIDE2D());
	FemNode2D Node1(new NODE2D());
	FemNode2D Node2(new NODE2D());
	FemNode2D Node3(new NODE2D());
	FemCell2D Cell(new CELL2D());

	Node1->X = Xhead;
	Node2->X = Xhead;
	Node3->X = Xhead;
	Node2->X[0] += this->MeshDist(Xhead);
	Node3->X[0] += this->MeshDist(Xhead) / 2.0;
	Node3->X[1] += this->MeshDist(Xhead) * sin(PI / 3.0);

	Cell->Nodes.push_back(Node1);
	Cell->Nodes.push_back(Node2);
	Cell->Nodes.push_back(Node3);
	Node1->Cells.push_back(Cell);
	Node2->Cells.push_back(Cell);
	Node3->Cells.push_back(Cell);
	this->Sides.push_back(Side1);
	this->Sides.push_back(Side2);
	this->Sides.push_back(Side3);
	this->Cells.push_back(Cell);
	this->Nodes.push_back(Node1);
	this->Nodes.push_back(Node2);
	this->Nodes.push_back(Node3);

	Side1->CreateSide(Node1, Node2, Cell);
	Side2->CreateSide(Node2, Node3, Cell);
	Side3->CreateSide(Node3, Node1, Cell);

	Node1->Nodes.push_back(Node2);
	Node2->Nodes.push_back(Node1);
	Node2->Nodes.push_back(Node3);
	Node3->Nodes.push_back(Node2);
	Node3->Nodes.push_back(Node1);
	Node1->Nodes.push_back(Node3);

	Side1->SideL = Side3;
	Side1->SideR = Side2;
	Side2->SideL = Side1;
	Side2->SideR = Side3;
	Side3->SideL = Side2;
	Side3->SideR = Side1;
}

//========================================================================================================================//
void MESH2D::Optimization(int Ndir, double eps, int Iter) {
	omp_set_num_threads(16);

	dVec2 force			= dVec2_zero;            //сила, действующая на узел
	dVec2 force_bound	= dVec2_zero;            //сила, действующая на граничный узел
	dVec2 bound			= dVec2_zero;
	dVec2 bound0		= dVec2_zero;
	dVec2 bound1		= dVec2_zero;
	bool sw0(false);
	bool sw1(false);
	double delt(0.0);           //отклонение расстояния между узлами от размера сетки size
	vector <dVec2> Beams;       //направления прилипания к границе
	double MinDistance(0.0);    //мимимальное расстояние до границы
	double CurDictance(0.0);    //текущее расчитанное расстояние до границы

	//формируем напрвления смещений (минимальный набор - 4)
	Beams.resize(Ndir);
	
	for (int i(0); i != Ndir; ++i) {
		Beams[i][0] = cos((double)i * 2.0 * PI / (double)Ndir);
		Beams[i][1] = sin((double)i * 2.0 * PI / (double)Ndir);
	}

	while(Iter--) {

		//прилипание к границе

		for(auto& n : Nodes)
			if (n->SwBound) {

				//определяем расстояние до границы методом дихотомии

				bool sw(false);
				MinDistance = INFINITY;
				dVec2 Aver(dVec2_zero);
				double Sum(0.0);

				//цикл по направлениям
				for (int i(0); i != Ndir; ++i) {
					bound1 = n->X;
					bound0 = n->X + Beams.at(i) * 2.0 * MeshDist(n->X);
					sw0 = Geometry(bound0);
					sw1 = Geometry(bound1);
					if (!(sw0 ^ sw1))  continue;
					if (sw0) { bound1 = bound0; bound0 = n->X; }

					do
					{
						bound = bound1 + bound0;
						bound /= 2.0;
						if (Geometry(bound)) bound1 = bound;
						else bound0 = bound;
					} 
					while (distance(bound0, bound1) > eps * MeshDist(n->X));

					bound1 = (bound1 + bound0) / 2.0;					
					Aver += bound1 / (EPSILON + distance(bound1, n->X));
					Sum += 1.0 / (EPSILON + distance(bound1, n->X));									
				}
		
				Aver /= Sum;
				n->Dir = Aver - n->X;
				n->X = n->X * (1 - eps) + Aver * eps;
			}

		//прилипание к границам на данной итерации завершено
		//приступаем к выравниванию сетки

#pragma omp parallel firstprivate(force, force_bound, delt)
		{

			#pragma omp for
			for (int i = 0; i < Nodes.size(); ++i) {

				FemNode2D n = Nodes[i];
				force = dVec2_zero;

				//цикл по соседним узлам
				for (auto& nn : n->Nodes) {
						delt = distance(n->X, nn->X) - (MeshDist(n->X));
						force -= normalize(n->X - nn->X) * delt;
					}

				//цикл по соседним ячейкам
				for (auto& cc : n->Cells) {
					cc->Xc = cc->GetBarycenter();
					delt = distance(n->X, cc->Xc) - (MeshDist(n->X) / sqrt(3.0));
					force -= normalize(n->X - cc->Xc) * delt;
				}

				//если узел не граничный
				if (!n->SwBound)
				{
					//смещаем координаты узла
					n->X += force * eps;
				}
				else {
					//если узел граничный
					//на узел дествует только проекция силы на границу
					force_bound = normalize(normal(n->Dir));
					force_bound *= dot(force, force_bound);
					n->X += force_bound * eps;
				}
			}
		}

	} //итерации	
}
//=============================================//
void MESH2D::Meshing(dVec2 Xhead) {

	omp_set_num_threads(16);
	AddTriangle(Xhead);
	vector <FemSide2D> Bound(Sides);
	Sides.clear();
	double CrossR, DotR;
	int index(0);

	while (true) {
		
		bool Sw(true);

		for (int i(index); i != Bound.size(); ++i)
			if (!Bound[i])
				index = i;
			else break;

		//цикл по внешней границе
		for (int i(index); i != Bound.size(); ++i)
			if (Bound[i]) {

				auto& s(Bound[i]);

				if (!s->SideR) {
					s = nullptr;
					continue;
				}
				
				//считаем расположение элементов границы
				CrossR = cross(s->Norm, s->SideR->Norm);
				DotR   = dot  (s->Norm, s->SideR->Norm);

				//если угол между элементами границами мал
				if ((CrossR <= 0.0) && (DotR <= 0.4)) {
				
					FemSide2D NewSide(new SIDE2D());
					FemCell2D NewCell(new CELL2D());
					NewCell->Nodes.push_back(s->NodeBegin);
					NewCell->Nodes.push_back(s->NodeEnd);
					NewCell->Nodes.push_back(s->SideR->NodeEnd);
					s->NodeBegin->Cells.push_back(NewCell);
					s->NodeEnd->Cells.push_back(NewCell);
					s->SideR->NodeEnd->Cells.push_back(NewCell);
					NewSide->CreateSide(s->NodeBegin, s->SideR->NodeEnd, NewCell);
					NewSide->SideR = s->SideR->SideR;
					NewSide->SideL = s->SideL;
					s->NodeBegin->Nodes.push_back(s->SideR->NodeEnd);
					s->SideR->NodeEnd->Nodes.push_back(s->NodeBegin);
					s->SideR->SideR->SideL = NewSide;
					s->SideL->SideR = NewSide;
					delete s->SideR;
					s->SideR->SideL = nullptr;
					s->SideR->SideR = nullptr;
					delete s;
					s = nullptr;		
					Cells.push_back(NewCell);
					Bound.push_back(NewSide);
					Sw = false;
					break;
				}
			}

		if (Sw)
			for (int i(index); i != Bound.size(); ++i)
				if (Bound[i]) {

					auto& s(Bound[i]);
					CrossR = cross(s->Norm, s->SideR->Norm);
					DotR   = dot  (s->Norm, s->SideR->Norm);

					if ((CrossR <= 0.0) && (DotR >= 0.4) || (CrossR >= 0.0)) {

						auto nx = s->Xc + s->Norm * MeshDist(s->Xc) * sin(PI / 3.0);
						if (!Geometry(nx)) continue;
						FemSide2D NewSide1(new SIDE2D());
						FemSide2D NewSide2(new SIDE2D());
						FemCell2D NewCell(new CELL2D());
						FemNode2D NewNode(new NODE2D());
						NewNode->X = nx;
						NewCell->Nodes.push_back(s->NodeBegin);
						NewCell->Nodes.push_back(s->NodeEnd);
						NewCell->Nodes.push_back(NewNode);
						s->NodeBegin->Cells.push_back(NewCell);
						s->NodeEnd->Cells.push_back(NewCell);
						NewNode->Cells.push_back(NewCell);
						NewSide1->CreateSide(s->NodeBegin, NewNode, NewCell);
						NewSide2->CreateSide(NewNode, s->NodeEnd, NewCell);
						s->NodeBegin->Nodes.push_back(NewNode);
						NewNode->Nodes.push_back(s->NodeBegin);
						s->NodeEnd->Nodes.push_back(NewNode);
						NewNode->Nodes.push_back(s->NodeEnd);
						NewSide1->SideR = NewSide2;
						NewSide1->SideL = s->SideL;
						NewSide2->SideR = s->SideR;
						NewSide2->SideL = NewSide1;
						s->SideL->SideR = NewSide1;
						s->SideR->SideL = NewSide2;
						delete s;
						s = nullptr;
						Bound.push_back(NewSide1);
						Bound.push_back(NewSide2);
						this->Nodes.push_back(NewNode);
						this->Cells.push_back(NewCell);
						Sw = false;
						break;
					}
				}

		if (Sw) break;
	}

	//определяем какие узлы лежат на границе

	for (auto& node : Nodes)
		node->SwBound = false;
		
	for (auto& s : Bound)
		if (s) {
			Sides.push_back(s);
			s->NodeBegin->SwBound = true;
			s->NodeEnd->SwBound = true;
		}
		
	for (int i(0); i != Nodes.size(); ++i)
		Nodes[i]->Num = i;

	for (int i(0); i != Cells.size(); ++i)
		this->Cells[i]->Num = i;

	//удаляем слабосвязные элементы 
	{
		vector <FemNode2D> VNODES;
		vector <FemSide2D> VSIDES;
		vector <FemCell2D> VCELLS;
		vector <FemNode2D> VNODES_DEL;
		vector <FemSide2D> VSIDES_DEL;
		vector <FemCell2D> VCELLS_DEL;

		//отмечаем удаляемые узлы
		for (auto& P : Nodes) {
			//если узел мало связан с сеткой
			if (P->Nodes.size() == 2) {
				//помечаем что его нужно удалить
				P->Num = -1;
				VNODES_DEL.push_back(P);

			}
			else VNODES.push_back(P);
			//иначе сохраняем указатели в вспогательный вектор нормальные узлы
		}

		//отмечаем удаляемые конечные элементы
		for (auto& C : Cells) {
			bool sw = true;
			for (auto& P : C->Nodes)
				if (P->Num == -1) {
					VCELLS_DEL.push_back(C);
					C->Num = -1;
					sw = false;
					break;
				}
			if (sw) VCELLS.push_back(C);
		}

		for (auto& P : VNODES_DEL) {
			FemSide2D NewSide(new SIDE2D());
			for (auto& c1 : P->Nodes[0]->Cells)
				for (auto& c2 : P->Nodes[1]->Cells)
					if ((c1 == c2) && (c1->Num != -1))
						NewSide->CreateSide(P->Nodes[0], P->Nodes[1], c1);
			VSIDES.push_back(NewSide);
		}

		//отмечаем удаляемые элементы границы расчетной области
		for (auto& S : Sides) {
			if (S->CellLink->Num == -1) {
				VSIDES_DEL.push_back(S);
				S->Len = -1.0;
			}
			else VSIDES.push_back(S);
		}

		this->Nodes = VNODES;
		this->Sides = VSIDES;
		this->Cells = VCELLS;

		for (int i(0); i != this->Nodes.size(); ++i)
			this->Nodes[i]->Num = i;

		for (int i(0); i != this->Cells.size(); ++i)
			this->Cells[i]->Num = i;

		//цикл по удаленным узлам
		for (auto& P : VNODES_DEL) {

			for (auto& PP : P->Nodes) {

				VNODES.clear();
				VCELLS.clear();

				//убираем ненужные указатели
				for (auto& PPP : PP->Nodes)
					if (PPP->Num != -1)
						VNODES.push_back(PPP);

				PP->Nodes.clear();
				PP->Nodes = VNODES;

				for (auto& C : PP->Cells)
					if (C->Num != -1)
						VCELLS.push_back(C);
				PP->Cells = VCELLS;
			}
		}
	}
		
	Optimization(24, 1E-2, 1000);
																			                                                 
	#pragma omp parallel
	{
		#pragma omp for
		for (int i(0); i<this->Cells.size(); ++i)
			Cells[i]->FormCalculate(this->Axis);
		#pragma omp for
		for (int i(0); i < this->Sides.size(); ++i)
			this->Sides[i]->FormCalculate(this->Axis);
	}

}
