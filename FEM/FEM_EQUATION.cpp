#pragma once
#include "FEM_EQUATION.h"
#include "FEM_RESULTS.h"

//============================================================//
void BOUNDARY_CONDITION::SetArea(bool(*Set)(dVec2 x), FemMesh2D Mesh) {

	//добавление в граничное условие элементы границы

	for (auto& s : Mesh->Sides)
		if (Set(s->NodeBegin->X) && Set(s->NodeEnd->X))
			Sides.push_back(s);
}
//============================================================//
void EQUATION::Intitialization(double(*InitCondition)(dVec2 X)) {

	for (auto& node : Mesh->Nodes) {
		auto value(InitCondition(node->X));
		auto num(node->Nodes.size() + 1);
		node->Value.push_back(value);
		node->MatrixFree.push_back(0.0);
		node->Matrix.push_back(vector<double>(num));
		auto& v = node->Matrix.back();
		fill(v.begin(), v.end(), 0.0);
	}
}
//==============================================================//
void EQUATION::ConstructMatrix() {

	//построение матрицы (постоянных значений)
	int index1(0);
	int index2(0);
	double value(0.0);

	//цикл по всем узлам
	for (auto& P : Mesh->Nodes) {

		//построение индексации разряженной матрицы
		P->SetInternalIndex();

		//цикл по структуре уравнения
		for(auto& term : TERMS) {

			//нестационарное слагаемое
			if (term.Type == TIME)
				for (int i(0); i != P->Cells.size(); ++i)
					for (int j(0); j != P->Cells[i]->Nodes.size(); ++j) {
						index1 = P->Index[i][j];
						index2 = P->IndexSelf[i];
						auto& cell = P->Cells[i];
						P->Matrix[Neq][index1] += cell->ShapeMoment2[index2][j] / term.Constant;
						P->MatrixPre[Neq][index1] -= cell->ShapeMoment2[index2][j] / term.Constant;
					}

			//оператор лапласа
			if (term.Type == LAPLASS)
				for (int i(0); i != P->Cells.size(); ++i)
					for (int j(0); j != P->Cells.at(i)->Nodes.size(); ++j) {
						index1 = P->Index[i][j];
						index2 = P->IndexSelf[i];
						auto& cell = P->Cells[i];
						value = term.Constant * cell->ShapeMoment0;
						P->Matrix[Neq][index1] -= value * cell->Form[1][index2] * cell->Form[1][j];
						P->Matrix[Neq][index1] -= value * cell->Form[2][index2] * cell->Form[2][j];
					}

			//само значение
			if (term.Type == THIS)
				for (int i(0); i != P->Cells.size(); ++i)                           
					for (int j(0); j != P->Cells[i]->Nodes.size(); ++j) {
						index1 = P->Index[i][j];
						index2 = P->IndexSelf[i];
						auto& cell = P->Cells[i];
						P->Matrix[Neq][index1] += cell->ShapeMoment2[index2][j] * term.Constant;
					}

		}
	}
}
//========================================================================================================================//
void EQUATION::UpdateMatrix() {

	//обновление свободных слагаемых
	//(непоятонные значения)

	//цикл по всем узлам

	for (auto& P : Mesh->Nodes) {

		P->MatrixFree[Neq] = 0.0;

		//цикл по структуре уравнения
		for(auto& term : TERMS) {

			//иточник (слагаемое без производных)
			if (term.Type == SOURCE) {

				//без множителей (количество множителей Len = 0)
				if (term.Len == 0)
					for (int i = 0; i < P->Cells.size(); ++i) {
						auto& cell = P->Cells[i];
						P->MatrixFree[Neq] += term.ConstCell(cell) *
							cell->ShapeMoment1[P->IndexSelf[i]];
					}

				//количество множителей 1
				if (term.Len == 1)
					for (int i(0); i != P->Cells.size(); ++i)
						for (int j(0); j != P->Cells[i]->Nodes.size(); ++j) {
							auto& cell = P->Cells[i];
							P->MatrixFree[Neq] += term.ConstCell(cell) *
								cell->ShapeMoment2[P->IndexSelf[i]][j] *
								cell->Nodes[j]->Value[term.param1[0]];

						}				
			}

			//слагаемое конвективного типа (при переходе к эйлеровым координатам)
			//как в уравнениях навье стокса для скоростей
			if (term.Type == CONVECTIVE) {

				if (term.Len == 0)
					for (int i(0); i != P->Cells.size(); ++i) {

						double sum(0.0);
						auto& cell = P->Cells[i];

						for (int j(0); j != cell->Nodes.size(); ++j)
							sum += (cell->ShapeMoment2[P->IndexSelf[i]][j]) *
								(cell->Nodes[j]->Value[term.param1[0]]);
						
						P->MatrixFree[Neq] += sum * cell->GetGradient(term.param1[1], term.param1[2]) *
								term.ConstCell(cell);
					}
			}

			//дивиргенция (как в уравнении неразрывности)
			if (term.Type == DIVERGENCE)
				if (term.Len == 0)					
					for (int i(0); i != P->Cells.size(); ++i) {
						
						double sum1(0.0);
						double sum2(0.0);
						auto& cell = P->Cells[i];
					
						for (int j(0); j != cell->Nodes.size(); ++j) {

							sum1 += (cell->ShapeMoment1[j])
								*(cell->Nodes[j]->Value[term.param1[0]]);
							
							sum2 += (cell->ShapeMoment1[j])
								*(cell->Nodes[j]->Value[term.param1[1]]);
						}
					
						P->MatrixFree[Neq] -= sum1 * cell->Form[1][P->IndexSelf[i]] * term.ConstCell(cell);
						P->MatrixFree[Neq] -= sum2 * cell->Form[2][P->IndexSelf[i]] * term.ConstCell(cell);
					}

			//производная по координате некоторой переменной
			if (term.Type == DIRIVATE)
				if (term.Len == 1)
					for (int i(0); i != P->Cells.size(); ++i) {
						auto& cell = P->Cells[i];
						P->MatrixFree[Neq] += (cell->ShapeMoment1[P->IndexSelf[i]])
							* cell->GetGradient(term.param1[0], term.param1[1]) * term.ConstCell(cell);
					}
	
		}
	}
	
	//учитываем граничные условия
	//цикл по граничным условиям
	for (auto& BC : this->BoundaryConditions)
		for(auto& S : BC.Sides)
			for (auto& term : BC.TERMS)
				for (int i(0); i != 3; ++i) {

					auto& P = S->CellLink->Nodes[i];

					if ((P != S->NodeBegin) && (P != S->NodeEnd)) 
						continue;

					if (term.Len == 0)
						P->MatrixFree[Neq] += term.ConstSide(S) * S->ShapeMoment1[i];
					
					if (term.Len == 1)
						for (int j(0); j != 3; j++) 
							P->MatrixFree[Neq] += term.ConstSide(S) * S->ShapeMoment2[i][j]
								* S->CellLink->Nodes[j]->Value[term.param1[0]];				

				}
}
//========================================================================================================================//
void EQUATION::InternalIteration() {

	for(auto& P : Mesh->Nodes) {
		
		double sum = 0.0;
		
		for (int i(0); i != P->Nodes.size(); ++i) {
			
			sum += P->Nodes[i]->Value[Neq] * P->Matrix[Neq][i];
			
			if (!Static)
				sum += P->Nodes[i]->ValuePre[Neq] * P->MatrixPre[Neq][i];
		}

		sum += P->Matrix.at(Neq)[7];
		
		if (!Static)
			sum += P->ValuePre[Neq] * P->MatrixPre[Neq][6];

		P->Value[Neq] = (1.0 - Relax) * P->Value[Neq] 
			- Relax * sum / (P->Matrix[Neq][6] + P->Matrix[Neq][8]);
	}
}
//===========================================================//
void EQUATION::NewTimeStep() {
	
	if (Static) return;

	for (auto& P : Mesh->Nodes)
		P->ValuePre[Neq] = P->Value[Neq];
}
//===========================================================//