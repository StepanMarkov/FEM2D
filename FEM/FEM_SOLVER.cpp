
#include "FEM_SOLVER.h"
//=======================================//
void SOLVER::Calculation() {

	//нумерация узлов
	for (int i = 0; i < this->MeshMain->Nodes.size(); i++)
		this->MeshMain->Nodes[i]->Num = i;

	//нумерация уравнений
	//построение матицы постоянных коэффициентов
	for (int i = 0; i < this->Equations.size(); i++) {
		this->Equations[i].Neq = i;
		this->Equations[i].ConstructMatrix();
	}

	//цикл по времени
	while (TimeSteps--) {

		int InterIter(MaxInterIter);
		
		//цикл по внутренним итерациям
		while (InterIter--) {

			for (auto& equation : Equations)
				equation.UpdateMatrix();
			for (auto& equation : Equations)
				equation.InternalIteration();

		}//завершение внутренних итераций

		//for (auto& equation : Equations)
		//	equation.NewTimeStep();

	}
}
//=======================================//


