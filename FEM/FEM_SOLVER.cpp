
#include "FEM_SOLVER.h"
//=======================================//
void SOLVER::Calculation() {

	//��������� �����
	for (int i = 0; i < this->MeshMain->Nodes.size(); i++)
		this->MeshMain->Nodes[i]->Num = i;

	//��������� ���������
	//���������� ������ ���������� �������������
	for (int i = 0; i < this->Equations.size(); i++) {
		this->Equations[i].Neq = i;
		this->Equations[i].ConstructMatrix();
	}

	//���� �� �������
	while (TimeSteps--) {

		int InterIter(MaxInterIter);
		
		//���� �� ���������� ���������
		while (InterIter--) {

			for (auto& equation : Equations)
				equation.UpdateMatrix();
			for (auto& equation : Equations)
				equation.InternalIteration();

		}//���������� ���������� ��������

		//for (auto& equation : Equations)
		//	equation.NewTimeStep();

	}
}
//=======================================//


