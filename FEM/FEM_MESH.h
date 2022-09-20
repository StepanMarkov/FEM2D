#pragma once

#include "FEM_BASE.h"

using namespace std;
using namespace MATH;
using namespace FEM;

namespace MESH
{
	struct MESH2D;
	typedef shared_ptr <MESH2D> FemMesh2D;

	//����� �����//

	struct  MESH2D {

		bool Axis;                      //���������� �������������� �� ������� ���������             
		bool(*Geometry)(dVec2 x);	    //��������� ��������� �������
		double(*MeshDist)(dVec2 x);     //������������� ������� ����� (� ������������)
		vector <FemNode2D> Nodes;       //��� ���� �������
		vector <FemSide2D> Sides;       //������� ��������� �������
		vector <FemCell2D> Cells;       //������ �������� �������
		vector <const char*> Names;		//�������� ����������

		//���������� ������������
		void AddTriangle(dVec2 Xhead);
		
		//���������� �����
		void Meshing(dVec2 Xhead);
		
		//����������� �����, ���������� ����� � �������
		//Ndir - �������� ���������� � �������
		//Eps - ��������, ������������ �������� ����������� (��������� ������������ ����������)
		//Iteration - ���������� �������� ��� �����������
		void Optimization(int Ndir, double Eps, int Iteration);
		void ReadFormatK(const char*);
	};
}
