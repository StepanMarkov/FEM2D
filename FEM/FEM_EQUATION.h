#pragma once
#include "FEM_MESH.h"

struct TERM;
struct EQUATION;
struct BOUNDARY_CONDITION;

using namespace std;
using namespace FEM;
using namespace MESH;

//================================================//
//����� ��� ���������� ��������� �������
//����� ��� ���������� ���������

struct TERM {

	short Type;                     //��� ����������
	short Len;                      //���������� ����������
	int param1[4];                  //���������� ������� � ���������   
	double Constant;                //����������� ������������������

	std::function<double(FemCell2D c)> ConstCell;  //��������� ���������� ��� ���������� ������
	std::function<double(FemSide2D c)> ConstSide;  //���������� ����������� ��� ����������� �������� �������
};
//================================================//
 
//����� ��������� �������//
struct BOUNDARY_CONDITION {
	vector <TERM> TERMS;								//��������� ���������� �������
	vector <FemSide2D> Sides;							//�������� �������, ��� ���������� ��������� �������
	void SetArea(bool(*Set)(dVec2 x), FemMesh2D Mesh);  //���������� ��������� �������
};

//================================================//
//����� ����������������� ���������//
struct EQUATION {
	
	FemMesh2D Mesh;                                                   //��������� �� �����
	vector <TERM> TERMS;                                              //��������� ���������
	vector <BOUNDARY_CONDITION> BoundaryConditions;                   //��������� �������
	double Relax;													  //����������� ������ ����������
	int Neq;														  //����� ����������������� ���������
	bool Static;													  //��������� �� �������� ��������� ������������

	void Intitialization(double(*InitCondition)(dVec2 X));            //�������������
	void ConstructMatrix();                                           //���������� �������
	void UpdateMatrix();                                              //������������� �������
	void InternalIteration();										  //���������� ��������
	void NewTimeStep();                                               //������� �� ����� ��������� ����
            
	enum TYPE_TERM { TIME, LAPLASS, THIS, SOURCE, CONVECTIVE, DIVERGENCE, DIRIVATE, LAPLASS_CONVECTIVE, TABLE};
};
//===========================================================================================================================//
