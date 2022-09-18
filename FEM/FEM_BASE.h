#pragma once

#include <memory>
#include <vector>
#include <functional>
#include "FEM_MATH.h"

using namespace std;
using namespace MATH;

namespace FEM
{
	struct NODE2D;                                        //����
	struct SIDE2D;                                        //������� ������� �������
	struct CELL2D;                                        //������
	//typedef shared_ptr <NODE2D> FemNode2D;
	//typedef shared_ptr <SIDE2D> FemSide2D;
	//typedef shared_ptr <CELL2D> FemCell2D;
	typedef NODE2D* FemNode2D;
	typedef SIDE2D* FemSide2D;
	typedef CELL2D* FemCell2D;

	//=================================//
    //����� ���� ��������� �����//
	struct NODE2D {

		int Num;		//���������� ��������� ����
		dVec2 X;		//���������� ����
		dVec2 Dir;      //��������������� ������

		vector <FemNode2D> Nodes;          //�������� ����
		vector <FemCell2D> Cells;          //�������� ������
		vector <double> Value;             //�������� ����������� ����������
		vector <double> ValuePre;          //�������� �� ���������� ��������� ����
		vector<vector <double>> Matrix;    //������������ �������
		vector<vector <double>> MatrixPre; //������������ ������� ��� �������� �� ���������� ��������� ����
		unsigned char IndexSelf[6];       //���������� ���������� (��� ���������� ����������� �������)
		unsigned char Index[6][3];        //[i][j], i - ����� �������� ������, j - ����� ����
		bool SwBool;	//��������������� ���������� ��� ������ �����
		bool SwBound;	//���������� ����� �� ���� �� �������

		void SetInternalIndex();			//����������
	};

	//======================================================================================================================================================================//
	//����� �������� ������� ��������� �����//

	struct SIDE2D {

		double Len;                          //����� �������� ������� (������� ������ ������� �����)
		dVec2 Norm;                          //������� � �������� �������
		dVec2 Xc;                            //�����
		FemNode2D NodeBegin;                 //������
		FemNode2D NodeEnd;                   //�����
		FemSide2D SideL;                     //����� � ��������� ������� �� ������� �������
		FemSide2D SideR;                     //����� � ��������� ������� ������ ������� �������
		FemCell2D CellLink;                  //����������, ������ �������� ����������� ��� �������
		double ShapeMoment1[3];              //������ ������� ������� ����� (�������������� �� �������)
		double ShapeMoment2[3][3];           //������ ������ ������� �����

		//������� �������
		void CreateSide(FemNode2D n1, FemNode2D n2, FemCell2D Cell);
		//������������� ������� ������� �����
		//bool axis - �������� �� ������� ��������� ��������������
		void FormCalculate(bool Axis);
		//�������������� ��������� ������� �� �������, Num - ����� ���������
		double IntegrateSide(const function<double(dVec2)>& Func, int Num);
	};


	//======================================================================================================================================================================//
	//����� ��������� ������//	
	struct CELL2D {

		int Num;                   //���������� ����� ���� ������
		double Area;               //������� �������� (�������������, ������������)
		dVec2 Xc;                  //����� ��������
		dVec2 Xmin;                //������ ����� ����������
		dVec2 Xmax;                //������� ������ ����������
		vector <FemNode2D> Nodes;  //�� ����� ����� �������
		//vector <FemSide2D> Sides;  //��������� �� �������� �������
		double Form[3][3];         //������������ �����
		double ShapeMoment0;       //������� ������ ������� �����
		double ShapeMoment1[3];    //������ ������ shape-������� (�������������� �� ��������)
		double ShapeMoment2[3][3]; //������ ������ shape-�������
		bool SwBound;			   //���������, �������� �� ���� ���������

		void FormCalculate(bool Axis);     //������������ ����������, ��������� �� �����
		bool PointInside(dVec2 px);        //����������, ������ �������� �� ���������� px
		dVec2 GetBarycenter();             //���������� ��������� ��������
		double GetGradient(int i, int k);  //���������� �������� ������� (�������� �� k)
		double IntegrateElem(const std::function<double(dVec2)>& Func, int Num);    //�������������� �� �������� ��������� �������
	};
}