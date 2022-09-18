#pragma once

#include <memory>
#include <vector>
#include <functional>
#include "FEM_MATH.h"

using namespace std;
using namespace MATH;

namespace FEM
{
	struct NODE2D;                                        //узел
	struct SIDE2D;                                        //элемент внешней границы
	struct CELL2D;                                        //ячейка
	//typedef shared_ptr <NODE2D> FemNode2D;
	//typedef shared_ptr <SIDE2D> FemSide2D;
	//typedef shared_ptr <CELL2D> FemCell2D;
	typedef NODE2D* FemNode2D;
	typedef SIDE2D* FemSide2D;
	typedef CELL2D* FemCell2D;

	//=================================//
    //класс узла расчетной сетки//
	struct NODE2D {

		int Num;		//глобальная нумерация узла
		dVec2 X;		//координата узла
		dVec2 Dir;      //вспомагательный вектор

		vector <FemNode2D> Nodes;          //соседние узлы
		vector <FemCell2D> Cells;          //соседние ячейки
		vector <double> Value;             //значения неизвестных переменных
		vector <double> ValuePre;          //значения на предыдущем временном слое
		vector<vector <double>> Matrix;    //коэффициенты матрицы
		vector<vector <double>> MatrixPre; //коэффициенты матрицы для значений на предыдущем временном слое
		unsigned char IndexSelf[6];       //определяет индексацию (для построения разреженной матрицы)
		unsigned char Index[6][3];        //[i][j], i - номер соседней ячейки, j - номер узла
		bool SwBool;	//вспомогательная переменная для обхода графа
		bool SwBound;	//определяет лежит ли узел на границе

		void SetInternalIndex();			//индексация
	};

	//======================================================================================================================================================================//
	//класс элемента границы расчетной сетки//

	struct SIDE2D {

		double Len;                          //длина элемента границы (нулевой момент функции формы)
		dVec2 Norm;                          //нормаль к элементу границы
		dVec2 Xc;                            //центр
		FemNode2D NodeBegin;                 //начало
		FemNode2D NodeEnd;                   //конец
		FemSide2D SideL;                     //связь с элементом границы по часовой стрелке
		FemSide2D SideR;                     //связь с элементом границы против часовой стрелке
		FemCell2D CellLink;                  //определяет, какому элементу принадлежит эта граница
		double ShapeMoment1[3];              //первый моменты функции формы (интегрирование по границе)
		double ShapeMoment2[3][3];           //второй момент функции формы

		//создает границу
		void CreateSide(FemNode2D n1, FemNode2D n2, FemCell2D Cell);
		//расчитывывает моменты функции формы
		//bool axis - является ли система координат цилиндрической
		void FormCalculate(bool Axis);
		//интегрирование некоторой функции по границе, Num - число разбиений
		double IntegrateSide(const function<double(dVec2)>& Func, int Num);
	};


	//======================================================================================================================================================================//
	//класс расчетной ячейки//	
	struct CELL2D {

		int Num;                   //глобальный номер узла ячейки
		double Area;               //площадь элемента (неприведенная, определитель)
		dVec2 Xc;                  //центр элемента
		dVec2 Xmin;                //нижняя левая координата
		dVec2 Xmax;                //верхняя правая координата
		vector <FemNode2D> Nodes;  //из каких узлов состоит
		//vector <FemSide2D> Sides;  //указатели на элементы границы
		double Form[3][3];         //коэффициенты формы
		double ShapeMoment0;       //нулевой момент функции формы
		double ShapeMoment1[3];    //первый момент shape-функции (интегрирование по элементу)
		double ShapeMoment2[3][3]; //второй момент shape-функции
		bool SwBound;			   //указывает, является ли узел граничным

		void FormCalculate(bool Axis);     //рассчитывает переменные, зависящие от формы
		bool PointInside(dVec2 px);        //определяет, внутри элемента ли координата px
		dVec2 GetBarycenter();             //возвращает барицентр элемента
		double GetGradient(int i, int k);  //возвращает градиент функции (проекция на k)
		double IntegrateElem(const std::function<double(dVec2)>& Func, int Num);    //интегрирование по элементу некоторой функции
	};
}