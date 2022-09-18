#pragma once
#include "FEM_MESH.h"

struct TERM;
struct EQUATION;
struct BOUNDARY_CONDITION;

using namespace std;
using namespace FEM;
using namespace MESH;

//================================================//
//класс для построения граничных условий
//класс для построения уравнений

struct TERM {

	short Type;                     //тип слагаемого
	short Len;                      //количество множителей
	int param1[4];                  //определяет индексы в слагаемом   
	double Constant;                //коэффициент пропорциональности

	std::function<double(FemCell2D c)> ConstCell;  //возращает коэффциент для конкретной ячейки
	std::function<double(FemSide2D c)> ConstSide;  //возвращает коэффициент для конкретного элемента границы
};
//================================================//
 
//КЛАСС ГРАНИЧНЫХ УСЛОВИЙ//
struct BOUNDARY_CONDITION {
	vector <TERM> TERMS;								//структура граничного условия
	vector <FemSide2D> Sides;							//элементы границы, где поставлено граничное условие
	void SetArea(bool(*Set)(dVec2 x), FemMesh2D Mesh);  //добавление элементов границы
};

//================================================//
//КЛАСС ДИФФЕРЕНЦИАЛЬНОГО УРАВНЕНИЯ//
struct EQUATION {
	
	FemMesh2D Mesh;                                                   //указатель на сетку
	vector <TERM> TERMS;                                              //структура уравнения
	vector <BOUNDARY_CONDITION> BoundaryConditions;                   //граничные условия
	double Relax;													  //коэффициент нижней релаксации
	int Neq;														  //номер дифференциального уравнения
	bool Static;													  //уравнение не содержит временных зависимостей

	void Intitialization(double(*InitCondition)(dVec2 X));            //инициализация
	void ConstructMatrix();                                           //построение матрицы
	void UpdateMatrix();                                              //корректировка матрицы
	void InternalIteration();										  //внутренняя итерация
	void NewTimeStep();                                               //переход на новый временной слой
            
	enum TYPE_TERM { TIME, LAPLASS, THIS, SOURCE, CONVECTIVE, DIVERGENCE, DIRIVATE, LAPLASS_CONVECTIVE, TABLE};
};
//===========================================================================================================================//
