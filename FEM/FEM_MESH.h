#pragma once

#include "FEM_BASE.h"

using namespace std;
using namespace MATH;
using namespace FEM;

namespace MESH
{
	struct MESH2D;
	typedef shared_ptr <MESH2D> FemMesh2D;

	//класс сетки//

	struct  MESH2D {

		bool Axis;                      //определяет цикиндрическая ли система координат             
		bool(*Geometry)(dVec2 x);	    //геометрия расчетной области
		double(*MeshDist)(dVec2 x);     //распределение размера сетки (в пространстве)
		vector <FemNode2D> Nodes;       //все узлы области
		vector <FemSide2D> Sides;       //граница расчетной облатси
		vector <FemCell2D> Cells;       //ячейки расчтной области
		vector <const char*> Names;		//названия переменных

		//добавление треугольника
		void AddTriangle(dVec2 Xhead);
		
		//построение сетки
		void Meshing(dVec2 Xhead);
		
		//сглаживание сетки, прилипание узлов к границы
		//Ndir - точность прилипания к границе
		//Eps - параметр, определяющий скорость оптимизации (наподобии коэффициента релаксации)
		//Iteration - количество итераций при оптимизации
		void Optimization(int Ndir, double Eps, int Iteration);
		void ReadFormatK(const char*);
	};
}
