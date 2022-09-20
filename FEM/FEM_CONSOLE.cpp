#pragma once
#include "FEM_SOLVER.h"
#include "FEM_RESULTS.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

FORCEINLINE bool Geom(dVec2 X)
{
	//if (X[0] * X[0] + X[1] * X[1] < 0.15 * 0.15) return false;
	if (X[0] < 0.0)  return false;
	if (X[0] > 0.3) return false;
	if (X[1] < 0.0)  return false;
	if (X[1] > 0.5) return false;
	return true;
}

double Init0(dVec2 X) {
	return 0.0;
}

bool LEFT(dVec2 X) {
	return ((X[0] >= -0.001) && (X[0] < 0.001));
}

bool RIGHT(dVec2 X) {
	return ((X[0] >= 0.3 - 0.001) && (X[0] < 0.3 + 0.001));
}

bool DOWN(dVec2 X) {
	return ((X[1] >= -0.001) && (X[1] < 0.001));
}

bool ALL(dVec2 X) {
	return true;
}

FORCEINLINE double MeshDist(dVec2 X)
{
	return 0.01;
}

const double koef = 1e+10;

int main() {

	SOLVER Solver;
	Solver.MeshMain->ReadFormatK("task_mesh_coarse.k");
	//ReadFormatK

	//TERM Term;
	//BOUNDARY_CONDITION Bound;
	//EQUATION Equation;
	//SOLVER Solver;
	//
	//Solver.MeshMain->Geometry = Geom;
	//Solver.MeshMain->Axis = false;
	//Solver.MeshMain->MeshDist = MeshDist;
	//Solver.MeshMain->Meshing({ 0.15,0.25 });
	//
	//
	//Equation.Mesh = Solver.MeshMain;
	//Equation.Static = true;
	//
	//Solver.MeshMain->Names.push_back("UX");
	//Equation.TERMS.clear();
	//Equation.BoundaryConditions.clear();
	//
	//Term.Type = EQUATION::LAPLASS;
	//Term.Constant = 1.0;
	//Equation.TERMS.push_back(Term);
	//
	//Term.Type = EQUATION::DIRIVATE;
	//Term.Len = 1;
	//Term.ConstCell = [&](FemCell2D c) {return 6.0; };
	//Term.param1[0] = 2;
	//Term.param1[1] = 0;
	//Equation.TERMS.push_back(Term);
	//
	//Bound.TERMS.clear();
	//Bound.Sides.clear();
	//Term.Len = 0;
	//Term.ConstSide = [](FemSide2D s) {return  -1.0 * koef; };
	//Bound.TERMS.push_back(Term);
	//Term.Len = 1;
	//Term.param1[0] = 0;
	//Term.ConstSide = [](FemSide2D s) {return  -koef; };
	//Bound.TERMS.push_back(Term);
	//Bound.SetArea(LEFT, Solver.MeshMain);
	//Equation.BoundaryConditions.push_back(Bound);
	//
	//Bound.TERMS.clear();
	//Bound.Sides.clear();
	//Term.Len = 0;
	//Term.ConstSide = [](FemSide2D s) {return  1.0 * koef; };
	//Bound.TERMS.push_back(Term);
	//Term.Len = 1;
	//Term.param1[0] = 0;
	//Term.ConstSide = [](FemSide2D s) {return  -koef; };
	//Bound.TERMS.push_back(Term);
	//Bound.SetArea(RIGHT, Solver.MeshMain);
	//Equation.BoundaryConditions.push_back(Bound);
	//
	//Equation.Intitialization(Init0);
	//Equation.Relax = 0.6;
	//Solver.Equations.push_back(Equation);
	////-----------------------------------------------------//
	//Solver.MeshMain->Names.push_back("UY");
	//Equation.TERMS.clear();
	//Equation.BoundaryConditions.clear();
	//
	//Term.Type = EQUATION::LAPLASS;
	//Term.Constant = 1.0;
	//Equation.TERMS.push_back(Term);
	//
	//Term.Type = EQUATION::DIRIVATE;
	//Term.Len = 1;
	//Term.ConstCell = [&](FemCell2D c) {return 6.0; };
	//Term.param1[0] = 2;
	//Term.param1[1] = 1;
	//Equation.TERMS.push_back(Term);
	//
	//Bound.TERMS.clear();
	//Bound.Sides.clear();
	//Term.Len = 0;
	//Term.ConstSide = [](FemSide2D s) {return  0.0 * koef; };
	//Bound.TERMS.push_back(Term);
	//Term.Len = 1;
	//Term.param1[0] = 1;
	//Term.ConstSide = [](FemSide2D s) {return  -koef; };
	//Bound.TERMS.push_back(Term);
	//Bound.SetArea(DOWN, Solver.MeshMain);
	//Equation.BoundaryConditions.push_back(Bound);
	//
	//Equation.Intitialization(Init0);
	//Equation.Relax = 0.6;
	//Solver.Equations.push_back(Equation);
	////----------------------------------------------------//
	//Solver.MeshMain->Names.push_back("TETA");
	//Equation.TERMS.clear();
	//Equation.BoundaryConditions.clear();
	//
	//Term.Type = EQUATION::LAPLASS;
	//Term.Constant = 1.0;
	//Equation.TERMS.push_back(Term);
	//
	//Bound.TERMS.clear();
	//Bound.Sides.clear();
	//Term.Len = 0;
	//Term.ConstSide = [](FemSide2D s) {
	//	return  (s->CellLink->GetGradient(0, 0) +
	//		s->CellLink->GetGradient(1, 1)) * koef;
	//};
	//Bound.TERMS.push_back(Term);
	//Term.Len = 1;
	//Term.param1[0] = 2;
	//Term.ConstSide = [](FemSide2D s) {return  -koef; };
	//Bound.TERMS.push_back(Term);
	//Bound.SetArea(ALL, Solver.MeshMain);
	//Equation.BoundaryConditions.push_back(Bound);
	//
	//Equation.Intitialization(Init0);
	//Equation.Relax = 0.5;
	//Solver.Equations.push_back(Equation);
	////----------------------------------------------------//
	//
	//cout << "Mesh" << endl;
	//cout << "Nodes " << Solver.MeshMain->Nodes.size() << endl;
	//cout << "Sides " << Solver.MeshMain->Sides.size() << endl;
	//cout << "Elems " << Solver.MeshMain->Cells.size() << endl;
	//
	//Solver.MaxInterIter = 20000;
	//Solver.TimeSteps = 1;
	//
	//time_t t;
	//t = -time(0);
	//Solver.Calculation();
	//t += time(0);
	//
	//cout << "Computational Time = " << t << " sec" << endl;
	//cout << " Done" << endl;
	//
	//WriteMeshVTK(Solver.MeshMain, "SOLUTION.vtk", 0);
	//system("PAUSE");
	//
	return 0;
}
