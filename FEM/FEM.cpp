#include "FEM.h"
#pragma once
#include "FEM_SOLVER.h"
#include "FEM_RESULTS.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

bool  LEFT(dVec2 X) { return ((X[0] >= -0.001) && (X[0] < 0.001));}
bool RIGHT(dVec2 X) { return ((X[0] >= 0.15 - 0.001) && (X[0] < 0.15 + 0.001));}
bool  DOWN(dVec2 X) { return ((X[1] >= -0.001) && (X[1] < 0.001)); }
bool   ALL(dVec2 X) { return true; }

const double koef = 1e+10;
const double e = 2E+11;
const double v = 0.3;
const double lambda = e * v / (1.0 + v) / (1 - 2.0 * v);
const double mu = e / 2.0 / (1.0 + v);
const double force = 1E+6;

int main() {

	TERM Term;
	BOUNDARY_CONDITION Bound;
	EQUATION Equation;
	SOLVER Solver;

	Solver.MeshMain->ReadFormatK("task_mesh_extrafine.k");
	Solver.MeshMain->Axis = false;
	Equation.Mesh = Solver.MeshMain;
	Equation.Static = true;
	
	//-----------------------------------------------------//
	Solver.MeshMain->Names.push_back("UX");
	Equation.TERMS.clear();
	Equation.BoundaryConditions.clear();
	
	Term.Type = EQUATION::LAPLASS;
	Term.Constant = 1.0;
	Equation.TERMS.push_back(Term);

	Term.Type = EQUATION::DIRIVATE;
	Term.ConstCell = [&](FemCell2D c) {return 1.0 + lambda / mu; };
	Term.param1[0] = 2;
	Term.param1[1] = 0;
	Equation.TERMS.push_back(Term);
	
	Bound.TERMS.clear();
	Bound.Sides.clear();
	Term.Len = 0;
	Term.ConstSide = [](FemSide2D s) {return  force / 2.0 / mu; };
	Bound.TERMS.push_back(Term);
	Term.Len = 1;
	Term.param1[0] = 2;
	Term.ConstSide = [](FemSide2D s) {return  -lambda / 2.0 / mu; };
	Bound.TERMS.push_back(Term);
	Bound.SetArea(RIGHT, Solver.MeshMain);
	Equation.BoundaryConditions.push_back(Bound);
	
	Bound.TERMS.clear();
	Bound.Sides.clear();
	Term.Len = 0;
	Term.ConstSide = [](FemSide2D s) {return  0.0 * koef; };
	Bound.TERMS.push_back(Term);
	Term.Len = 1;
	Term.param1[0] = 0;
	Term.ConstSide = [](FemSide2D s) {return  -koef; };
	Bound.TERMS.push_back(Term);
	Bound.SetArea(LEFT, Solver.MeshMain);
	Equation.BoundaryConditions.push_back(Bound);
	
	Equation.Intitialization();
	Equation.Relax = 0.5;
	Solver.Equations.push_back(Equation);
	//-----------------------------------------------------//
	Solver.MeshMain->Names.push_back("UY");
	Equation.TERMS.clear();
	Equation.BoundaryConditions.clear();
	
	Term.Type = EQUATION::LAPLASS;
	Term.Constant = 1.0;
	Equation.TERMS.push_back(Term);
	
	Term.Type = EQUATION::DIRIVATE;
	Term.ConstCell = [&](FemCell2D c) {return 1.0 + lambda / mu; };
	Term.param1[0] = 2;
	Term.param1[1] = 1;
	Equation.TERMS.push_back(Term);
	
	Bound.TERMS.clear();
	Bound.Sides.clear();
	Term.Len = 0;
	Term.ConstSide = [](FemSide2D s) {return  0.0 * koef; };
	Bound.TERMS.push_back(Term);
	Term.Len = 1;
	Term.param1[0] = 1;
	Term.ConstSide = [](FemSide2D s) {return  -koef; };
	Bound.TERMS.push_back(Term);
	Bound.SetArea(DOWN, Solver.MeshMain);
	Equation.BoundaryConditions.push_back(Bound);
	
	Equation.Intitialization();
	Equation.Relax = 0.5;
	Solver.Equations.push_back(Equation);
	//----------------------------------------------------//
	Solver.MeshMain->Names.push_back("TETA");
	Equation.TERMS.clear();
	Equation.BoundaryConditions.clear();
	
	Term.Type = EQUATION::LAPLASS;
	Term.Constant = 1.0;
	Equation.TERMS.push_back(Term);
	
	Bound.TERMS.clear();
	Bound.Sides.clear();
	Term.Len = 0;
	Term.ConstSide = [](FemSide2D s) {
		return  (s->CellLink->GetGradient(0, 0) +
			s->CellLink->GetGradient(1, 1)) * koef;
	};
	Bound.TERMS.push_back(Term);
	Term.Len = 1;
	Term.param1[0] = 2;
	Term.ConstSide = [](FemSide2D s) {return  -koef; };
	Bound.TERMS.push_back(Term);
	Bound.SetArea(ALL, Solver.MeshMain);
	Equation.BoundaryConditions.push_back(Bound);
	
	Equation.Intitialization();
	Equation.Relax = 0.5;
	Solver.Equations.push_back(Equation);
	////----------------------------------------------------//
	Solver.MeshMain->Names.push_back("SXX");
	Equation.TERMS.clear();
	Equation.BoundaryConditions.clear();
	Term.Type = EQUATION::THIS;
	Term.Constant = 1.0;
	Equation.TERMS.push_back(Term);
	
	Term.Type = EQUATION::SOURCE;
	Term.Len = 1;
	Term.param1[0] = 2;
	Term.ConstCell = [](FemCell2D c) {return  -lambda; };
	Equation.TERMS.push_back(Term);

	Term.Type = EQUATION::DIRIVATE;
	Term.param1[0] = 0;
	Term.param1[1] = 0;
	Term.ConstCell = [](FemCell2D c) {return  -2*mu; };
	Equation.Intitialization();
	Equation.Relax = 1.0;
	Solver.Equations.push_back(Equation);
	
	Solver.MaxInterIter = 5000;
	Solver.TimeSteps = 1;
	Solver.Calculation();
	WriteMeshVTK(Solver.MeshMain, "SOLUTION.vtk", 0);
	system("PAUSE");

	return 0;
}
