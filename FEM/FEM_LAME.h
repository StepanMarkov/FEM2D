#include "FEM_SOLVER.h"

double MeshSize(dVec2 X) {
	return 0.0001;
	//return min(0.001+(X.length2() - 0.03 * 0.03), 0.015);
}

bool LEFT(dVec2 X) { return ((X[0] >= 0.0 - MeshSize(X) * 0.1) && (X[0] < 0.0 + MeshSize(X) * 0.1)); }
bool RIGHT(dVec2 X) { return ((X[0] >= 0.15 - MeshSize(X) * 0.1) && (X[0] < 0.15 + MeshSize(X) * 0.1)); }
bool  DOWN(dVec2 X) { return ((X[1] >= -MeshSize(X) * 0.1) && (X[1] < MeshSize(X) * 0.1)); }
bool  UP(dVec2 X) { return ((X[1] >= 0.25-MeshSize(X) * 0.1) && (X[1] < 0.25 + MeshSize(X) * 0.1)); }
bool NOTRIGHT(dVec2 X) { return !RIGHT(X); }
bool NOTDOWN(dVec2 X) { return !DOWN(X); };
bool NOT_LEFT_RIGHT(dVec2 X) { return (!LEFT(X)) && (!RIGHT(X)); }
bool NOT_LEFT_DOWN(dVec2 X) { return (!LEFT(X)) && (!DOWN(X)); }
bool OTHER(dVec2 X) { return (!LEFT(X)) && (!RIGHT(X)) && (!DOWN(X)) && (!UP(X)); }
bool   ALL(dVec2 X) { return true; }

const double koef = 1e+14;
const double e = 2E+11;
const double v = 0.3;
const double lambda = e * v / (1.0 + v) / (1 - 2.0 * v);
const double mu = e / 2.0 / (1.0 + v);
const double force = 1E+6;

void SetLameEquations(SOLVER& Solver) {

	TERM Term;
	BOUNDARY_CONDITION Bound;
	EQUATION Equation;

	Equation.Mesh = Solver.MeshMain;
	Equation.Static = true;

	Solver.MeshMain->Names.push_back("UX");
	Equation.TERMS.clear();
	Equation.BoundaryConditions.clear();

	Term.Type = EQUATION::LAPLASS;
	Term.Constant = 1.0;
	Equation.TERMS.push_back(Term);

	Term.Type = EQUATION::DIRIVATE;
	Term.ConstCell = [&](FemCell2D c) {return 1.0/(1.0-2.0*v); };
	Term.param1[0] = 2;
	Term.param1[1] = 1;
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

	Bound.TERMS.clear();
	Bound.Sides.clear();
	Term.Len = 1;
	Term.param1[0] = 2;
	Term.ConstSide = [](FemSide2D s) {return  -s->Norm[0] * lambda / 2.0 / mu; };
	Bound.TERMS.push_back(Term);

	Term.Len = 0;
	Term.ConstSide = [](FemSide2D s) {return s->Norm[1] * s->CellLink->GetGradient(0, 2); };
	Bound.TERMS.push_back(Term);

	Bound.SetArea(OTHER, Solver.MeshMain);
	Bound.SetArea(UP, Solver.MeshMain);
	Equation.BoundaryConditions.push_back(Bound);


	Equation.Intitialization();
	Equation.Relax = 0.7;
	Solver.Equations.push_back(Equation);

	Solver.MeshMain->Names.push_back("UY");
	Equation.TERMS.clear();
	Equation.BoundaryConditions.clear();

	Term.Type = EQUATION::LAPLASS;
	Term.Constant = 1.0;
	Equation.TERMS.push_back(Term);

	Term.Type = EQUATION::DIRIVATE;
	Term.ConstCell = [&](FemCell2D c) {return  1.0 / (1.0 - 2.0 * v); };
	Term.param1[0] = 2;
	Term.param1[1] = 2;
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

	Bound.TERMS.clear();
	Bound.Sides.clear();
	Term.Len = 1;
	Term.param1[0] = 2;
	Term.ConstSide = [](FemSide2D s) {return  -s->Norm[1] * lambda / 2.0 / mu; };
	Bound.TERMS.push_back(Term);

	Term.Len = 0;
	Term.ConstSide = [](FemSide2D s) {return s->Norm[0] * s->CellLink->GetGradient(1, 1); };
	Bound.TERMS.push_back(Term);

	Bound.SetArea(RIGHT, Solver.MeshMain);
	Bound.SetArea(OTHER, Solver.MeshMain);
	Bound.SetArea(UP, Solver.MeshMain);
	Equation.BoundaryConditions.push_back(Bound);

	Equation.Intitialization();
	Equation.Relax = 0.7;
	Solver.Equations.push_back(Equation);

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
		return  (s->CellLink->GetGradient(0, 1) +
			s->CellLink->GetGradient(1, 2)) * koef;
	};
	Bound.TERMS.push_back(Term);
	Term.Len = 1;
	Term.param1[0] = 2;
	Term.ConstSide = [](FemSide2D s) {return  -koef; };
	Bound.TERMS.push_back(Term);
	Bound.SetArea(ALL, Solver.MeshMain);
	Equation.BoundaryConditions.push_back(Bound);

	Equation.Intitialization();
	Equation.Relax = 0.7;
	Solver.Equations.push_back(Equation);

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
	Term.param1[1] = 1;
	Term.ConstCell = [](FemCell2D c) {return  -2.0 * mu; };
	Equation.TERMS.push_back(Term);
	Equation.Intitialization();
	Equation.Relax = 1.0;
	Solver.Equations.push_back(Equation);

}