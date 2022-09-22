#include "FEM.h"
#pragma once
#include "FEM_SOLVER.h"
#include "FEM_RESULTS.h"
#include "FEM_LAME.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

bool Geom(dVec2 X) {
	if (X.length2() <= 0.03 * 0.03) return false;
	if (X[0] > 0.15) return false;
	if (X[1] < 0.0)  return false;
	if (X[1] > 0.25) return false;
	return true;
}



int main() {

	SOLVER Solver;
	Solver.MeshMain->ReadFormatK("task_mesh_fine.k");
	Solver.MeshMain->Axis = false;
	SetLameEquations(Solver);
	Solver.MaxInterIter = 20000;
	Solver.TimeSteps = 1;
	Solver.Calculation();
	WriteMeshVTK(Solver.MeshMain, "SOLUTION.vtk");
	system("PAUSE");
	return 0;
}
