#include "FEM.h"
#pragma once
#include "FEM_SOLVER.h"
#include "FEM_RESULTS.h"
#include "FEM_LAME.h"
#include <iostream>
#include <fstream>
#include <ctime>

using namespace std;

int main() {

	SOLVER Solver;
	Solver.MeshMain->ReadFormatK("task_mesh_fine.k");
	Solver.MeshMain->Axis = false;
	SetLameEquations(Solver);
	Solver.MaxInterIter = 2000;
	Solver.TimeSteps = 1;
	Solver.Calculation();
	WriteMeshVTK(Solver.MeshMain, "SOLUTION.vtk");
	system("PAUSE");
	return 0;
}
