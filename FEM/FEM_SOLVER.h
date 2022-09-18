#pragma once
#include "FEM_EQUATION.h"

struct SOLVER;

//класс решателя
struct SOLVER {

	unsigned int MaxInterIter;
	unsigned int TimeSteps;
	FemMesh2D MeshMain;
	vector <EQUATION> Equations;

	void Calculation();
	
	SOLVER() {FemMesh2D M(new MESH2D()); MeshMain = M;}
	~SOLVER() {}
};