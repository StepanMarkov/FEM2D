#pragma once

#include "FEM_SOLVER.h"
#include <fstream>

using namespace FEM;
using namespace MESH;

void WriteMeshVTK(FemMesh2D M, const char *Name, int NumStep);