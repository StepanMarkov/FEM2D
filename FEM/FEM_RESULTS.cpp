#include "FEM_RESULTS.h"
#include <string>

void WriteMeshVTK(FemMesh2D M, const char *Name, int NumStep) {

	ofstream Out = ofstream(ofstream(std::to_string(NumStep) + "SOLUTION.vtk"));
	
	Out << "# vtk DataFile Version 2.0" << endl;
	Out << "Unstructured Grid Example" << endl;
	Out << "ASCII" << endl;
	Out << "DATASET UNSTRUCTURED_GRID" << endl;
	Out << "POINTS " << M->Nodes.size() << ' ' << "float" << endl;
	
	for (auto& N : M->Nodes)
		Out << N->X[0] << ' ' << N->X[1] << ' ' << 0.0 << endl;
	
	int Size = M->Cells.size();
	
	for (auto& C : M->Cells)
		Size += C->Nodes.size();
	
	Out << "CELLS " << M->Cells.size() << ' ' << Size << endl;
	
	for(auto& C : M->Cells) {
		
		Out << C->Nodes.size() << ' ';
		for (auto& N : C->Nodes)
			Out << N->Num << ' ';
		Out << endl;
	}
	
	Out << "CELL_TYPES " << M->Cells.size() << endl;

	for (int i(0); i != M->Cells.size(); ++i)
		Out << "5" << endl;
	
	Out << "POINT_DATA " << M->Nodes.size() << endl;
	
	for (int i(0); i != M->Names.size(); i++) {	
		Out << "SCALARS " << M->Names[i] << " float 1" << endl;
		Out << "LOOKUP_TABLE default" << endl;
	
		for (auto& N : M->Nodes)
			Out << N->Value[i] << endl;
	}
	
	Out.close();
}

