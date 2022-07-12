#include "Tools/Potential/AbstractPotential.h"
#include "Tools/Potential/OneDConstantCellPotential.h"
#include "Tools/Potential/TwoDConstantCellPotential.h"
#include "Tools/Potential/ThreeDConstantCellPotential.h"
#include "Tools/Potential/HardBoxPyramidQuantumDotThreeDConstantCellPotential.h"
#include "Tools/Potential/DotEmbeddedWellThreeDConstantCellPotential.h"

#include <iostream>
#include <fstream>
#include <math.h>

using std::cout;
using std::endl;
using std::ofstream;
using std::ifstream;

int main(int argc, char** argv)
{  

  DotEmbeddedWellThreeDConstantCellPotential* potential = new DotEmbeddedWellThreeDConstantCellPotential(13, 11, 9, 1, 1, 2, 5, 4, 2);  
  potential->ConstructPotential (1.1, 0.7);
  potential->SavePotential("DotWellPotential.txt");  

  /*
  HardBoxPyramidQuantumDotThreeDConstantCellPotential* potential = new HardBoxPyramidQuantumDotThreeDConstantCellPotential(50, 50, 30, 6, 10, 4, 20, 5);  
  potential->LoadDiagram("Diagram.txt");
  potential->ConstructPotential(0.1, 0.295, -0.001, 0.0245, 0.0245, -0.001, 0.9, 2.64, false, "Parameter.txt");
  potential->SavePotentialWithConstantField("PyramidDotPotential.txt"); 
  //Potential->SaveDiagram("Diagrambis.txt");
  */

  return 0;
}
