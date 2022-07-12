#ifndef _LEFTORTHONORMALIZEALGORITHM_H 
#define _LEFTORTHONORMALIZEALGORITHM_H 

#include "ComplexMPOperatorOBC.h"
#include "Architecture/AbstractArchitecture.h"
#include "LanczosAlgorithm/LanczosManager.h"
#include "Matrix/RealDiagonalMatrix.h"

class LeftOrthonormalizeAlgorithm 
{

 protected:
  
  unsigned int PhysicalDimension;
  ComplexMatrix * Mps;
  ComplexMatrix * MpsInLeftForm;
  ComplexMatrix CenterMatrix;
  ComplexMatrix OldCenterMatrix;
  double Accuracy;  
  double Eigenvalue;
  
 public:
  
  LeftOrthonormalizeAlgorithm (ComplexMatrix * mps, unsigned int physicalDimension, double accuracy, ComplexMatrix & initialGuess);
  ~LeftOrthonormalizeAlgorithm();
  void RunAlgorithm();
  
 protected:
  
};

#endif
