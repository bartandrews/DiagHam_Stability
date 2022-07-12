#ifndef _DMRGINFINITESIZECOMPLEXMAINTASK_H 
#define _DMRGINFINITESIZECOMPLEXMAINTASK_H

#include "ComplexMPOperatorOBC.h"
#include "Architecture/AbstractArchitecture.h"
#include "LanczosAlgorithm/LanczosManager.h"
#include "Matrix/RealDiagonalMatrix.h"

class DMRGInfiniteSizeComplexMainTask 
{

 protected:
  unsigned int PhysicalDimension;
  AbstractArchitecture * Architecture;
  LanczosManager *  AlgorithmManager;
  AbstractMPOperatorOBC * MPOperator;
  RealDiagonalMatrix PreviousSingularValues;
  int MaximumBondDimension;
  double PreviousEnergy;
  
 public:
  
  DMRGInfiniteSizeComplexMainTask( unsigned int physicalDimension,AbstractMPOperatorOBC * mPOperator, int MaximumBondDimension, AbstractArchitecture * architecture, LanczosManager* lanczos);
  virtual ~DMRGInfiniteSizeComplexMainTask();
  void RunAlgorithm();
  
 protected:

//  void InitializeLatticeUsingIDMRGAndStatePrediction();
  void TwoSiteOptimizationUsingLanczosLanczosAlgorithm (ComplexMPSSite * leftSite , ComplexMPSSite * rightSite, RealDiagonalMatrix & singularValues);
  //  void TwoSiteOptimizationUsingLanczosLanczosAlgorithmAndStatePrediction ( ComplexMPSSite * leftSite , ComplexMPSSite * rightSite, RealDiagonalMatrix & singularValues, ComplexVector * statePredict);
  double ComputeFidelity (ComplexMPSSite * leftSite, RealDiagonalMatrix & singularValues);

};

#endif
