#ifndef _DMRGFINITESIZEREALOBCMAINTASK_H 
#define _DMRGFINITESIZEREALOBCMAINTASK_H

#include "RealMPOperatorOBC.h"
#include "Architecture/AbstractArchitecture.h"
#include "LanczosAlgorithm/LanczosManager.h"

class   RealDiagonalMatrix;

class DMRGFiniteSizeRealOBCMainTask 
{
 protected:
  const int NbrSites;
  
  RealMPSSite * LatticeSite;
  AbstractArchitecture * Architecture;
  LanczosManager *  AlgorithmManager;
  AbstractMPOperatorOBC * MPOperator;
  RealDiagonalMatrix PreviousSingularValues;
  int NbrSweep;
  int MaximumBondDimension;
  double PreviousEnergy;
 public:
  
  DMRGFiniteSizeRealOBCMainTask(RealMPSSite * latticeSite, AbstractMPOperatorOBC * mPOperator, int nbrSites, int NbrSweep,int MaximumBondDimension,  AbstractArchitecture * architecture, LanczosManager* lanczos);
  virtual ~DMRGFiniteSizeRealOBCMainTask();
  void RunAlgorithm();
  
 protected:
  void InitializeLattice();
  void InitializeLatticeUsingIDMRG();
  void InitializeLatticeUsingIDMRGAndStatePrediction();
  void OptimizeUsingLanczosLanczosAlgorithm (int siteIndex);
  void TwoSiteOptimizationUsingLanczosLanczosAlgorithm (RealMPSSite * leftSite , RealMPSSite * rightSite, RealDiagonalMatrix & singularValues);
  void TwoSiteOptimizationUsingLanczosLanczosAlgorithmAndStatePrediction ( RealMPSSite * leftSite , RealMPSSite * rightSite, RealDiagonalMatrix & singularValues, RealVector * statePredict);
};

#endif
