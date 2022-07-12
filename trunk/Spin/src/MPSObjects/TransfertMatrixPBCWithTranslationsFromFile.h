#ifndef _TRANSFERTMATRIXPBCWITHTRANSLAIIONSFROMFILE_H
#define _TRANSFERTMATRIXPBCWITHTRANSLAIIONSFROMFILE_H



#include "Hamiltonian/AbstractHamiltonian.h"
#include "MPSObjects/AbstractTransfertMatrixPBC.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractSpinChainWithTranslations.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"


class TransfertMatrixPBCWithTranslationsFromFile : public AbstractTransfertMatrixPBC
{
 protected:

  // momentum along the x direction
  int XMomentum;
  // periodicity in the x direction
  int MaxXMomentum;

  //array containing all the phase factors that are needed when computing matrix elements
  Complex* ExponentialFactors;

 public:

  TransfertMatrixPBCWithTranslationsFromFile ();
  TransfertMatrixPBCWithTranslationsFromFile (MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture);
  ~TransfertMatrixPBCWithTranslationsFromFile ();  

  

  // multiply a vector by the current hamiltonian and store result in another vector
  // low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // return value = reference on vectorwhere result has been stored
  //  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination);
  //  virtual ComplexVector& LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination);

  // multiply a vector by the current hamiltonian and store result in another vector
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);    

  virtual void SetHilbertSpace(AbstractHilbertSpace * hilbertSpace);

  // evaluate all exponential factors
  //   
  virtual void EvaluateExponentialFactors();
    
};

#endif

