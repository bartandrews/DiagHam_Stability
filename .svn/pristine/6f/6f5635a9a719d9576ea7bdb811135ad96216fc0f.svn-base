#ifndef _COMPLEXPEPSTRANSFERTMATRIXPBCWITHTRANSLAIIONS_H
#define _COMPLEXPEPSTRANSFERTMATRIXPBCWITHTRANSLAIIONS_H

#include "MPSObjects/ComplexPEPSTransfertMatrixPBC.h"
#include "MPSObjects/AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractDoubledSpinChainWithTranslations.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class AbstractMPSSite;

class ComplexPEPSTransfertMatrixPBCWithTranslations : public  ComplexPEPSTransfertMatrixPBC 
{
 protected:

  // momentum along the x direction
  int XMomentum;
  // periodicity in the x direction
  int MaxXMomentum;

 public:

  ComplexPEPSTransfertMatrixPBCWithTranslations ();
  ComplexPEPSTransfertMatrixPBCWithTranslations (MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture);
  ~ComplexPEPSTransfertMatrixPBCWithTranslations ();  
  
  // multiply a vector by the current hamiltonian for a given range of indices 
  // and store result in another vector, low level function (no architecture optimization)
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  //  virtual RealVector& LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent);
  virtual ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent);
//  virtual ComplexVector* LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent);

  virtual void SetHilbertSpace(AbstractHilbertSpace * hilbertSpace);
  
 protected:

    
};


#endif

