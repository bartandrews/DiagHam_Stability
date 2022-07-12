#ifndef _COMPLEXPEPSABSUBLATTICETRANSFERTMATRIXPBCWITHTRANSLAIIONS_H
#define _COMPLEXPEPSABSUBLATTICETRANSFERTMATRIXPBCWITHTRANSLAIIONS_H

#include "MPSObjects/ComplexPEPSTransfertMatrixPBC.h"
#include "MPSObjects/AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractDoubledSpinChainWithTranslations.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class AbstractMPSSite;

class ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations : public  AbstractTransfertMatrixPBC
{
 protected:

  // momentum along the x direction
  int XMomentum;
  // periodicity in the x direction
  int MaxXMomentum;

  Complex **** ValuesNonZeroTensorElementTopLeft;
  int *** NbrNonZeroTensorElementTopLeft;
  int **** IndiceRightNonZeroTensorElementTopLeft;
  int **** IndiceBottomNonZeroTensorElementTopLeft;

  unsigned long * NewHilbertSpace1;
  unsigned long * NewHilbertSpace2;
  unsigned long * OldHilbertSpace;
  ComplexVector * TmpVector1;
  ComplexVector * TmpVector2;
  ComplexVector * StartVector;
  ComplexVector * EndVector;
  RealDiagonalMatrix * BoundaryMatrix; 
  
 public:

   ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations ();
   ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations (MultiColumnASCIIFile & tensorElementsFileA,MultiColumnASCIIFile & tensorElementsFileB ,AbstractArchitecture * architecture);
   ~ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations ();  
  
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
  inline  int GetBondDimension() const {return this->MPOBondDimension; }; 

  
 protected:

  inline void GetBraAndKetIndexFromCommonIndex(unsigned int communIndex,unsigned int & braIndex ,unsigned int & ketIndex);
  inline unsigned int GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex );

  void LowLevelAddMultiplyOnLastSite(int topValue, int sublatticeIndex);
  void LowLevelAddMultiplyOnAnySite(int position, int sublatticeIndex);
  void LowLevelAddMultiplyOnFirstSite(int topIndice, int sublatticeIndex);
    
  virtual inline long GetNewIndexFromOldIndex(unsigned long oldIndex, int oldPhysicalSpin, int newPhysicalSpin, int oldVirtualSpin, int newVirtualSpin, int position);
  virtual void InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile, int sublatticeIndex);

    
};

inline void ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::GetBraAndKetIndexFromCommonIndex(unsigned int communIndex,unsigned int & braIndex ,unsigned int & ketIndex)
{
  braIndex = communIndex%this->MPOBondDimension;
  ketIndex = communIndex/this->MPOBondDimension;
}


inline unsigned int ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex )
{
  return ketIndex  * this->MPOBondDimension + braIndex;
}


inline long ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::GetNewIndexFromOldIndex(unsigned long oldIndex, int oldPhysicalSpin, int newPhysicalSpin, int oldVirtualSpin, int newVirtualSpin, int position)
{ 
  return  oldIndex + (newPhysicalSpin - oldPhysicalSpin) *this->PowerD[position] + newVirtualSpin - oldVirtualSpin;
}


#endif

