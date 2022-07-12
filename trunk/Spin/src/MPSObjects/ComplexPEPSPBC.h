#ifndef _COMPLEXPEPSPBC_H
#define _COMPLEXPEPSPBC_H

#include "MPSObjects/AbstractTransfertMatrixPBC.h"
#include "MPSObjects/AbstractMPSSite.h"
#include "Tensor/Tensor3.h"
#include "HilbertSpace/AbstractDoubledSpinChain.h"
#include "Architecture/AbstractArchitecture.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

class AbstractMPSSite;

class ComplexPEPSPBC 
{
 protected:
  
  unsigned int NbrNonZeroElements;
  unsigned int PhysicalDimension;
  unsigned int MPOBondDimension;
  
  AbstractArchitecture * Architecture;
  
  Complex **** ValuesNonZeroTensorElementTopLeft;
  int *** NbrNonZeroTensorElementTopLeft;
  int **** IndiceRightNonZeroTensorElementTopLeft;
  int **** IndiceBottomNonZeroTensorElementTopLeft;
  
 public:

  ComplexPEPSPBC ();
  ComplexPEPSPBC(MultiColumnASCIIFile & tensorElementsFile, AbstractArchitecture * architecture = 0);
  ~ComplexPEPSPBC ();  
    
  virtual void PrintTensorElements();
  
  inline int GetBondDimension() const {return this->MPOBondDimension; }; 

  virtual ComplexVector ComputeFockSpaceRepresentationOfAPEPS (int lx, int lylogtwo, ComplexDiagonalMatrix virtualSymmetryHorizontal, ComplexDiagonalMatrix virtualSymmetryVertical, bool horizontalFlag, bool verticalFlag);
  virtual ComplexVector ComputeFockSpaceRepresentationOfAPEPSSzConstraint (int lx, int lylogtwo, int sz, ComplexDiagonalMatrix virtualSymmetryHorizontal, ComplexDiagonalMatrix virtualSymmetryVertical, bool horizontalFlag, bool verticalFlag);
  
  virtual ComplexVector ComputeFockSpaceRepresentationOfAPEPSSzConstraint (int lx, int lylogtwo, int sz,  ComplexDiagonalMatrix virtualSymmetryHorizontal, ComplexDiagonalMatrix virtualSymmetryVertical, ComplexVector LeftVector, ComplexVector RightVector, bool horizontalFlag, bool verticalFlag);
  virtual ComplexVector ComputeFockSpaceRepresentationOfAPEPS (int lx, int lylogtwo, ComplexDiagonalMatrix virtualSymmetryHorizontal,  ComplexDiagonalMatrix virtualSymmetryVertical, ComplexVector LeftVector, ComplexVector RightVector, bool horizontalFlag, bool verticalFlag);

  virtual  void ComputeBlockTensor ();

 protected:
  
  inline void GetBraAndKetIndexFromCommonIndex(unsigned int communIndex,unsigned int & braIndex ,unsigned int & ketIndex);
  inline unsigned int GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex );

/*  void LowLevelAddMultiplyOnLastSite(int topValue);
  void LowLevelAddMultiplyOnAnySite(int position);
  void LowLevelAddMultiplyOnFirstSite(int topIndice);*/
    
  virtual void InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile);
};


inline void ComplexPEPSPBC::GetBraAndKetIndexFromCommonIndex(unsigned int communIndex,unsigned int & braIndex ,unsigned int & ketIndex)
{
  braIndex = communIndex%this->MPOBondDimension;
  ketIndex = communIndex/this->MPOBondDimension;
}


inline unsigned int  ComplexPEPSPBC::GetCommonIndexFromBraAndKetIndices(unsigned int braIndex, unsigned int ketIndex )
{
  return ketIndex  * this->MPOBondDimension + braIndex;
}



#endif

