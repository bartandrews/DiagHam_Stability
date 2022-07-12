#include "ComplexPEPSTransfertMatrixPBCWithTranslations.h"

#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

ComplexPEPSTransfertMatrixPBCWithTranslations::ComplexPEPSTransfertMatrixPBCWithTranslations ()
{
}

ComplexPEPSTransfertMatrixPBCWithTranslations::ComplexPEPSTransfertMatrixPBCWithTranslations(MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture): ComplexPEPSTransfertMatrixPBC (tensorElementsFile,architecture)
{
}

ComplexPEPSTransfertMatrixPBCWithTranslations::~ComplexPEPSTransfertMatrixPBCWithTranslations()
{
}

void ComplexPEPSTransfertMatrixPBCWithTranslations::SetHilbertSpace(AbstractHilbertSpace * hilbertSpace)
{
  this->HilbertSpace = (AbstractSpinChain *) hilbertSpace;
  this->XMomentum= ((AbstractDoubledSpinChainWithTranslations *) this->HilbertSpace)->GetMomentum();
  
  if (  this->ChainLength != this->HilbertSpace->GetSpinChainLength() )
    {
      delete this->TmpVector1; 
      delete this->TmpVector2;
      delete this->EndVector;
      delete this->StartVector;
      this->ChainLength = this->HilbertSpace->GetSpinChainLength();
      this->MaxXMomentum= this->HilbertSpace->GetSpinChainLength();
      
      delete [] this->PowerD;
      
      this->PowerD = new int[this->ChainLength+2];
      this->PowerD[0] = 1;
      for(int i =1; i <=this->ChainLength+1; i++)
	this->PowerD[i] = this->PowerD[i-1] * (this->MPOBondDimension * this->MPOBondDimension);
      
      this->TmpVector1 = new ComplexVector (this->PowerD[this->ChainLength+1],true);
      this->TmpVector2 = new ComplexVector (this->PowerD[this->ChainLength+1],true);
      this->EndVector =  new ComplexVector (this->PowerD[this->ChainLength],true);
      this->StartVector =  new ComplexVector (this->PowerD[this->ChainLength],true);
    }
  unsigned long MemoryCost =  (2*this->PowerD[this->ChainLength+1] + 2*this->PowerD[this->ChainLength])*sizeof(Complex);
  cout <<"Memory Cost " <<MemoryCost<<endl;
}


ComplexVector& ComplexPEPSTransfertMatrixPBCWithTranslations::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  ComplexVector * TmpPointorVector;
  this->EndVector->ClearVector();
  ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->ConvertToGeneralSpaceWithMomentum(vSource,(*this->StartVector));
  for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
    {	  
      this->LowLevelAddMultiplyOnFirstSite(IndiceTop);
      for(int Position = 1; Position <this->ChainLength - 1;Position++)
	{
	  this->LowLevelAddMultiplyOnAnySite(Position);
	  TmpPointorVector = this->TmpVector1;
	  this->TmpVector1 = this->TmpVector2;
	  this->TmpVector2 = TmpPointorVector;
	  this->TmpVector2->ClearVector();
	}
      this->LowLevelAddMultiplyOnLastSite(IndiceTop);
      this->TmpVector1->ClearVector();
    }
  ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->AddConvertFromGeneralSpaceWithMomentum((*this->EndVector),vDestination);
  return vDestination;
}
