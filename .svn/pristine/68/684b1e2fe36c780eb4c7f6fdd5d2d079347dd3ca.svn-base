#include "TransfertMatrixPBCWithTranslationsFromFile.h"
#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

TransfertMatrixPBCWithTranslationsFromFile::TransfertMatrixPBCWithTranslationsFromFile()
{
  this->ExponentialFactors=0;
}



TransfertMatrixPBCWithTranslationsFromFile::TransfertMatrixPBCWithTranslationsFromFile (MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture): AbstractTransfertMatrixPBC(tensorElementsFile,architecture)
{
  this->ExponentialFactors=0;
}

TransfertMatrixPBCWithTranslationsFromFile::~TransfertMatrixPBCWithTranslationsFromFile()
{
  delete[] this->ExponentialFactors;
}

RealVector& TransfertMatrixPBCWithTranslationsFromFile::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  cout <<"Using wrong function RealVector& TransfertMatrixPBCWithTranslationsFromFile::LowLevelMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)"<<endl;
  return vDestination;
}


ComplexVector& TransfertMatrixPBCWithTranslationsFromFile::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector & vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent +  nbrComponent;
  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
      this->HilbertSpace->GetBosonicOccupation(i,this->TemporaryArray);
      
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension;  IndiceTop++)
	{	  
	  Dim += this->EvaluateNbrResultingState(IndiceTop,this->ChainLength-1,IndiceTop);
	}

      double * ResultingCoefficient = new double [Dim];
      unsigned long * ResultingIndex = new unsigned long[Dim];
      unsigned long Tmp=0;
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension;  IndiceTop++)
	{	  
	  Tmp=this->GenerateResultingStateAndCoefficient(IndiceTop,this->ChainLength-1,IndiceTop,ResultingCoefficient,ResultingIndex,Tmp);
	}
      int NbrTranslation;
      for (int p = 0; p < Dim; p++)
	{

	  int DestinationState = ((AbstractSpinChainWithTranslations *) this->HilbertSpace)->FindStateIndex(((AbstractSpinChainWithTranslations *) this->HilbertSpace)->FindCanonicalForm(ResultingIndex[p], NbrTranslation));
	  if (DestinationState < this->HilbertSpace->GetHilbertSpaceDimension() ) 
	    {
	      vDestination[DestinationState] += ResultingCoefficient[p]*vSource[i] * this->ExponentialFactors[NbrTranslation] * ((AbstractSpinChainWithTranslations *) this->HilbertSpace)->GetRescalingFactor(i,DestinationState);
	    }
	}
     delete [] ResultingCoefficient; delete [] ResultingIndex;
    }
  return vDestination;
}

// evaluate all exponential factors
//   

void TransfertMatrixPBCWithTranslationsFromFile::EvaluateExponentialFactors()
{
  this->ExponentialFactors = new Complex[this->MaxXMomentum];
  for (int i = 0; i < this->MaxXMomentum; ++i)
    { 
      this->ExponentialFactors[i] = Phase(2.0 * M_PI * ((this->XMomentum * ((double) i) / ((double) this->MaxXMomentum))));
    }
}

void TransfertMatrixPBCWithTranslationsFromFile::SetHilbertSpace(AbstractHilbertSpace * hilbertSpace)
{
  this->HilbertSpace = (AbstractSpinChain *) hilbertSpace;
  this->XMomentum= ((AbstractSpinChainWithTranslations *) this->HilbertSpace)->GetMomentum();
  this->ChainLength = this->HilbertSpace->GetSpinChainLength();
  this->MaxXMomentum= this->ChainLength;
  delete [] this->ExponentialFactors;
  this->EvaluateExponentialFactors();
}
