#include "ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations.h"

#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;

ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations ()
{
}

ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations(MultiColumnASCIIFile & tensorElementsFileA, MultiColumnASCIIFile & tensorElementsFileB, AbstractArchitecture * architecture)
{
  this->NbrNonZeroTensorElementTopLeft = new int ** [2];
  this->IndiceBottomNonZeroTensorElementTopLeft = new int *** [2];
  this->IndiceRightNonZeroTensorElementTopLeft = new int *** [2];
  this->ValuesNonZeroTensorElementTopLeft = new Complex *** [2];
  this->InitializeTensorsElements(tensorElementsFileA,0);
  this->InitializeTensorsElements(tensorElementsFileB,1);
  this->Architecture = architecture;
  this->TemporaryArray = 0;
  this->TmpVector1 = 0;
  this->TmpVector2 = 0;
  this->StartVector = 0;
  this->EndVector = 0;
  this->ChainLength = 0;
  this->BoundaryMatrix = 0; 
}

ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::~ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations()
{
  delete this->BoundaryMatrix;
  for(int SublatticeIndex = 0; SublatticeIndex <2;  SublatticeIndex++)
    {
      for (int i = 0; i < this->MPOBondDimension *  this->MPOBondDimension; i++)
	{
	  for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
	    {
	      delete [] this->IndiceBottomNonZeroTensorElementTopLeft[SublatticeIndex][i][j];
	      delete [] this->IndiceRightNonZeroTensorElementTopLeft[SublatticeIndex][i][j];
	      delete [] this->ValuesNonZeroTensorElementTopLeft[SublatticeIndex][i][j];
	    }
	  delete [] this->IndiceBottomNonZeroTensorElementTopLeft[SublatticeIndex][i];
	  delete [] this->IndiceRightNonZeroTensorElementTopLeft[SublatticeIndex][i];
	  delete [] this->ValuesNonZeroTensorElementTopLeft[SublatticeIndex][i];
	  delete [] this->NbrNonZeroTensorElementTopLeft[SublatticeIndex][i];
	}
      delete [] this->IndiceBottomNonZeroTensorElementTopLeft[SublatticeIndex];
      delete [] this->IndiceRightNonZeroTensorElementTopLeft[SublatticeIndex];
      delete [] this->ValuesNonZeroTensorElementTopLeft[SublatticeIndex];
      delete [] this->NbrNonZeroTensorElementTopLeft[SublatticeIndex];
      
    }
  delete [] this->IndiceBottomNonZeroTensorElementTopLeft;
  delete [] this->IndiceRightNonZeroTensorElementTopLeft;
  delete [] this->ValuesNonZeroTensorElementTopLeft;
  delete [] this->NbrNonZeroTensorElementTopLeft;
  delete this->TmpVector1; delete  this->TmpVector2; delete  this->StartVector; delete this->EndVector;
}

void ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::SetHilbertSpace(AbstractHilbertSpace * hilbertSpace)
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
      this->MaxXMomentum= this->HilbertSpace->GetSpinChainLength()/2;   
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

ComplexVector& ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::LowLevelAddMultiply(ComplexVector & vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  ComplexVector * TmpPointorVector;
  this->EndVector->ClearVector();
  this->StartVector->ClearVector();
  ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->ConvertToGeneralSpaceWithMomentum(vSource,(*this->StartVector));

  for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
    {	 
      this->LowLevelAddMultiplyOnFirstSite(IndiceTop,0);
      for(int Position = 1; Position <this->ChainLength - 1;Position++)
	{
	  this->LowLevelAddMultiplyOnAnySite(Position,Position%2);
	  TmpPointorVector = this->TmpVector1;
	  this->TmpVector1 = this->TmpVector2;
	  this->TmpVector2 = TmpPointorVector;
	  this->TmpVector2->ClearVector();
	}
      this->LowLevelAddMultiplyOnLastSite(IndiceTop,1);
      this->TmpVector1->ClearVector();
    }
  this->StartVector->Copy(*(this->EndVector));
  this->EndVector->ClearVector();
  for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
    {	 
      this->LowLevelAddMultiplyOnFirstSite(IndiceTop,1);
      for(int Position = 1; Position <this->ChainLength - 1;Position++)
	{
	  this->LowLevelAddMultiplyOnAnySite(Position,(Position+1)%2);
	  TmpPointorVector = this->TmpVector1;
	  this->TmpVector1 = this->TmpVector2;
	  this->TmpVector2 = TmpPointorVector;
	  this->TmpVector2->ClearVector();
	}
      this->LowLevelAddMultiplyOnLastSite(IndiceTop,0);
      this->TmpVector1->ClearVector();
    }
  
  ((AbstractDoubledSpinChainWithTranslations * )this->HilbertSpace)->AddConvertFromGeneralSpaceWithMomentum((*this->EndVector),vDestination);
  return vDestination;
}


void ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile, int sublatticeIndex)
{
  this->NbrNonZeroElements = tensorElementsFile.GetNbrLines();
  int* IndexVertical  = tensorElementsFile.GetAsIntegerArray (0);
  int* IndexLeft = tensorElementsFile.GetAsIntegerArray (1);
  int* IndexUp = tensorElementsFile.GetAsIntegerArray (2);
  int* IndexRight = tensorElementsFile.GetAsIntegerArray (3);
  int* IndexDown  = tensorElementsFile.GetAsIntegerArray (4);
  Complex * ElementsValues = tensorElementsFile.GetAsComplexArray (5);

  int TmpPhysicalDimension = 0;
  int TmpMPODimension = 0;
  int VerticalDimension = 0;
  unsigned long MemoryCost = 0;  
  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      if (IndexVertical[i] >  VerticalDimension)
	{
	  VerticalDimension = IndexVertical[i];
	}
      if (IndexDown[i] > TmpPhysicalDimension)
	TmpPhysicalDimension = IndexDown[i];
      if (IndexLeft[i] > TmpMPODimension)
	TmpMPODimension = IndexLeft[i];
    }
  
  this->PhysicalDimension = TmpPhysicalDimension+1;
  this->MPOBondDimension =  TmpMPODimension+1;
  VerticalDimension++;
  
  this->NbrNonZeroTensorElementTopLeft[sublatticeIndex] = new int * [this->MPOBondDimension*this->MPOBondDimension];
  this->IndiceBottomNonZeroTensorElementTopLeft[sublatticeIndex] = new int ** [this->MPOBondDimension*this->MPOBondDimension];
  this->IndiceRightNonZeroTensorElementTopLeft[sublatticeIndex] = new int ** [this->MPOBondDimension*this->MPOBondDimension];
  this->ValuesNonZeroTensorElementTopLeft[sublatticeIndex] = new Complex ** [this->MPOBondDimension*this->MPOBondDimension];
  
  MemoryCost+= (sizeof(int*) + sizeof(Complex**) + 2*sizeof(int**))*this->MPOBondDimension*this->MPOBondDimension;

  for (int i = 0; i < this->MPOBondDimension * this->MPOBondDimension; i++)
    {
      this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i] = new int [this->PhysicalDimension*this->PhysicalDimension];
      this->IndiceBottomNonZeroTensorElementTopLeft[sublatticeIndex][i] = new int * [this->PhysicalDimension*this->PhysicalDimension];
      this->IndiceRightNonZeroTensorElementTopLeft[sublatticeIndex][i] = new int * [this->PhysicalDimension*this->PhysicalDimension];
      this->ValuesNonZeroTensorElementTopLeft[sublatticeIndex][i] = new Complex * [this->PhysicalDimension*this->PhysicalDimension];
      MemoryCost+=this->PhysicalDimension*this->PhysicalDimension*(2*sizeof(int*)+sizeof(double*)+sizeof(int));
      for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
	{
	  this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i][j] = 0;
	}
    }

  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      for(int j = 0 ; j < this->NbrNonZeroElements; j++)
	{
	  if (IndexVertical[i]==IndexVertical[j]) 
	    this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][GetCommonIndexFromBraAndKetIndices(IndexUp[i],IndexUp[j])][GetCommonIndexFromBraAndKetIndices(IndexLeft[i],IndexLeft[j])]++;
	} 
    }
  
  for (int i = 0; i < this->MPOBondDimension * this->MPOBondDimension; i++)
    {
      for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
	{
	  this->IndiceBottomNonZeroTensorElementTopLeft[sublatticeIndex][i][j] = new int [this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i][j]];
	  this->IndiceRightNonZeroTensorElementTopLeft[sublatticeIndex][i][j] =  new int [this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i][j]];
	  this->ValuesNonZeroTensorElementTopLeft[sublatticeIndex][i][j] = new Complex [this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i][j]];
	  MemoryCost+= (2*sizeof(int)+ sizeof(Complex))*this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i][j];
	  this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i][j]=0;
	}
    }
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      for(int j = 0 ; j < this->NbrNonZeroElements; j++)
	{
	  if (IndexVertical[i]==IndexVertical[j]) 
	    {
	      int IndiceUp = GetCommonIndexFromBraAndKetIndices(IndexUp[i],IndexUp[j]);
	      int IndiceLeft = GetCommonIndexFromBraAndKetIndices(IndexLeft[i],IndexLeft[j]);
	      this->IndiceBottomNonZeroTensorElementTopLeft[sublatticeIndex][IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][IndiceUp][IndiceLeft]] = GetCommonIndexFromBraAndKetIndices(IndexDown[i],IndexDown[j]);
	      this->IndiceRightNonZeroTensorElementTopLeft[sublatticeIndex][IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][IndiceUp][IndiceLeft]] = GetCommonIndexFromBraAndKetIndices(IndexRight[i],IndexRight[j]);
	      this->ValuesNonZeroTensorElementTopLeft[sublatticeIndex][IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][IndiceUp][IndiceLeft]] =  ElementsValues[i]*Conj(ElementsValues[j]);
	      this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][IndiceUp][IndiceLeft]++;
	    }
	}
    }
  
  cout <<" Physical Dimension = " <<  this->PhysicalDimension<<endl;
  cout <<" MPO Dimension = " <<  this->MPOBondDimension <<endl;
  cout <<"Memory Cost = "<<MemoryCost <<endl;
  delete [] IndexVertical;
  delete [] IndexDown;
  delete [] IndexUp;
  delete [] IndexLeft;
  delete [] IndexRight;
  delete [] ElementsValues;
}

void ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::LowLevelAddMultiplyOnFirstSite(int topIndice, int sublatticeIndex)
{
  for(int i = 0; i< this->StartVector->GetVectorDimension();i++) 
    {
      if ( Norm((*this->StartVector)[i]) > 1e-13 )
	{
	  for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][topIndice][i%this->PowerD[1]]; p++)
	    {
	      (*this->TmpVector1)[this->GetNewIndexFromOldIndex(this->PowerD[1]*i,i%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[sublatticeIndex][topIndice][i%this->PowerD[1]][p],0,this->IndiceBottomNonZeroTensorElementTopLeft[sublatticeIndex][topIndice][i%this->PowerD[1]][p],1)]+=this->ValuesNonZeroTensorElementTopLeft[sublatticeIndex][topIndice][i%this->PowerD[1]][p]*  (*this->StartVector)[i];
	    }
	}
    }
}


void ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::LowLevelAddMultiplyOnAnySite(int position, int sublatticeIndex)
{
  for(int i = 0; i< this->TmpVector1->GetVectorDimension();i++) 
    {
      if ( Norm((*this->TmpVector1)[i]) > 1e-13 )
	{
	  for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]] ; p++)
	    {
	      (*this->TmpVector2)[this->GetNewIndexFromOldIndex(i, (i/this->PowerD[position+1])%this->PowerD[1], this->IndiceRightNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]][p], i%this->PowerD[1], this->IndiceBottomNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]][p],position+1)] += this->ValuesNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][(i/this->PowerD[position+1])%this->PowerD[1]][p] *  (*this->TmpVector1)[i];
	    }
	}
    }
}
  

void ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations::LowLevelAddMultiplyOnLastSite(int topValue, int sublatticeIndex)
{
  for(int i = 0; i< this->TmpVector1->GetVectorDimension();i++) 
    {
      if ( Norm((*this->TmpVector1)[i]) > 1e-13 )
	{
	  for (int p = 0; p < this->NbrNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][i/this->PowerD[this->ChainLength]] ; p++)
	    {
	      if (this->IndiceBottomNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][i/this->PowerD[this->ChainLength]][p] == topValue)
		{
		  (*this->EndVector)[this->GetNewIndexFromOldIndex(i,i/this->PowerD[this->ChainLength], this->IndiceRightNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][i/this->PowerD[this->ChainLength]][p],i%this->PowerD[1],0,this->ChainLength)/this->PowerD[1]] +=  (*this->TmpVector1)[i] * this->ValuesNonZeroTensorElementTopLeft[sublatticeIndex][i%this->PowerD[1]][i/this->PowerD[this->ChainLength]][p];
		  
		}
	    }
	}
    }
}


