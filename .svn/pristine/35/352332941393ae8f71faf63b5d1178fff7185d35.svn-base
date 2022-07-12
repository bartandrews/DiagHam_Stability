#include "MPSObjects/AbstractPEPSTransfertMatrixPBC.h"
#include "HilbertSpace/AbstractDoubledSpinChain.h"
#include <iostream>

using std::cout;
using std::endl;

AbstractPEPSTransfertMatrixPBC::AbstractPEPSTransfertMatrixPBC ()
{
}

AbstractPEPSTransfertMatrixPBC::AbstractPEPSTransfertMatrixPBC(MultiColumnASCIIFile & tensorElementsFile,AbstractArchitecture * architecture)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->Architecture = architecture;
  this->TemporaryArray = 0;
}


AbstractPEPSTransfertMatrixPBC::~AbstractPEPSTransfertMatrixPBC()
{
}




RealVector& AbstractPEPSTransfertMatrixPBC::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent +  nbrComponent;
  int * IndiceLeftBra = new int[this->ChainLength];
  int * IndiceLeftKet = new int[this->ChainLength];
  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
//      ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetBosonicOccupation(i,IndiceLeftBra,IndiceLeftKet);
      for (int k =0; k < this->ChainLength;k++)
	{
	  this->TemporaryArray[k]= this->GetCommonIndexFromBraAndKetIndices(IndiceLeftBra[k],IndiceLeftKet[k]);
	}
      
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
	{	  
	  Dim += this->EvaluateNbrResultingState(IndiceTop,this->ChainLength-1,IndiceTop);
	}

      double * ResultingCoefficient = new double [Dim];
      unsigned long * ResultingIndexBra = new unsigned long[Dim];
      unsigned long * ResultingIndexKet = new unsigned long[Dim];
      unsigned long Tmp=0;
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension ;  IndiceTop++)
	{	  
	  Tmp=this->GenerateResultingStateAndCoefficient(IndiceTop,this->ChainLength-1,IndiceTop,ResultingCoefficient,ResultingIndexBra,ResultingIndexKet,Tmp);
	}
      for (int p = 0; p < Dim; p++)
	{
	  vDestination[((AbstractDoubledSpinChain * )this->HilbertSpace)->FindStateIndex(ResultingIndexBra[p],ResultingIndexKet[p])] += ResultingCoefficient[p]*vSource[i];
	}
      delete [] ResultingCoefficient; delete [] ResultingIndexBra;  delete [] ResultingIndexKet; 
    }
  delete [] IndiceLeftBra;
  delete [] IndiceLeftKet;
  return vDestination;
}


void AbstractPEPSTransfertMatrixPBC::InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile)
{
  this->NbrNonZeroElements = tensorElementsFile.GetNbrLines();
  int* IndexVertical  = tensorElementsFile.GetAsIntegerArray (0);
  int* IndexLeft = tensorElementsFile.GetAsIntegerArray (1);
  int* IndexUp = tensorElementsFile.GetAsIntegerArray (2);
  int* IndexRight = tensorElementsFile.GetAsIntegerArray (3);
  int* IndexDown  = tensorElementsFile.GetAsIntegerArray (4);
  double * ElementsValues = tensorElementsFile.GetAsDoubleArray (5);

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
  
  this->NbrNonZeroTensorElementTopLeft = new int * [this->MPOBondDimension*this->MPOBondDimension];
  this->IndiceBottomNonZeroTensorElementTopLeft = new int ** [this->MPOBondDimension*this->MPOBondDimension];
  this->IndiceRightNonZeroTensorElementTopLeft = new int ** [this->MPOBondDimension*this->MPOBondDimension];
  this->ValuesNonZeroTensorElementTopLeft = new double ** [this->MPOBondDimension*this->MPOBondDimension];

  MemoryCost+= (sizeof(int*) + sizeof(double**) + 2*sizeof(int**))*this->MPOBondDimension*this->MPOBondDimension;

  for (int i = 0; i < this->MPOBondDimension * this->MPOBondDimension; i++)
    {
      this->NbrNonZeroTensorElementTopLeft[i] = new int [this->PhysicalDimension*this->PhysicalDimension];
      this->IndiceBottomNonZeroTensorElementTopLeft[i] = new int * [this->PhysicalDimension*this->PhysicalDimension];
      this->IndiceRightNonZeroTensorElementTopLeft[i] = new int * [this->PhysicalDimension*this->PhysicalDimension];
      this->ValuesNonZeroTensorElementTopLeft[i] = new double * [this->PhysicalDimension*this->PhysicalDimension];
      MemoryCost+=this->PhysicalDimension*this->PhysicalDimension*(2*sizeof(int*)+sizeof(double*)+sizeof(int));
      for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
	{
	  this->NbrNonZeroTensorElementTopLeft[i][j] = 0;
	}
    }

  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      for(int j = 0 ; j < this->NbrNonZeroElements; j++)
	{
	  if (IndexVertical[i]==IndexVertical[j]) 
	    this->NbrNonZeroTensorElementTopLeft[GetCommonIndexFromBraAndKetIndices(IndexUp[i],IndexUp[j])][GetCommonIndexFromBraAndKetIndices(IndexLeft[i],IndexLeft[j])]++;
	} 
    }
  
  for (int i = 0; i < this->MPOBondDimension * this->MPOBondDimension; i++)
    {
      for(int j = 0; j < this->PhysicalDimension * this->PhysicalDimension; j++)
	{
	  this->IndiceBottomNonZeroTensorElementTopLeft[i][j] = new int [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  this->IndiceRightNonZeroTensorElementTopLeft[i][j] =  new int [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  this->ValuesNonZeroTensorElementTopLeft[i][j] = new double [this->NbrNonZeroTensorElementTopLeft[i][j]];
	  MemoryCost+= (2*sizeof(int)+ sizeof(double))*this->NbrNonZeroTensorElementTopLeft[i][j];
	  this->NbrNonZeroTensorElementTopLeft[i][j]=0;
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
	      this->IndiceBottomNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]] = GetCommonIndexFromBraAndKetIndices(IndexDown[i],IndexDown[j]);
	      this->IndiceRightNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]] = GetCommonIndexFromBraAndKetIndices(IndexRight[i],IndexRight[j]);
	      this->ValuesNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft][this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]] =  ElementsValues[i]*ElementsValues[j];
	      this->NbrNonZeroTensorElementTopLeft[IndiceUp][IndiceLeft]++;
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


int AbstractPEPSTransfertMatrixPBC::GenerateResultingStateAndCoefficient(int indiceTop, int chainSize, int lastIndice, double * coefArray, unsigned long * stateArrayBra, unsigned long * stateArrayKet, unsigned long pos)
{
  unsigned int IndexBra, IndexKet;
  if(chainSize==0)
    {
      for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]]; i++)
	{
	  if(this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]][i] ==  lastIndice ) 
	    {
	      this->GetBraAndKetIndexFromCommonIndex(this->IndiceRightNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]][i], IndexBra, IndexKet);
//	      stateArrayBra[pos] = ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateBra(IndexBra,0);
	//      stateArrayKet[pos] = ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateKet(IndexKet,0);
	      coefArray[pos]  = this->ValuesNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[0]][i];
	      pos++;
	    }
	}
      return pos;
    }
  int TmpPos = pos;
  
  for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]]; i++)
    {
      TmpPos =this->GenerateResultingStateAndCoefficient(this->IndiceBottomNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]][i],chainSize-1, lastIndice, coefArray, stateArrayBra, stateArrayKet, pos);
      for(;pos <TmpPos;++pos)
	{
	  this->GetBraAndKetIndexFromCommonIndex(this->IndiceRightNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]][i], IndexBra, IndexKet);
//	  stateArrayBra[pos] |= ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateBra(IndexBra,chainSize);
//	  stateArrayKet[pos] |= ((AbstractDoubledSpinChain * )this->HilbertSpace)->EncodeSiteStateKet(IndexKet,chainSize);
//	  coefArray[pos] *= this->ValuesNonZeroTensorElementTopLeft[indiceTop][this->TemporaryArray[chainSize]][i];
	}
    }
  return pos;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector * AbstractPEPSTransfertMatrixPBC::LowLevelMultipleAddMultiply(RealVector * vSources, RealVector* vDestinations, int nbrVectors, int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent +  nbrComponent;
  int * IndiceLeftBra = new int[this->ChainLength];
  int * IndiceLeftKet = new int[this->ChainLength];
  for (int i = firstComponent ; i  <  LastComponent; i++)
    {
      int Dim = 0;
//    ((AbstractDoubledSpinChain * )this->HilbertSpace)->GetBosonicOccupation(i,IndiceLeftBra,IndiceLeftKet);
      for (int k =0; k < this->ChainLength;k++)
	{
	  this->TemporaryArray[k]= this->GetCommonIndexFromBraAndKetIndices(IndiceLeftBra[k],IndiceLeftKet[k]);
	}
      
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension;  IndiceTop++)
	{	  
	  Dim += this->EvaluateNbrResultingState(IndiceTop,this->ChainLength-1,IndiceTop);
	}

      double * ResultingCoefficient = new double [Dim];
      unsigned long * ResultingIndexBra = new unsigned long[Dim];
      unsigned long * ResultingIndexKet = new unsigned long[Dim];
      unsigned long Tmp=0;
      for (int  IndiceTop = 0 ;  IndiceTop < this->MPOBondDimension *  this->MPOBondDimension ;  IndiceTop++)
	{	  
	  Tmp=this->GenerateResultingStateAndCoefficient(IndiceTop,this->ChainLength-1,IndiceTop,ResultingCoefficient,ResultingIndexBra,ResultingIndexKet,Tmp);
	}
      for (int p = 0; p < Dim; p++)
	{
	  int Index = ((AbstractDoubledSpinChain * )this->HilbertSpace)->FindStateIndex(ResultingIndexBra[p],ResultingIndexKet[p]);
	  for(int t=0 ; t < nbrVectors; t++)
	    vDestinations[t][Index] += ResultingCoefficient[p]*vSources[t][i];
	}
      delete [] ResultingCoefficient; delete [] ResultingIndexBra;  delete [] ResultingIndexKet; 
    }
  delete [] IndiceLeftBra;
  delete [] IndiceLeftKet;
  return vDestinations;
}


void  AbstractPEPSTransfertMatrixPBC::PrintTensorElements()
{
  cout <<"#Tensor IndiceLeft IndiceTop  IndicexRight IndiceBottom Values" <<endl;
  for(int IndiceLeft=0; IndiceLeft < this->PhysicalDimension * this->PhysicalDimension; IndiceLeft++)
    {
      for(int IndiceTop =0; IndiceTop < this->MPOBondDimension *this->MPOBondDimension ; IndiceTop++)
	{
	  for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft]; i++)
	    {
	      cout <<IndiceLeft<< " "<< IndiceTop  <<" "<< this->IndiceRightNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft][i]<< " "<< this->IndiceBottomNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft][i]<<" "<<this->ValuesNonZeroTensorElementTopLeft[IndiceTop][IndiceLeft][i]<<endl;
	    }
	}
    }
}
