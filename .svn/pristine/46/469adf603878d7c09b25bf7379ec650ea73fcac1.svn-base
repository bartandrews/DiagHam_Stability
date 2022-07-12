#include "ComplexPEPSPBC.h"

#include <iostream>
#include <sys/time.h>
#include "GeneralTools/ArrayTools.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainFull.h"

#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"


using std::cout;
using std::endl;

ComplexPEPSPBC::ComplexPEPSPBC ()
{
}


ComplexPEPSPBC::ComplexPEPSPBC(MultiColumnASCIIFile & tensorElementsFile, AbstractArchitecture * architecture)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->Architecture = architecture;
}



ComplexPEPSPBC::~ComplexPEPSPBC()
{
  for(int k = 0; k < this->PhysicalDimension; k++)
    {
      for (int i = 0; i < this->MPOBondDimension ; i++)
	{
	  for(int j = 0; j < this->MPOBondDimension ; j++)
	    {
	      delete [] this->IndiceBottomNonZeroTensorElementTopLeft[k][i][j];
	      delete [] this->IndiceRightNonZeroTensorElementTopLeft[k][i][j];
	      delete [] this->ValuesNonZeroTensorElementTopLeft[k][i][j];
	    }
	  delete [] this->IndiceBottomNonZeroTensorElementTopLeft[k][i];
	  delete [] this->IndiceRightNonZeroTensorElementTopLeft[k][i];
	  delete [] this->ValuesNonZeroTensorElementTopLeft[k][i];
	  delete [] this->NbrNonZeroTensorElementTopLeft[k][i];
	}
      delete [] this->IndiceBottomNonZeroTensorElementTopLeft[k];
      delete [] this->IndiceRightNonZeroTensorElementTopLeft[k];
      delete [] this->ValuesNonZeroTensorElementTopLeft[k];
      delete [] this->NbrNonZeroTensorElementTopLeft[k];
    }
}


void ComplexPEPSPBC::InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile)
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
  unsigned long MemoryCost = 0;  
  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      if (IndexVertical[i] >  TmpPhysicalDimension)
	{
	  TmpPhysicalDimension = IndexVertical[i];
	}
      if (IndexLeft[i] > TmpMPODimension)
	TmpMPODimension = IndexLeft[i];
    }
  
  this->PhysicalDimension = TmpPhysicalDimension+1;
  this->MPOBondDimension =  TmpMPODimension+1;
  
  this->NbrNonZeroTensorElementTopLeft = new int ** [this->PhysicalDimension];
  this->IndiceBottomNonZeroTensorElementTopLeft = new int *** [this->PhysicalDimension];
  this->IndiceRightNonZeroTensorElementTopLeft = new int *** [this->PhysicalDimension];
  this->ValuesNonZeroTensorElementTopLeft = new Complex *** [this->PhysicalDimension];
  
//  MemoryCost+= (sizeof(int*) + sizeof(Complex**) + 2*sizeof(int**))*this->MPOBondDimension*this->MPOBondDimension;

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      this->NbrNonZeroTensorElementTopLeft[i] = new int * [this->MPOBondDimension];
      this->IndiceBottomNonZeroTensorElementTopLeft[i] = new int ** [this->MPOBondDimension];
      this->IndiceRightNonZeroTensorElementTopLeft[i] = new int ** [this->MPOBondDimension];
      this->ValuesNonZeroTensorElementTopLeft[i] = new Complex ** [this->MPOBondDimension];
      
      for (int p = 0; p < this->MPOBondDimension ; p++)
	{
	  this->NbrNonZeroTensorElementTopLeft[i][p] = new int [this->MPOBondDimension];
	  this->IndiceBottomNonZeroTensorElementTopLeft[i][p] = new int * [this->MPOBondDimension];
	  this->IndiceRightNonZeroTensorElementTopLeft[i][p] = new int * [this->MPOBondDimension];
	  this->ValuesNonZeroTensorElementTopLeft[i][p] = new Complex * [this->MPOBondDimension];
	  //      MemoryCost+=this->PhysicalDimension*this->PhysicalDimension*(2*sizeof(int*)+sizeof(double*)+sizeof(int));
	  for(int j = 0; j <  this->MPOBondDimension; j++)
	    {
	      this->NbrNonZeroTensorElementTopLeft[i][p][j] = 0;
	    }
	}
    }
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      this->NbrNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]]++;
    }
  
  for (int k = 0; k <  this->PhysicalDimension; k++)
    {  
      for (int i = 0; i < this->MPOBondDimension; i++)
	{
	  for(int j = 0; j < this->MPOBondDimension; j++)
	    {
	      this->IndiceBottomNonZeroTensorElementTopLeft[k][i][j] = new int [this->NbrNonZeroTensorElementTopLeft[k][i][j]];
	      this->IndiceRightNonZeroTensorElementTopLeft[k][i][j] =  new int [this->NbrNonZeroTensorElementTopLeft[k][i][j]];
	      this->ValuesNonZeroTensorElementTopLeft[k][i][j] = new Complex [this->NbrNonZeroTensorElementTopLeft[k][i][j]];
//	      MemoryCost+= (2*sizeof(int)+ sizeof(Complex))*this->NbrNonZeroTensorElementTopLeft[i][j];
	      this->NbrNonZeroTensorElementTopLeft[k][i][j]=0;
	    }
	}
    }
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      this->IndiceBottomNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]][this->NbrNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]]] = IndexDown[i];
      this->IndiceRightNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]][this->NbrNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]]] = IndexRight[i];
      this->ValuesNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]][this->NbrNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]]] =  ElementsValues[i];
      this->NbrNonZeroTensorElementTopLeft[IndexVertical[i]][IndexUp[i]][IndexLeft[i]]++;
    }
  
  cout <<" Physical Dimension = " <<  this->PhysicalDimension<<endl;
  cout <<" MPO Dimension = " <<  this->MPOBondDimension <<endl;
//  cout <<"Memory Cost = "<<MemoryCost <<endl;
  delete [] IndexVertical;
  delete [] IndexDown;
  delete [] IndexUp;
  delete [] IndexLeft;
  delete [] IndexRight;
  delete [] ElementsValues;
}

void ComplexPEPSPBC::PrintTensorElements()
{
  cout <<"#Tensor PhysicalIndice IndiceLeft IndiceTop  IndicexRight IndiceBottom Values" <<endl;
  for(int IndicePhysical=0; IndicePhysical < this->PhysicalDimension; IndicePhysical++)
    {
      for(int IndiceLeft=0; IndiceLeft <  this->MPOBondDimension ; IndiceLeft++)
	{
	  for(int IndiceTop =0; IndiceTop < this->MPOBondDimension ; IndiceTop++)
	    {
	      for (int i = 0; i < this->NbrNonZeroTensorElementTopLeft[IndicePhysical][IndiceTop][IndiceLeft]; i++)
		{
		  cout <<IndicePhysical<<" "<<IndiceLeft<< " "<< IndiceTop  <<" "<< this->IndiceRightNonZeroTensorElementTopLeft[IndicePhysical][IndiceTop][IndiceLeft][i]<< " "<< this->IndiceBottomNonZeroTensorElementTopLeft[IndicePhysical][IndiceTop][IndiceLeft][i]<<" "<<this->ValuesNonZeroTensorElementTopLeft[IndicePhysical][IndiceTop][IndiceLeft][i]<<endl;
		}
	    }
	}
    }
}

ComplexVector ComplexPEPSPBC::ComputeFockSpaceRepresentationOfAPEPS (int lx, int lylogtwo, ComplexDiagonalMatrix virtualSymmetryHorizontal, ComplexDiagonalMatrix virtualSymmetryVertical, bool horizontalFlag, bool verticalFlag)
{
  int Ly = 1;
  for(int i = 0 ; i <lylogtwo ; i++)
    {
      Ly*=2;
    }
  int Lx = lx;
  
  Spin1_2ChainFull Space (Lx*Ly); 
  ComplexVector PEPS (Space.GetHilbertSpaceDimension(), true);  

  int NbrColumnMatrix=this->PhysicalDimension*this->PhysicalDimension; 
  int DimColumMatrix = this->MPOBondDimension*  this->MPOBondDimension;
  
  ComplexMatrix *** ColumnMatrix = new ComplexMatrix ** [NbrColumnMatrix];
  for(int i = 0; i < NbrColumnMatrix ; i++)
    {
      ColumnMatrix[i] = new ComplexMatrix * [this->MPOBondDimension];
      for(int j = 0 ; j < this->MPOBondDimension ; j++)
	{
	  ColumnMatrix[i][j] = new ComplexMatrix[this->MPOBondDimension];
	  for(int k = 0 ; k < this->MPOBondDimension ; k++)
	    {
	      ColumnMatrix[i][j][k] = ComplexMatrix(DimColumMatrix,DimColumMatrix,true) ;
	    }
	}
    }
  
  for(int i = 0; i < NbrColumnMatrix; i++)
    {
      for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumMatrix;  NewIndiceLeft++)
	{
	  for(int TopIndice = 0 ;  TopIndice <   this->MPOBondDimension ;  TopIndice++)
	    {
	      for (int NonZeroElementUp =0 ; NonZeroElementUp < this->NbrNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension];  NonZeroElementUp++ )
		{
		  int NewRightUpTensorIndex = this->IndiceRightNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  int NewBottomtUpTensorIndex = this->IndiceBottomNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  Complex TmpValues = this->ValuesNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft%this->MPOBondDimension][NonZeroElementUp];
		  for (int NonZeroElementDown =0 ;NonZeroElementDown <  this->NbrNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension];  NonZeroElementDown++ )
		    { 
		      ColumnMatrix[i][TopIndice][this->IndiceBottomNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]].AddToMatrixElement(NewIndiceLeft, this->IndiceRightNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]  *   this->MPOBondDimension + NewRightUpTensorIndex, this->ValuesNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension][NonZeroElementDown] * TmpValues);
		    }
		}
	    }
	}
    }

  if (Ly == 2 ) 
    {
      int NbrColumnTwoMatrix=NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix;
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      
      Complex Tmp; 
      Complex Tmp3 = 1.0; 
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      if (horizontalFlag)
			virtualSymmetryHorizontal.GetMatrixElement(TopIndice,TopIndice,Tmp3);
		      ColumnMatrix[i][TopIndice][TopIndice].GetMatrixElement(NewIndiceLeft, NewIndiceRight ,Tmp); 
		      ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp);
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();    
	    }
	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex Tmp,Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int NewIndiceLeft = 0 ; NewIndiceLeft < DimColumTwoMatrix; NewIndiceLeft++)
	    {	  
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p = 0 ; p < Ly ; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();    
	    }
	}
      return PEPS;
    }
  else
    {
      int NbrColumnTwoMatrix=NbrColumnMatrix*NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix*DimColumMatrix;
      
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i =0; i < NbrColumnTwoMatrix ;i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      Complex Tmp,Tmp2;
      Complex Tmp3 = 1.0;
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      for(int MiddleIndice = 0;  MiddleIndice <this->MPOBondDimension ;   MiddleIndice++)
			{
			  ColumnMatrix[i% NbrColumnMatrix][TopIndice][MiddleIndice].GetMatrixElement(NewIndiceLeft%DimColumMatrix, NewIndiceRight%DimColumMatrix  ,Tmp); 
			  ColumnMatrix[i/NbrColumnMatrix][MiddleIndice][TopIndice].GetMatrixElement(NewIndiceLeft/DimColumMatrix, NewIndiceRight/DimColumMatrix  ,Tmp2); 
			  if (horizontalFlag)
			    virtualSymmetryHorizontal.GetMatrixElement(MiddleIndice,MiddleIndice,Tmp3);
			  ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp*Tmp2*Tmp3);
			}
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0; i < Space.GetHilbertSpaceDimension(); i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();    
	    }
	  return  PEPS;
	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex Tmp,Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {	  
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p=0; p < Ly; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();     
	    }
	  return  PEPS;
	}
    }
}


ComplexVector ComplexPEPSPBC::ComputeFockSpaceRepresentationOfAPEPSSzConstraint (int lx, int lylogtwo, int sz, ComplexDiagonalMatrix virtualSymmetryHorizontal, ComplexDiagonalMatrix virtualSymmetryVertical, bool horizontalFlag, bool verticalFlag)
{
  int Ly = 1;
  for(int i = 0 ; i <lylogtwo ; i++)
    {
      Ly*=2;
    }
  int Lx = lx;
  
  Spin1_2ChainNew Space (Lx*Ly,sz, 100000); 
  ComplexVector PEPS (Space.GetHilbertSpaceDimension(), true);  

  int NbrColumnMatrix=this->PhysicalDimension*this->PhysicalDimension; 
  int DimColumMatrix = this->MPOBondDimension*  this->MPOBondDimension;
  
  ComplexMatrix *** ColumnMatrix = new ComplexMatrix ** [NbrColumnMatrix];
  for(int i = 0; i < NbrColumnMatrix ; i++)
    {
      ColumnMatrix[i] = new ComplexMatrix * [this->MPOBondDimension];
      for(int j = 0 ; j < this->MPOBondDimension ; j++)
	{
	  ColumnMatrix[i][j] = new ComplexMatrix[this->MPOBondDimension];
	  for(int k = 0 ; k < this->MPOBondDimension ; k++)
	    {
	      ColumnMatrix[i][j][k] = ComplexMatrix(DimColumMatrix,DimColumMatrix,true) ;
	    }
	}
    }
  
  for(int i = 0; i < NbrColumnMatrix; i++)
    {
      for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumMatrix;  NewIndiceLeft++)
	{
	  for(int TopIndice = 0 ;  TopIndice <   this->MPOBondDimension ;  TopIndice++)
	    {
	      for (int NonZeroElementUp = 0 ; NonZeroElementUp < this->NbrNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension];  NonZeroElementUp++ )
		{
		  int NewRightUpTensorIndex = this->IndiceRightNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  int NewBottomtUpTensorIndex = this->IndiceBottomNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  Complex TmpValues = this->ValuesNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft%this->MPOBondDimension][NonZeroElementUp];
		  for (int NonZeroElementDown =0 ;NonZeroElementDown <  this->NbrNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension];  NonZeroElementDown++ )
		    { 
		      ColumnMatrix[i][TopIndice][this->IndiceBottomNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]].AddToMatrixElement(NewIndiceLeft, this->IndiceRightNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]  *   this->MPOBondDimension + NewRightUpTensorIndex, this->ValuesNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension][NonZeroElementDown] * TmpValues);
		    }
		}
	    }
	}
    }

  if (Ly == 2 ) 
    {
      int NbrColumnTwoMatrix = NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix;
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      
      Complex Tmp; 
      Complex Tmp3 = 1.0; 
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      ColumnMatrix[i][TopIndice][TopIndice].GetMatrixElement(NewIndiceLeft, NewIndiceRight ,Tmp); 
		      if (horizontalFlag)
			{
			  virtualSymmetryHorizontal.GetMatrixElement(TopIndice,TopIndice,Tmp3);
			  Tmp*=Tmp3;
			}
		      ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp);
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();    
	    }
	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex  Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int NewIndiceLeft = 0 ; NewIndiceLeft < DimColumTwoMatrix ; NewIndiceLeft++)
	    {
	      Tmp = 1.0; 
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p = 0 ; p < Ly ; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();    
	    }
	}
      return PEPS;
    }
  else
    {
      int NbrColumnTwoMatrix=NbrColumnMatrix*NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix*DimColumMatrix;
      
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i =0; i < NbrColumnTwoMatrix ;i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      Complex Tmp,Tmp2;
      Complex Tmp3 = 1.0;
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      for(int MiddleIndice = 0;  MiddleIndice <this->MPOBondDimension ;   MiddleIndice++)
			{
			  ColumnMatrix[i% NbrColumnMatrix][TopIndice][MiddleIndice].GetMatrixElement(NewIndiceLeft%DimColumMatrix, NewIndiceRight%DimColumMatrix  ,Tmp); 
			  ColumnMatrix[i/NbrColumnMatrix][MiddleIndice][TopIndice].GetMatrixElement(NewIndiceLeft/DimColumMatrix, NewIndiceRight/DimColumMatrix  ,Tmp2); 
			  if (horizontalFlag)
			    virtualSymmetryHorizontal.GetMatrixElement(MiddleIndice,MiddleIndice,Tmp3);
			  ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp*Tmp2*Tmp3);
			}
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0; i < Space.GetHilbertSpaceDimension(); i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();    
	    }
	  return  PEPS;
	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex Tmp,Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {	  
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p=0; p < Ly; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p = 1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      PEPS[i]= TmpResult.ComplexTr();     
	    }
	  return  PEPS;
	}
    }
}

// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

void ComplexPEPSPBC::ComputeBlockTensor ()
{
  int NbrColumnMatrix = this->PhysicalDimension*this->PhysicalDimension; 
  int DimColumMatrix = this->MPOBondDimension*  this->MPOBondDimension;

  ComplexMatrix *** ColumnMatrix = new ComplexMatrix ** [NbrColumnMatrix];
  for(int i =0; i <  NbrColumnMatrix ;i++)
    {
      ColumnMatrix[i] = new ComplexMatrix * [this->MPOBondDimension];
      for(int j=  0 ; j <this->MPOBondDimension;j++)
	{
	  ColumnMatrix[i][j] = new ComplexMatrix[this->MPOBondDimension];
	  for(int k = 0; k < this->MPOBondDimension ; k++)
	    {
	      ColumnMatrix[i][j][k] = ComplexMatrix(DimColumMatrix,DimColumMatrix,true) ;
	    }
	}
    }
  
  for(int i = 0; i < NbrColumnMatrix; i++)
    {
      for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumMatrix;  NewIndiceLeft++)
	{
	  for(int TopIndice = 0 ;  TopIndice <   this->MPOBondDimension ;  TopIndice++)
	    {
	      for (int NonZeroElementUp =0 ; NonZeroElementUp < this->NbrNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension];  NonZeroElementUp++ )
		{
		  int NewRightUpTensorIndex = this->IndiceRightNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  int NewBottomtUpTensorIndex = this->IndiceBottomNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  Complex TmpValues = this->ValuesNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft%this->MPOBondDimension][NonZeroElementUp];
		  for (int NonZeroElementDown =0 ;NonZeroElementDown <  this->NbrNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension];  NonZeroElementDown++ )
		    { 
		      ColumnMatrix[i][TopIndice][this->IndiceBottomNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]].AddToMatrixElement(NewIndiceLeft, this->IndiceRightNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]  *   this->MPOBondDimension + NewRightUpTensorIndex, this->ValuesNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension][NonZeroElementDown] * TmpValues);
		    }
		}
	    }
	}
    }

  int BlockedTensorPhysicalDimension = NbrColumnMatrix*NbrColumnMatrix; 
  int BlockedTensorVirtualDimension  = DimColumMatrix;
  
  ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [BlockedTensorPhysicalDimension];
  for(int i =0; i < BlockedTensorPhysicalDimension ;i++)
    {
      ColumnTwoMatrix[i] = ComplexMatrix(BlockedTensorVirtualDimension* BlockedTensorVirtualDimension,BlockedTensorVirtualDimension* BlockedTensorVirtualDimension,true);
    }
  Complex Tmp,Tmp2;  
  for(int i =0; i < BlockedTensorPhysicalDimension ;i++)
    {
      for(int TopIndice = 0; TopIndice < BlockedTensorVirtualDimension;  TopIndice++)
	{
	  for(int DownIndice = 0; DownIndice <BlockedTensorVirtualDimension;  DownIndice++)
	    {
	      ComplexMatrix TmpMatrix = ColumnMatrix[i%NbrColumnMatrix][TopIndice%this->MPOBondDimension][DownIndice%this->MPOBondDimension] * ColumnMatrix[i/NbrColumnMatrix][TopIndice/this->MPOBondDimension][DownIndice/this->MPOBondDimension]; 
	      for(int  Left = 0; Left < BlockedTensorVirtualDimension;   Left++)
		{
		  for(int  Right = 0; Right < BlockedTensorVirtualDimension;   Right++)
		    {
		      TmpMatrix.GetMatrixElement(Left,Right,Tmp); 
		      ColumnTwoMatrix[i].AddToMatrixElement(Left*BlockedTensorVirtualDimension+ TopIndice, Right * BlockedTensorVirtualDimension + DownIndice,Tmp);
		    }
		}
	    }
	}
    }
  
  for(int i =0; i < BlockedTensorPhysicalDimension ;i++)
    {
      for(int L=0; L <  BlockedTensorVirtualDimension; L++)
	{
	  for(int U=0; U <  BlockedTensorVirtualDimension; U++)
	    {
	      for(int R=0; R <  BlockedTensorVirtualDimension; R++)
		{
		  for(int D=0; D <  BlockedTensorVirtualDimension; D++)
		    {
		      ColumnTwoMatrix[i].GetMatrixElement(L*BlockedTensorVirtualDimension+U,R * BlockedTensorVirtualDimension+D,Tmp);
		      if (Norm(Tmp) !=0.0 )
		      cout <<i<<" "<<L << " "<<U << " "<<R<<" " << D<<" " <<Tmp<<endl;
		    }
		}
	    }
	}
    }
}



// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ComplexPEPSPBC::ComputeFockSpaceRepresentationOfAPEPSSzConstraint (int lx, int lylogtwo, int sz,  ComplexDiagonalMatrix virtualSymmetryHorizontal, ComplexDiagonalMatrix virtualSymmetryVertical, ComplexVector LeftVector, ComplexVector RightVector, bool horizontalFlag, bool verticalFlag)
{
//  this->PrintTensorElements();
  int Ly = 1;
  for(int i = 0 ; i <lylogtwo ; i++)
    {
      Ly*=2;
    }
  int Lx = lx;
  
  Spin1_2ChainNew Space (Lx*Ly,sz, 100000); 
  ComplexVector PEPS (Space.GetHilbertSpaceDimension(), true);  

  int NbrColumnMatrix=this->PhysicalDimension*this->PhysicalDimension; 
  int DimColumMatrix = this->MPOBondDimension*  this->MPOBondDimension;
  
  ComplexMatrix *** ColumnMatrix = new ComplexMatrix ** [NbrColumnMatrix];
  for(int i = 0; i < NbrColumnMatrix ; i++)
    {
      ColumnMatrix[i] = new ComplexMatrix * [this->MPOBondDimension];
      for(int j = 0 ; j < this->MPOBondDimension ; j++)
	{
	  ColumnMatrix[i][j] = new ComplexMatrix[this->MPOBondDimension];
	  for(int k = 0 ; k < this->MPOBondDimension ; k++)
	    {
	      ColumnMatrix[i][j][k] = ComplexMatrix(DimColumMatrix,DimColumMatrix,true) ;
	    }
	}
    }
  
  for(int i = 0; i < NbrColumnMatrix; i++)
    {
      for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumMatrix;  NewIndiceLeft++)
	{
	  for(int TopIndice = 0 ;  TopIndice <   this->MPOBondDimension ;  TopIndice++)
	    {
	      for (int NonZeroElementUp =0 ; NonZeroElementUp < this->NbrNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension];  NonZeroElementUp++ )
		{
		  int NewRightUpTensorIndex = this->IndiceRightNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  int NewBottomtUpTensorIndex = this->IndiceBottomNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  Complex TmpValues = this->ValuesNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft%this->MPOBondDimension][NonZeroElementUp];
		  for (int NonZeroElementDown = 0 ; NonZeroElementDown <  this->NbrNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension] ; NonZeroElementDown++ )
		    { 
		      ColumnMatrix[i][TopIndice][this->IndiceBottomNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]].AddToMatrixElement(NewIndiceLeft, this->IndiceRightNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]  *   this->MPOBondDimension + NewRightUpTensorIndex, this->ValuesNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension][NonZeroElementDown] * TmpValues);
		    }
		}
	    }
	}
    }

  if (Ly == 2 ) 
    {
      int NbrColumnTwoMatrix=NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix;
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      
      Complex Tmp; 
      Complex Tmp3 = 1.0; 
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      if (horizontalFlag)
			virtualSymmetryHorizontal.GetMatrixElement(TopIndice,TopIndice,Tmp3);
		      ColumnMatrix[i][TopIndice][TopIndice].GetMatrixElement(NewIndiceLeft, NewIndiceRight ,Tmp); 
		      ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp);
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      
	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		for(int p = 0; p < RightVector.GetVectorDimension(); p++)
		  {
		    TmpResult.GetMatrixElement(t,p,Tmp);
		    PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
		  }
	    }
	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex Tmp,Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int NewIndiceLeft = 0 ; NewIndiceLeft < DimColumTwoMatrix; NewIndiceLeft++)
	    {	  
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p = 0 ; p < Ly ; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      
	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		for(int p = 0; p < RightVector.GetVectorDimension(); p++)
		  {
		    TmpResult.GetMatrixElement(t,p,Tmp);
		    PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
		  }
	    }
	}
      return PEPS;
    }
  else
    {
      int NbrColumnTwoMatrix=NbrColumnMatrix*NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix*DimColumMatrix;
      
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i =0; i < NbrColumnTwoMatrix ;i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      Complex Tmp,Tmp2;
      Complex Tmp3 = 1.0;
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      for(int MiddleIndice = 0;  MiddleIndice <this->MPOBondDimension ;   MiddleIndice++)
			{
			  ColumnMatrix[i% NbrColumnMatrix][TopIndice][MiddleIndice].GetMatrixElement(NewIndiceLeft%DimColumMatrix, NewIndiceRight%DimColumMatrix  ,Tmp); 
			  ColumnMatrix[i/NbrColumnMatrix][MiddleIndice][TopIndice].GetMatrixElement(NewIndiceLeft/DimColumMatrix, NewIndiceRight/DimColumMatrix  ,Tmp2); 
			  if (horizontalFlag)
			    virtualSymmetryHorizontal.GetMatrixElement(MiddleIndice,MiddleIndice,Tmp3);
			  ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp*Tmp2*Tmp3);
			}
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0; i < Space.GetHilbertSpaceDimension(); i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		{
		  if( Norm(LeftVector[t]) > 1e-7  ) 
		    {
		      for(int p = 0; p < RightVector.GetVectorDimension(); p++)
			{
			  if( Norm(RightVector[p]) > 1e-7  ) 
			    {
			      TmpResult.GetMatrixElement(t,p,Tmp);
			      PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
			    }
			}
		      }
		}
	    }
	  return  PEPS;
	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex Tmp,Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {	  
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p=0; p < Ly; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      
	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		for(int p = 0; p < RightVector.GetVectorDimension(); p++)
		  {
		    TmpResult.GetMatrixElement(t,p,Tmp);
		    PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
		  }
	    }
	  return  PEPS;
	}
    }
}




// convert a state defined in the real space basis into a state in the (Kx,Ky) basis
//
// state = reference on the state to convert
// space = pointer to the Hilbert space where state is defined
// return value = state in the (Kx,Ky) basis

ComplexVector ComplexPEPSPBC::ComputeFockSpaceRepresentationOfAPEPS (int lx, int lylogtwo, ComplexDiagonalMatrix virtualSymmetryHorizontal,  ComplexDiagonalMatrix virtualSymmetryVertical, ComplexVector LeftVector, ComplexVector RightVector, bool horizontalFlag, bool verticalFlag)
{
  int Ly = 1;
  for(int i = 0 ; i <lylogtwo ; i++)
    {
      Ly*=2;
    }
  
  int Lx = lx;
  
  Spin1_2ChainFull Space (Lx*Ly);   
  ComplexVector PEPS (Space.GetHilbertSpaceDimension(), true);
  
  int NbrColumnMatrix = this->PhysicalDimension*this->PhysicalDimension; 
  int DimColumMatrix = this->MPOBondDimension*  this->MPOBondDimension;
  
  ComplexMatrix *** ColumnMatrix = new ComplexMatrix ** [NbrColumnMatrix];
  for(int i = 0; i < NbrColumnMatrix ; i++)
    {
      ColumnMatrix[i] = new ComplexMatrix * [this->MPOBondDimension];
      for(int j = 0 ; j < this->MPOBondDimension ; j++)
	{
	  ColumnMatrix[i][j] = new ComplexMatrix[this->MPOBondDimension];
	  for(int k = 0 ; k < this->MPOBondDimension ; k++)
	    {
	      ColumnMatrix[i][j][k] = ComplexMatrix(DimColumMatrix,DimColumMatrix,true) ;
	    }
	}
    }
  Complex Tmp;  
  for(int i = 0; i < NbrColumnMatrix; i++)
    {
      for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumMatrix;  NewIndiceLeft++)
	{
	  for(int TopIndice = 0 ;  TopIndice <   this->MPOBondDimension ;  TopIndice++)
	    {
	      for (int NonZeroElementUp =0 ; NonZeroElementUp < this->NbrNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension];  NonZeroElementUp++ )
		{
		  int NewRightUpTensorIndex = this->IndiceRightNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  int NewBottomtUpTensorIndex = this->IndiceBottomNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft% this->MPOBondDimension][NonZeroElementUp];
		  Complex TmpValues = this->ValuesNonZeroTensorElementTopLeft[i%this->PhysicalDimension][TopIndice][NewIndiceLeft%this->MPOBondDimension][NonZeroElementUp];
		  for (int NonZeroElementDown =0 ;NonZeroElementDown <  this->NbrNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension];  NonZeroElementDown++ )
		    { 
		      ColumnMatrix[i][TopIndice][this->IndiceBottomNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]].AddToMatrixElement(NewIndiceLeft, this->IndiceRightNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/ this->MPOBondDimension][NonZeroElementDown]  *   this->MPOBondDimension + NewRightUpTensorIndex, this->ValuesNonZeroTensorElementTopLeft[i/this->PhysicalDimension][NewBottomtUpTensorIndex][NewIndiceLeft/this->MPOBondDimension][NonZeroElementDown] * TmpValues);
		    }
		}
	    }
	}
    }
  
  if (Ly == 2 ) 
    {
      int NbrColumnTwoMatrix=NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix;
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      
      Complex Tmp; 
      Complex Tmp3 = 1.0; 
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      if (horizontalFlag)
			virtualSymmetryHorizontal.GetMatrixElement(TopIndice,TopIndice,Tmp3);
		      ColumnMatrix[i][TopIndice][TopIndice].GetMatrixElement(NewIndiceLeft, NewIndiceRight ,Tmp); 
		      ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp);
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI = i;
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult = TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}

	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		for(int p = 0; p < RightVector.GetVectorDimension(); p++)
		  {
		    TmpResult.GetMatrixElement(t,p,Tmp);
		    PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
		  }
	    }

	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex Tmp,Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int NewIndiceLeft = 0 ; NewIndiceLeft < DimColumTwoMatrix; NewIndiceLeft++)
	    {	  
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p = 0 ; p < Ly ; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p = 1; p < Lx ; p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      
	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		for(int p = 0; p < RightVector.GetVectorDimension(); p++)
		  {
		    TmpResult.GetMatrixElement(t,p,Tmp);
		    PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
		  }
	    }
	}
      return PEPS;
    }
  else
    {
      int NbrColumnTwoMatrix = NbrColumnMatrix*NbrColumnMatrix; 
      int DimColumTwoMatrix = DimColumMatrix*DimColumMatrix;
      
      ComplexMatrix * ColumnTwoMatrix = new ComplexMatrix [NbrColumnTwoMatrix];
      for(int i =0; i < NbrColumnTwoMatrix ;i++)
	{
	  ColumnTwoMatrix[i] = ComplexMatrix(DimColumTwoMatrix,DimColumTwoMatrix,true);
	}
      Complex Tmp,Tmp2;
      Complex Tmp3 = 1.0;
      for(int i = 0 ; i < NbrColumnTwoMatrix ; i++)
	{
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {
	      for(int  NewIndiceRight = 0; NewIndiceRight < DimColumTwoMatrix;  NewIndiceRight++)
		{
		  for(int TopIndice = 0; TopIndice <this->MPOBondDimension ;  TopIndice++)
		    {
		      for(int MiddleIndice = 0;  MiddleIndice <this->MPOBondDimension ;   MiddleIndice++)
			{
			  ColumnMatrix[i% NbrColumnMatrix][TopIndice][MiddleIndice].GetMatrixElement(NewIndiceLeft%DimColumMatrix, NewIndiceRight%DimColumMatrix,Tmp); 
			  ColumnMatrix[i/NbrColumnMatrix][MiddleIndice][TopIndice].GetMatrixElement(NewIndiceLeft/DimColumMatrix, NewIndiceRight/DimColumMatrix,Tmp2); 
			  if (horizontalFlag)
			    virtualSymmetryHorizontal.GetMatrixElement(MiddleIndice,MiddleIndice,Tmp3);
			  ColumnTwoMatrix[i].AddToMatrixElement(NewIndiceLeft, NewIndiceRight,Tmp*Tmp2*Tmp3);
			}
		    }
		}
	    }
	}
      
      if (verticalFlag == false)
	{
	  for(int i = 0; i < Space.GetHilbertSpaceDimension(); i++)
	    {
	      int TmpI =  Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		for(int p = 0; p < RightVector.GetVectorDimension(); p++)
		  {
		    TmpResult.GetMatrixElement(t,p,Tmp);
		    PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
		  }
	    }
	  return  PEPS;
	}
      else
	{
	  ComplexMatrix * ColumnTwoMatrixWithAZ = new ComplexMatrix [NbrColumnTwoMatrix];
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i].Copy(ColumnTwoMatrix[i]);
	    }
	  
	  Complex Tmp,Tmp2;
	  ComplexDiagonalMatrix ZFactor (DimColumTwoMatrix,true);
	  for(int  NewIndiceLeft = 0; NewIndiceLeft < DimColumTwoMatrix;  NewIndiceLeft++)
	    {	  
	      int TmpNewIndiceLeft = NewIndiceLeft;
	      for(int p=0; p < Ly; p++)
		{
		  Tmp *= virtualSymmetryVertical[TmpNewIndiceLeft% this->MPOBondDimension];
		  TmpNewIndiceLeft/= this->MPOBondDimension;
		}
	      ZFactor[NewIndiceLeft] = Tmp;
	    }
	  
	  for(int i =0; i < NbrColumnTwoMatrix ;i++)
	    {
	      ColumnTwoMatrixWithAZ[i]= ZFactor * ColumnTwoMatrix[i];
	    }
	  
	  for(int i = 0 ; i < Space.GetHilbertSpaceDimension() ; i++)
	    {
	      int TmpI = Space.StateDescription[i];
	      ComplexMatrix TmpResult = ColumnTwoMatrixWithAZ[TmpI%NbrColumnTwoMatrix];
	      for(int p =1; p <Lx;p++)
		{
		  TmpI/=NbrColumnTwoMatrix;
		  TmpResult= TmpResult*ColumnTwoMatrix[TmpI%NbrColumnTwoMatrix];
		}
	      for(int t = 0; t < LeftVector.GetVectorDimension(); t++)
		for(int p = 0; p < RightVector.GetVectorDimension(); p++)
		  {
		    TmpResult.GetMatrixElement(t,p,Tmp);
		    PEPS[i]+= Conj(LeftVector[t])* RightVector[p] * Tmp;    
		  }
	    }
	  return  PEPS;
	}
    }
}
