#include "MPOPeratorSixVertexModelTransferMatrixSquare.h"

using std::cout;
using std::endl;



MPOPeratorSixVertexModelTransferMatrixSquare::MPOPeratorSixVertexModelTransferMatrixSquare()
{
}



MPOPeratorSixVertexModelTransferMatrixSquare::MPOPeratorSixVertexModelTransferMatrixSquare( AbstractArchitecture * architecture)
{
  this->PhysicalDimension = 2;
  this->MPOBondDimension = 4; 
  this->InitializeTensorsElements();
  this->Architecture = architecture;
  this->LeftVector = new double[this->MPOBondDimension];
  this->RightVector = new double[this->MPOBondDimension];
  this->RightVector[0] = 0;  this->RightVector[1] = 1.0;  this->RightVector[2] = -1.0; this->RightVector[3] = 0;
  this->LeftVector[0] = 0;  this->LeftVector[1] = -1.0;  this->LeftVector[2] = 1.0; this->LeftVector[3] = 0;
}



MPOPeratorSixVertexModelTransferMatrixSquare::~MPOPeratorSixVertexModelTransferMatrixSquare()
{
  delete [] LeftVector;
  delete [] RightVector;
  delete [] ElementsValues;
  delete [] IndexValues;     
}


void MPOPeratorSixVertexModelTransferMatrixSquare::InitializeTensorsElements()
{
  double T[this->PhysicalDimension][this->PhysicalDimension][this->PhysicalDimension][this->PhysicalDimension];
  
  for (int j = 0; j < this->PhysicalDimension ; j++)
    {
      for (int k=0; k < this->PhysicalDimension; k++)
	{
	  for (int l=0; l <  this->PhysicalDimension; l++)
	    {
	      for (int m = 0;  m <  this->PhysicalDimension; m++)
		T[j][k][l][m] = 0;
	    }
	}
    }
  
  T[1][1][1][1] = 0.5;
  T[0][0][0][0] = 0.5;
  T[1][0][1][0] = -0.5;
  T[0][1][0][1] = -0.5;
  T[1][1][0][0] = 1;
  T[0][0][1][1] = 1;
  
  this->NbrNonZeroElements = 0;
  for (int j = 0; j <  this->MPOBondDimension; j++)
    {
      for (int k = 0; k < this->PhysicalDimension; k++)
	{
	  for (int l = 0; l< this->MPOBondDimension ; l++)
	    {
	      for (int m=0; m< this->PhysicalDimension; m++)
		{
		  double Tmp = 0.0;
		  for (int p=0; p < this->PhysicalDimension; p++)
		    {
		      Tmp += T[j% this->PhysicalDimension][k][l% this->PhysicalDimension][p] * T[j/ this->PhysicalDimension][p][l/ this->PhysicalDimension][m];
		    }
		  if (Tmp != 0.0)
		    {
		      this->NbrNonZeroElements++;
		    }
		}
	    }
	}
    }
  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  this->ElementsValues = new double [this->NbrNonZeroElements];
  this->IndexValues = new unsigned int[this->NbrNonZeroElements];

  this->NbrNonZeroElements=0;  
  for (int j = 0; j<  this->MPOBondDimension; j++)
    {
      for (int k = 0; k < this->PhysicalDimension; k++)
	{
	  for (int l = 0; l< this->MPOBondDimension ; l++)
	    {
	      for (int m = 0;  m < this->PhysicalDimension; m++)
		{
		  double Tmp = 0.0;
		  for (int p=0; p < this->PhysicalDimension; p++)
		    {
		      Tmp += T[j% this->PhysicalDimension][k][l% this->PhysicalDimension][p] * T[j/ this->PhysicalDimension][p][l/ this->PhysicalDimension][m];
		    }
		  if (Tmp != 0.0)
		    {
		      this->ElementsValues[this->NbrNonZeroElements]=Tmp;
		      this->IndexValues[this->NbrNonZeroElements] = GetTensorIndexFromAllIndices(m, k, l,j);
		      this->NbrNonZeroElements++;
		    }
		}
	    }
	}
    }   
}


void MPOPeratorSixVertexModelTransferMatrixSquare::PrintTensorElements()
{
  cout <<"#Tensor index indexDown indexUp indexLeft indexRight Check Index Values" <<endl;
  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
      int Tmp = GetTensorIndexFromAllIndices( MPOIndiceDown,  MPOIndiceUp,  MPOIndiceLeft,  MPOIndiceRight);
      cout << this->IndexValues[i] <<" "<<MPOIndiceDown<< " "<< MPOIndiceUp<< " "<< MPOIndiceLeft<< " "<<MPOIndiceRight<<" " <<Tmp<<" "<< this->ElementsValues[i]<<endl;
    }
}
