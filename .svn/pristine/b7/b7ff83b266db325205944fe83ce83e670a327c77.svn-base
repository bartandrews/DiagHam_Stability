#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Matrix/ComplexMatrix.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

using std::cout;
using std::endl;

double SigmaZ(int i, int j);
double SigmaX(int i, int j);
Complex SigmaY(int i, int j);

int InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile, ComplexMatrix *& tensorMatrix, int & TmpMPODimension);

int main(int argc, char** argv)
{
  cout.precision(14);
  
  // some running options and help
  OptionManager Manager ("EDTransfertMatrixMPOGivenByInputFile" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  
  Manager += SystemGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\n', "tensor-file", "name of the file containing the eigenstate to be displayed");
  (*SystemGroup) += new BooleanOption ('c', "complex", "use complex version of the code");

  (*SystemGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type GenericOverlap -h" << endl;
      return -1;
    }

  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  bool ComplexFlag = Manager.GetBoolean("complex");
  MultiColumnASCIIFile TensorsElementsDefinition;
  if (TensorsElementsDefinition.Parse(Manager.GetString("tensor-file")) == false)
    {
      TensorsElementsDefinition.DumpErrors(cout) << endl;
      return -1;
    } 
  
  if (ComplexFlag)
    {
      if (TensorsElementsDefinition.GetNbrColumns() != 6)
	{
	  cout <<" The tensor file should have 6 columnns"<<endl;
	}
    }
  else
    {
      if (TensorsElementsDefinition.GetNbrColumns() != 5)
	{
	  cout <<" The tensor file should have 5 columnns"<<endl;
	}
    }

  ComplexMatrix *  TensorMatrix = 0;
  int TmpMPODimension =0;
  int PhysicalDimension = InitializeTensorsElements(TensorsElementsDefinition, TensorMatrix, TmpMPODimension) ;

/*
  ComplexMatrix *  ResultSigmaZ = new ComplexMatrix[PhysicalDimension];
  ComplexMatrix *  ResultSigmaX = new ComplexMatrix[PhysicalDimension];
  ComplexMatrix *  ResultSigmaY = new ComplexMatrix[PhysicalDimension];

  for (int i =0; i <PhysicalDimension; i++)
    {
      ResultSigmaZ[i]= ComplexMatrix( TensorMatrix[0].GetNbrRow(), TensorMatrix[0].GetNbrRow(),true );
      ResultSigmaX[i]= ComplexMatrix( TensorMatrix[0].GetNbrRow(), TensorMatrix[0].GetNbrRow(),true ); 
      ResultSigmaY[i]= ComplexMatrix( TensorMatrix[0].GetNbrRow(), TensorMatrix[0].GetNbrRow(),true );
    }

  for (int i =0; i <PhysicalDimension; i++)
    {
      for (int j =0; j < TensorMatrix[0].GetNbrRow(); j++)
	{
	  for (int k =0; k < TensorMatrix[0].GetNbrRow(); k++)
	    {
	      Complex Tmp;
	      for( int primei = 0; primei <PhysicalDimension ;  primei++)
		{
		  TensorMatrix[primei].GetMatrixElement(j,k,Tmp);
		  ResultSigmaZ[i].AddToMatrixElement(j,k, SigmaZ(i,primei)*Tmp);  
		  ResultSigmaX[i].AddToMatrixElement(j,k, SigmaX(i,primei)*Tmp);
		  ResultSigmaY[i].AddToMatrixElement(j,k, SigmaY(i,primei)*Tmp);
		}

	      int oldLeft = j % TmpMPODimension;
	      for( int lPrime = 0; lPrime < TmpMPODimension ;  lPrime++)
		{
		  int NewJ =  j -  oldLeft +  lPrime;
		  TensorMatrix[i].GetMatrixElement(NewJ,k,Tmp);
		  ResultSigmaZ[i].AddToMatrixElement(j,k, -SigmaZ(oldLeft,lPrime)*Tmp); 
		  ResultSigmaX[i].AddToMatrixElement(j,k, -SigmaX(oldLeft,lPrime)*Tmp); 
		  ResultSigmaY[i].AddToMatrixElement(j,k, -Conj(SigmaY(oldLeft,lPrime))*Tmp); 
		}
	     
	      int oldUp = j / TmpMPODimension;
	      for( int uPrime = 0; uPrime < TmpMPODimension ;  uPrime++)
		{
		  int NewJ =  j + (uPrime -  oldUp)* TmpMPODimension;
		  TensorMatrix[i].GetMatrixElement(NewJ,k,Tmp);
		  ResultSigmaZ[i].AddToMatrixElement(j,k, SigmaZ(oldUp,uPrime)*Tmp); 
		  ResultSigmaX[i].AddToMatrixElement(j,k, SigmaX(oldUp,uPrime)*Tmp); 
		  ResultSigmaY[i].AddToMatrixElement(j,k, SigmaY(oldUp,uPrime)*Tmp); 
		}

	      int oldRight = k % TmpMPODimension;
	      for( int rPrime = 0; rPrime < TmpMPODimension ;  rPrime++)
		{
		  int NewK =  k + (rPrime -  oldRight);
		  TensorMatrix[i].GetMatrixElement(j,NewK,Tmp);
		  ResultSigmaZ[i].AddToMatrixElement(j,k, SigmaZ(oldRight,rPrime)*Tmp); 
		  ResultSigmaX[i].AddToMatrixElement(j,k, SigmaX(oldRight,rPrime)*Tmp); 
		  ResultSigmaY[i].AddToMatrixElement(j,k, SigmaY(oldRight,rPrime)*Tmp); 
		}

	      int oldDown = k / TmpMPODimension;
	      for( int dPrime = 0; dPrime < TmpMPODimension ;  dPrime++)
		{
		  int NewK =  k + (dPrime -  oldDown)* TmpMPODimension;
		  TensorMatrix[i].GetMatrixElement(j,NewK,Tmp);
		  ResultSigmaZ[i].AddToMatrixElement(j,k, -SigmaZ(oldDown,dPrime)*Tmp); 
		  ResultSigmaX[i].AddToMatrixElement(j,k, -SigmaX(oldDown,dPrime)*Tmp);  
		  ResultSigmaY[i].AddToMatrixElement(j,k, -Conj(SigmaY(oldDown,dPrime))*Tmp); 
		}
	      
	      
	    }
	}
    }
    
    for (int i =0; i <PhysicalDimension; i++)
    {
    cout <<  ResultSigmaZ[i]<<endl;
    cout <<  ResultSigmaX[i]<<endl;
    cout <<  ResultSigmaY[i]<<endl;
    }
*/

  ComplexMatrix *  ResultGaussLaw = new ComplexMatrix[PhysicalDimension];
  for (int i =0; i <PhysicalDimension; i++)
    {
      ResultGaussLaw[i]= ComplexMatrix( TensorMatrix[0].GetNbrRow(), TensorMatrix[0].GetNbrRow(),true );
    }
  
  double ChargeMatrix[3];
  ChargeMatrix[0]=1;
  ChargeMatrix[1]=1;
  ChargeMatrix[2]=-3;
  for (int i =0; i <PhysicalDimension; i++)
    {
      Complex Tmp;      
      for (int j =0; j < TensorMatrix[0].GetNbrRow(); j++)
	{
	  for (int k =0; k < TensorMatrix[0].GetNbrRow(); k++)
	    {
	      TensorMatrix[i].GetMatrixElement(j,k,Tmp);
	      int Left = j % TmpMPODimension;
	      int Up = j / TmpMPODimension;
	      int Right = k % TmpMPODimension;
	      int Down = k / TmpMPODimension;
	      ResultGaussLaw[i].AddToMatrixElement(j,k, Tmp*(ChargeMatrix[Left] +ChargeMatrix[Up] + ChargeMatrix[Right] +  ChargeMatrix[Down]));
	    }
	}
      cout << ResultGaussLaw[i]<<endl;
    }
}

double SigmaZ(int i, int j)
{
  if (i==j )
    {
      if ( i ==0)
	return 1.0;
      if ( i ==1)
	return -1.0;
    }
  return 0.0;
};


double SigmaX (int i, int j)
{
  if ((i ==0)&& (j==1))
    return 1.0;
  if  ((i ==1)&& (j==0 )) 
    return 1.0;
  return 0.0;
};

Complex SigmaY(int i, int j)
{
  if ( (i ==0)&& (j==1 ))
    return Complex (0.0,1.0);
  if ( (i ==1)&& (j==0 ) )
    return Complex (0.0,-1.0);
  return Complex (0.0,0.0);
};


int InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile, ComplexMatrix *& tensorMatrix, int & TmpMPODimension)
{
  int NbrNonZeroElements = tensorElementsFile.GetNbrLines();
  int* IndexVertical  = tensorElementsFile.GetAsIntegerArray (0);
  int* IndexLeft = tensorElementsFile.GetAsIntegerArray (1);
  int* IndexUp = tensorElementsFile.GetAsIntegerArray (2);
  int* IndexRight = tensorElementsFile.GetAsIntegerArray (3);
  int* IndexDown  = tensorElementsFile.GetAsIntegerArray (4);

  Complex * ElementsValues = tensorElementsFile.GetAsComplexArray (5);
  int VerticalDimension = 0;


  TmpMPODimension = 0;
  for(int i = 0 ; i < NbrNonZeroElements; i++)
    {
      if (IndexVertical[i] >  VerticalDimension)
	{
	  VerticalDimension = IndexVertical[i];
	}
      if (IndexLeft[i] > TmpMPODimension)
	TmpMPODimension = IndexLeft[i];
    }
  
  TmpMPODimension++;
  VerticalDimension++;
  
  tensorMatrix = new ComplexMatrix [VerticalDimension];
  
  for (int i =0; i <VerticalDimension; i++)
    {
      tensorMatrix[i] = ComplexMatrix(TmpMPODimension* TmpMPODimension,TmpMPODimension* TmpMPODimension,true );
    }
  
  cout <<"VerticalDimension = "<<VerticalDimension<<endl;
  cout <<"TmpMPODimension = " <<  TmpMPODimension<<endl;
  

  for(int i = 0 ; i < NbrNonZeroElements; i++)
    {
      tensorMatrix[IndexVertical[i]].SetMatrixElement(IndexLeft[i]+ TmpMPODimension *IndexUp[i],IndexRight[i] + TmpMPODimension *IndexDown[i],ElementsValues[i]);
    }
  
  delete [] IndexVertical;
  delete [] IndexDown;
  delete [] IndexUp;
  delete [] IndexLeft;
  delete [] IndexRight;
  delete [] ElementsValues;
  return VerticalDimension;
};
