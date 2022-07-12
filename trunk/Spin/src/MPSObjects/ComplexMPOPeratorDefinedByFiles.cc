#include "ComplexMPOPeratorDefinedByFiles.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

using std::cout;
using std::endl;
 


ComplexMPOPeratorDefinedByFiles::ComplexMPOPeratorDefinedByFiles()
{
}



ComplexMPOPeratorDefinedByFiles::ComplexMPOPeratorDefinedByFiles(MultiColumnASCIIFile & tensorElementsFile, MultiColumnASCIIFile & boundaryVectorsFile, AbstractArchitecture * architecture)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->InitializeBoundaryVectors(boundaryVectorsFile);
  this->Architecture= architecture;
}



ComplexMPOPeratorDefinedByFiles::~ComplexMPOPeratorDefinedByFiles()
{
  delete [] LeftVector;
  delete [] RightVector;
  delete [] ElementsValues;
  delete [] IndexValues;    
}


void ComplexMPOPeratorDefinedByFiles::InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile)
{
  this->NbrNonZeroElements = tensorElementsFile.GetNbrLines();
  int* IndexLeft = tensorElementsFile.GetAsIntegerArray (0);  
  int* IndexUp = tensorElementsFile.GetAsIntegerArray (1);
  int* IndexRight = tensorElementsFile.GetAsIntegerArray (2);
  int* IndexDown  = tensorElementsFile.GetAsIntegerArray (3);
  this->ElementsValues = tensorElementsFile.GetAsComplexArray (4);
  int TmpPhysicalDimension = 0;
  int TmpMPODimension = 0;
  
  
  cout <<"Nbr Non zero Elements = "<<  this->NbrNonZeroElements<<endl;
  this->IndexValues = new unsigned int[this->NbrNonZeroElements];
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      if (IndexDown[i] > TmpPhysicalDimension)
	TmpPhysicalDimension = IndexDown[i];
      if (IndexLeft[i] > TmpMPODimension)
	TmpMPODimension = IndexLeft[i];
    }
  
  this->PhysicalDimension = TmpPhysicalDimension+1;
  this->MPOBondDimension =  TmpMPODimension+1;
  
  for(int i = 0 ; i < this->NbrNonZeroElements; i++)
    {
      this->IndexValues[i] = this->GetTensorIndexFromAllIndices(IndexDown[i],IndexUp[i], IndexLeft[i] ,IndexRight[i]);
    }
  cout <<" Physical Dimension = " <<  this->PhysicalDimension<<endl;;
  cout <<" MPO Dimension = " <<  this->MPOBondDimension <<endl;;
  delete [] IndexDown;
  delete [] IndexUp;
  delete [] IndexLeft;
  delete [] IndexRight;
}


void ComplexMPOPeratorDefinedByFiles::InitializeBoundaryVectors(MultiColumnASCIIFile & boundaryVectorsFile)
{
  if ( boundaryVectorsFile.GetNbrLines() !=   this->MPOBondDimension)
  {
    cout <<"error size of the boundary vectors are not compatible with the one" <<endl;
  }
  else
  {
    if ( boundaryVectorsFile.GetNbrColumns() == 2  )
      {
	this->LeftVector = boundaryVectorsFile.GetAsComplexArray (0); 
	this->RightVector = boundaryVectorsFile.GetAsComplexArray (1);
      }
    else
      {
	if ( boundaryVectorsFile.GetNbrColumns() == 1  )
	  {
	    this->LeftVector = boundaryVectorsFile.GetAsComplexArray (0); 
	    this->RightVector = boundaryVectorsFile.GetAsComplexArray (0);
	  }
      }
  }
}


