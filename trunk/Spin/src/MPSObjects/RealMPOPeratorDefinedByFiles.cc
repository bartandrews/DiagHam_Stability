#include "RealMPOPeratorDefinedByFiles.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

using std::cout;
using std::endl;
 


RealMPOPeratorDefinedByFiles::RealMPOPeratorDefinedByFiles()
{
}



RealMPOPeratorDefinedByFiles::RealMPOPeratorDefinedByFiles(MultiColumnASCIIFile & tensorElementsFile, MultiColumnASCIIFile & boundaryVectorsFile,AbstractArchitecture * architecture)
{
  this->InitializeTensorsElements(tensorElementsFile);
  this->InitializeBoundaryVectors(boundaryVectorsFile);
  this->Architecture = architecture;
}



RealMPOPeratorDefinedByFiles::~RealMPOPeratorDefinedByFiles()
{
  delete [] LeftVector;
  delete [] RightVector;
  delete [] ElementsValues;
  delete [] IndexValues;    
}


void RealMPOPeratorDefinedByFiles::InitializeTensorsElements(MultiColumnASCIIFile & tensorElementsFile)
{

 this->NbrNonZeroElements = tensorElementsFile.GetNbrLines();

 int* IndexDown  = tensorElementsFile.GetAsIntegerArray (0);
 int* IndexUp = tensorElementsFile.GetAsIntegerArray (1);
 int* IndexLeft = tensorElementsFile.GetAsIntegerArray (2);
 int* IndexRight = tensorElementsFile.GetAsIntegerArray (3);
 this->ElementsValues = tensorElementsFile.GetAsDoubleArray (4);
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


void RealMPOPeratorDefinedByFiles::InitializeBoundaryVectors(MultiColumnASCIIFile & boundaryVectorsFile)
{
  if ( boundaryVectorsFile.GetNbrLines() !=   this->MPOBondDimension)
  {
    cout <<"error size of the boundary vectors are not compatible with the one" <<endl;
  }
  else
  {
    if ( boundaryVectorsFile.GetNbrColumns() == 2  )
{
     this->LeftVector = boundaryVectorsFile.GetAsDoubleArray (0); 
     this->RightVector = boundaryVectorsFile.GetAsDoubleArray (1);
}
else
{
   if ( boundaryVectorsFile.GetNbrColumns() == 1  )
{
     this->LeftVector = boundaryVectorsFile.GetAsDoubleArray (0); 
     this->RightVector = boundaryVectorsFile.GetAsDoubleArray (0);
}
}
  }

}


