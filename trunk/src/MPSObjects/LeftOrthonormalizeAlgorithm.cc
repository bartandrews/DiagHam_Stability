#include <iostream>
#include <sys/time.h>

#include "Matrix/ComplexMatrix.h"
#include "LeftOrthonormalizeAlgorithm.h"
#include "CompletelyPositiveMap.h"
#include "LanczosAlgorithm/ComplexArnoldiCPMapsAlgorithm.h"

using std::cout;
using std::endl;

LeftOrthonormalizeAlgorithm::LeftOrthonormalizeAlgorithm (ComplexMatrix * mps, unsigned int physicalDimension, double accuracy, ComplexMatrix & initialGuess)
{
  this->Mps = mps;
  this->PhysicalDimension = physicalDimension;
  this->Accuracy = accuracy;
  this->CenterMatrix = initialGuess;
  this->MpsInLeftForm = new ComplexMatrix[this->PhysicalDimension];
  for(int i = 0; i < this->PhysicalDimension; i++)
    this->MpsInLeftForm[i] = ComplexMatrix(this->Mps[0].GetNbrRow(),this->Mps[0].GetNbrColumn() ,true);
}


LeftOrthonormalizeAlgorithm::~LeftOrthonormalizeAlgorithm()
{
/*  for(int i = 0; i < this->PhysicalDimension; i++)
    delete this->MpsInLeftForm[i];*/
  delete [] this->MpsInLeftForm;
}

void LeftOrthonormalizeAlgorithm::RunAlgorithm()
{
  ComplexMatrix TmpUnitary (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
  ComplexMatrix TmpUpperTriangular (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
  this->CenterMatrix.QRDecompositionFromLapack (TmpUnitary,  TmpUpperTriangular);
  TmpUpperTriangular/= TmpUpperTriangular.FrobeniusNorm();
  
  for(int i = 0; i < PhysicalDimension; i++)
    {
      this->MpsInLeftForm[i] = TmpUpperTriangular * this->Mps[i];
      cout <<this->MpsInLeftForm[i]<<endl;
    }
  
  this->OldCenterMatrix = TmpUpperTriangular;

  ComplexMatrix ACenter (this->MpsInLeftForm[0].GetNbrRow() * PhysicalDimension, this->MpsInLeftForm[0].GetNbrColumn(),true);
//  cout <<this->MpsInLeftForm[0].GetNbrRow() * PhysicalDimension<<" "<<  this->MpsInLeftForm[0].GetNbrColumn() <<endl;
  
  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j <  ACenter.GetNbrColumn() ; j++)
	{
//	  cout << this->MpsInLeftForm[i/this->MpsInLeftForm[0].GetNbrRow()].GetMatrixElement(i%this->MpsInLeftForm[0].GetNbrRow(),j) <<endl;
	  ACenter.SetMatrixElement(i,j, this->MpsInLeftForm[i/this->MpsInLeftForm[0].GetNbrRow()].GetMatrixElement(i%this->MpsInLeftForm[0].GetNbrRow(),j));
	}
    }
  
  cout << ACenter<<endl;

  ComplexMatrix TmpUnitary2 (ACenter.GetNbrRow(), ACenter.GetNbrRow(), true);
  
  ACenter.QRDecompositionFromLapack (TmpUnitary2,this->CenterMatrix);
  
  cout <<"After   ACenter.QRDecompositionFromLapack (TmpUnitary2, this->CenterMatrix)"<<endl;
  
  cout << this->CenterMatrix<<endl;
  cout <<TmpUnitary2<<endl;

  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j < this->MpsInLeftForm[0].GetNbrColumn() ; j++)
	{
	  this->MpsInLeftForm[i/this->MpsInLeftForm[0].GetNbrRow()].SetMatrixElement(i%this->MpsInLeftForm[0].GetNbrRow(),j , TmpUnitary2.GetMatrixElement(i,j));
	}
    }
  
  for(int i = 0; i < PhysicalDimension; i++)
    {
      cout <<this->MpsInLeftForm[i]<<endl;
    }
  
  this->Eigenvalue = this->CenterMatrix.FrobeniusNorm();
  this->CenterMatrix/=  this->Eigenvalue;
  double ActualAccuracy = (this->CenterMatrix - this->OldCenterMatrix).FrobeniusNorm();
  cout <<"Actual Accuracy = " << ActualAccuracy<<endl;
  while( ActualAccuracy > this->Accuracy )
    {
      CompletelyPositiveMap Map (this->PhysicalDimension, this->Mps, this->MpsInLeftForm);
      ComplexArnoldiCPMapsAlgorithm Arnoldi (&Map,  ActualAccuracy*0.1, 1, 100, true);  
      Arnoldi.InitializeLanczosAlgorithm(this->CenterMatrix);
      Arnoldi.RunLanczosAlgorithm(1);
      while (Arnoldi.TestConvergence() == false)
	{
	  timeval TotalStartingTime;
	  timeval TotalEndingTime;
	  gettimeofday (&(TotalStartingTime), 0);
	  Arnoldi.RunLanczosAlgorithm(1);
	  gettimeofday (&(TotalEndingTime), 0);
	  double Dt = (double) ((TotalEndingTime.tv_sec - TotalStartingTime.tv_sec) + 
				((TotalEndingTime.tv_usec - TotalStartingTime.tv_usec) / 1000000.0));		      
	  cout << "iteration done in " << Dt << "s" << endl;
	}
      
      this->CenterMatrix = Arnoldi.GetGroundState();
      
      ComplexMatrix TmpUnitary (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
      ComplexMatrix TmpUpperTriangular (this->CenterMatrix.GetNbrRow(),this->CenterMatrix.GetNbrRow(),true);
      this->CenterMatrix.QRDecompositionFromLapack (TmpUnitary,  TmpUpperTriangular);
      TmpUpperTriangular/= TmpUpperTriangular.FrobeniusNorm();
      
      for(int i = 0; i < PhysicalDimension; i++)
	this->MpsInLeftForm[i] = TmpUpperTriangular * this->Mps[i];
      
      this->OldCenterMatrix = TmpUpperTriangular;
      
  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j <  ACenter.GetNbrColumn() ; j++)
	{
	  ACenter.SetMatrixElement(i,j, this->MpsInLeftForm[i/this->MpsInLeftForm[i].GetNbrRow()].GetMatrixElement(i%this->MpsInLeftForm[i].GetNbrRow(),j));
	}
    }
      
  ComplexMatrix TmpUnitary2 ( ACenter.GetNbrRow(),  ACenter.GetNbrColumn());
  
  ACenter.QRDecompositionFromLapack (TmpUnitary2, this->CenterMatrix);
  
  for(int i = 0; i < ACenter.GetNbrRow() ; i++)
    {
      for(int j = 0; j <  ACenter.GetNbrColumn() ; j++)
	{
	  this->MpsInLeftForm[i/this->MpsInLeftForm[i].GetNbrRow()].SetMatrixElement(i%this->MpsInLeftForm[i].GetNbrRow(),j , TmpUnitary2.GetMatrixElement(i,j));
	}
    }
  
      
      this->Eigenvalue = this->CenterMatrix.FrobeniusNorm();
      this->CenterMatrix/=  this->Eigenvalue;
      ActualAccuracy = (this->CenterMatrix - this->OldCenterMatrix).FrobeniusNorm();
    }
}
