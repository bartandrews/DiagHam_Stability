#include "DMRGFiniteSizeRealOBCMainTask.h"
#include "RealMPSSite.h"
#include <iostream>
#include <sys/time.h>
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/MonoProcessorArchitecture.h"
#include "LanczosAlgorithm/LanczosManager.h"
#include "LanczosAlgorithm/AbstractLanczosAlgorithm.h"
#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

using std::cout;
using std::endl;


DMRGFiniteSizeRealOBCMainTask::DMRGFiniteSizeRealOBCMainTask(RealMPSSite * latticeSite, AbstractMPOperatorOBC * mPOperator, int nbrSites, int nbrSweep, int maximumBondDimension,  AbstractArchitecture * architecture, LanczosManager* lanczos) :  NbrSites(nbrSites)
{
  this->LatticeSite = latticeSite;
  this->MPOperator = mPOperator;
  this->NbrSweep = nbrSweep;
  this->MaximumBondDimension = maximumBondDimension;
  this->Architecture = architecture;
  this->PreviousEnergy = 1e50;
  this->AlgorithmManager = lanczos;
}


DMRGFiniteSizeRealOBCMainTask::~DMRGFiniteSizeRealOBCMainTask(){ }

void DMRGFiniteSizeRealOBCMainTask::RunAlgorithm()
{
//  this->InitializeLattice();
//  this->InitializeLatticeUsingIDMRG();
  this->InitializeLatticeUsingIDMRGAndStatePrediction();

  this->MPOperator->SetDMRGFlag(false);
  for(int CurrentSweep = 0; CurrentSweep <  this->NbrSweep; CurrentSweep++)
    {
      cout <<"Sweep : "<<CurrentSweep<<endl;
      
      for (int i = 0; i <  this->NbrSites - 1; i++)
	{
	  cout <<"From left to right "<<" site "<<i <<endl;
          this->MPOperator->SetSite(&this->LatticeSite[i]);
	  this->OptimizeUsingLanczosLanczosAlgorithm (i);
	  this->LatticeSite[i].BringMInLeftCanonicalFormCareful();  
	}
     for (int i =  this->NbrSites - 1; i >  0; i--)
	{
	  cout <<"From right to left,"<<" site "<<i <<endl;
	  this->MPOperator->SetSite(&this->LatticeSite[i]);
	  this->OptimizeUsingLanczosLanczosAlgorithm (i);
	  this->LatticeSite[i].BringMInRightCanonicalFormCareful();
	}
    }
}




void DMRGFiniteSizeRealOBCMainTask::InitializeLattice()
{
  this->MPOperator->SetDMRGFlag(false);
  for (int i = 0 ; i < NbrSites ; i++) 
    {
      LatticeSite[i].InitializeWithRandomMatrices();
    }

  for (int i = NbrSites - 1 ; i>0 ; i--) 
    {
      LatticeSite[i].BringMInRightCanonicalFormCareful(); 
    }
}



void DMRGFiniteSizeRealOBCMainTask::InitializeLatticeUsingIDMRG()
{
  this->MPOperator->SetDMRGFlag(true);
  RealDiagonalMatrix SingularValues;

  for (int i = 0 ; i < NbrSites/2 ; i++) 
    {
       cout <<" Initializazion Site " << i <<endl;
       this->MPOperator->SetSiteLeftAndRight(&(this->LatticeSite[i]),&this->LatticeSite[NbrSites-1-i]);
       this->TwoSiteOptimizationUsingLanczosLanczosAlgorithm (&LatticeSite[i] ,&this->LatticeSite[NbrSites-1-i],SingularValues);
       cout <<"End Initializazion Site " << i <<endl;
    }



  RealMatrix * TmpM = LatticeSite[NbrSites/2 - 1].GetM();
  for(int i = 0; i <  this->MPOperator->GetPhysicalDimension() ; i++)
  {
    TmpM[i] = TmpM[i]*SingularValues;
  }

  for (int i = NbrSites/2 - 1 ; i>0 ; i--) 
    {
      LatticeSite[i].BringMInRightCanonicalFormCareful(); 
    }
}


void DMRGFiniteSizeRealOBCMainTask::InitializeLatticeUsingIDMRGAndStatePrediction()
{
  this->MPOperator->SetDMRGFlag(true);
  RealDiagonalMatrix SingularValues;
  RealVector * StatePredict = 0;
  int HalfNumberSite = NbrSites/2;
  this->MPOperator->SetSiteLeftAndRight(&(this->LatticeSite[0]),&this->LatticeSite[NbrSites-1]);
  this->TwoSiteOptimizationUsingLanczosLanczosAlgorithm(&LatticeSite[0] ,&this->LatticeSite[NbrSites-1],SingularValues);
  this->PreviousSingularValues = SingularValues;
  for (int i = 1 ; i < HalfNumberSite -1; i++) 
    {
       cout <<" Initializazion Site " << i <<endl;
       this->MPOperator->SetSiteLeftAndRight(&(this->LatticeSite[i]),&this->LatticeSite[NbrSites-1-i]);

       this->TwoSiteOptimizationUsingLanczosLanczosAlgorithmAndStatePrediction(&LatticeSite[i] ,&this->LatticeSite[NbrSites-1-i],SingularValues,StatePredict);
       StatePredict = LatticeSite[i].StatePrediction(&this->LatticeSite[NbrSites-1-i], SingularValues, this->PreviousSingularValues);
       this->PreviousSingularValues = SingularValues;
       cout <<"End Initializazion Site " << i <<endl;
    }

  this->MPOperator->SetSiteLeftAndRight(&(this->LatticeSite[HalfNumberSite-1]),&this->LatticeSite[HalfNumberSite]);
  this->TwoSiteOptimizationUsingLanczosLanczosAlgorithmAndStatePrediction(&LatticeSite[HalfNumberSite-1] ,&this->LatticeSite[HalfNumberSite],SingularValues,StatePredict);

  RealMatrix * TmpM = LatticeSite[NbrSites/2 - 1].GetM();
  for(int i = 0; i <  this->MPOperator->GetPhysicalDimension() ; i++)
  {
    TmpM[i] = TmpM[i]*SingularValues;
  }

  for (int i = NbrSites/2 - 1 ; i>0 ; i--) 
    {
      LatticeSite[i].BringMInRightCanonicalFormCareful(); 
    }
}



void DMRGFiniteSizeRealOBCMainTask::OptimizeUsingLanczosLanczosAlgorithm (int siteIndex)
{
  if (this->MPOperator->GetHilbertSpaceDimension() < 500 )
    {
    RealSymmetricMatrix HRep (this->MPOperator->GetHilbertSpaceDimension(), true);
    this->MPOperator->GetHamiltonian(HRep);

    if (this->MPOperator->GetHilbertSpaceDimension() > 1)
     {
#ifdef __LAPACK__
      RealDiagonalMatrix TmpDiag (this->MPOperator->GetHilbertSpaceDimension());
      RealMatrix Q(this->MPOperator->GetHilbertSpaceDimension(), this->MPOperator->GetHilbertSpaceDimension());
      HRep.LapackDiagonalize(TmpDiag, Q);
      cout <<"Highest energy = " << TmpDiag[0]<<" change = " <<  (( TmpDiag[0] - this->PreviousEnergy) /this->PreviousEnergy)<< endl;
      this->PreviousEnergy = TmpDiag[0];
      this->LatticeSite[siteIndex].UpdateFromVector(&Q[0]);


#endif

}
}
else
{
  RealVector * TmpVector = 0;
  AbstractArchitecture * MonoProc = new MonoProcessorArchitecture;
  AbstractLanczosAlgorithm* LanczosAlgorithm = this->AlgorithmManager->GetLanczosAlgorithm(MonoProc, true);
  this->LatticeSite[siteIndex].GetMatrixInVectorForm(TmpVector);
  LanczosAlgorithm->InitializeLanczosAlgorithm(*TmpVector);
  double GroundStateEnergy;
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = 0;
  LanczosAlgorithm->SetHamiltonian(this->MPOperator);
  cout << "Run Lanczos Algorithm" << endl; 
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  timeval TotalCurrentTime;
  double Dt;
  gettimeofday (&(TotalStartingTime), 0);
  int StartTimeSecond = TotalStartingTime.tv_sec;
  LanczosAlgorithm->RunLanczosAlgorithm(3);
  CurrentNbrIterLanczos = 4;
  RealTriDiagonalSymmetricMatrix TmpMatrix;  
  int CurrentTimeSecond = TotalCurrentTime.tv_sec;
  if ((LanczosAlgorithm->TestConvergence() == true))
	{
	  TmpMatrix.Copy(LanczosAlgorithm->GetDiagonalizedMatrix());
	  TmpMatrix.SortMatrixUpOrder();
          Lowest = TmpMatrix.DiagonalElement(0);
	}
  while ((LanczosAlgorithm->TestConvergence() == false)&&( CurrentNbrIterLanczos < 2000 ))
  {
     ++CurrentNbrIterLanczos;
     LanczosAlgorithm->RunLanczosAlgorithm(1);
     TmpMatrix.Copy(LanczosAlgorithm->GetDiagonalizedMatrix());
     TmpMatrix.SortMatrixUpOrder();
     Lowest = TmpMatrix.DiagonalElement(0);
     Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
     PreviousLowest = Lowest; 
     cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << " ";
     gettimeofday (&(TotalEndingTime), 0);
     CurrentTimeSecond = TotalEndingTime.tv_sec;
     Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
		((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
	      cout << "(" << Dt << " s for step " << CurrentNbrIterLanczos << ")"<< endl;;
     TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
     TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
   }
      GroundStateEnergy = Lowest;
      cout << endl;
      cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	   << CurrentNbrIterLanczos << endl;
      RealVector GroundState = LanczosAlgorithm->GetGroundState();
      cout <<"Highest energy = " << GroundStateEnergy<<" change = " <<  ((GroundStateEnergy - this->PreviousEnergy) /this->PreviousEnergy)<< endl;
      this->PreviousEnergy = GroundStateEnergy;
      this->LatticeSite[siteIndex].UpdateFromVector(&GroundState);
}

}


void DMRGFiniteSizeRealOBCMainTask::TwoSiteOptimizationUsingLanczosLanczosAlgorithm ( RealMPSSite * leftSite , RealMPSSite * rightSite, RealDiagonalMatrix & singularValues)
{
  if (this->MPOperator->GetTwoSitesHilbertSpaceDimension() < 500 )
    {
    int Dimension = this->MPOperator->GetTwoSitesHilbertSpaceDimension();
    cout <<"Dimension = " <<Dimension<<endl;
    RealSymmetricMatrix HRep (Dimension, true);
    this->MPOperator->GetTwoSitesHamiltonian(HRep);
    if (Dimension > 1)
     {
#ifdef __LAPACK__
      RealDiagonalMatrix TmpDiag (Dimension);
      RealMatrix Q(Dimension,Dimension);
      HRep.LapackDiagonalize(TmpDiag, Q);
      cout <<"Highest energy = " << TmpDiag[0]<<" change = " <<  (( TmpDiag[0] - this->PreviousEnergy) /this->PreviousEnergy)<< endl;
      leftSite->SymmetricUpdateOfTwoSites(rightSite, &Q[0],singularValues);
     this->PreviousEnergy = TmpDiag[0];
#endif

}
}
else
{
  RealVector * TmpVector = 0;
  AbstractArchitecture * MonoProc = new MonoProcessorArchitecture;
  AbstractLanczosAlgorithm* LanczosAlgorithm = this->AlgorithmManager->GetLanczosAlgorithm(MonoProc, true);
//  this->LatticeSite[siteIndex].GetMatrixInVectorForm(TmpVector);
  LanczosAlgorithm->SetHamiltonian(this->MPOperator);
  LanczosAlgorithm->InitializeLanczosAlgorithm();
  double GroundStateEnergy;
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = 0;
  LanczosAlgorithm->SetHamiltonian(this->MPOperator);
  cout << "Run Lanczos Algorithm" << endl;
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  timeval TotalCurrentTime;
  double Dt;
  gettimeofday (&(TotalStartingTime), 0);
  int StartTimeSecond = TotalStartingTime.tv_sec;
  LanczosAlgorithm->RunLanczosAlgorithm(3);
  CurrentNbrIterLanczos = 4;
  RealTriDiagonalSymmetricMatrix TmpMatrix;  
  int CurrentTimeSecond = TotalCurrentTime.tv_sec;

  if ((LanczosAlgorithm->TestConvergence() == true))
  {
      TmpMatrix.Copy(LanczosAlgorithm->GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(0);
  }

  while ((LanczosAlgorithm->TestConvergence() == false)&&( CurrentNbrIterLanczos < 2000 ))
  {
     ++CurrentNbrIterLanczos;
     LanczosAlgorithm->RunLanczosAlgorithm(1);
     TmpMatrix.Copy(LanczosAlgorithm->GetDiagonalizedMatrix());
     TmpMatrix.SortMatrixUpOrder();
     Lowest = TmpMatrix.DiagonalElement(0);
     Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
     PreviousLowest = Lowest; 
     cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << " ";
     gettimeofday (&(TotalEndingTime), 0);
     CurrentTimeSecond = TotalEndingTime.tv_sec;
     Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
		((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
	      cout << "(" << Dt << " s for step " << CurrentNbrIterLanczos << ")"<< endl;;
     TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
     TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
   }
      GroundStateEnergy = Lowest;
      cout << endl;
      cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	   << CurrentNbrIterLanczos << endl;
      RealVector GroundState = LanczosAlgorithm->GetGroundState();
      cout <<"Highest energy = " << GroundStateEnergy<<" change = " <<  ((GroundStateEnergy - this->PreviousEnergy) /this->PreviousEnergy)<< endl;
      this->PreviousEnergy = GroundStateEnergy;
      leftSite->SymmetricUpdateOfTwoSites(rightSite, &GroundState,singularValues);
}
}



void DMRGFiniteSizeRealOBCMainTask::TwoSiteOptimizationUsingLanczosLanczosAlgorithmAndStatePrediction ( RealMPSSite * leftSite , RealMPSSite * rightSite, RealDiagonalMatrix & singularValues, RealVector * statePredict)
{
  if (this->MPOperator->GetTwoSitesHilbertSpaceDimension() < 500 )
    {
    int Dimension = this->MPOperator->GetTwoSitesHilbertSpaceDimension();
    cout <<"Dimension = " <<Dimension<<endl;
    RealSymmetricMatrix HRep (Dimension, true);
    this->MPOperator->GetTwoSitesHamiltonian(HRep);
    if (Dimension > 1)
     {
#ifdef __LAPACK__
      RealDiagonalMatrix TmpDiag (Dimension);
      RealMatrix Q(Dimension,Dimension);
      HRep.LapackDiagonalize(TmpDiag, Q);
      cout <<"Highest energy = " << TmpDiag[0]<<" change = " <<  (( TmpDiag[0] - this->PreviousEnergy) /this->PreviousEnergy)<< endl;
     cout <<" Before      this->LatticeSite[0].SymmetricUpdateOfTwoSites(leftSite,rightSite, &Q[0],singularValues)"<<endl;
     leftSite->SymmetricUpdateOfTwoSites(rightSite, &Q[0],singularValues);
     this->PreviousEnergy = TmpDiag[0];
     cout <<" After      this->LatticeSite[0].SymmetricUpdateOfTwoSites(leftSite,rightSite, &Q[0],singularValues)"<<endl;
#endif

}
}
else
{
  AbstractArchitecture * MonoProc = new MonoProcessorArchitecture;
  AbstractLanczosAlgorithm* LanczosAlgorithm = this->AlgorithmManager->GetLanczosAlgorithm(MonoProc, true);
  LanczosAlgorithm->SetHamiltonian(this->MPOperator);
  LanczosAlgorithm->InitializeLanczosAlgorithm(*statePredict);
  double GroundStateEnergy;
  double Precision = 1.0;
  double PreviousLowest = 1e50;
  double Lowest = PreviousLowest;
  int CurrentNbrIterLanczos = 0;
  LanczosAlgorithm->SetHamiltonian(this->MPOperator);
  cout << "Run Lanczos Algorithm" << endl;
  timeval TotalStartingTime;
  timeval TotalEndingTime;
  timeval TotalCurrentTime;
  double Dt;
  gettimeofday (&(TotalStartingTime), 0);
  int StartTimeSecond = TotalStartingTime.tv_sec;
  LanczosAlgorithm->RunLanczosAlgorithm(3);
  CurrentNbrIterLanczos = 4;
  RealTriDiagonalSymmetricMatrix TmpMatrix;  
  int CurrentTimeSecond = TotalCurrentTime.tv_sec;

  if ((LanczosAlgorithm->TestConvergence() == true))
  {
      TmpMatrix.Copy(LanczosAlgorithm->GetDiagonalizedMatrix());
      TmpMatrix.SortMatrixUpOrder();
      Lowest = TmpMatrix.DiagonalElement(0);
  }

  while ((LanczosAlgorithm->TestConvergence() == false)&&( CurrentNbrIterLanczos < 2000 ))
  {
     ++CurrentNbrIterLanczos;
     LanczosAlgorithm->RunLanczosAlgorithm(1);
     TmpMatrix.Copy(LanczosAlgorithm->GetDiagonalizedMatrix());
     TmpMatrix.SortMatrixUpOrder();
     Lowest = TmpMatrix.DiagonalElement(0);
     Precision = fabs((PreviousLowest - Lowest) / PreviousLowest);
     PreviousLowest = Lowest; 
     cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << " ";
     gettimeofday (&(TotalEndingTime), 0);
     CurrentTimeSecond = TotalEndingTime.tv_sec;
     Dt = (double) (TotalEndingTime.tv_sec - TotalCurrentTime.tv_sec) + 
		((TotalEndingTime.tv_usec - TotalCurrentTime.tv_usec) / 1000000.0);		      
	      cout << "(" << Dt << " s for step " << CurrentNbrIterLanczos << ")"<< endl;;
     TotalCurrentTime.tv_usec = TotalEndingTime.tv_usec;
     TotalCurrentTime.tv_sec = TotalEndingTime.tv_sec;
   }
      GroundStateEnergy = Lowest;
      cout << endl;
      cout << (TmpMatrix.DiagonalElement(0)) << " " << Lowest << " " << Precision << "  Nbr of iterations = " 
	   << CurrentNbrIterLanczos << endl;
      RealVector GroundState = LanczosAlgorithm->GetGroundState();
      cout <<"Highest energy = " << GroundStateEnergy<<" change = " <<  ((GroundStateEnergy - this->PreviousEnergy) /this->PreviousEnergy)<< endl;
      this->PreviousEnergy = GroundStateEnergy;
      leftSite->SymmetricUpdateOfTwoSites(rightSite, &GroundState,singularValues);
}
}
