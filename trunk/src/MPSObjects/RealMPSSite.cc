#include <cmath>
#include "RealMPSSite.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "AbstractMPOperatorOBC.h"
#include "GeneralTools/GarbageFlag.h"

using std::endl;
using std::abs;
using std::cout;

RealMPSSite::RealMPSSite()
{
  this->SiteOnLeft = 0;
  this->SiteOnRight = 0;
  this->M = 0;
  this->L = 0;
  this->R = 0;
}

RealMPSSite::RealMPSSite(unsigned int physicalDimension, RealMPSSite * siteOnLeft, RealMPSSite * siteOnRight , unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator) : AbstractMPSSite ( physicalDimension, bondDimension, mPOperator)
{
  this->SiteOnLeft = siteOnLeft;
  this->SiteOnRight = siteOnRight;
  this->Flag.Initialize();
  this->M = new RealMatrix[this->PhysicalDimension];
  this->L = 0;
  this->R = 0;
}

RealMPSSite::~RealMPSSite()
{
  delete [] this->M;
  delete this->L;
  delete this->R;
}



// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

RealMPSSite & RealMPSSite::operator = (const RealMPSSite & site) 
{
  AbstractMPSSite::operator=(site);
  this->SiteOnLeft = site.SiteOnLeft;
  this->SiteOnRight = site.SiteOnRight;
  delete this->M;
  this->M = new RealMatrix[this->PhysicalDimension];
  for(int i= 0; i <this->PhysicalDimension; i++)
    this->M[i] = site.M[i];
  this->L = site.L;
  this->R = site.R;
  return *this;
}


void RealMPSSite::InitializeLeft(RealMatrix * newA)
{
  delete this->M;
  this->M = newA;
  this->BondDimensionLeft = this->M[0].GetNbrRow();
  this->BondDimensionRight = this->M[0].GetNbrColumn();
  delete this->L;
  this->L = new Tensor3<double>(this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
  this->L->PrintTensor();
}


void RealMPSSite::UpdateFromVector(RealVector * psi)
{
  delete [] this->M;
  this->M = new RealMatrix [this->PhysicalDimension];
  
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      this->M[i] = RealMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    }

  for (int physicalIndice = 0; physicalIndice < this->PhysicalDimension; physicalIndice++)
    {
      for (int LeftIndice = 0; LeftIndice < this->BondDimensionLeft; LeftIndice++)
	{
	  for (int RightIndice = 0; RightIndice < this->BondDimensionRight; RightIndice++)
	    {
	      this->M[physicalIndice](LeftIndice,RightIndice) = (*psi)[this->GetVectorOneSiteIndice(LeftIndice,RightIndice,physicalIndice)]; 
	    }
	}
    }
}


void RealMPSSite::InitializeRight(RealMatrix * newB)
{
  delete this->M;
  this->M = newB;
  this->BondDimensionLeft = this->M[0].GetNbrRow();
  this->BondDimensionRight = this->M[0].GetNbrColumn();
  delete this->R;
  this->R = new Tensor3<double>(this->BondDimensionLeft,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionLeft,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeR(*this->R);
  this->R->PrintTensor();
}

bool RealMPSSite::CheckLeftNormalization()
{
  RealMatrix Result(this->BondDimensionRight,this->BondDimensionRight,true);
  for (int i = 0 ; i< this->PhysicalDimension ; i++)
    {
      RealMatrix Tmp = this->M[i].DuplicateAndTranspose();
      Result += Tmp * (this->M[i]);
    }
  for(int i = 0 ; i < this->BondDimensionRight; i++)
    {
      if ( abs(Result(i,i)-1) > 1e-13  )
	  {
	    return false;
	  }
      for(int j = i+1 ; j < this->BondDimensionRight ; j++)
	{
	  if (abs(Result(i,j)) > 1e-13)
	    {
	      return false;
	    }
	}
    }
  return true;
}


bool RealMPSSite::CheckRightNormalization()
{
  RealMatrix Result(this->BondDimensionLeft,this->BondDimensionLeft,true);
  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      RealMatrix Tmp = this->M[i].DuplicateAndTranspose();
      Result += (this->M[i])* Tmp;
    }
  for(int i = 0 ; i < this->BondDimensionLeft; i++)
    {
      if ( abs(Result(i,i)-1) > 1e-13  )
	{
	  return false;
	}
      for(int j = i+1;j< this->BondDimensionLeft; j++)
	{
	  if (abs(Result(i,j)) > 1e-13)
	    {
	      return false;
	    }
	}
      
    }
  return true;
}

//can be used only if all matrices on the right sites are right-normalized
void RealMPSSite::BringMInRightCanonicalForm()
{
// cout << "BringMInRightCanonicalForm() for site" << this->SitePosition<<endl;
 RealMatrix TmpMatrix (this->BondDimensionLeft,this->BondDimensionRight * this->PhysicalDimension, true);
  
  for(int i = 0 ; i < this->BondDimensionLeft ; i++)
    {
      for(int p = 0 ; p < this->PhysicalDimension ; p++)
	{
	  for(int j = 0; j < this->BondDimensionRight; j++)
	    { 
	      TmpMatrix(i, this->PhysicalDimension * j + p) = this->M[p](i,j);
	    }
	}
    }
//  cout << TmpMatrix<<endl;
  RealMatrix U,V;
  RealDiagonalMatrix SingularValues;

  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
//  cout <<"SingularValues  = "<< SingularValues<<endl;
//  cout <<U <<endl;
// cout <<V <<endl;
  U = U * SingularValues;
  double Entropy = 0.0;
  for(int i = 0; i < SingularValues.GetNbrRow(); i++)
{
   SingularValues[i]*=SingularValues[i];
   Entropy -= SingularValues[i]*log(SingularValues[i]);
}

//  cout <<U<<endl;;
 cout <<"Entropy = "<< Entropy<<" " << SingularValues[SingularValues.GetNbrRow()-1]<< endl;
  for(int i =0 ; i <  this->PhysicalDimension; i++)
    {
      for(int j = 0 ; j <  this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = V(j,this->PhysicalDimension * k + i);
	    }
	}

       ((RealMPSSite * )this->SiteOnLeft)->M[i] = ((RealMPSSite * )this->SiteOnLeft)->M[i]*U;
       }
  delete this->R;
  this->R = new Tensor3<double> (this->BondDimensionLeft,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionLeft,true) ;
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeR(*this->R);
}



// can be used only if all matrices on the left sites are left-normalized

void RealMPSSite::BringMInLeftCanonicalForm()
{
  RealMatrix TmpMatrix (this->PhysicalDimension*this->BondDimensionLeft,this->BondDimensionRight, true);
  
  for(int i = 0; i < this->BondDimensionLeft; i++)
    {
      for(int p = 0; p <   this->PhysicalDimension ; p++)
	{
	  for(int j = 0; j < this->BondDimensionRight; j++)
	    { 
	      TmpMatrix(this->PhysicalDimension*i + p,j) = this->M[p](i,j);
	    }
	}
    }
  
  RealMatrix U,V;
  RealDiagonalMatrix SingularValues;
//  cout <<"check singular value decomposition in void MPSSite::BringMInLeftCanonicalForm()"<<endl;
//  cout <<  TmpMatrix<<endl;
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
//  cout << SingularValues<<endl;
//  cout <<U<<endl;
//  cout <<V<<endl;
//  V.Transpose();
  V = SingularValues*V;
//  cout <<V<<endl;
double Entropy=0.0;
  for(int i = 0; i < SingularValues.GetNbrRow(); i++)
{
   SingularValues[i]*=SingularValues[i];
   Entropy -= SingularValues[i]*log(SingularValues[i]);
}
//  cout <<U<<endl;;
 cout <<"Entropy = "<< Entropy<<" " << SingularValues[SingularValues.GetNbrRow()-1]<< endl;

  for(int i = 0 ; i <  this->PhysicalDimension; i++)
    {
      for(int j = 0 ; j <  this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = U(this->PhysicalDimension*j + i,k);
	    }
	}
      ((RealMPSSite * )this->SiteOnRight)->M[i] =  V * ((RealMPSSite * )this->SiteOnRight)->M[i];
    }

  delete this->L;
  this->L = new Tensor3<double> (this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
}


void RealMPSSite::GetMatrixInVectorForm(RealVector *& resultInvector)
{
  delete resultInvector;
  resultInvector = new RealVector((long int) this->PhysicalDimension*this->BondDimensionLeft*this->BondDimensionRight,true);

  for (int PhysicalIndice = 0 ; PhysicalIndice < this->PhysicalDimension ; PhysicalIndice++)
    {
      for (int LeftIndice = 0 ; LeftIndice <  this->BondDimensionLeft ; LeftIndice++)
	{
	  for (int RightIndice = 0 ; RightIndice <  this->BondDimensionRight ; RightIndice++)
	    {
	      (*resultInvector)[this->GetVectorOneSiteIndice(LeftIndice,RightIndice,PhysicalIndice)] = this->M[PhysicalIndice](LeftIndice,RightIndice); 
	    }
	}
    }
}


void RealMPSSite::InitializeWithRandomMatrices()
{
   for (int i = 0; i < this->PhysicalDimension; i++)
   {
    this->M[i] = RealMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    for (int j = 0; j < this->BondDimensionLeft ; j++)
     for (int k = 0; k < this->BondDimensionRight ; k++)
      this->M[i](j,k) = ((double) rand() / (RAND_MAX) - 0.5);
   }
}



void RealMPSSite::ComputeDensityMatrixLeft()
{
   RealSymmetricMatrix TmpMatrix (this->PhysicalDimension*this->BondDimensionLeft, true);

   for (int i = 0; i < this->PhysicalDimension; i++)
   {
    for (int j = 0; j < this->BondDimensionLeft ; j++)
    {
   	for (int k = 0; k < this->PhysicalDimension; k++)
       {
           for (int l = 0; l < this->BondDimensionLeft ; l++)
	   {
  
           for (int t = 0; t < this->BondDimensionRight ; t++)
	   {
              TmpMatrix.AddToMatrixElement(j*this->PhysicalDimension+i,l*this->PhysicalDimension+k,this->M[i](j,t)*this->M[k](l,t));
           }
    }
   }
  }
  }
     RealDiagonalMatrix TmpDiag (TmpMatrix.GetNbrRow());
#ifdef __LAPACK__     
     TmpMatrix.LapackDiagonalize(TmpDiag);
#else
     TmpMatrix.Diagonalize(TmpDiag);
#endif
     TmpDiag.SortMatrixDownOrder();
}


void RealMPSSite::ComputeDensityMatrixRight()
{
   RealSymmetricMatrix TmpMatrix (this->PhysicalDimension*this->BondDimensionRight, true);

   for (int i = 0; i < this->PhysicalDimension; i++)
   {
    for (int j = 0; j < this->BondDimensionRight ; j++)
    {
   	for (int k = 0; k < this->PhysicalDimension; k++)
       {
           for (int l = 0; l < this->BondDimensionRight ; l++)
	   {
  
           for (int t = 0; t < this->BondDimensionLeft ; t++)
	   {
              TmpMatrix.AddToMatrixElement(j*this->PhysicalDimension+i,l*this->PhysicalDimension+k,this->M[i](t,j)*this->M[k](t,l));
           }
   }
   }
  }
 }
     RealDiagonalMatrix TmpDiag (TmpMatrix.GetNbrRow());
#ifdef __LAPACK__     
     TmpMatrix.LapackDiagonalize(TmpDiag);
#else
     TmpMatrix.Diagonalize(TmpDiag);
#endif
     TmpDiag.SortMatrixDownOrder();
}


void RealMPSSite::SymmetricUpdateOfTwoSites(RealMPSSite * rightSite, RealVector * psi, RealDiagonalMatrix & SingularValues)
{
  RealMatrix TmpMatrix (this->BondDimensionLeft *this->PhysicalDimension , rightSite->BondDimensionRight * this->PhysicalDimension, true);

  // Index = LinearizedPhysicalIndice + SquarePhysicalDimension*(LeftA + RightA*BondDimensionLeft))

  for(int LeftIndice = 0; LeftIndice < this->BondDimensionLeft; LeftIndice++)
    {
      for(int LeftPhysicalDimension = 0; LeftPhysicalDimension <   this->PhysicalDimension ; LeftPhysicalDimension++)
	{
	  for(int RightIndice = 0; RightIndice < rightSite->BondDimensionRight; RightIndice++)
	  { 
      for(int RightPhysicalDimension = 0; RightPhysicalDimension <   this->PhysicalDimension ; RightPhysicalDimension++)
	{
	      TmpMatrix(this->PhysicalDimension*LeftIndice + LeftPhysicalDimension , this->PhysicalDimension*RightIndice + RightPhysicalDimension) = (*psi)[this->GetVectorTwoSiteIndice(LeftIndice,RightIndice, LeftPhysicalDimension +this->PhysicalDimension*RightPhysicalDimension)];
	    }
	}
    }
  }

  RealMatrix U,V;
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
  RealMatrix Tmp = U * SingularValues * V;


  double Entropy=0.0;
  unsigned int KeptStates = 0;


  for(int i = 0; i < SingularValues.GetNbrRow(); i++)
  {
   if (SingularValues[i] > 1e-20)
       KeptStates++;
  }

  if ( KeptStates >  this->MaxBondDimension)
       KeptStates = this->MaxBondDimension;
 double * KeptSingularValues = new double[KeptStates];

for(int i = 0; i < KeptStates; i++)
  {
   KeptSingularValues[i] = SingularValues[i];
   SingularValues[i]*=SingularValues[i];
   Entropy -= SingularValues[i]*log(SingularValues[i]);
  }
 double RejectedWeight = 0.0;
 for(int i = KeptStates; i < SingularValues.GetNbrRow(); i++)
  {
   RejectedWeight +=SingularValues[i];
  }

 SingularValues = RealDiagonalMatrix(KeptSingularValues,KeptStates);
 cout <<"Entropy = "<< Entropy<<" " <<  RejectedWeight<< endl;

 delete [] this->M;
 delete [] rightSite->M;
 delete this->L;
 delete rightSite->R;

 this->M = new RealMatrix [this->PhysicalDimension];
 rightSite->M = new RealMatrix [this->PhysicalDimension];  


  this->SetRightDimension(KeptStates);
  rightSite->SetLeftDimension(KeptStates);
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      this->M[i] = RealMatrix(this->BondDimensionLeft,KeptStates, true);
      rightSite->M[i] = RealMatrix(KeptStates, rightSite->BondDimensionRight, true);
    } 

  for(int  LeftPhysicalDimension = 0 ;  LeftPhysicalDimension <  this->PhysicalDimension;  LeftPhysicalDimension++)
    {
      for(int j = 0 ; j < this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  KeptStates ; k++)
	    {
	      this->M[LeftPhysicalDimension](j,k) = U(this->PhysicalDimension*j +  LeftPhysicalDimension,k);
	    }
	}
    }

  for(int  RightPhysicalDimension = 0 ;  RightPhysicalDimension <  this->PhysicalDimension;  RightPhysicalDimension++)
    {
      for(int j = 0 ; j < rightSite->BondDimensionRight; j++)
	{
	  for(int k = 0 ; k <  KeptStates ; k++)
	    {
	      rightSite->M[RightPhysicalDimension](k,j) = V(k,this->PhysicalDimension*j +  RightPhysicalDimension);
	    }
	}
    }

  this->L = new Tensor3<double> (this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);

  rightSite->R = new Tensor3<double> (rightSite->BondDimensionLeft,rightSite->OperatorToBeMinimized->GetMPODimension(),rightSite->BondDimensionLeft,true);
  rightSite->OperatorToBeMinimized->SetSite(rightSite);
  rightSite->OperatorToBeMinimized->ComputeR(*rightSite->R);
}



RealVector *  RealMPSSite::StatePrediction(RealMPSSite * rightSite, RealDiagonalMatrix & SingularValues, RealDiagonalMatrix & OldSingularValues)
{
 RealMatrix * TmpA = new RealMatrix [this->PhysicalDimension];
 RealMatrix * TmpB = new RealMatrix [this->PhysicalDimension];  

 for(int i = 0; i < this->PhysicalDimension; i++)
 {
     TmpA[i] = this->M[i] * SingularValues;
     TmpB[i] = SingularValues * rightSite->M[i] / OldSingularValues;
}

 RealVector * PredictedPsi = new RealVector((long) this->BondDimensionRight*rightSite->BondDimensionLeft*this->PhysicalDimension*this->PhysicalDimension ,true);

  for(int LeftIndice = 0; LeftIndice < this->BondDimensionRight; LeftIndice++)
    {
      for(int LeftPhysicalDimension = 0; LeftPhysicalDimension <   this->PhysicalDimension ; LeftPhysicalDimension++)
	{
	  for(int RightIndice = 0; RightIndice < rightSite->BondDimensionLeft; RightIndice++)
	  { 
      for(int RightPhysicalDimension = 0; RightPhysicalDimension <   this->PhysicalDimension ; RightPhysicalDimension++)
	{

      for(int k = 0; k <   OldSingularValues.GetNbrColumn() ; k++)
	{

//             (*PredictedPsi)[(long) LeftPhysicalDimension+ this->PhysicalDimension*RightPhysicalDimension + this->PhysicalDimension *this->PhysicalDimension * (LeftIndice + j * this->BondDimensionRight)] =   TmpB[LeftPhysicalDimension](LeftIndice,k) *  TmpA[RightPhysicalDimension](k ,RightIndice);
             (*PredictedPsi)[this->SiteOnRight->GetVectorTwoSiteIndice(LeftIndice,RightIndice, LeftPhysicalDimension +this->PhysicalDimension*RightPhysicalDimension)] =   TmpB[LeftPhysicalDimension](LeftIndice,k) *  TmpA[RightPhysicalDimension](k ,RightIndice);
	 }
	}
    }
  }
}

  return PredictedPsi;
}
