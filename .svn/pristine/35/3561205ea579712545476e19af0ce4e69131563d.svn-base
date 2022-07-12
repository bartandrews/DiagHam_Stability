
#include <cmath>
#include "AbstractMPSSite.h"

#include "AbstractMPOperatorOBC.h"
#include "GeneralTools/GarbageFlag.h"

using std::endl;
using std::abs;
using std::cout;

AbstractMPSSite::AbstractMPSSite()
{
  this->PhysicalDimension = 0;
}

AbstractMPSSite::AbstractMPSSite(unsigned int physicalDimension, unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator)
{
  this->PhysicalDimension = physicalDimension;
  this->SquarePhysicalDimension =  this->PhysicalDimension *  this->PhysicalDimension;
  this->OperatorToBeMinimized = mPOperator;
  this->Flag.Initialize();
  this->MaxBondDimension = bondDimension;
}

AbstractMPSSite::~AbstractMPSSite()
{
}



// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

AbstractMPSSite & AbstractMPSSite::operator = (const AbstractMPSSite & site) 
{
  this->PhysicalDimension = site.PhysicalDimension;
  this->SquarePhysicalDimension = site.SquarePhysicalDimension;
  this->OperatorToBeMinimized = site.OperatorToBeMinimized;
  this->BondDimensionLeft = site.BondDimensionLeft;
  this->BondDimensionRight = site.BondDimensionRight;
  this->MaxBondDimension = site.MaxBondDimension;
  
  return *this;
}

/*

void MPSSite::InitializeLeft(RealMatrix * newA)
{
  delete this->M;
  this->M = newA;
  this->BondDimensionLeft = this->M[0].GetNbrRow();
  this->BondDimensionRight = this->M[0].GetNbrColumn();
  delete this->L;
  this->L = new Tensor3<double>(this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
//  cout <<this->BondDimensionLeft<<" "<< this->BondDimensionRight<<" "<<this->PhysicalDimension<<endl;
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
  this->L->PrintTensor();
}
*/

/*

void MPSSite::UpdateFromVector(RealVector * psi)
{
  delete [] this->M;
  this->M = new RealMatrix [this->PhysicalDimension];
  
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      this->M[i] = RealMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    }

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      for (int j = 0; j < this->BondDimensionLeft; j++)
	{
	  for (int k = 0; k < this->BondDimensionRight; k++)
	    {
	      this->M[i](j,k) = (*psi)[(long int)this->BondDimensionRight*(this->BondDimensionLeft*i+j) + k]; 
	    }
	}
    }
}
*/

/*
void MPSSite::InitializeRight(RealMatrix * newB)
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
*/


bool AbstractMPSSite::CheckLeftNormalization()
{
  cout <<"error using low-level AbstractMPSSite::CheckLeftNormalization()" <<endl;
  return true;
}


bool AbstractMPSSite::CheckRightNormalization()
{
  cout <<"error using low-level AbstractMPSSite::CheckRightNormalization()" <<endl;
  return false;
}

//can be used only if all matrices on the right sites are right-normalized
void AbstractMPSSite::BringMInRightCanonicalForm()
{
  cout <<"error using low-level AbstractMPSSite::BringMInRightCanonicalForm()" <<endl;
}



// can be used only if all matrices on the right sites are right-normalized
// check the result

void AbstractMPSSite::BringMInRightCanonicalFormCareful()
{
  this->BringMInRightCanonicalForm();
if(this->CheckRightNormalization())
  ;
  else 
    cout <<"Right Normalization issue in BringMinRightCanonicalFormCareful()"<<endl;
}


// can be used only if all matrices on the left sites are left-normalized

void AbstractMPSSite::BringMInLeftCanonicalForm()
{
 cout <<"error using low-level AbstractMPSSite::BringMInLeftCanonicalForm()" <<endl;
}


// can be used only if all matrices on the left sites are left-normalized
// check the result

void AbstractMPSSite::BringMInLeftCanonicalFormCareful()
{
  this->BringMInLeftCanonicalForm();
  if(this->CheckLeftNormalization())
    ;
  else 
    cout <<"Left Normalization issue in BringMinLeftCanonicalFormCareful()"<<endl;
}


/*
void MPSSite::GetMatrixInVectorForm(RealVector *& resultInvector)
{
  delete resultInvector;
  resultInvector = new RealVector((long int) this->PhysicalDimension*this->BondDimensionLeft*this->BondDimensionRight,true);

  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      for (int j = 0 ; j <  this->BondDimensionLeft ; j++)
	{
	  for (int k = 0 ; k <  this->BondDimensionRight ; k++)
	    {
	      (*resultInvector)[(long int)this->BondDimensionRight*(this->BondDimensionLeft*i+j) + k] = this->M[i](j,k); 
	    }
	}
    }
}
*/

/*
void MPSSite::InitializeWithRandomMatrices()
{
//   cout <<" start initialization i = "<<this->SitePosition <<endl;;
//   cout <<this->BondDimensionLeft<< " " << this->BondDimensionRight<<endl;
   for (int i = 0; i < this->PhysicalDimension; i++)
   {
    this->M[i] = RealMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    for (int j = 0; j < this->BondDimensionLeft ; j++)
     for (int k = 0; k < this->BondDimensionRight ; k++)
      this->M[i](j,k) = ((double) rand() / (RAND_MAX) - 0.5);
   }
}
*/


/*
void MPSSite::ComputeDensityMatrixLeft()
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
     TmpMatrix.LapackDiagonalize(TmpDiag);
     TmpDiag.SortMatrixDownOrder();
}


void MPSSite::ComputeDensityMatrixRight()
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
     TmpMatrix.LapackDiagonalize(TmpDiag);
     TmpDiag.SortMatrixDownOrder();

}


void MPSSite::SymmetricUpdateOfTwoSites(MPSSite * leftSite , MPSSite * rightSite, RealVector * psi, RealDiagonalMatrix & SingularValues)
{
  RealMatrix TmpMatrix (leftSite->BondDimensionLeft *this->PhysicalDimension , rightSite->BondDimensionRight * this->PhysicalDimension, true);
  // Index = LinearizedPhysicalIndice + SquarePhysicalDimension*(LeftA + RightA*BondDimensionLeft))

  for(int i = 0; i < leftSite->BondDimensionLeft; i++)
    {
      for(int LeftPhysicalDimension = 0; LeftPhysicalDimension <   this->PhysicalDimension ; LeftPhysicalDimension++)
	{
	  for(int j = 0; j < rightSite->BondDimensionRight; j++)
	  { 
      for(int RightPhysicalDimension = 0; RightPhysicalDimension <   this->PhysicalDimension ; RightPhysicalDimension++)
	{
	      TmpMatrix(this->PhysicalDimension*i + LeftPhysicalDimension , this->PhysicalDimension*j + RightPhysicalDimension) = (*psi)[(long) LeftPhysicalDimension+ this->PhysicalDimension*RightPhysicalDimension + this->PhysicalDimension *this->PhysicalDimension * (i + j * leftSite->BondDimensionLeft)];
	    }
	}
    }
  }

  RealMatrix U,V;
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
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

 delete [] leftSite->M;
 delete [] rightSite->M;
 delete leftSite->L;
 delete rightSite->R;

 leftSite->M = new RealMatrix [this->PhysicalDimension];
 rightSite->M = new RealMatrix [this->PhysicalDimension];  


  leftSite->SetRightDimension(KeptStates);
  rightSite->SetLeftDimension(KeptStates);
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      leftSite->M[i] = RealMatrix(leftSite->BondDimensionLeft,KeptStates, true);
      rightSite->M[i] = RealMatrix(KeptStates, rightSite->BondDimensionRight, true);
    } 

  for(int  LeftPhysicalDimension = 0 ;  LeftPhysicalDimension <  this->PhysicalDimension;  LeftPhysicalDimension++)
    {
      for(int j = 0 ; j < leftSite->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  KeptStates ; k++)
	    {
	      leftSite->M[LeftPhysicalDimension](j,k) = U(this->PhysicalDimension*j +  LeftPhysicalDimension,k);
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
  leftSite->L = new Tensor3<double> (leftSite->BondDimensionRight,leftSite->OperatorToBeMinimized->GetMPODimension(),leftSite->BondDimensionRight,true);
  leftSite->OperatorToBeMinimized->SetSite(leftSite);
  leftSite->OperatorToBeMinimized->ComputeL(*leftSite->L);


  rightSite->R = new Tensor3<double> (rightSite->BondDimensionLeft,rightSite->OperatorToBeMinimized->GetMPODimension(),rightSite->BondDimensionLeft,true);
  rightSite->OperatorToBeMinimized->SetSite(rightSite);
  rightSite->OperatorToBeMinimized->ComputeR(*rightSite->R);
}
*/
