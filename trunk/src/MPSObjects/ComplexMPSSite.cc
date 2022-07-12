#include <cmath>
#include "ComplexMPSSite.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "AbstractMPOperatorOBC.h"
#include "GeneralTools/GarbageFlag.h"

using std::endl;
using std::abs;
using std::cout;

ComplexMPSSite::ComplexMPSSite()
{
  this->SiteOnLeft = 0;
  this->SiteOnRight = 0;
  this->M = 0;
  this->L = 0;
  this->R = 0;
}

ComplexMPSSite::ComplexMPSSite(unsigned int physicalDimension, ComplexMPSSite * siteOnLeft, ComplexMPSSite * siteOnRight , unsigned int bondDimension, AbstractMPOperatorOBC * mPOperator) : AbstractMPSSite (physicalDimension, bondDimension, mPOperator)
{
  this->SiteOnLeft = siteOnLeft;
  this->SiteOnRight = siteOnRight;
  this->Flag.Initialize();
  this->M = new ComplexMatrix[this->PhysicalDimension];
  this->L = 0;
  this->R = 0;
}


ComplexMPSSite::~ComplexMPSSite()
{
  delete [] this->M;
  delete this->L;
  delete this->R;
}



// assignement (without duplicating datas)
//
// M = matrix to copy
// return value = reference on modified matrix

ComplexMPSSite & ComplexMPSSite::operator = (const ComplexMPSSite & site) 
{
  AbstractMPSSite::operator=(site);
  this->SiteOnLeft = site.SiteOnLeft;
  this->SiteOnRight = site.SiteOnRight;
  delete this->M;
  this->M = new ComplexMatrix[this->PhysicalDimension];
  for(int i= 0; i <this->PhysicalDimension; i++)
    this->M[i] = site.M[i];
  this->L = site.L;
  this->R = site.R;
  return *this;
}


void ComplexMPSSite::InitializeLeft(ComplexMatrix * newA)
{
  delete this->M;
  this->M = newA;
  this->BondDimensionLeft = this->M[0].GetNbrRow();
  this->BondDimensionRight = this->M[0].GetNbrColumn();
  delete this->L;
  this->L = new Tensor3<Complex>(this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
  this->L->PrintTensor();
}


void ComplexMPSSite::UpdateFromVector(ComplexVector * psi)
{
  delete [] this->M;
  this->M = new ComplexMatrix [this->PhysicalDimension];
  
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      this->M[i] = ComplexMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
    }

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      for (int j = 0; j < this->BondDimensionLeft; j++)
	{
	  for (int k = 0; k < this->BondDimensionRight; k++)
	    {
	      this->M[i].SetMatrixElement(j,k, (*psi)[(long int)this->BondDimensionRight*(this->BondDimensionLeft*i+j) + k]); 
	    }
	}
    }
}


void ComplexMPSSite::InitializeRight(ComplexMatrix * newB)
{
  delete this->M;
  this->M = newB;
  this->BondDimensionLeft = this->M[0].GetNbrRow();
  this->BondDimensionRight = this->M[0].GetNbrColumn();
  delete this->R;
  this->R = new Tensor3<Complex>(this->BondDimensionLeft,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionLeft,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeR(*this->R);
  this->R->PrintTensor();
}

bool ComplexMPSSite::CheckLeftNormalization()
{
  ComplexMatrix Result(this->BondDimensionRight,this->BondDimensionRight,true);
  for (int i = 0 ; i< this->PhysicalDimension ; i++)
    {
      ComplexMatrix Tmp = this->M[i].GetAdjoint(); //CHECK!!
      Result += Tmp * (this->M[i]);
    }
  for(int i = 0 ; i < this->BondDimensionRight; i++)
    {
      if ( abs(Norm(Result.GetMatrixElement(i,i)) - 1.0) > 1e-13  )
	  {
	    return false;
	  }
      for(int j = i+1 ; j < this->BondDimensionRight ; j++)
	{
	  if ( abs(Norm(Result.GetMatrixElement(i,j))) > 1e-13)
	    {
	      return false;
	    }
	}
    }
  return true;
}


bool ComplexMPSSite::CheckRightNormalization()
{
  ComplexMatrix Result(this->BondDimensionLeft,this->BondDimensionLeft,true);
  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      ComplexMatrix Tmp = this->M[i].GetAdjoint(); //CHECK
      Result += (this->M[i])* Tmp;
    }
  for(int i = 0 ; i < this->BondDimensionLeft; i++)
    {
      if ( abs(Norm(Result.GetMatrixElement(i,i)) - 1.0) > 1e-13  )
	{
	  return false;
	}
      for(int j = i+1;j< this->BondDimensionLeft; j++)
	{
  	  if ( abs(Norm(Result.GetMatrixElement(i,j))) > 1e-13)
	    {
	      return false;
	    }
	}
      
    }
  return true;
}

//can be used only if all matrices on the right sites are right-normalized
void ComplexMPSSite::BringMInRightCanonicalForm()
{
  
  ComplexMatrix TmpMatrix (this->BondDimensionLeft,this->BondDimensionRight * this->PhysicalDimension, true);
  
  for(int i = 0 ; i < this->BondDimensionLeft ; i++)
    {
      for(int p = 0 ; p < this->PhysicalDimension ; p++)
	{
	  for(int j = 0; j < this->BondDimensionRight; j++)
	    { 
	      TmpMatrix.SetMatrixElement(i, this->PhysicalDimension * j + p,this->M[p].GetMatrixElement(i,j));
	    }
	}
    }

  //  cout << TmpMatrix<<endl;
  ComplexMatrix U,V;
  RealDiagonalMatrix SingularValues;
  
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
  //  cout <<"SingularValues  = "<< SingularValues<<endl;
  //  cout <<"U = "<<U <<endl;
  //  cout <<"V = "<<V <<endl;
  U = U * SingularValues;
  //  cout <<" U * SingularValues = "<<U <<endl;
  double Entropy = 0.0;
  for(int i = 0; i < SingularValues.GetNbrRow(); i++)
    {
      SingularValues[i]*=SingularValues[i];
      Entropy -= SingularValues[i]*log(SingularValues[i]);
    }
  
  
  cout <<"Entropy = "<< Entropy<<" " << SingularValues[SingularValues.GetNbrRow()-1]<< endl;
  for(int i =0 ; i <  this->PhysicalDimension; i++)
    {
      for(int j = 0 ; j <  this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  this->BondDimensionRight; k++)
	    {
	      this->M[i].SetMatrixElement(j,k, V.GetMatrixElement(j,this->PhysicalDimension * k + i));
	    }
	}
      
      (((ComplexMPSSite * )this->SiteOnLeft)->M[i]) = (((ComplexMPSSite * )this->SiteOnLeft)->M[i]) *U;
    }
  delete this->R;
  this->R = new Tensor3<Complex> (this->BondDimensionLeft,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionLeft,true) ;
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeR(*this->R);
}



// can be used only if all matrices on the left sites are left-normalized

void ComplexMPSSite::BringMInLeftCanonicalForm()
{
  ComplexMatrix TmpMatrix (this->PhysicalDimension*this->BondDimensionLeft,this->BondDimensionRight, true);
  
  for(int i = 0; i < this->BondDimensionLeft; i++)
    {
      for(int p = 0; p <   this->PhysicalDimension ; p++)
	{
	  for(int j = 0; j < this->BondDimensionRight; j++)
	    { 
	      TmpMatrix.SetMatrixElement(this->PhysicalDimension*i + p,j,this->M[p].GetMatrixElement(i,j));
	    }
	}
    }
  
  ComplexMatrix U,V;
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
	      this->M[i].SetMatrixElement(j,k,U.GetMatrixElement(this->PhysicalDimension*j + i,k));
	    }
	}
      ((ComplexMPSSite * )this->SiteOnRight)->M[i] =  V * ((ComplexMPSSite * )this->SiteOnRight)->M[i];
    }

  delete this->L;
  this->L = new Tensor3<Complex> (this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
}


void ComplexMPSSite::GetMatrixInVectorForm(ComplexVector *& resultInvector)
{
//  cout <<"in GetMatrixInVectorForm(ComplexVector *& resultInvector)"<<endl;
  delete resultInvector;

  resultInvector = new ComplexVector((long int) this->PhysicalDimension*this->BondDimensionLeft*this->BondDimensionRight,true);

  for (int i = 0 ; i < this->PhysicalDimension ; i++)
    {
      for (int j = 0 ; j <  this->BondDimensionLeft ; j++)
	{
	  for (int k = 0 ; k <  this->BondDimensionRight ; k++)
	    {
	      (*resultInvector)[(long int)this->BondDimensionRight*(this->BondDimensionLeft*i+j) + k] = this->M[i].GetMatrixElement(j,k); 
	    }
	}
    }
}


void ComplexMPSSite::InitializeWithRandomMatrices()
{
   for (int i = 0; i < this->PhysicalDimension; i++)
     {
       this->M[i] = ComplexMatrix(this->BondDimensionLeft,this->BondDimensionRight, true);
       for (int j = 0; j < this->BondDimensionLeft ; j++)
	 for (int k = 0; k < this->BondDimensionRight ; k++)
	   {
	     this->M[i].SetMatrixElement(j,k,((double) rand() / (RAND_MAX) - 0.5),((double) rand() / (RAND_MAX) - 0.5));
	   }
     }
}



void ComplexMPSSite::ComputeDensityMatrixLeft()
{
  HermitianMatrix TmpMatrix (this->PhysicalDimension*this->BondDimensionLeft, true);
  
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
		      TmpMatrix.AddToMatrixElement(j*this->PhysicalDimension+i,l*this->PhysicalDimension+k, Conj(this->M[i](j,t)) * this->M[k](l,t));
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


void ComplexMPSSite::ComputeDensityMatrixRight()
{
   HermitianMatrix TmpMatrix (this->PhysicalDimension*this->BondDimensionRight, true);

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
		       TmpMatrix.AddToMatrixElement(j*this->PhysicalDimension+i,l*this->PhysicalDimension+k,Conj(this->M[i](t,j))*this->M[k](t,l));
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


void ComplexMPSSite::SymmetricUpdateOfTwoSites(ComplexMPSSite * rightSite, ComplexVector * psi, RealDiagonalMatrix & SingularValues)
{
  ComplexMatrix TmpMatrix (this->BondDimensionLeft *this->PhysicalDimension , rightSite->BondDimensionRight * this->PhysicalDimension, true);
  
  for(int i = 0; i < this->BondDimensionLeft; i++)
    {
      for(int LeftPhysicalDimension = 0; LeftPhysicalDimension <   this->PhysicalDimension ; LeftPhysicalDimension++)
	{
	  for(int j = 0; j < rightSite->BondDimensionRight; j++)
	    { 
	      for(int RightPhysicalDimension = 0; RightPhysicalDimension <   this->PhysicalDimension ; RightPhysicalDimension++)
		{
		  TmpMatrix.SetMatrixElement(this->PhysicalDimension*i + LeftPhysicalDimension , this->PhysicalDimension*j + RightPhysicalDimension , (*psi)[(long) LeftPhysicalDimension+ this->PhysicalDimension*RightPhysicalDimension + this->PhysicalDimension *this->PhysicalDimension * (i + j * this->BondDimensionLeft)]);
		}
	    }
	}
    }
  
  ComplexMatrix U,V;
  TmpMatrix.SingularValueDecomposition(U,SingularValues,V,false);
  double Entropy=0.0;
  unsigned int KeptStates = 0;
  double TotalWeigth = 0.0;
  for(int i = 0; i < SingularValues.GetNbrRow(); i++)
    {
      TotalWeigth+=(SingularValues[i]*SingularValues[i]);
      cout <<SingularValues[i]<<" ";
      if (SingularValues[i] > 1e-20)
	KeptStates++;
    }
  cout <<endl;  
  cout <<"Total weigth = "<<  TotalWeigth <<endl;
  if ( KeptStates >  this->MaxBondDimension)
    KeptStates = this->MaxBondDimension;
  double * KeptSingularValues = new double[KeptStates];
  
  for(int i = 0; i < KeptStates; i++)
    {
      KeptSingularValues[i] = SingularValues[i];
      SingularValues[i] *= SingularValues[i];
      Entropy -= SingularValues[i]*log(SingularValues[i]);
    }
  double RejectedWeight = 0.0;
  for(int i = KeptStates; i < SingularValues.GetNbrRow(); i++)
    {
      SingularValues[i] *= SingularValues[i];
      cout <<SingularValues[i]<<" ";
      RejectedWeight +=SingularValues[i];
    }
  cout <<endl;
  SingularValues = RealDiagonalMatrix(KeptSingularValues,KeptStates);
  cout <<"Entropy = "<< Entropy<<" " <<RejectedWeight<< endl;
  
  delete [] this->M;
  delete [] rightSite->M;
  delete this->L;
  delete rightSite->R;
  
  this->M = new ComplexMatrix [this->PhysicalDimension];
  rightSite->M = new ComplexMatrix [this->PhysicalDimension];  
  
  
  this->SetRightDimension(KeptStates);
  rightSite->SetLeftDimension(KeptStates);
  for(int i = 0; i < this->PhysicalDimension; i++)
    {
      this->M[i] = ComplexMatrix(this->BondDimensionLeft,KeptStates, true);
      rightSite->M[i] = ComplexMatrix(KeptStates, rightSite->BondDimensionRight, true);
    } 
  
  for(int  LeftPhysicalDimension = 0 ;  LeftPhysicalDimension <  this->PhysicalDimension;  LeftPhysicalDimension++)
    {
      for(int j = 0 ; j < this->BondDimensionLeft; j++)
	{
	  for(int k = 0 ; k <  KeptStates ; k++)
	    {
              Complex Tmp = U.GetMatrixElement(this->PhysicalDimension*j +  LeftPhysicalDimension,k);
	      //	      this->M[LeftPhysicalDimension].SetMatrixElement(j,k, U.GetMatrixElement(this->PhysicalDimension*j +  LeftPhysicalDimension,k)); 
	      this->M[LeftPhysicalDimension].SetMatrixElement(j,k,Tmp); 
	      
	    }
	}
    }
  
  for(int  RightPhysicalDimension = 0 ;  RightPhysicalDimension <  this->PhysicalDimension;  RightPhysicalDimension++)
    {
      for(int j = 0 ; j < rightSite->BondDimensionRight; j++)
	{
	  for(int k = 0 ; k <  KeptStates ; k++)
	    {
              Complex Tmp = V.GetMatrixElement(k,this->PhysicalDimension*j +  RightPhysicalDimension);	
	      rightSite->M[RightPhysicalDimension].SetMatrixElement(k,j, Tmp);
	      //	      rightSite->M[RightPhysicalDimension].SetMatrixElement(k,j, V.GetMatrixElement(k,this->PhysicalDimension*j +  RightPhysicalDimension));
	    }
	}
    }
  
  if( this->CheckLeftNormalization() == false)
    cout <<"left noralisation issue invoid ComplexMPSSite::SymmetricUpdateOfTwoSites(ComplexMPSSite * leftSite , ComplexMPSSite * rightSite, ComplexVector * psi, RealDiagonalMatrix & SingularValues)" <<endl;
  
  
  if( rightSite->CheckRightNormalization() == false)
    cout <<"right normalisation issue invoid ComplexMPSSite::SymmetricUpdateOfTwoSites(ComplexMPSSite * leftSite , ComplexMPSSite * rightSite, ComplexVector * psi, RealDiagonalMatrix & SingularValues)" <<endl;
  
  this->L = new Tensor3<Complex> (this->BondDimensionRight,this->OperatorToBeMinimized->GetMPODimension(),this->BondDimensionRight,true);
  this->OperatorToBeMinimized->SetSite(this);
  this->OperatorToBeMinimized->ComputeL(*this->L);
  
  rightSite->R = new Tensor3<Complex> (rightSite->BondDimensionLeft,rightSite->OperatorToBeMinimized->GetMPODimension(),rightSite->BondDimensionLeft,true);
  rightSite->OperatorToBeMinimized->SetSite(rightSite);
  rightSite->OperatorToBeMinimized->ComputeR(*rightSite->R);
}




ComplexVector *  ComplexMPSSite::StatePrediction(ComplexMPSSite * rightSite, RealDiagonalMatrix & SingularValues, RealDiagonalMatrix & OldSingularValues)
{
 ComplexMatrix * TmpA = new ComplexMatrix [this->PhysicalDimension];
 ComplexMatrix * TmpB = new ComplexMatrix [this->PhysicalDimension];  

 for(int i = 0; i < this->PhysicalDimension; i++)
 {
     TmpA[i] = this->M[i] * SingularValues;
     TmpB[i] =  ((SingularValues * rightSite->M[i]) / OldSingularValues);
}
 
 ComplexVector * PredictedPsi = new ComplexVector((long) this->BondDimensionRight*rightSite->BondDimensionLeft*this->PhysicalDimension*this->PhysicalDimension ,true);

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
             (*PredictedPsi)[this->SiteOnRight->GetVectorTwoSiteIndice(LeftIndice,RightIndice, LeftPhysicalDimension +this->PhysicalDimension*RightPhysicalDimension)] =   TmpB[LeftPhysicalDimension].GetMatrixElement(LeftIndice,k) *  TmpA[RightPhysicalDimension].GetMatrixElement(k ,RightIndice);
	 }
	}
    }
  }
}
  return PredictedPsi;
}
