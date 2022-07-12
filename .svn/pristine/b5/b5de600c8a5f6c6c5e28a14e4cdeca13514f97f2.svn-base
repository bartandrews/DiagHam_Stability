#include "RealMPOperatorOBC.h"
#include "Tensor/Tensor3.h"
#include "Matrix/RealSymmetricMatrix.h"
#include <iostream>
#include <sys/time.h>
#include "Architecture/ArchitectureOperation/TensorMatrixContractionOperation.h"
#include "Architecture/ArchitectureOperation/MPOApplyOnTensorOperation.h"
#include "Architecture/ArchitectureOperation/TensorVectorContractionOperation.h"

using std::cout;
using std::endl;

RealMPOperatorOBC::RealMPOperatorOBC()
{
  this->ElementsValues = 0;
  this->LeftVector = 0;
  this->RightVector = 0;
}

RealMPOperatorOBC::~RealMPOperatorOBC()
{
}

void RealMPOperatorOBC::ComputeLCore(Tensor3<double> & L)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  
  Tensor3<double> & LeftL = ((RealMPSSite *) this->Site)->GetPreviousL();
  
  Tensor3<double> * B =  new  Tensor3<double> [this->PhysicalDimension];
  RealMatrix * M = ((RealMPSSite *) this->Site)->GetM();

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
     B[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
     for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
     {
     for(int LeftB = 0; LeftB < this->MPOBondDimension; LeftB++)
     {
      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
              double & Tmp =  B[i](RightA,LeftB,LeftC);
	      for(int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
		{
		   Tmp +=  LeftL(LeftA,LeftB,LeftC)* M[i](LeftA,RightA);
		}
	    }
           }
	}
    }

  Tensor3<double> * A = new Tensor3<double>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
  A[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
  }

  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
     this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
     for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
     {
      for (int RightA = 0; RightA < BondDimensionRight; RightA++)
	{
	      A[MPOIndiceUp](RightA, MPOIndiceRight,LeftC) +=  B[MPOIndiceDown](RightA,MPOIndiceLeft,LeftC) * this->ElementsValues[i];
	    }
	}
    }

  delete [] B;
  

  for (int RightC = 0; RightC <  BondDimensionRight;  RightC++)
  {
	  for (int  RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	    { 
             for (int RightA = 0; RightA < BondDimensionRight; RightA++)
             {
              double & Tmp = L(RightA,RightB,RightC);
              for (int i = 0; i < this->PhysicalDimension; i++)
               {
	      for (int LeftC = 0;  LeftC < BondDimensionLeft; LeftC++)
		{ 
		   Tmp += M[i](LeftC,RightC) * A[i](RightA,RightB,LeftC);
		   }
		}
	    }
	}
    } 


  delete [] A;
}


void RealMPOperatorOBC::ComputeLCoreBis(Tensor3<double> & L)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  
  Tensor3<double> & LeftL = ((RealMPSSite *) this->Site)->GetPreviousL();  
  Tensor3<double> * B =  new  Tensor3<double> [this->PhysicalDimension];
  RealMatrix * M = ((RealMPSSite *) this->Site)->GetM();

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
     B[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
     TensorMatrixContractionOperation<double,0,0> Operation (&LeftL, &B[i],&M[i],0,0);
     Operation.SetIndicesRange(0,BondDimensionRight);
     Operation.ApplyOperation(this->Architecture);
    }

  
  Tensor3<double> * A = new Tensor3<double>[this->PhysicalDimension];

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
  A[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
    }

  MPOApplyOnTensorOperation<double> Operation (A, B,this,false);
  Operation.SetIndicesRange(0,BondDimensionRight);
  Operation.ApplyOperation(this->Architecture);
  delete [] B;
  

    for (int i = 0; i < this->PhysicalDimension; i++)
    {
      TensorMatrixContractionOperation<double,0,0> Operation (&A[i], &L,&M[i],2,0);
      Operation.SetIndicesRange(0,BondDimensionRight);
      Operation.ApplyOperation(this->Architecture);
   }

  delete [] A;
}


void RealMPOperatorOBC::ComputeRCore(Tensor3<double> & R)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  Tensor3<double> & RightR = ((RealMPSSite *) this->Site)->GetNextR();
  
  Tensor3<double> * B =  new  Tensor3<double>  [this->PhysicalDimension];
  RealMatrix * M = ((RealMPSSite *) this->Site)->GetM();
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
      for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
	{
	  for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	    {
	      for(int RightC = 0; RightC < BondDimensionRight; RightC++)
		{
		  double & Tmp =  B[i](LeftA,RightB,RightC);
		  for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		    {
		      Tmp +=  RightR(RightA,RightB,RightC) * M[i](LeftA,RightA);
		    }
		}
	    }
	}
    }
  
  
  Tensor3<double> * A = new Tensor3<double> [this->PhysicalDimension];
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }
  
  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
      
      for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
	  for (int RightC = 0;  RightC < BondDimensionRight;  RightC++)
	    {
	      A[MPOIndiceUp](LeftA, MPOIndiceLeft,RightC) +=  B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) * this->ElementsValues[i];
	    }
	}
    } 
  
  delete [] B;
  
  for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
    {
      for (int LeftC = 0; LeftC <  BondDimensionLeft;  LeftC++)
	{
	  for (int  LeftB = 0;  LeftB < this->MPOBondDimension;  LeftB++)
	    { 
              double & Tmp =  R(LeftA,LeftB,LeftC);
	      for (int i = 0; i < this->PhysicalDimension; i++)
		{
		  for (int RightC = 0;  RightC < BondDimensionRight; RightC++)
		    { 
		      Tmp += M[i](LeftC,RightC) * A[i](LeftA,LeftB,RightC);
		    }
		}
	    }
	}
    }
  
  delete [] A;
}



void RealMPOperatorOBC::ComputeRCoreBis(Tensor3<double> & R)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  Tensor3<double> & RightR = ((RealMPSSite *) this->Site)->GetNextR();
  
  Tensor3<double> * B =  new  Tensor3<double>  [this->PhysicalDimension];
  RealMatrix * M = ((RealMPSSite *) this->Site)->GetM();

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
      TensorMatrixContractionOperation<double,0,0> Operation (&RightR, &B[i],&M[i],0,1);
      Operation.SetIndicesRange(0,BondDimensionRight);
      Operation.ApplyOperation(this->Architecture);
   }

  Tensor3<double> * A = new Tensor3<double> [this->PhysicalDimension];
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }

  MPOApplyOnTensorOperation<double> Operation (A, B,this,true);
  Operation.SetIndicesRange(0,BondDimensionLeft);
  Operation.ApplyOperation(this->Architecture);
  delete [] B;

  for (int i = 0; i < this->PhysicalDimension; i++)
{
  TensorMatrixContractionOperation<double,0,0> Operation (&A[i], &R,&M[i],2,1);
  Operation.SetIndicesRange(0,BondDimensionLeft);
  Operation.ApplyOperation(this->Architecture);
}
  delete [] A;
}




void RealMPOperatorOBC::MPOApplyOnTensorOnTheRightCore(Tensor3<double> * result, Tensor3<double> * source, 
				       int firstComponent, int nbrComponent)
{
  cout <<"in the parallelized part"<< firstComponent << " " << nbrComponent<<endl;
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
   for (int LeftA =  firstComponent;  LeftA < firstComponent+nbrComponent ;  LeftA++)
	 {
 for (int i =  0; i < this->NbrNonZeroElements; i++)
       {

       this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
     
	  for (int RightC = 0;  RightC < BondDimensionRight;  RightC++)
	    {
	       result[MPOIndiceUp](LeftA, MPOIndiceLeft,RightC) +=   source[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) * this->ElementsValues[i];
	    }
	}
   } 
}


void RealMPOperatorOBC::MPOApplyOnTensorOnTheLeftCore(Tensor3<double> * result, Tensor3<double> * source, 
				       int firstComponent, int nbrComponent)
{
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
   unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
       for (int RightA =  firstComponent;  RightA < firstComponent+nbrComponent ;  RightA++)
	 {
 for (int i =  0; i < this->NbrNonZeroElements; i++)
       {

       this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
     
	  for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	    {
	       result[MPOIndiceUp](RightA, MPOIndiceRight,LeftC) +=   source[MPOIndiceDown](RightA,MPOIndiceLeft,LeftC) * this->ElementsValues[i];
	    }
	}
   } 
}




void RealMPOperatorOBC::LowLevelMultiplyCoreFirst(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent)
{
  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  int MaxLinearizedIndex = firstComponent + nbrComponent;
  int TmpIndice,TmpIndice2 ;
  
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
  for (int LinearizedIndexMiddleTop = firstComponent ; LinearizedIndexMiddleTop < MaxLinearizedIndex ; LinearizedIndexMiddleTop++)
  { 
     TmpIndice = LinearizedIndexMiddleTop *BondDimensionRight;
     TmpIndice2  = LinearizedIndexMiddleTop * BondDimensionLeft;
     for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
              double & Tmp =  result[i][TmpIndice2+LeftA];
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		  Tmp +=   (*source)[TmpIndice+RightA] * vSource[this->Site->GetVectorOneSiteIndice(LeftA, RightA, i)];
		}
	    }
         }
}
}


void RealMPOperatorOBC::LowLevelMultiplyCoreSecond(Tensor3<double> * leftTensor, Tensor3<double> * source , RealVector & vDestination, int firstComponent, int nbrComponent)
{
  long int LastComponent = firstComponent +  nbrComponent;
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  int LeftC, RightC, PhysicalIndice;

/*  for (int i = 0; i < this->PhysicalDimension; i++)
      {
cout <<" source [" << i << "]"<<endl;
       source[i].PrintTensor();
       }
   cout <<" leftTensor" <<endl;
       leftTensor->PrintTensor();*/

  for(long int Index =  firstComponent; Index < LastComponent ;Index++)
  {
    double & Tmp = vDestination[Index];
    this->Site->DecodeVectorOneSiteIndice(Index, LeftC, RightC,PhysicalIndice);
  for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
       for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
          Tmp +=  (*leftTensor)(LeftA,RightB, LeftC) *  source[PhysicalIndice](LeftA, RightB, RightC);
        }
    }
  }
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& RealMPOperatorOBC::LowLevelMultiplyCore(RealVector& vSource, RealVector& vDestination, 
				       int firstComponent, int nbrComponent)
{

  int BondDimensionRight = this->Site->GetBondDimensionRight();
  int BondDimensionLeft = this->Site->GetBondDimensionLeft();
  Tensor3<double> & RightR = ((RealMPSSite *) this->Site)->GetNextR();
  Tensor3<double> & LeftL = ((RealMPSSite *) this->Site)->GetPreviousL();
 
  Tensor3<double> * B =  new  Tensor3<double>  [this->PhysicalDimension];
  int TmpIndice,TmpIndice2 ;
  int MaxLinearizedIndex = BondDimensionRight * this->MPOBondDimension;

  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
   }
 
  TensorVectorContractionOperation<double> Operation(&RightR, B,  &vSource,this,true);
  Operation.SetIndicesRange(0,MaxLinearizedIndex);
  Operation.ApplyOperation(this->Architecture);

/* for (int i = 0; i < this->PhysicalDimension; i++)
    {
     for (int LinearizedIndexMiddleTop = 0 ; LinearizedIndexMiddleTop < MaxLinearizedIndex ; LinearizedIndexMiddleTop++)
     { 
      TmpIndice = LinearizedIndexMiddleTop *BondDimensionRight;
      TmpIndice2  = LinearizedIndexMiddleTop * BondDimensionLeft;
      for (int LeftA = 0; LeftA < BondDimensionLeft; LeftA++)
	{
              double & Tmp =  B[i][TmpIndice2+LeftA];
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{
		  Tmp +=  RightR[TmpIndice+RightA] * vSource[this->Site->GetVectorOneSiteIndice(LeftA, RightA, i)];
		}
	    }
         }
   }*/

  Tensor3<double> * A =  new  Tensor3<double>  [this->PhysicalDimension];
  for (int i = 0; i < this->PhysicalDimension; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }

//  MPOApplyOnTensorOperation<double> Operation (B, A,this,true);
//  Operation.SetIndicesRange(0,BondDimensionLeft);
//  Operation.ApplyOperation(this->Architecture);


  unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;

 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {

      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);

     
	  for (int RightC = 0;  RightC < BondDimensionRight;  RightC++)
	    {
       for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	 {
	      A[MPOIndiceUp](LeftA, MPOIndiceLeft,RightC) +=  B[MPOIndiceDown](LeftA,MPOIndiceRight,RightC) * this->ElementsValues[i];
	    }
	}
   } 

  delete [] B;

  TensorVectorContractionOperation<double> Operation1(&LeftL, A,  &vDestination,this,false);
  Operation1.SetIndicesRange(0,this->GetHilbertSpaceDimension());
  Operation1.ApplyOperation(this->Architecture);

/*
  int LeftC, RightC, PhysicalIndice;
  for(long int Index =  firstComponent; Index <  LastComponent ;Index++)
  {
    double & Tmp = vDestination[Index];
    this->Site->DecodeVectorOneSiteIndice(Index, LeftC, RightC,PhysicalIndice);
  for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
       for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	{
          Tmp +=  LeftL(LeftA,RightB, LeftC) * A[PhysicalIndice](LeftA, RightB, RightC);
        }
    }
  }*/

 delete [] A;
 return vDestination;
}




void RealMPOperatorOBC::LowLevelMultiplyCoreTwoSitesFirst(Tensor3<double> * result, Tensor3<double> * source , RealVector & vSource, int firstComponent, int nbrComponent)
{
  int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;

  int BondDimensionRight = this->SiteRight->GetBondDimensionRight();
  int BondDimensionLeft = this->SiteLeft->GetBondDimensionLeft();
  int MaxLinearizedIndex = firstComponent + nbrComponent;
  int TmpIndice,TmpIndice2 ;
  
  for (int i = 0; i < SquarePhysicalDimension; i++)
    {
  for (int LinearizedIndexMiddleTop = firstComponent ; LinearizedIndexMiddleTop < MaxLinearizedIndex ; LinearizedIndexMiddleTop++)
  { 
     TmpIndice = LinearizedIndexMiddleTop *BondDimensionRight;
     TmpIndice2  = LinearizedIndexMiddleTop * BondDimensionLeft;
      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	{
              double & Tmp =  result[i][RightA+TmpIndice];
           for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	    {
		  Tmp +=   (*source)[TmpIndice2 + LeftA] * vSource[this->SiteLeft->GetVectorTwoSiteIndice(LeftA,RightA,i)];
	    }
         }
}
}
}


void RealMPOperatorOBC::LowLevelMultiplyCoreTwoSitesSecond(Tensor3<double> * rightTensor, Tensor3<double> * source , RealVector & vDestination, int firstComponent, int nbrComponent)
{
  long int LastComponent = firstComponent +  nbrComponent;
  int BondDimensionRight = this->SiteRight->GetBondDimensionRight();
  int LeftC, RightC, PhysicalIndice;

  for(long int Index =  firstComponent; Index < LastComponent ;Index++)
  {
    double & Tmp = vDestination[Index];
    this->Site->DecodeVectorTwoSiteIndice(Index, LeftC, RightC,PhysicalIndice);
      for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
       for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	{
          Tmp +=  (*rightTensor)(RightA,RightB,RightC) *  source[PhysicalIndice](RightA,RightB,LeftC);
        }
    }
  }
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& RealMPOperatorOBC::LowLevelMultiplyTwoSitesCore(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent)
{
//  int TmpIndice,TmpIndice2 ;
  int BondDimensionRight = this->SiteRight->GetBondDimensionRight();
  int BondDimensionLeft = this->SiteLeft->GetBondDimensionLeft();
  int MaxLinearizedIndex = BondDimensionRight * this->MPOBondDimension;
  Tensor3<double> & RightR = ((RealMPSSite *) this->SiteRight)->GetNextR();
  Tensor3<double> & LeftL = ((RealMPSSite *) this->SiteLeft)->GetPreviousL();

  int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;

  Tensor3<double> * B =  new  Tensor3<double>  [SquarePhysicalDimension];

  for (int LinearizedPhysicalIndice= 0 ; LinearizedPhysicalIndice < SquarePhysicalDimension; LinearizedPhysicalIndice++)
  {
      B[LinearizedPhysicalIndice] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
  }


/*
  for (int LinearizedIndexMiddleTop = 0 ; LinearizedIndexMiddleTop < MaxLinearizedIndex ; LinearizedIndexMiddleTop++)
    { 
      TmpIndice = LinearizedIndexMiddleTop * BondDimensionRight;
      TmpIndice2  = LinearizedIndexMiddleTop * BondDimensionLeft;

      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	{
          double & Tmp = B[LinearizedPhysicalIndice][RightA+TmpIndice];

           for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	    {
		  Tmp +=  LeftL[TmpIndice2 + LeftA] * vSource[this->SiteLeft->GetVectorTwoSiteIndice(LeftA,RightA,LinearizedPhysicalIndice)];
	    }
	    }
     }
   } 
*/

  TensorVectorContractionOperation<double> Operation(&LeftL, B,  &vSource,this,true,true);
  Operation.SetIndicesRange(0,MaxLinearizedIndex);
  Operation.ApplyOperation(this->Architecture);



  Tensor3<double> * A =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      A[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
    }
 
   unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceMiddle,MPOIndiceRight;

   for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceMiddle);
    
 for(int PhysicalIndiceRight= 0; PhysicalIndiceRight <this->PhysicalDimension ; PhysicalIndiceRight++)
     {
      for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	{
	  for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	    {
              A[MPOIndiceUp +  this->PhysicalDimension *  PhysicalIndiceRight](RightA, MPOIndiceMiddle,LeftC) += B[MPOIndiceDown +  this->PhysicalDimension *  PhysicalIndiceRight](RightA,MPOIndiceLeft,LeftC) * this->ElementsValues[i];
	    }
	}
    }
   }



  delete [] B;

  Tensor3<double> * C =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      C[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
    }

 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
     this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceMiddle,  MPOIndiceRight);
     for(int PhysicalIndiceLeft = 0; PhysicalIndiceLeft <this->PhysicalDimension ; PhysicalIndiceLeft++)
     {
      for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	{
	  for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	    {
              C[PhysicalIndiceLeft +  this->PhysicalDimension *   MPOIndiceUp](RightA, MPOIndiceRight,LeftC) += A[PhysicalIndiceLeft +  this->PhysicalDimension *  MPOIndiceDown](RightA,MPOIndiceMiddle,LeftC) * this->ElementsValues[i];
	    }
	}
    }
   }

 delete [] A;

  TensorVectorContractionOperation<double> Operation1(&RightR, C, &vDestination,this,false,true);
  Operation1.SetIndicesRange(0,this->GetTwoSitesHilbertSpaceDimension());
  Operation1.ApplyOperation(this->Architecture);

/*
  int LeftC, RightC, PhysicalIndice;
  for(int Index =  firstComponent; Index < LastComponent ;Index++)
  {
    double & Tmp = vDestination[Index];
    this->Site->DecodeVectorTwoSiteIndice(Index, LeftC, RightC,PhysicalIndice);
       for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	{
      for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
            Tmp +=  RightR(RightA,RightB,RightC) *  C[PhysicalIndice](RightA,RightB,LeftC);
        }
    }
 }*/

 delete [] C;

 return vDestination;
}


// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix
RealSymmetricMatrix& RealMPOperatorOBC::GetTwoSitesHamiltonian (RealSymmetricMatrix & M)
{
  RealVector TmpV1 (this->GetTwoSitesHilbertSpaceDimension(), true);
  RealVector TmpV2 (this->GetTwoSitesHilbertSpaceDimension(), true);
  for (int i = 0; i < this->GetTwoSitesHilbertSpaceDimension(); i++)
    {
      TmpV1[i] = 1.0;
      this->AbstractMPOperatorOBC::LowLevelMultiplyTwoSites(TmpV1, TmpV2);
      if (this->LeftHamiltonianVectorMultiplicationFlag == false)
	{
	  for (int j = i; j < this->GetTwoSitesHilbertSpaceDimension(); j++)
	    {
	      M.SetMatrixElement(i, j, TmpV2[j]);
	    }
	}
      else
	{
	  for (int j = i; j < this->GetHilbertSpaceDimension(); j++)
	    {
	      M.SetMatrixElement(j, i, TmpV2[j]);
	    }
	}
      TmpV1[i] = 0.0;
    }
  return M;  
}



void RealMPOperatorOBC::ComputeL(Tensor3<double> & L)
{
  if (this->Site->GetBondDimensionLeft() == 1)
    {
      int BondDimensionRight = this->Site->GetBondDimensionRight();
      RealMatrix * M =  ((RealMPSSite * )this->Site)->GetM();
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
          this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
	  
	  for (int RightC = 0;RightC < this->Site->GetBondDimensionRight() ; RightC++ )
	    {
	      for (int RightA = 0;RightA < this->Site->GetBondDimensionRight() ; RightA++ )
		{
		  
		  L(RightA, MPOIndiceRight,RightC) +=  M[MPOIndiceUp](0,RightC) * this->ElementsValues[i] * this->LeftVector[MPOIndiceLeft] * M[MPOIndiceDown](0,RightA);
		}
	    }
	}
    }
  else
    {
  this->ComputeLCore(L);
//  this->ComputeLCoreBis(L);
}

}


void RealMPOperatorOBC::ComputeR(Tensor3<double> & R)
{
  if (this->Site->GetBondDimensionRight() == 1)
    {
      RealMatrix * M = ((RealMPSSite *)this->Site)->GetM();
      int BondDimensionLeft = this->Site->GetBondDimensionLeft();
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
	  this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
	  
	  for (int LeftC = 0;LeftC < this->Site->GetBondDimensionLeft() ; LeftC++ )
	    {
	      for (int LeftA = 0; LeftA < this->Site->GetBondDimensionLeft() ; LeftA++ )
		{
		  R(LeftA, MPOIndiceLeft,LeftC) +=  M[MPOIndiceUp](LeftC,0) * this->ElementsValues[i]*this->RightVector[MPOIndiceRight] * M[MPOIndiceDown](LeftA,0);
		}
	    }
	}
    }
  else
    this->ComputeRCore(R);
}



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& RealMPOperatorOBC::LowLevelMultiplyOneSite(RealVector& vSource, RealVector& vDestination, 
						       int firstComponent, int nbrComponent)
{
  vDestination.ClearVector();
  if (this->Site->GetBondDimensionLeft() == 1)
    {
      int BondDimensionRight = this->Site->GetBondDimensionRight(); 
      Tensor3<double> & RightR = ((RealMPSSite * )this->Site)->GetNextR();
      Tensor3<double> * B = new Tensor3<double>[this->PhysicalDimension];
      
      for (int i = 0; i < this->PhysicalDimension; i++)
	{
	  B[i] = Tensor3<double>(this->MPOBondDimension,BondDimensionRight,1,true);
	  for (int RightB = 0; RightB < this->MPOBondDimension ; RightB++)
	    {
	      for(int RightC = 0; RightC < BondDimensionRight; RightC++)
		{
		  for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		    {
		      B[i](RightB,RightC,0) +=  RightR(RightA,RightB,RightC) * vSource[(long int)BondDimensionRight*i + RightA]; 
		    }
		}
	    }
	}
      
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
	  this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
	  for (int NewRight = 0;  NewRight < BondDimensionRight;  NewRight++)
	    {
	      vDestination[(long int)BondDimensionRight*MPOIndiceUp + NewRight] +=  B[MPOIndiceDown](MPOIndiceRight,NewRight,0) *  this->ElementsValues[i] * this->LeftVector[MPOIndiceLeft];
	    }
	}
      
      delete [] B;
      return  vDestination;
    }

  if (this->Site->GetBondDimensionRight() == 1)
    {
      int BondDimensionLeft = this->Site->GetBondDimensionLeft(); 
      Tensor3<double> & LeftL = ((RealMPSSite * )this->Site)->GetPreviousL();
      Tensor3<double> * B = new Tensor3<double>[this->PhysicalDimension];
      
      for (int i = 0; i < this->PhysicalDimension; i++)
	{
	  B[i] = Tensor3<double>(this->MPOBondDimension,BondDimensionLeft,1,true);
	  
	  for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
	    {
	      for (int LeftB = 0; LeftB < this->MPOBondDimension ; LeftB++)
		{
		  for(int LeftA  = 0;  LeftA < BondDimensionLeft;  LeftA++)
		    {
		      B[i](LeftB,LeftC,0) +=  LeftL(LeftA,LeftB,LeftC) * vSource[(long int)BondDimensionLeft*i + LeftA];
		    }
		}
	    }
	}
      
      unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
      for (int i = 0; i < this->NbrNonZeroElements; i++)
	{
	  this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);      
	  
	  for (int NewLeft = 0;  NewLeft < BondDimensionLeft;  NewLeft++)
	    {
              vDestination[(long int)BondDimensionLeft*MPOIndiceUp + NewLeft] +=  B[MPOIndiceDown](MPOIndiceLeft,NewLeft,0) *  this->ElementsValues[i] * this->RightVector[MPOIndiceRight];
            }
	}
      
      delete [] B;
      return vDestination;
    }
  
  return this->LowLevelMultiplyCore(vSource,vDestination,firstComponent,nbrComponent);
}
 
// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& RealMPOperatorOBC::LowLevelMultiplyTwoSites(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent)
{

 vDestination.ClearVector();

 if ((this->SiteLeft->GetBondDimensionLeft() == 1) && (this->SiteRight->GetBondDimensionRight() ==  1 ) )
   {
     
     int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;
     
     Tensor3<double> * B =  new  Tensor3<double>  [SquarePhysicalDimension];
     for (int i = 0; i < SquarePhysicalDimension; i++)
       {
	 B[i] = Tensor3<double>(this->MPOBondDimension,1,1,true);
       }
     
     unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceRight;
     for(int  PhysicalIndiceRight = 0;  PhysicalIndiceRight <this->PhysicalDimension ; PhysicalIndiceRight++)
       {
	 for (int i = 0; i < this->NbrNonZeroElements; i++)
	   {
	     this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight); 
	     B[MPOIndiceUp +  this->PhysicalDimension *  PhysicalIndiceRight](MPOIndiceRight,0,0) +=  this->ElementsValues[i] * this->LeftVector[MPOIndiceLeft] * vSource[ (long int) MPOIndiceDown +  this->PhysicalDimension *  PhysicalIndiceRight];
	   }
	 
       }
     for(int PhysicalIndiceLeft= 0;   PhysicalIndiceLeft <this->PhysicalDimension ; PhysicalIndiceLeft++)
       {
	 for (int i = 0; i < this->NbrNonZeroElements; i++)
	   {
	     this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceRight);
	     vDestination[(long) PhysicalIndiceLeft + this->PhysicalDimension *MPOIndiceUp] += this->ElementsValues[i] * this->RightVector[MPOIndiceRight] *  B[PhysicalIndiceLeft +  this->PhysicalDimension *  MPOIndiceDown](MPOIndiceLeft,0,0);
	   }
       }
     
     delete [] B;
     return vDestination;
     
   }
 
 return this->LowLevelMultiplyTwoSitesCore(vSource,vDestination,firstComponent,nbrComponent);
}


void  RealMPOperatorOBC::PrintTensorElements()
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



// multiply a vector by the current hamiltonian for a given range of indices 
// and store result in another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& RealMPOperatorOBC::LowLevelMultiplyTwoSitesCoreBis(RealVector& vSource, RealVector& vDestination,  int firstComponent, int nbrComponent)
{

/*  timeval TotalStartingTime;
  timeval TotalEndingTime;
  gettimeofday (&TotalStartingTime, 0)
;
*/
  int BondDimensionRight = this->SiteRight->GetBondDimensionRight();
  int BondDimensionLeft = this->SiteLeft->GetBondDimensionLeft();

  Tensor3<double> & RightR = ((RealMPSSite *) this->SiteRight)->GetNextR();
  Tensor3<double> & LeftL = ((RealMPSSite *) this->SiteLeft)->GetPreviousL();

  int SquarePhysicalDimension = this->PhysicalDimension * this->PhysicalDimension;

  Tensor3<double> * B =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension; i++)
    {
      B[i] = Tensor3<double>(BondDimensionRight,this->MPOBondDimension,BondDimensionLeft,true);
    }

//  int LinearizedPhysicalIndice = PhysicalIndiceLeft +  this->PhysicalDimension *  PhysicalIndiceRight;

  for (int LinearizedPhysicalIndice= 0 ; LinearizedPhysicalIndice < SquarePhysicalDimension; LinearizedPhysicalIndice++)
  {


  for (int LeftB = 0; LeftB < this->MPOBondDimension ; LeftB++)
    {
	  for(int LeftC = 0; LeftC < BondDimensionLeft; LeftC++)
	    {
	      for(int RightA = 0;  RightA < BondDimensionRight;  RightA++)
		{

           for (int LeftA = 0;  LeftA < BondDimensionLeft;  LeftA++)
	    {
		  B[LinearizedPhysicalIndice](RightA,LeftB,LeftC) +=  LeftL(LeftA,LeftB,LeftC) * vSource[this->SiteLeft->GetVectorTwoSiteIndice(LeftA,RightA, LinearizedPhysicalIndice)];
		}
	    }
	}
    }


}



  Tensor3<double> * A =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      A[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }
 
unsigned int MPOIndiceDown,MPOIndiceLeft,MPOIndiceUp,MPOIndiceMiddle,MPOIndiceRight;

 for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceLeft,  MPOIndiceMiddle);
    
 for(int PhysicalIndiceRight= 0; PhysicalIndiceRight <this->PhysicalDimension ; PhysicalIndiceRight++)
     {
      for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	{
	  for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	    {
              A[MPOIndiceUp +  this->PhysicalDimension *  PhysicalIndiceRight](LeftC, MPOIndiceMiddle,RightA) += B[MPOIndiceDown +  this->PhysicalDimension *  PhysicalIndiceRight](RightA,MPOIndiceLeft,LeftC) * this->ElementsValues[i];
	    }
	}
    }
   }



  delete [] B;

  Tensor3<double> * C =  new  Tensor3<double>  [SquarePhysicalDimension];
  for (int i = 0; i < SquarePhysicalDimension ; i++)
    {
      C[i] = Tensor3<double>(BondDimensionLeft,this->MPOBondDimension,BondDimensionRight,true);
    }

  for (int i = 0; i < this->NbrNonZeroElements; i++)
    {
      this->GetAllIndicesFromTensorIndex(this->IndexValues[i], MPOIndiceDown, MPOIndiceUp, MPOIndiceMiddle,  MPOIndiceRight);
     for(int PhysicalIndiceLeft= 0; PhysicalIndiceLeft <this->PhysicalDimension ; PhysicalIndiceLeft++)
     {
      for (int LeftC = 0;  LeftC < BondDimensionLeft;  LeftC++)
	{
	  for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	    {
              C[PhysicalIndiceLeft +  this->PhysicalDimension *   MPOIndiceUp](LeftC, MPOIndiceRight,RightA) += A[PhysicalIndiceLeft +  this->PhysicalDimension *  MPOIndiceDown](LeftC,MPOIndiceMiddle,RightA) * this->ElementsValues[i];
	    }
	}
    }
   } 

 delete [] A;

 // Index = LinearizedPhysicalIndice + SquarePhysicalDimension*(LeftA + RightA*BondDimensionLeft))

  int LastComponent = firstComponent + nbrComponent;

/*  gettimeofday (&TotalEndingTime, 0);
  double  Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
  cout <<"First Part " << Dt << "s" << endl;
  gettimeofday (&TotalStartingTime, 0);*/
  for(int Index =  firstComponent; Index < LastComponent ;Index++)
  {
       for (int RightA = 0;  RightA < BondDimensionRight;  RightA++)
	{
      for (int RightB = 0;  RightB < this->MPOBondDimension;  RightB++)
	{
            vDestination[Index] +=  RightR(RightA,RightB,Index/(SquarePhysicalDimension*BondDimensionLeft)) *  C[Index%SquarePhysicalDimension](Index/SquarePhysicalDimension%BondDimensionLeft , RightB,RightA);
        }
    }
 }

/*
gettimeofday (&TotalEndingTime, 0);
Dt = (((double) (TotalEndingTime.tv_sec - TotalStartingTime.tv_sec)) +
 		(((double) (TotalEndingTime.tv_usec - TotalStartingTime.tv_usec)) / 1000000.0));
 cout <<"Second Part " << Dt << "s" << endl;*/

 delete [] C;

 return vDestination;
}
