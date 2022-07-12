////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of abstract quantum Hall hamiltonian associated            //
//     to particles on a sphere with spin and n-body interaction terms        //
//                                                                            //
//                        last modification : 26/08/2008                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Hamiltonian/ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;
using std::ofstream;
using std::ifstream;
using std::ios;


// destructor
//

ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::~ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian()
{
  for (int k=0; k <= this->MaxNBody; ++k)
    {
      if (this->NBodyFlags[k])
	{
	  for (int i = 0; i < this->NbrSpinSectors[k]; ++i)
	    {
	      for (int m = 0; m < this->NbrNBodySpinMomentumSectorSum[k][i]; ++m)
		delete[] this->NBodyInteractionFactors[k][i][m];
	      if (this->NbrNBodySpinMomentumSectorSum[k][i] != 0)
		{
		  for (int m = 0; m < this->NbrNBodySpinMomentumSectorSum[k][i]; ++m)
		    {
		      delete[] this->NBodySpinMomentumSectorIndicesPerSum[k][i][m];
		    }
		}
	      delete[] this->NBodyInteractionFactors[k][i];
	      delete[] this->NBodySpinMomentumSectorIndicesPerSum[k][i];
	      delete[] this->NbrNBodySpinMomentumSectorIndicesPerSum[k][i];
	      delete[] this->SpinIndices[k][i];
	    }
	  delete [] this->NbrNBodySpinMomentumSectorSum[k];
	  delete [] this->NBodyInteractionFactors[k];
	  delete [] this->NbrNBodySpinMomentumSectorIndicesPerSum[k];
	  delete [] this->NBodySpinMomentumSectorIndicesPerSum[k];
	  delete [] this->SpinIndices[k];
	  delete [] this->SpinIndicesShort[k];
	}
    }
  delete [] this->NBodyFlags;
  delete [] this->NBodySign;
  delete [] this->SpinIndices;
  delete [] this->SpinIndicesShort;
  delete [] this->NbrSpinSectors;
  delete [] this->NbrNBodySpinMomentumSectorSum;
  delete [] this->NBodyInteractionFactors;
  delete [] this->NBodySpinMomentumSectorIndicesPerSum;
  delete [] this->NbrNBodySpinMomentumSectorIndicesPerSum;
  
}


// symmetrize interaction factors to enable hermitian matrix multiplication
// return = true upon success
bool ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::HermitianSymmetrizeInteractionFactors()
{
  if (HermitianSymmetryFlag)
    return true;

  if (this->Particles->HaveOrder()==false)
    {
      cout << "Hamiltonian tried to use hermitian symmetry, but this is not implemented in HilbertSpace!"<<endl;
      HermitianSymmetryFlag=false;
      return false;
    }

  cout << "Using hermitian symmetry not implemented, yet."<<endl;

  /* 
  int *M = new int[2];
  int *N = new int[2];

  // single particle terms
  if (NbrHoppingTerms>0)
    {
      int TmpNbrHoppingTerms = 0;
      int *Flags = new int[this->NbrHoppingTerms];
      for (int j = 0; j < NbrHoppingTerms; ++j) 
	{
	  M[0] = this->KineticQi[j];
	  N[0] = this->KineticQf[j];
	  Flags[j] = this->Particles->CheckOrder(M, N, 1);
	  // cout << "M="<<M[0]<<", N="<<N[0]<<", order: "<<Flags[j]<<" element: "<<HoppingTerms[j]<<endl;
	  if (Flags[j]>0)
	    ++TmpNbrHoppingTerms;
	  else if (Flags[j]==0)
	    {
	      ++TmpNbrHoppingTerms;
	      HoppingTerms[j]*=0.5;
	    }
	}
      Complex *TmpHoppingTerms = new Complex[TmpNbrHoppingTerms];
      int *TmpQi = new int[TmpNbrHoppingTerms];
      int *TmpQf = new int[TmpNbrHoppingTerms];
      int Pos=0;
      for (int j = 0; j < this->NbrHoppingTerms; ++j)
	if (Flags[j]>=0)
	  {
	    TmpHoppingTerms[Pos]=this->HoppingTerms[j];
	    TmpQi[Pos]=this->KineticQi[j];
	    TmpQf[Pos]=this->KineticQf[j];
	    ++Pos;
	  }
      delete [] this->HoppingTerms;
      delete [] this->KineticQi;
      delete [] this->KineticQf;
      this->HoppingTerms = TmpHoppingTerms;
      this->KineticQi = TmpQi;
      this->KineticQf = TmpQf; 
      this->NbrHoppingTerms = TmpNbrHoppingTerms;
    }
  
  if (this->NbrQ12Indices == 0)
    {
      if (NbrInteractionFactors>0)
	{
	  int TmpNbrInteractionFactors = 0;
	  int *Flags = new int[NbrInteractionFactors];
	  for (int j = 0; j < NbrInteractionFactors; ++j) 
	    {
	      M[0] = this->Q1Value[j];
	      M[1] = this->Q2Value[j];
	      N[0] = this->Q3Value[j];
	      N[1] = this->Q4Value[j];
	      Flags[j] = this->Particles->CheckOrder (M, N, 2);
	      cout << "Flag("<<this->Q1Value[j]<<", "<<this->Q2Value[j]<<", "<<this->Q3Value[j]<<", "<<this->Q4Value[j]<<")="<<Flags[j]<<endl;
	      if (Flags[j]>0)
		++TmpNbrInteractionFactors;
	      else
		{
		  if (Flags[j]==0)
		    {
		      ++TmpNbrInteractionFactors;
		      this->InteractionFactors[j]*=0.5; // diagonal term: make up for double counting
		    }
		  else
		    {
		      cout << "Discarding element "<<this->Q1Value[j]<<", "<<this->Q2Value[j]<<", "<<this->Q3Value[j]<<", "<<this->Q4Value[j]<<endl;
		    }
		}
	    }
	  Complex* TmpInteractionFactors = new Complex[TmpNbrInteractionFactors];
	  int* TmpQ1Value = new int[TmpNbrInteractionFactors];
	  int* TmpQ2Value = new int[TmpNbrInteractionFactors];
	  int* TmpQ3Value = new int[TmpNbrInteractionFactors];
	  int* TmpQ4Value = new int[TmpNbrInteractionFactors];
	  int Pos=0;
	  for (int j = 0; j < NbrInteractionFactors; ++j)
	    {
	      if (Flags[j]>=0)
		{
		  TmpInteractionFactors[Pos]=InteractionFactors[j];
		  TmpQ1Value[Pos]=Q1Value[j];
		  TmpQ2Value[Pos]=Q2Value[j];
		  TmpQ3Value[Pos]=Q3Value[j];
		  TmpQ4Value[Pos]=Q4Value[j];
		  ++Pos;
		}
	    }
	  delete [] InteractionFactors;
	  delete [] Q1Value;
	  delete [] Q2Value;
	  delete [] Q3Value;
	  delete [] Q4Value;
	  this->InteractionFactors = TmpInteractionFactors;
	  this->NbrInteractionFactors = TmpNbrInteractionFactors;
	  this->Q1Value = TmpQ1Value;
	  this->Q2Value = TmpQ2Value;
	  this->Q3Value = TmpQ3Value;
	  this->Q4Value = TmpQ4Value;
	  delete [] Flags;
	}
    }
  else
    {
      int OldNbrQ34Values;
      int* OldQ3PerQ12;
      int* OldQ4PerQ12;
      int TmpNbrQ12Values = 0;
      int* Q12Flags = new int[this->NbrQ12Indices];
      int TmpNbrQ34Values;
      int* TmpQ3PerQ12;
      int* TmpQ4PerQ12;
      int* Q34Flags;
      // quick 5count of interaction factors
      int OldNbrInteractionFactors=0;
      for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	OldNbrInteractionFactors+=this->NbrQ34Values[q12];      
      int TmpNbrInteractionFactors=0;
      Complex *TmpInteractionFactors=new Complex[OldNbrInteractionFactors];
      int Pos=0;
      for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	{
	  M[0]=this->Q1Value[q12];
	  M[1]=this->Q2Value[q12];
	  OldNbrQ34Values = this->NbrQ34Values[q12];
	  OldQ3PerQ12 = this->Q3PerQ12[q12];
	  OldQ4PerQ12 = this->Q4PerQ12[q12];
	  TmpNbrQ34Values = 0;
	  Q34Flags = new int[OldNbrQ34Values];
	  for (int q34 = 0; q34 < OldNbrQ34Values; ++q34)
	    {
	      N[0]=OldQ3PerQ12[q34];
	      N[1]=OldQ4PerQ12[q34];
	      Q34Flags[q34] = this->Particles->CheckOrder(M, N, 2);
	      //cout << "Flag("<<M[0]<<", "<<M[1]<<", "<<N[0]<<", "<<N[1]<<")="<<Q34Flags[q34]<<endl;
	      if (Q34Flags[q34]>0)
		{
		  ++TmpNbrQ34Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=this->InteractionFactors[Pos];
		}
	      else if (Q34Flags[q34]==0)
		{
		  ++TmpNbrQ34Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=0.5*this->InteractionFactors[Pos];
		}
// 	      else
// 		{
// 		  cout << "Discarding element "<<M[0]<<", "<<M[1]<<", "<<N[0]<<", "<<N[1]<<endl;
// 		}
	      
	      ++Pos;
	    }
	  if (TmpNbrQ34Values>0)
	    {
	      //cout << "Q1="<<M[0]<<", Q2="<<M[1]<<": ";
	      ++TmpNbrQ12Values;
	      Q12Flags[q12]=1;
	      TmpQ3PerQ12 = new int[TmpNbrQ34Values];
	      TmpQ4PerQ12 = new int[TmpNbrQ34Values];
	      int Pos2=0;
	      for (int q34 = 0; q34 < OldNbrQ34Values; ++q34)
		if (Q34Flags[q34]>=0)
		  {
		    TmpQ3PerQ12[Pos2]=OldQ3PerQ12[q34];
		    TmpQ4PerQ12[Pos2]=OldQ4PerQ12[q34];
		    //cout << " Q3: " << TmpQ3PerQ12[Pos2];
		    //cout << " Q4: " << TmpQ4PerQ12[Pos2];
		    Pos2++;
		  }
	      //cout << endl;
	      delete [] OldQ3PerQ12;
	      delete [] OldQ4PerQ12;
	      this->Q3PerQ12[q12] = TmpQ3PerQ12;
	      this->Q4PerQ12[q12] = TmpQ4PerQ12;
	      this->NbrQ34Values[q12] = TmpNbrQ34Values;
 	    }
	  else
	    {
	      Q12Flags[q12]=-1;
	      delete [] OldQ3PerQ12;
	      delete [] OldQ4PerQ12;
	    }
	}
      if (this->NbrQ12Indices!=TmpNbrQ12Values)
	{
	  int *NewQ1Value=new int[TmpNbrQ12Values];
	  int *NewQ2Value=new int[TmpNbrQ12Values];
	  int **NewQ3PerQ12=new int*[TmpNbrQ12Values];
	  int **NewQ4PerQ12=new int*[TmpNbrQ12Values];
	  int *NewNbrQ34Values=new int[TmpNbrQ12Values];
	  Pos = 0;
	  for (int q12 = 0; q12 < this->NbrQ12Indices; ++q12)
	    if (Q12Flags[q12]>0)
	      {
		NewQ1Value[Pos]=this->Q1Value[q12];
		NewQ2Value[Pos]=this->Q2Value[q12];
		NewQ3PerQ12[Pos]=this->Q3PerQ12[q12];
		NewQ4PerQ12[Pos]=this->Q4PerQ12[q12];
		NewNbrQ34Values[Pos]=this->NbrQ34Values[q12];
		++Pos;
	      }
	  delete [] this->Q1Value;
	  delete [] this->Q2Value;
	  delete [] this->Q3PerQ12;
	  delete [] this->Q4PerQ12;
	  delete [] this->NbrQ34Values;
	  this->Q1Value=NewQ1Value;
	  this->Q2Value=NewQ2Value;
	  this->Q3PerQ12=NewQ3PerQ12;
	  this->Q4PerQ12=NewQ4PerQ12;
	  this->NbrQ34Values=NewNbrQ34Values;
	  this->NbrQ12Indices=TmpNbrQ12Values;
	}
      // reduce size of table InteractionFactors to match new size, and copy contents
      delete [] this->InteractionFactors;
      this->InteractionFactors = new Complex[TmpNbrInteractionFactors];
      for (int i=0; i<TmpNbrInteractionFactors; ++i)
	this->InteractionFactors[i]=TmpInteractionFactors[i];
      delete [] TmpInteractionFactors;
    }

  // diagonal terms are always the same... so we're done

  delete [] M;
  delete [] N;
  this->HermitianSymmetryFlag=true;
  
  return true;*/
  return false;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
											   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  Complex Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      if (this->FullTwoBodyFlag == true)
	for (int i = firstComponent; i < LastComponent; ++i)
	  this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination) ;
      for (int k = 3; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    for (int i = firstComponent; i < LastComponent; ++i)
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
	  }
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSource, vDestination);
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      Coefficient = vSource[k];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[k++] += this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  int l =  PosMod + firstComponent + this->PrecalculationShift;
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      Coefficient = vSource[l];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      vDestination[l] += this->HamiltonianShift * Coefficient;
	      l += this->FastMultiplicationStep;
	      ++Pos;
	    }
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  for (l = 0; l < this->FastMultiplicationStep; ++l)
	    if (PosMod != l)
	      {	
		    
		if (this->FullTwoBodyFlag == true)
		  for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		    this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination);
		for (int k = 3; k <= this->MaxNBody; ++k)
		  if (this->NBodyFlags[k] == true)
		    {
		      for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSource, vDestination, k) ;
		    }
		this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSource, vDestination);
	      }
	  delete TmpParticles;
	}
    }
//   if (this->S2Hamiltonian != 0)
//     this->S2Hamiltonian->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
												int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  Complex Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      Complex* Coefficient2 = new Complex [nbrVectors];
      ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
      if (this->FullTwoBodyFlag == true)
	{
	  for (int i = firstComponent; i < LastComponent; ++i)
	    this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2) ;
	}
      for (int k = 3; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  {
	    for (int i = firstComponent; i < LastComponent; ++i)
	      this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2);
	  }	    
      this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent, LastComponent, 1, vSources, vDestinations, nbrVectors);
      delete TmpParticles;
      delete[] Coefficient2;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  int j;
	  int Pos;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  Coefficient2[l] = vSources[l][k];
		  vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
	  int* TmpIndexArray;
	  Complex* TmpCoefficientArray; 
	  Complex* Coefficient2 = new Complex [nbrVectors];
	  int j;
	  int Pos2;
	  int TmpNbrInteraction;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  int Pos = firstComponent / this->FastMultiplicationStep; 
	  int PosMod = firstComponent % this->FastMultiplicationStep;
	  if (PosMod != 0)
	    {
	      ++Pos;
	      PosMod = this->FastMultiplicationStep - PosMod;
	    }
	  int l =  PosMod + firstComponent + this->PrecalculationShift;
	  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
	    {
	      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
	      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
	      for (int k = 0; k < nbrVectors; ++k)
		{
		  Coefficient2[k] = vSources[k][l];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos2 = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int k = 0; k < nbrVectors; ++k)
		    {
		      vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
		    }
		}
	      l += this->FastMultiplicationStep;
	      ++Pos;
	    }
	  firstComponent += this->PrecalculationShift;
	  LastComponent += this->PrecalculationShift;
	  for (l = 0; l < this->FastMultiplicationStep; ++l)
	    if (PosMod != l)
	      {	
		if (this->FullTwoBodyFlag == true)
		  {
		    for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
		      this->EvaluateMNTwoBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, Coefficient2) ;
		  }
		for (int k = 3; k <= this->MaxNBody; ++k)
		  if (this->NBodyFlags[k] == true)
		    {
		      for (int i = firstComponent + l; i < LastComponent; i += this->FastMultiplicationStep)
			this->EvaluateMNNBodyAddMultiplyComponent(TmpParticles, i, vSources, vDestinations, nbrVectors, k, Coefficient2) ;
		    }
		this->EvaluateMNOneBodyAddMultiplyComponent(TmpParticles, firstComponent + l, LastComponent, this->FastMultiplicationStep, vSources, vDestinations, nbrVectors);
	      }
	  delete[] Coefficient2;
	  delete TmpParticles;
	}
    }
//   if (this->S2Hamiltonian != 0)
//     this->S2Hamiltonian->LowLevelMultipleAddMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
  return vDestinations;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// lastComponent  = index of the last component that has to be precalcualted
// return value = number of non-zero matrix element

long ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int lastComponent)
{
  long Memory = 0;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();
  int LastComponent = lastComponent + firstComponent;

  if ((this->OneBodyInteractionFactorsdowndown != 0) || (this->OneBodyInteractionFactorsupup != 0)||(this->OneBodyInteractionFactorsupdown != 0))
    EvaluateMNOneBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory); 
  if (this->FullTwoBodyFlag == true)
    EvaluateMNTwoBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, Memory); 
  for (int k = 3; k <= this->MaxNBody; ++k)
    if (this->NBodyFlags[k] == true)
      EvaluateMNNBodyFastMultiplicationMemoryComponent(TmpParticles, firstComponent, LastComponent, k, Memory);
  delete TmpParticles;
  
  return Memory;
}


// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that have to be precalcualted
//
void ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{
  int LastComponent = nbrComponent + firstComponent;

  int* TmpIndexArray;
  Complex* TmpCoefficientArray;
  long ColumnIndex;
  ParticleOnSphereWithSpin* TmpParticles = (ParticleOnSphereWithSpin*) this->Particles->Clone();

  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  long Pos = firstComponent / this->FastMultiplicationStep; 
  int PosMod = firstComponent % this->FastMultiplicationStep;
  if (PosMod != 0)
    {
      ++Pos;
      PosMod = this->FastMultiplicationStep - PosMod;
    }
  for (int i = PosMod + firstComponent; i < LastComponent; i += this->FastMultiplicationStep)
    {
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      ColumnIndex = 0l;

      if ((this->OneBodyInteractionFactorsdowndown != 0) || (this->OneBodyInteractionFactorsupup != 0)||(this->OneBodyInteractionFactorsupdown != 0))
	this->EvaluateMNOneBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);
      
      if (this->FullTwoBodyFlag == true)
	this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, TmpIndexArray, TmpCoefficientArray, ColumnIndex);

      for (int k = 3; k <= this->MaxNBody; ++k)
	if (this->NBodyFlags[k] == true)
	  this->EvaluateMNNBodyFastMultiplicationComponent(TmpParticles, i, k, TmpIndexArray, TmpCoefficientArray, ColumnIndex);

      ++Pos;
    }
  delete TmpParticles;
}



// get all indices needed to characterize a completly skew symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::GetAllSkewSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
											  int**& sortedIndicesPerSum)
{
  long** BinomialCoefficients = GetBinomialCoefficients(nbrValues);
  long NbrElements = BinomialCoefficients[nbrValues][nbrIndices];
  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Min = nbrIndices - 1;
  int Max;
  int Step;
  int Pos = 0;
  for (int i = nbrValues - 1; i >= Min; --i)
    {
      Step = BinomialCoefficients[i][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      Min = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1] - 1;
	  for (; Max >= Min; --Max)
	    {
	      Step = BinomialCoefficients[Max][Min];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }

  int MaxSum = ((nbrValues - 1) * nbrIndices) - (((nbrIndices - 1) * (nbrIndices)) / 2);
  int MinSum = (nbrIndices * (nbrIndices - 1)) / 2;
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    ++nbrSortedIndicesPerSum[Sum[i]];
  long* TmpPos = new long [MaxSum + 1];
  for (int i = MinSum; i <= MaxSum; ++i)
    {
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      nbrSortedIndicesPerSum[i] = 0;
      TmpPos[i] = 0l;      
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = Indices[i][j];
	  ++TmpPos[Pos];
	}
      ++nbrSortedIndicesPerSum[Pos];
      delete[] Indices[i];
    }
  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] BinomialCoefficients[i];
  delete[] BinomialCoefficients;
  return NbrElements;
}

// get all indices needed to characterize a completly symmetric tensor, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndices = number of indices 
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::GetAllSymmetricIndices (int nbrValues, int nbrIndices, int*& nbrSortedIndicesPerSum, 
										     int**& sortedIndicesPerSum,
										     double**& sortedIndicesPerSumSymmetryFactor)
{
  long** DimensionSymmetricGroup;
  if (nbrValues >= nbrIndices)
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrValues);
  else
    DimensionSymmetricGroup = GetIrreducibleRepresentationDimensionSymmetricGroup(nbrIndices);
  long NbrElements = DimensionSymmetricGroup[nbrValues][nbrIndices];

  int** Indices = new int* [NbrElements];
  int* Sum = new int [NbrElements];
  int Max;
  int Step;
  int Pos = 0;
  int TmpNbrIndices;
  for (int i = nbrValues - 1; i >= 0; --i)
    {
      Step = DimensionSymmetricGroup[i + 1][nbrIndices - 1];
      for (int j = 0; j < Step; ++j)
	{
	  Indices[Pos] = new int [nbrIndices];
	  Indices[Pos][0] = i;
	  Sum[Pos] = i;
	  ++Pos;
	}
    }
  for (int i = 1; i < nbrIndices; ++i)
    {
      int Pos = 0;
      TmpNbrIndices = nbrIndices - i - 1;
      while (Pos < NbrElements)
	{
	  Max = Indices[Pos][i - 1];
	  Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	  for (int j = 0; j < Step; ++j)
	    {
	      Indices[Pos][i] = Max;
	      Sum[Pos] += Max;
	      ++Pos;
	    }
	  --Max;
	  for (; Max >= 0; --Max)
	    {
	      Step = DimensionSymmetricGroup[Max + 1][TmpNbrIndices];
	      for (int j = 0; j < Step; ++j)
		{
		  Indices[Pos][i] = Max;
		  Sum[Pos] += Max;
		  ++Pos;
		}
	    }
	}
    }

  int MaxSum = (nbrValues - 1) * nbrIndices;
  long* TmpPos = new long [MaxSum + 1];
  nbrSortedIndicesPerSum = new int [MaxSum + 1];
  sortedIndicesPerSum = new int* [MaxSum + 1];
  sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
  for (int i = 0; i <= MaxSum; ++i)
    nbrSortedIndicesPerSum[i] = 0;
  for (int i = 0; i < NbrElements; ++i)
    {
      ++nbrSortedIndicesPerSum[Sum[i]];
    }
  for (int i = 0; i <= MaxSum; ++i)
    {
      TmpPos[i] = 0l;
      sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * nbrIndices];
      sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
      nbrSortedIndicesPerSum[i] = 0;
    }
  for (int i = 0; i < NbrElements; ++i)
    {   
      Pos = Sum[i];
      Max = nbrSortedIndicesPerSum[Pos];
      int* TmpIndices = Indices[i];
      for (int j = 0; j < nbrIndices; ++j)
	{
	  sortedIndicesPerSum[Pos][TmpPos[Pos]] = TmpIndices[j];
	  ++TmpPos[Pos];
	}
      double& SymmetryFactor = sortedIndicesPerSumSymmetryFactor[Pos][Max];
      SymmetryFactor = 1.0;
      for (int j = 1; j < nbrIndices; ++j)
	{
	  int TmpSymmetryFactor = 1;
	  while ((j < nbrIndices) && (TmpIndices[j - 1] == TmpIndices[j]))
	    {
	      ++TmpSymmetryFactor;
	      ++j;
	    }
	  if (TmpSymmetryFactor != 1)
	    for (int k = 2; k <= TmpSymmetryFactor; ++k)
	      SymmetryFactor *= (double) k;
	}
      delete[] TmpIndices;
      SymmetryFactor = 1.0 / SymmetryFactor;
      ++nbrSortedIndicesPerSum[Pos];
    }

  delete[] TmpPos;
  delete[] Sum;
  delete[]Indices;
  for (int i = 0; i <= nbrValues; ++i)
    delete[] DimensionSymmetricGroup[i];
  delete[] DimensionSymmetricGroup;
  return NbrElements;
}

// get all indices needed to characterize a  tensor made of two completly skew symmetric  sets of indices, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
// nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long  ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::GetAllTwoSetSkewSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, 
												int**& sortedIndicesPerSum)
{
  int* TmpNbrSortedIndicesPerSumUp;
  int** TmpSortedIndicesPerSumUp;
  this->GetAllSkewSymmetricIndices(nbrValues, nbrIndicesUp, TmpNbrSortedIndicesPerSumUp, TmpSortedIndicesPerSumUp);
  int MaxSumUp = ((nbrValues - 1) * nbrIndicesUp) - (((nbrIndicesUp - 1) * nbrIndicesUp) / 2);
  int MinSumUp = (nbrIndicesUp * (nbrIndicesUp - 1)) / 2;
  long NbrElements = 0l;
  if (nbrIndicesDown == 1)
    {
      int MaxSum = MaxSumUp + nbrValues - 1;
      int MinSum = MinSumUp;
      nbrSortedIndicesPerSum = new int [MaxSum + 1];
      sortedIndicesPerSum = new int* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      nbrSortedIndicesPerSum[MinSum] = TmpNbrSortedIndicesPerSumUp[MinSumUp];
      sortedIndicesPerSum[MinSum] = new int [nbrSortedIndicesPerSum[MinSum] * (nbrIndicesUp + 1)];     
      int Lim = nbrSortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSum = sortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MinSum];
      int Pos = 0;
      int Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	  TmpSortedIndicesPerSum[Pos2++] = 0;
	}
      nbrSortedIndicesPerSum[MaxSum] = TmpNbrSortedIndicesPerSumUp[MaxSumUp];
      sortedIndicesPerSum[MaxSum] = new int [nbrSortedIndicesPerSum[MaxSum] * (nbrIndicesUp + 1)];     
      Lim = nbrSortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSum = sortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MaxSumUp];
      Pos = 0;
      Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	  TmpSortedIndicesPerSum[Pos2++] = nbrValues - 1;
	}
      for (int i = MinSum + 1; i < MaxSum; ++i)
	{
	  int TmpMinSumUp = i - nbrValues + 1;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + 1)];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  for (int l = 0; l < nbrIndicesUp; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
		  TmpSortedIndicesPerSum[Pos2++] = i - j;
		}
	    }
	  
	}
    }
  else
    {
      int* TmpNbrSortedIndicesPerSumDown;
      int** TmpSortedIndicesPerSumDown;
      this->GetAllSkewSymmetricIndices(nbrValues, nbrIndicesDown, TmpNbrSortedIndicesPerSumDown, TmpSortedIndicesPerSumDown);
      int MaxSumDown = ((nbrValues - 1) * nbrIndicesDown) - (((nbrIndicesDown - 1) * nbrIndicesDown) / 2);
      int MinSumDown = (nbrIndicesDown * (nbrIndicesDown - 1)) / 2;
      int MaxSum = MaxSumUp + MaxSumDown;
      int MinSum = MinSumUp + MinSumDown;
      sortedIndicesPerSum = new int* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      for (int i = MinSum; i <= MaxSum; ++i)
	{
	  int TmpMinSumUp = i - MaxSumDown;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i - MinSumDown;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j] * TmpNbrSortedIndicesPerSumDown[i - j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + nbrIndicesDown)];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  int Pos3 = 0;
		  int* TmpSortedIndicesPerSumDown2 = TmpSortedIndicesPerSumDown[i - j];
		  int Lim2 = TmpNbrSortedIndicesPerSumDown[i - j];
		  for (int m = 0; m < Lim2; ++m)
		    {
		      for (int l = 0; l < nbrIndicesUp; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos + l];
		      for (int l = 0; l < nbrIndicesDown; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumDown2[Pos3++];		      
		    }
		  Pos += nbrIndicesUp;
		}
	    }
	}
      for (long i = MinSumDown; i <= MaxSumDown; ++i)
	delete[] TmpSortedIndicesPerSumDown[i];
      delete[] TmpNbrSortedIndicesPerSumDown;
      delete[] TmpSortedIndicesPerSumDown;
    }
  for (long i = MinSumUp; i <= MaxSumUp; ++i)
    delete[] TmpSortedIndicesPerSumUp[i];
  delete[] TmpNbrSortedIndicesPerSumUp;
  delete[] TmpSortedIndicesPerSumUp;
  return NbrElements;
}

// get all indices needed to characterize a  tensor made of two completly symmetric  sets of indices, sorted by the sum of the indices
//
// nbrValues = number of different values an index can have
// nbrIndicesUp = number of indices for the first set of indices (i.e. spin up)
// nbrIndicesDown = number of indices for the first set of indices (i.e. spin down), warning nbrIndicesDown should lower of equal to nbrIndicesUp
// nbrSortedIndicesPerSum = reference on a array where the number of group of indices per each index sum value is stored
// sortedIndicesPerSum = reference on a array where group of indices are stored (first array dimension corresponding to sum of the indices)
// sortedIndicesPerSumSymmetryFactor = reference on a array where symmetry factor (aka inverse of the product of the factorial of the number 
//                                      of time each index appears) are stored (first array dimension corresponding to sum of the indices)
// return value = total number of index groups

long ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::GetAllTwoSetSymmetricIndices (int nbrValues, int nbrIndicesUp, int nbrIndicesDown, int*& nbrSortedIndicesPerSum, int**& sortedIndicesPerSum,
											   double**& sortedIndicesPerSumSymmetryFactor)
{
  int* TmpNbrSortedIndicesPerSumUp;
  int** TmpSortedIndicesPerSumUp;
  double** TmpSortedIndicesPerSumSymmetryFactorUp;
  this->GetAllSymmetricIndices(nbrValues, nbrIndicesUp, TmpNbrSortedIndicesPerSumUp, TmpSortedIndicesPerSumUp,
			       TmpSortedIndicesPerSumSymmetryFactorUp);
  int MaxSumUp = (nbrValues - 1) * nbrIndicesUp;
  int MinSumUp = 0;
  long NbrElements = 0l;
  if (nbrIndicesDown == 1)
    {
      int MaxSum = MaxSumUp + nbrValues - 1;
      int MinSum = MinSumUp;
      nbrSortedIndicesPerSum = new int [MaxSum + 1];
      sortedIndicesPerSum = new int* [MaxSum + 1];
      sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      nbrSortedIndicesPerSum[MinSum] = TmpNbrSortedIndicesPerSumUp[MinSumUp];
      sortedIndicesPerSum[MinSum] = new int [nbrSortedIndicesPerSum[MinSum] * (nbrIndicesUp + 1)];     
      sortedIndicesPerSumSymmetryFactor[MinSum] = new double [nbrSortedIndicesPerSum[MinSum]];
      int Lim = nbrSortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSum = sortedIndicesPerSum[MinSum];
      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MinSum];
      double* TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[MinSum];
      double* TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[MinSum];
	int Pos = 0;
      int Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    {
	      TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	    }
	  TmpSortedIndicesPerSum[Pos2++] = 0;
	  TmpSortedIndicesPerSumSymmetryFactor[i] = TmpSortedIndicesPerSumSymmetryFactorUp2[i];
	}
      nbrSortedIndicesPerSum[MaxSum] = TmpNbrSortedIndicesPerSumUp[MaxSumUp];
      sortedIndicesPerSum[MaxSum] = new int [nbrSortedIndicesPerSum[MaxSum] * (nbrIndicesUp + 1)];     
      sortedIndicesPerSumSymmetryFactor[MaxSum] = new double [nbrSortedIndicesPerSum[MaxSum]];
      Lim = nbrSortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSum = sortedIndicesPerSum[MaxSum];
      TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[MaxSumUp];
      TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[MaxSum];
      TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[MaxSumUp];
      Pos = 0;
      Pos2 = 0;
      for (int i = 0; i < Lim; ++i)
	{
	  for (int j = 0; j < nbrIndicesUp; ++j)
	    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
	  TmpSortedIndicesPerSum[Pos2++] = nbrValues - 1;
	  TmpSortedIndicesPerSumSymmetryFactor[i] = TmpSortedIndicesPerSumSymmetryFactorUp2[i];
	}
      for (int i = MinSum + 1; i < MaxSum; ++i)
	{
	  int TmpMinSumUp = i - nbrValues + 1;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + 1)];
	  sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  double* TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[i];
	  int Pos3 = 0;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      double* TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  for (int l = 0; l < nbrIndicesUp; ++l)
		    TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos++];
		  TmpSortedIndicesPerSum[Pos2++] = i - j;
		  TmpSortedIndicesPerSumSymmetryFactor[Pos3++] = TmpSortedIndicesPerSumSymmetryFactorUp2[k];
		}
	    }
	}
    }
  else
    {
      int* TmpNbrSortedIndicesPerSumDown;
      int** TmpSortedIndicesPerSumDown;
      double** TmpSortedIndicesPerSumSymmetryFactorDown;
      this->GetAllSymmetricIndices(nbrValues, nbrIndicesDown, TmpNbrSortedIndicesPerSumDown, TmpSortedIndicesPerSumDown,
				   TmpSortedIndicesPerSumSymmetryFactorDown);
      int MaxSumDown = (nbrValues - 1) * nbrIndicesDown;
      int MinSumDown = 0;
      int MaxSum = MaxSumUp + MaxSumDown;
      int MinSum = MinSumUp + MinSumDown;
      sortedIndicesPerSum = new int* [MaxSum + 1];
      sortedIndicesPerSumSymmetryFactor = new double* [MaxSum + 1];
      for (int i = 0; i <= MaxSum; ++i)
	nbrSortedIndicesPerSum[i] = 0;
      for (int i = MinSum; i <= MaxSum; ++i)
	{
	  int TmpMinSumUp = i - MaxSumDown;
	  if (TmpMinSumUp < MinSumUp)
	    TmpMinSumUp = MinSumUp;
	  int TmpMaxSumUp = i - MinSumDown;
	  if (TmpMaxSumUp > MaxSumUp)
	    TmpMaxSumUp = MaxSumUp;
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)
	    nbrSortedIndicesPerSum[i] += TmpNbrSortedIndicesPerSumUp[j] * TmpNbrSortedIndicesPerSumDown[i - j];
	  NbrElements += (long) nbrSortedIndicesPerSum[i];
	  sortedIndicesPerSum[i] = new int [nbrSortedIndicesPerSum[i] * (nbrIndicesUp + nbrIndicesDown)];
	  sortedIndicesPerSumSymmetryFactor[i] = new double [nbrSortedIndicesPerSum[i]];
	  int Pos2 = 0;
	  int* TmpSortedIndicesPerSum = sortedIndicesPerSum[i];
	  double* TmpSortedIndicesPerSumSymmetryFactor = sortedIndicesPerSumSymmetryFactor[i];
	  for (int j = TmpMinSumUp; j <= TmpMaxSumUp; ++j)	    
	    {
	      int* TmpSortedIndicesPerSumUp2 = TmpSortedIndicesPerSumUp[j];
	      double* TmpSortedIndicesPerSumSymmetryFactorUp2 = TmpSortedIndicesPerSumSymmetryFactorUp[j];
	      int Lim = TmpNbrSortedIndicesPerSumUp[j];
	      int Pos = 0;
	      int Pos4 = 0;
	      for (int k = 0; k < Lim; ++k)
		{
		  int Pos3 = 0;
		  int* TmpSortedIndicesPerSumDown2 = TmpSortedIndicesPerSumDown[i - j];
		  double* TmpSortedIndicesPerSumSymmetryFactorDown2 = TmpSortedIndicesPerSumSymmetryFactorDown[i - j];
		  int Lim2 = TmpNbrSortedIndicesPerSumDown[i - j];
		  for (int m = 0; m < Lim2; ++m)
		    {
		      for (int l = 0; l < nbrIndicesUp; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumUp2[Pos + l];
		      for (int l = 0; l < nbrIndicesDown; ++l)
			TmpSortedIndicesPerSum[Pos2++] = TmpSortedIndicesPerSumDown2[Pos3++];		      
		      TmpSortedIndicesPerSumSymmetryFactor[Pos4] = (TmpSortedIndicesPerSumSymmetryFactorUp2[k] * 
								    TmpSortedIndicesPerSumSymmetryFactorDown2[m]);
		      ++Pos4;
		    }
		  Pos += nbrIndicesUp;
		}
	    }
	}
      for (long i = MinSumDown; i <= MaxSumDown; ++i)
	delete[] TmpSortedIndicesPerSumDown[i];
      delete[] TmpNbrSortedIndicesPerSumDown;
      delete[] TmpSortedIndicesPerSumDown;
    }
  for (long i = MinSumUp; i <= MaxSumUp; ++i)
    delete[] TmpSortedIndicesPerSumUp[i];
  delete[] TmpNbrSortedIndicesPerSumUp;
  delete[] TmpSortedIndicesPerSumUp;
  return NbrElements;
}


// compute the permutation array for the interaction indices
// 
// permutations = arrays where permuted indices will be stored
// permutationSign = array with the sign of each permutation (initialized only when dealing with fermionic statistics)
// nbodyIndex = order of interaction to consider
// return value = number of permutations

int ParticleOnLatticeWithSpinChernInsulatorNBodyHamiltonian::ComputePermutations(int**& permutations, double*& permutationSign, int nbodyIndex)
{
  int NbrPermutations = 1;
  for (int i = 1; i <= nbodyIndex; ++i)
    NbrPermutations *= i;
  permutations = new int*[NbrPermutations];

  permutations[0] = new int [nbodyIndex];
  for (int i = 0; i < nbodyIndex; ++i)
    permutations[0][i] = i;

  if (this->Particles->GetParticleStatistic() == ParticleOnSphere::FermionicStatistic)
    {
      permutationSign = new double[NbrPermutations]; 
      permutationSign[0] = 1.0;
      double TmpSign = 1.0;
      for (int i = 1; i < NbrPermutations; ++i)
	{
	  permutations[i] = new int [nbodyIndex];
	  for (int j = 0; j < nbodyIndex; ++j)
	    permutations[i][j] = permutations[i - 1][j];
	  int* TmpArrayPerm = permutations[i];
	  int Pos1 = nbodyIndex - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = nbodyIndex - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  int TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  TmpSign *= -1.0;
	  Pos2 = nbodyIndex - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	      TmpSign *= -1.0;
	    }
	  permutationSign[i] = TmpSign;
	}
    }
  else
    {
      for (int i = 1; i < NbrPermutations; ++i)
	{
	  permutations[i] = new int [nbodyIndex];
	  for (int j = 0; j < nbodyIndex; ++j)
	    permutations[i][j] = permutations[i - 1][j];
	  int* TmpArrayPerm = permutations[i];
	  int Pos1 = nbodyIndex - 1;
	  while (TmpArrayPerm[Pos1 - 1] >= TmpArrayPerm[Pos1])
	    --Pos1;
	  --Pos1;
	  int Pos2 = nbodyIndex - 1;      
	  while (TmpArrayPerm[Pos2] <= TmpArrayPerm[Pos1])
	    --Pos2;
	  int TmpIndex = TmpArrayPerm[Pos1];
	  TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	  TmpArrayPerm[Pos2] = TmpIndex;
	  Pos2 = nbodyIndex - 1;   
	  Pos1++;
	  while (Pos1 < Pos2)
	    {
	      TmpIndex = TmpArrayPerm[Pos1];
	      TmpArrayPerm[Pos1] = TmpArrayPerm[Pos2];
	      TmpArrayPerm[Pos2] = TmpIndex;
	      ++Pos1;
	      --Pos2;
	    }
	}
    }
  return NbrPermutations;
}
