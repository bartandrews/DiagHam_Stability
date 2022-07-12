////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of quatum Hall hamiltonian associated              //
//                          to particles on a sphere with                     //
//                                                                            //
//                        last modification : 24/03/2003                      //
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
#include "Hamiltonian/AbstractQHEOnSphereHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <cstring>
#include <sys/time.h>
#include <fstream>


using std::ofstream;
using std::ifstream;
using std::ios;
using std::cout;
using std::endl;
using std::ostream;


// default constructor
//
AbstractQHEOnSphereHamiltonian::AbstractQHEOnSphereHamiltonian()
{
  this->NbrM12Indices = 0;
  this->NbrInteractionPerComponent=0;
  this->LoadBalancingArray=0;
  this->NbrBalancedTasks=0;
  this->FastMultiplicationStep=0;
  this->HermitianSymmetryFlag=false;
}

// destructor
//

AbstractQHEOnSphereHamiltonian::~AbstractQHEOnSphereHamiltonian()
{
  if (this->NbrM12Indices != 0)
    {
      for (int i = 0; i < this->NbrM12Indices; ++i)
	delete[] this->M3Values[i];
      delete[] this->M3Values;
      delete[] this->NbrM3Values;
    }
  if (this->FastMultiplicationFlag == true)
    {
      if (this->DiskStorageFlag == false)
	{
	  long MinIndex;
	  long MaxIndex;
	  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
	  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
	  int ReducedDim = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if ((ReducedDim * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
	    ++ReducedDim;
	  for (int i = 0; i < ReducedDim; ++i)
	    {
	      delete[] this->InteractionPerComponentIndex[i];
	      delete[] this->InteractionPerComponentCoefficient[i];
	    }
	  delete[] this->InteractionPerComponentIndex;
	  delete[] this->InteractionPerComponentCoefficient;
	}
      else
	{
	  remove (this->DiskStorageFileName);
	  delete[] this->DiskStorageFileName;
	}
      delete[] this->NbrInteractionPerComponent;
      this->FastMultiplicationFlag = false;
    }
  if (this->LoadBalancingArray!=NULL)
    delete [] this->LoadBalancingArray;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void AbstractQHEOnSphereHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->InteractionFactors;
  this->Particles = (ParticleOnSphere*) hilbertSpace;
  this->EvaluateInteractionFactors();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* AbstractQHEOnSphereHamiltonian::GetHilbertSpace ()
{
  return this->Particles;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int AbstractQHEOnSphereHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Particles->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void AbstractQHEOnSphereHamiltonian::ShiftHamiltonian (double shift)
{
  this->HamiltonianShift = shift;
}

// ask if Hamiltonian implements hermitian symmetry operations
//
bool AbstractQHEOnSphereHamiltonian::IsHermitian()
{
  return HermitianSymmetryFlag;
}

// ask if Hamiltonian implements conjugate methods
//
bool AbstractQHEOnSphereHamiltonian::IsConjugate()
{
  return false;
}

// symmetrize interaction factors to enable hermitian matrix multiplication
// return = true upon success
bool AbstractQHEOnSphereHamiltonian::HermitianSymmetrizeInteractionFactors()
{
  if (HermitianSymmetryFlag)
    return true;

  if (this->Particles->HaveOrder()==false)
    {
      cout << "Hamiltonian tried to use hermitian symmetry, but this is not implemented in HilbertSpace!"<<endl;
      HermitianSymmetryFlag=false;
      return false;
    }

  cout << "Using hermitian symmetry"<<endl;
  
  int *M = new int[2];
  int *N = new int[2];
  if (this->NbrM12Indices == 0)
    {
      int TmpNbrInteractionFactors = 0;
      int *Flags = new int[NbrInteractionFactors];
      for (int j = 0; j < NbrInteractionFactors; ++j) 
	{
	  M[0] = this->M1Value[j];
	  M[1] = this->M2Value[j];
	  N[0] = this->M3Value[j];
	  N[1] = this->M1Value[j] + this->M2Value[j] - this->M3Value[j];
	  Flags[j] = this->Particles->CheckOrder (M, N, 2);
	  if (Flags[j]>0)
	    ++TmpNbrInteractionFactors;
	  else if (Flags[j]==0)
	    {
	      ++TmpNbrInteractionFactors;
	      this->InteractionFactors[j]*=0.5; // diagonal term: make up for double counting
	    }
	}
      double* TmpInteractionFactors = new double[TmpNbrInteractionFactors];
      int* TmpM1Value = new int[TmpNbrInteractionFactors];
      int* TmpM2Value = new int[TmpNbrInteractionFactors];
      int* TmpM3Value = new int[TmpNbrInteractionFactors];
      int Pos=0;
      for (int j = 0; j < NbrInteractionFactors; ++j)
	{
	  if (Flags[j]>=0)
	    {
	      TmpInteractionFactors[Pos]=InteractionFactors[j];
	      TmpM1Value[Pos]=M1Value[j];
	      TmpM2Value[Pos]=M2Value[j];
	      TmpM3Value[Pos]=M3Value[j];
	      ++Pos;
	    }
	}
      delete [] InteractionFactors;
      delete [] M1Value;
      delete [] M2Value;
      delete [] M3Value;
      this->InteractionFactors = TmpInteractionFactors;
      this->NbrInteractionFactors = TmpNbrInteractionFactors;
      this->M1Value = TmpM1Value;
      this->M2Value = TmpM2Value;
      this->M3Value = TmpM3Value;
      delete [] Flags;
    }
  else
    {
      int SumIndices;
      int OldNbrM3Values;
      int* OldM3Values;
      int TmpNbrM12Values = 0;
      int* M12Flags = new int[this->NbrM12Indices];
      int TmpNbrM3Values;
      int* TmpM3Values;
      int* M3Flags;
      // quick 5count of interaction factors
      int OldNbrInteractionFactors=0;
      for (int m1 = 0; m1 < this->NbrM12Indices; ++m1)
	OldNbrInteractionFactors+=this->NbrM3Values[m1];
      int TmpNbrInteractionFactors=0;
      double *TmpInteractionFactors=new double[OldNbrInteractionFactors];
      int Pos=0;
      for (int m1 = 0; m1 < this->NbrM12Indices; ++m1)
	{
	  M[0]=this->M1Value[m1];
	  M[1]=this->M2Value[m1];
	  SumIndices = this->M1Value[m1] + this->M2Value[m1];
	  OldNbrM3Values = this->NbrM3Values[m1];
	  OldM3Values = this->M3Values[m1];
	  TmpNbrM3Values = 0;
	  M3Flags = new int[OldNbrM3Values];
	  for (int m3 = 0; m3 < OldNbrM3Values; ++m3)
	    {
	      N[0]=OldM3Values[m3];
	      N[1]=SumIndices - OldM3Values[m3];
	      M3Flags[m3] = this->Particles->CheckOrder(M, N, 2);
	      // cout << "Check Order ="<< M3Flags[m3]<<endl;
	      if (M3Flags[m3]>0)
		{
		  ++TmpNbrM3Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=this->InteractionFactors[Pos];
		}
	      else if (M3Flags[m3]==0)
		{
		  ++TmpNbrM3Values;
		  TmpInteractionFactors[TmpNbrInteractionFactors++]=0.5*this->InteractionFactors[Pos];
		}
	      ++Pos;
	    }
	  if (TmpNbrM3Values>0)
	    {
	      //cout << "M1="<<M[0]<<", M2="<<M[1]<<": ";
	      ++TmpNbrM12Values;
	      M12Flags[m1]=1;
	      TmpM3Values = new int[TmpNbrM3Values];
	      int Pos2=0;
	      for (int m3 = 0; m3 < OldNbrM3Values; ++m3)
		if (M3Flags[m3]>=0)
		  {
		    TmpM3Values[Pos2]=OldM3Values[m3];
		    //    cout << " " << TmpM3Values[Pos2];
		    Pos2++;
		  }
	      //cout << endl;
	      delete [] OldM3Values;
	      this->M3Values[m1] = TmpM3Values;
	      this->NbrM3Values[m1] = TmpNbrM3Values;
 	    }
	  else
	    {
	      M12Flags[m1]=-1;
	      delete [] OldM3Values;
	    }
	}
      if (this->NbrM12Indices!=TmpNbrM12Values)
	{
	  int **NewM3Values=new int*[TmpNbrM12Values];
	  int *NewNbrM3Values=new int[TmpNbrM12Values];
	  Pos = 0;
	  for (int m1 = 0; m1 < this->NbrM12Indices; ++m1)
	    if (M12Flags[m1]>0)
	      {
		NewM3Values[Pos]=this->M3Values[m1];
		NewNbrM3Values[Pos]=this->NbrM3Values[m1];
		++Pos;
	      }
	  delete [] this->M3Values;
	  delete [] this->NbrM3Values;
	  this->M3Values=NewM3Values;
	  this->NbrM3Values=NewNbrM3Values;
	}
      // reduce size of table InteractionFactors to match new size, and copy contents
      delete [] this->InteractionFactors;
      this->InteractionFactors = new double[TmpNbrInteractionFactors];
      for (int i=0; i<TmpNbrInteractionFactors; ++i)
	this->InteractionFactors[i]=TmpInteractionFactors[i];
      delete [] TmpInteractionFactors;

//       Pos=0;
//       for (int m1=0; m1<this->NbrM12Indices; ++m1)
// 	for (int m3=0; m3<this->NbrM3Values[m1]; ++m3)
// 	  cout << "New Interaction Factor ["<<this->M1Value[m1]<<", "<<this->M2Value[m1]<<", "<<this->M3Values[m1][m3]<<"]="
// 	       << this->InteractionFactors[Pos]<<endl;
    }

  if (this->OneBodyTermFlag == true)
    {
      int TmpNbrOneBodyInteractionFactors=0;
      int *Flags = new int[this->NbrOneBodyInteractionFactors];
      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	{
	  M[0] = this->OneBodyMValues[j];
	  N[0] = this->OneBodyNValues[j];
	  Flags[j] = this->Particles->CheckOrder(M, N, 1);
	  if (Flags[j]>0)
	    ++TmpNbrOneBodyInteractionFactors;
	  else if (Flags[j]==0)
	    {
	      ++TmpNbrOneBodyInteractionFactors;
	      OneBodyInteractionFactors[j]*=0.5;
	    }
	}
      double *TmpOneBodyInteractionFactors = new double[TmpNbrOneBodyInteractionFactors];
      int *TmpOneBodyMValues = new int[TmpNbrOneBodyInteractionFactors];
      int *TmpOneBodyNValues = new int[TmpNbrOneBodyInteractionFactors];
      int Pos=0;
      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	if (Flags[j]>=0)
          { 
	    TmpOneBodyInteractionFactors[Pos]=this->OneBodyInteractionFactors[j];
	    TmpOneBodyMValues[Pos]=this->OneBodyMValues[j];
	    TmpOneBodyNValues[Pos]=this->OneBodyNValues[j];
            Pos++;
          }
      delete [] this->OneBodyInteractionFactors;
      delete [] this->OneBodyMValues;      
      delete [] this->OneBodyNValues;

      this->OneBodyInteractionFactors = TmpOneBodyInteractionFactors;
      this->OneBodyMValues = TmpOneBodyMValues;
      this->OneBodyNValues = TmpOneBodyNValues;
      this->NbrOneBodyInteractionFactors = TmpNbrOneBodyInteractionFactors;
    }

  delete [] M;
  delete [] N;
  this->HermitianSymmetryFlag=true;
  return true;
}



// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnSphereHamiltonian::MatrixElement (RealVector& V1, RealVector& V2) 
{
  double x = 0.0;
  int dim = this->Particles->GetHilbertSpaceDimension();
  for (int i = 0; i < dim; i++)
    {
    }
  return Complex(x);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex AbstractQHEOnSphereHamiltonian::MatrixElement (ComplexVector& V1, ComplexVector& V2) 
{
  return Complex();
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      //cout << "AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply, FastMultiplicationFlag == false"<<endl;
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->NbrM12Indices == 0)
	{
	  for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		}
	    }
	  m1 = this->M1Value[ReducedNbrInteractionFactors];
	  m2 = this->M2Value[ReducedNbrInteractionFactors];
	  m3 = this->M3Value[ReducedNbrInteractionFactors];
	  TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
      else
	{
	  double Coefficient2;
	  int SumIndices;
	  int TmpNbrM3Values;
	  int* TmpM3Values;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ReducedNbrInteractionFactors = 0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      Coefficient *= vSource[i];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      TmpM3Values = this->M3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	    }
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}

      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      TmpInteraction = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];		  
		}
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
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
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  if (this->L2Operator != 0)
    this->L2Operator->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  Complex Coefficient1;
  if (this->FastMultiplicationFlag == false)
    {
      //cout << "AbstractQHEOnSphereHamiltonian::LowLevelAddMultiply, FastMultiplicationFlag == false"<<endl;
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->NbrM12Indices == 0)
	{
	  for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		}
	    }
	  m1 = this->M1Value[ReducedNbrInteractionFactors];
	  m2 = this->M2Value[ReducedNbrInteractionFactors];
	  m3 = this->M3Value[ReducedNbrInteractionFactors];
	  TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
      else
	{
	  double Coefficient2;
	  int SumIndices;
	  int TmpNbrM3Values;
	  int* TmpM3Values;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ReducedNbrInteractionFactors = 0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient1 = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient1 != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      Coefficient1 = Coefficient1 * vSource[i];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      TmpM3Values = this->M3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    vDestination[Index] += Coefficient1 * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	    }
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}

      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      TmpInteraction = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    vDestination[Index] += Coefficient * TmpInteraction * vSource[i];		  
		}
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  Complex Coefficient;
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
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  if (this->L2Operator != 0)
  {
    cout << "Warning, using method LowLevelAddMultiply, which is not implemented for L2 operators" << endl;
    this->L2Operator->LowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  }
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
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
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	if (this->NbrM12Indices == 0)
	  {
	    for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	      {
		m1 = this->M1Value[j];
		m2 = this->M2Value[j];
		m3 = this->M3Value[j];
		TmpInteraction = this->InteractionFactors[j];
		m4 = m1 + m2 - m3;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		    if (Index < Dim)
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		  }
	      }
	    m1 = this->M1Value[ReducedNbrInteractionFactors];
	    m2 = this->M2Value[ReducedNbrInteractionFactors];
	    m3 = this->M3Value[ReducedNbrInteractionFactors];
	    TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	    m4 = m1 + m2 - m3;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		if (Index < Dim)
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		vDestination[i] += this->HamiltonianShift * vSource[i];
	      }
	  }
	else
	  {
	    double Coefficient2;
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		  {
		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		    if (Coefficient != 0.0)
		      {
			SumIndices = this->M1Value[m1] + this->M2Value[m1];
			Coefficient *= vSource[i];
			TmpNbrM3Values = this->NbrM3Values[m1];
			TmpM3Values = this->M3Values[m1];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			    if (Index < Dim)			
			      vDestination[Index] += Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		  }
		vDestination[i] += this->HamiltonianShift * vSource[i];
	      }
	    
	  }
	if (this->OneBodyTermFlag == true)
	  {
	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	      {
		m1 = this->OneBodyMValues[j];
		m2 = this->OneBodyNValues[j];
		TmpInteraction = this->OneBodyInteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		    if (Index < Dim)
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];		  
		  }
	      }
	  }
      }

  delete TmpParticles;
  return vDestination;
}




// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiplyPartialFastMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  Complex Coefficient1;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
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
      Coefficient1 = vSource[l];
      for (j = 0; j < TmpNbrInteraction; ++j)
	vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient1;
      vDestination[l] += this->HamiltonianShift * Coefficient1;
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {
	if (this->NbrM12Indices == 0)
	  {
	    for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	      {
		m1 = this->M1Value[j];
		m2 = this->M2Value[j];
		m3 = this->M3Value[j];
		TmpInteraction = this->InteractionFactors[j];
		m4 = m1 + m2 - m3;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		    if (Index < Dim)
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		  }
	      }
	    m1 = this->M1Value[ReducedNbrInteractionFactors];
	    m2 = this->M2Value[ReducedNbrInteractionFactors];
	    m3 = this->M3Value[ReducedNbrInteractionFactors];
	    TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	    m4 = m1 + m2 - m3;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		if (Index < Dim)
		  vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		vDestination[i] += this->HamiltonianShift * vSource[i];
	      }
	  }
	else
	  {
	    double Coefficient2;
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		  {
		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		    if (Coefficient != 0.0)
		      {
			SumIndices = this->M1Value[m1] + this->M2Value[m1];
			Coefficient1 *= Coefficient * vSource[i];
			TmpNbrM3Values = this->NbrM3Values[m1];
			TmpM3Values = this->M3Values[m1];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			    if (Index < Dim)			
			      vDestination[Index] += Coefficient1 * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		  }
		vDestination[i] += this->HamiltonianShift * vSource[i];
	      }
	    
	  }
	if (this->OneBodyTermFlag == true)
	  {
	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	      {
		m1 = this->OneBodyMValues[j];
		m2 = this->OneBodyNValues[j];
		TmpInteraction = this->OneBodyInteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		    if (Index < Dim)
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];		  
		  }
	      }
	  }
      }

  delete TmpParticles;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  double Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int l = 0; l < this->BufferSize; ++l)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int i = 0; i < Lim; ++i)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& AbstractQHEOnSphereHamiltonian::LowLevelAddMultiplyDiskStorage(ComplexVector& vSource, ComplexVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  Complex Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int l = 0; l < this->BufferSize; ++l)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int i = 0; i < Lim; ++i)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  return vDestination;
}


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::LowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->NbrM12Indices == 0)
	{
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += Coefficient * vSources[l][i];
		    }
		}
	    }
	}
      else
	{
	  double Coefficient2;
	  int SumIndices;
	  int TmpNbrM3Values;
	  int* TmpM3Values;
	  double* TmpCoefficients = new double[nbrVectors];
	  int ReducedNbrInteractionFactors;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ReducedNbrInteractionFactors = 0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      TmpM3Values = this->M3Values[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    for (int l = 0; l < nbrVectors; ++l)
			      vDestinations[l][Index] += TmpCoefficients[l] * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	    }
	  delete[] TmpCoefficients;
	}
      for (int l = 0; l < nbrVectors; ++l)
	{
	  RealVector& TmpSourceVector = vSources[l];
	  RealVector& TmpDestinationVector = vDestinations[l];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      TmpInteraction = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += Coefficient * vSources[l][i];
		    }
		}
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  double* Coefficient2 = new double [nbrVectors];
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
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
		    vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
   }
  if (this->L2Operator != 0)
    for (int l = 0; l < nbrVectors; ++l)
      this->L2Operator->LowLevelAddMultiply(vSources[l], vDestinations[l], firstComponent, nbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::LowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  double* Coefficient2 = new double [nbrVectors];
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos2;
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
	  vDestinations[k][l] += this->HamiltonianShift * Coefficient2[k];
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
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {	
	if (this->NbrM12Indices == 0)
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][Index] += Coefficient * vSources[l][i];
		    }
		}
	    }
	else
	  {
	    double Coefficient2;
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    double* TmpCoefficients = new double[nbrVectors];
	    int ReducedNbrInteractionFactors;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		  {
		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		    if (Coefficient != 0.0)
		      {
			SumIndices = this->M1Value[m1] + this->M2Value[m1];
			TmpNbrM3Values = this->NbrM3Values[m1];
			TmpM3Values = this->M3Values[m1];
			for (int l = 0; l < nbrVectors; ++l)
			  TmpCoefficients[l] = Coefficient * vSources[l][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			    if (Index < Dim)
			      for (int l = 0; l < nbrVectors; ++l)
				vDestinations[l][Index] += TmpCoefficients[l] * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		  }
	      }
	    delete[] TmpCoefficients;
	  }
	for (int l = 0; l < nbrVectors; ++l)
	  {
	    RealVector& TmpSourceVector = vSources[l];
	    RealVector& TmpDestinationVector = vDestinations[l];
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	  }
	if (this->OneBodyTermFlag == true)
	  {
	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	      {
		m1 = this->OneBodyMValues[j];
		m2 = this->OneBodyNValues[j];
		TmpInteraction = this->OneBodyInteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		    if (Index < Dim)
		      {
			Coefficient *= TmpInteraction;
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][Index] += Coefficient * vSources[l][i];
		      }
		  }
	      }
	  }
      }
  delete[] Coefficient2;
  delete TmpParticles;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::LowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  double Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* Coefficient2 = new double [nbrVectors];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  int Pos;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int m = 0; m < this->BufferSize; ++m)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      Coefficient2[l] = vSources[l][k];
	      vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
	    }
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int m = 0; m < Lim; ++m)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      Coefficient2[l] = vSources[l][k];
	      vDestinations[l][k] += this->HamiltonianShift * Coefficient2[l];
	    }
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		    }
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  delete[] Coefficient2;
  return vDestinations;
}

// low level classes applying conjugate Hamiltonian


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::ConjugateLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									 int firstComponent, int nbrComponent)
{
  //cout << "AbstractQHEOnSphereHamiltonian::ConjugateLowLevelAddMultiply"<<endl;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      //cout << "AbstractQHEOnSphereHamiltonian::ConjugateLowLevelAddMultiply, FastMultiplicationFlag == false"<<endl;
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->NbrM12Indices == 0)
	{
	  for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		}
	    }
	  m1 = this->M1Value[ReducedNbrInteractionFactors];
	  m2 = this->M2Value[ReducedNbrInteractionFactors];
	  m3 = this->M3Value[ReducedNbrInteractionFactors];
	  TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
      else
	{
	  double Coefficient2;
	  double TmpSum;
	  int SumIndices;
	  int TmpNbrM3Values;
	  int* TmpM3Values; 
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ReducedNbrInteractionFactors = 0;
	      TmpSum=0.0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      // previously at this point:
		      // Coefficient *= vSource[i];
		      // to optimize memory access at this point, would require implementation of
		      // a method AdAd, which returns the intermediate state, and AA which returns an index
		      // however, this also entails a different storage of the matrix elements, so for now, we avoid doing this.
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      TmpM3Values = this->M3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)			
			    TmpSum += vSource[Index] * Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	      vDestination[i] += TmpSum;
	    }
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}

      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      TmpInteraction = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		}
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
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
	      double TmpSum=0.0;
	      for (j = 0; j < TmpNbrInteraction; ++j)
		TmpSum +=  TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
	      vDestination[k] += TmpSum + this->HamiltonianShift * vSource[k];
	      k++;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->ConjugateLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->ConjugateLowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  if (this->L2Operator != 0)
    this->L2Operator->ConjugateLowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::ConjugateLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  //cout << "AbstractQHEOnSphereHamiltonian::ConjugateLowLevelAddMultiplyPartialFastMultiply"<<endl;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
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
      double TmpSum=0.0;
      for (j = 0; j < TmpNbrInteraction; ++j)
	TmpSum +=  TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
      vDestination[l] += TmpSum + this->HamiltonianShift * vSource[l];
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {		
	if (this->NbrM12Indices == 0)
	  {
	    for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	      {
		m1 = this->M1Value[j];
		m2 = this->M2Value[j];
		m3 = this->M3Value[j];
		TmpInteraction = this->InteractionFactors[j];
		m4 = m1 + m2 - m3;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		    if (Index < Dim)
		      vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		  }
	      }
	    m1 = this->M1Value[ReducedNbrInteractionFactors];
	    m2 = this->M2Value[ReducedNbrInteractionFactors];
	    m3 = this->M3Value[ReducedNbrInteractionFactors];
	    TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	    m4 = m1 + m2 - m3;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		if (Index < Dim)
		  vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		vDestination[i] += this->HamiltonianShift * vSource[i];
	      }
	  }
	else
	  {
	    double Coefficient2;
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		double TmpSum=0.0;
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		  {
		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		    if (Coefficient != 0.0)
		      {
			SumIndices = this->M1Value[m1] + this->M2Value[m1];
			TmpNbrM3Values = this->NbrM3Values[m1];
			TmpM3Values = this->M3Values[m1];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			    if (Index < Dim)			
			      TmpSum += vSource[Index] * Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		  }
		vDestination[i] += TmpSum + this->HamiltonianShift * vSource[i];
	      }
	    
	  }
	if (this->OneBodyTermFlag == true)
	  {
	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	      {
		m1 = this->OneBodyMValues[j];
		m2 = this->OneBodyNValues[j];
		TmpInteraction = this->OneBodyInteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		    if (Index < Dim)
		      vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		  }
	      }
	  }
      }

  delete TmpParticles;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::ConjugateLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int l = 0; l < this->BufferSize; ++l)
	{
	  double TmpSum=0.0;
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		TmpSum +=  TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += TmpSum + this->HamiltonianShift * vSource[k];
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int i = 0; i < Lim; ++i)
	{
	  double TmpSum = 0.0;
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		TmpSum +=  TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += TmpSum + this->HamiltonianShift *  vSource[k];
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  return vDestination;
}


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::ConjugateLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
									int firstComponent, int nbrComponent)
{
  //cout << "AbstractQHEOnSphereHamiltonian::ConjugateLowLevelMultipleAddMultiply"<<endl;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->NbrM12Indices == 0)
	{
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += Coefficient * vSources[l][Index];
		    }
		}
	    }
	}
      else
	{
	  double Coefficient2;
	  int SumIndices;
	  int TmpNbrM3Values;
	  int* TmpM3Values;
	  double* TmpSum = new double[nbrVectors];
	  int ReducedNbrInteractionFactors;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ReducedNbrInteractionFactors = 0;
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSum[l] = 0.0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      TmpM3Values = this->M3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    {
			      Coefficient2 *= this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient;
			      for (int l = 0; l < nbrVectors; ++l)
				TmpSum[l] += Coefficient2 * vSources[l][Index];
			      ++ReducedNbrInteractionFactors;
			    }
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i]+=TmpSum[l];
	    }
	  delete[] TmpSum;
	}
      for (int l = 0; l < nbrVectors; ++l)
	{
	  RealVector& TmpSourceVector = vSources[l];
	  RealVector& TmpDestinationVector = vDestinations[l];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      TmpInteraction = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += Coefficient * vSources[l][Index];
		    }
		}
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  double* TmpSum = new double [nbrVectors];
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
	  int j;
	  int Index;
	  int TmpNbrInteraction;
	  int k = firstComponent;
	  firstComponent -= this->PrecalculationShift;
	  LastComponent -= this->PrecalculationShift;
	  cout << "before " << firstComponent << " " << LastComponent << "  ";
	  for (int i = firstComponent; i < (firstComponent + 10); ++i)
	    {
	      cout << (i + this->PrecalculationShift) << " : ";
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  cout << vDestinations[l][i + this->PrecalculationShift] << " ";
		}
	      cout << "| ";
	    }
	  cout << endl;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSum[l] = 0.0;
	      TmpNbrInteraction = this->NbrInteractionPerComponent[i];
	      TmpIndexArray = this->InteractionPerComponentIndex[i];
	      TmpCoefficientArray = this->InteractionPerComponentCoefficient[i];
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Index = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    TmpSum[l] +=  Coefficient * vSources[l][Index];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * vSources[l][k];
		++k;
	    }
	  delete[] TmpSum;
	  cout << "after " << k << "  ; ";
	  for (int i = firstComponent; i < (firstComponent + 10); ++i)
	    {
	      cout << (i + this->PrecalculationShift) << " : ";
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  cout << vDestinations[l][i + this->PrecalculationShift] << " ";
		}
	      cout << "| ";
	    }
	  cout << endl;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->LowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->LowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
   }
  if (this->L2Operator != 0)
    for (int l = 0; l < nbrVectors; ++l)
      this->L2Operator->LowLevelAddMultiply(vSources[l], vDestinations[l], firstComponent, nbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  //cout << "AbstractQHEOnSphereHamiltonian::ConjugateLowLevelMultipleAddMultiplyPartialFastMultiply"<<endl;
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  double *TmpSum = new double[nbrVectors];
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  // double* Coefficient2 = new double [nbrVectors];
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos2;
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
      for (int n = 0; n < nbrVectors; ++n)
	TmpSum[n] = 0.0;
      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    TmpSum[k] += Coefficient  * vSources[k][Pos2];
	}
      for (int k = 0; k < nbrVectors; ++k)
	vDestinations[k][l] += TmpSum[k] + this->HamiltonianShift * vSources[k][l];
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {	
	if (this->NbrM12Indices == 0)
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			vDestinations[l][i] += Coefficient * vSources[l][Index];
		    }
		}
	    }
	else
	  {
	    double Coefficient2;
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    int ReducedNbrInteractionFactors;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		ReducedNbrInteractionFactors = 0;
		for (int l = 0; l < nbrVectors; ++l)
		  TmpSum[l] = 0.0; 
		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		  {
		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		    if (Coefficient != 0.0)
		      {
			SumIndices = this->M1Value[m1] + this->M2Value[m1];
			TmpNbrM3Values = this->NbrM3Values[m1];
			TmpM3Values = this->M3Values[m1];
			
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			    if (Index < Dim)
			      {
				Coefficient2 *= this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient;
				for (int l = 0; l < nbrVectors; ++l)
				  TmpSum[l] += Coefficient2 * vSources[l][Index];
			      }
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		  }
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][i]+=TmpSum[l];
	      }
	    delete[] TmpSum;
	  }
	for (int l = 0; l < nbrVectors; ++l)
	  {
	    RealVector& TmpSourceVector = vSources[l];
	    RealVector& TmpDestinationVector = vDestinations[l];
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	  }
	if (this->OneBodyTermFlag == true)
	  {
	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	      {
		m1 = this->OneBodyMValues[j];
		m2 = this->OneBodyNValues[j];
		TmpInteraction = this->OneBodyInteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		    if (Index < Dim)
		      {
			Coefficient *= TmpInteraction;
			for (int l = 0; l < nbrVectors; ++l)
			  vDestinations[l][i] += Coefficient * vSources[l][Index];
		      }
		  }
	      }
	  }
      }
  delete[] TmpSum;
  delete TmpParticles;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::ConjugateLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  double Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  double *TmpSum = new double[nbrVectors]; 
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  int Pos;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;


  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int m = 0; m < this->BufferSize; ++m)
	{
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l]=0.0;
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    TmpSum[l] +=  Coefficient * vSources[l][Pos];
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * vSources[l][k];
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int m = 0; m < Lim; ++m)
	{
	  for (int l = 0; l < nbrVectors; ++l)
	    TmpSum[l]=0.0;
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    TmpSum[l] +=  Coefficient * vSources[l][Pos];
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * vSources[l][k];
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] TmpSum;
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  return vDestinations;
}

// Methods applying both the Hamiltonian and its Hermitian conjugate at the same time

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::HermitianLowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
								int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->NbrM12Indices == 0)
	{
	  for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      vDestination[Index] += Coefficient * vSource[i];
		      vDestination[i] += Coefficient * vSource[Index];
		    }
		}
	    }
	  m1 = this->M1Value[ReducedNbrInteractionFactors];
	  m2 = this->M2Value[ReducedNbrInteractionFactors];
	  m3 = this->M3Value[ReducedNbrInteractionFactors];
	  TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	  m4 = m1 + m2 - m3;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < Dim)
		{
		  Coefficient *= TmpInteraction;
		  vDestination[Index] += Coefficient * vSource[i];
		  vDestination[i] += Coefficient * vSource[Index];
		}
	      vDestination[i] += this->HamiltonianShift * vSource[i];
	    }
	}
      else
	{
	  double TmpSum;
	  double Coefficient2;
	  double Coefficient3;
	  int SumIndices;
	  int TmpNbrM3Values;
	  int* TmpM3Values;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      ReducedNbrInteractionFactors = 0;
	      TmpSum=0.0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      Coefficient3 = Coefficient*vSource[i];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      TmpM3Values = this->M3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < Dim)
			    {
			      TmpSum += vSource[Index] * Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			      vDestination[Index] += this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2 * Coefficient3;
			    }
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	      vDestination[i] += TmpSum;
	    }
	  for (int i = firstComponent; i < LastComponent; ++i)
	    vDestination[i] += this->HamiltonianShift * vSource[i];
	}

      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      TmpInteraction = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    {
		      vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
		      vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		      //TmpParticles->PrintState(cout, i);
                      //cout << " ;  " << m1 << " " << m2 << " --> "; TmpParticles->PrintState(cout, Index); cout << " Coeff= "<<Coefficient << " "<< TmpInteraction<<endl;  
                      //exit(2);
		    }
		}
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
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
	      double TmpSum=0.0;
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  TmpSum +=  TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
		}
	      vDestination[k++] += TmpSum + this->HamiltonianShift * Coefficient;
	    }
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelAddMultiplyPartialFastMultiply(vSource, vDestination, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelAddMultiplyDiskStorage(vSource, vDestination, firstComponent, nbrComponent);
	    }
	}
    }
  if (this->L2Operator != 0)
    this->L2Operator->HermitianLowLevelAddMultiply(vSource, vDestination, firstComponent, nbrComponent);
  return vDestination;
}


// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::HermitianLowLevelAddMultiplyPartialFastMultiply(RealVector& vSource, RealVector& vDestination, 
										   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  double TmpSum;
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
      TmpSum=0.0;
      TmpNbrInteraction = this->NbrInteractionPerComponent[Pos];
      TmpIndexArray = this->InteractionPerComponentIndex[Pos];
      TmpCoefficientArray = this->InteractionPerComponentCoefficient[Pos];
      Coefficient = vSource[l];
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  // xxx segfault here: reading non-initialized
	  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
	  TmpSum +=  TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
	}
      vDestination[l] += TmpSum + this->HamiltonianShift * Coefficient;
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  int ReducedNbrInteractionFactors = this->NbrInteractionFactors - 1;  
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {		
	if (this->NbrM12Indices == 0)
	  {
	    for (int j = 0; j < ReducedNbrInteractionFactors; ++j) 
	      {
		m1 = this->M1Value[j];
		m2 = this->M2Value[j];
		m3 = this->M3Value[j];
		TmpInteraction = this->InteractionFactors[j];
		m4 = m1 + m2 - m3;
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		    if (Index < Dim)
		      {
			Coefficient *= TmpInteraction;
			vDestination[Index] += Coefficient * vSource[i];
			vDestination[i] += Coefficient * vSource[Index];
		      }
		  }
	      }
	    m1 = this->M1Value[ReducedNbrInteractionFactors];
	    m2 = this->M2Value[ReducedNbrInteractionFactors];
	    m3 = this->M3Value[ReducedNbrInteractionFactors];
	    TmpInteraction = this->InteractionFactors[ReducedNbrInteractionFactors];
	    m4 = m1 + m2 - m3;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		if (Index < Dim)
		  {
		    Coefficient *= TmpInteraction;
		    vDestination[Index] += Coefficient * vSource[i];
		    vDestination[i] += Coefficient * vSource[Index] + this->HamiltonianShift * vSource[i];
		  }
		else
		  vDestination[i] += this->HamiltonianShift * vSource[i];
	      }
	  }
	else
	  {
	    double Coefficient2;
	    double Coefficient3;
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		TmpSum=0.0;
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		  {
		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		    if (Coefficient != 0.0)
		      {
			SumIndices = this->M1Value[m1] + this->M2Value[m1];
			Coefficient3 = Coefficient * vSource[i];
			TmpNbrM3Values = this->NbrM3Values[m1];
			TmpM3Values = this->M3Values[m1];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			    if (Index < Dim)
			      {
				vDestination[Index] += this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2 * Coefficient3;
				TmpSum += this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient * Coefficient2 * vSource[Index];
			      }
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		  }
		vDestination[i] += TmpSum + this->HamiltonianShift * vSource[i];
	      }
	    
	  }
	if (this->OneBodyTermFlag == true)
	  {
	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	      {
		m1 = this->OneBodyMValues[j];
		m2 = this->OneBodyNValues[j];
		TmpInteraction = this->OneBodyInteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		    if (Index < Dim)
		      {
			vDestination[Index] += Coefficient * TmpInteraction * vSource[i];
			vDestination[i] += Coefficient * TmpInteraction * vSource[Index];
		      }
		  }
	      }
	  }
      }

  delete TmpParticles;
  return vDestination;
}

// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
// using disk storage option
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& AbstractQHEOnSphereHamiltonian::HermitianLowLevelAddMultiplyDiskStorage(RealVector& vSource, RealVector& vDestination, 
									   int firstComponent, int nbrComponent)
{
  double Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  double TmpSum;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int l = 0; l < this->BufferSize; ++l)
	{
	  TmpSum=0.0;
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  TmpSum += TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += TmpSum + this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int i = 0; i < Lim; ++i)
	{
	  TmpSum=0.0;
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  Coefficient = vSource[k];
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  vDestination[TmpIndexArray[j]] +=  TmpCoefficientArray[j] * Coefficient;
		  TmpSum += TmpCoefficientArray[j] * vSource[TmpIndexArray[j]];
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  vDestination[k] += TmpSum + this->HamiltonianShift * Coefficient;
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  return vDestination;
}


// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::HermitianLowLevelMultipleAddMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										 int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  if (this->FastMultiplicationFlag == false)
    {
      int Index;
      int m1;
      int m2;
      int m3;
      int m4;
      double TmpInteraction;
      ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
      if (this->NbrM12Indices == 0)
	{
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			{
			  vDestinations[l][Index] += Coefficient * vSources[l][i];
			  vDestinations[l][i] += Coefficient * vSources[l][Index];
			}
		    }
		}
	    }
	}
      else
	{
	  double Coefficient2;
	  int SumIndices;
	  int TmpNbrM3Values;
	  int* TmpM3Values;
	  double* TmpCoefficients = new double[nbrVectors];
	  double* TmpSum = new double[nbrVectors];
	  int ReducedNbrInteractionFactors;
	  for (int i = firstComponent; i < LastComponent; ++i)
	    {
	      for (int l = 0; l < nbrVectors; ++l)
		TmpSum[l] = 0.0;
	      ReducedNbrInteractionFactors = 0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      TmpM3Values = this->M3Values[m1];
		      for (int l = 0; l < nbrVectors; ++l)
			TmpCoefficients[l] = Coefficient * vSources[l][i];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index <= Dim)
			    {
			      Coefficient2 *= this->InteractionFactors[ReducedNbrInteractionFactors];
			      for (int l = 0; l < nbrVectors; ++l)
				{
				  vDestinations[l][Index] += TmpCoefficients[l] * Coefficient2;
				  TmpSum[l] += Coefficient * Coefficient2 * vSources[l][Index];
				}
			    }
			  ++ReducedNbrInteractionFactors;
			}
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}
	      for (int l = 0; l < nbrVectors; ++l)
		vDestinations[l][i] += TmpSum[l];
	    }
	  delete[] TmpCoefficients;
	  delete[] TmpSum;
	}
      for (int l = 0; l < nbrVectors; ++l)
	{
	  RealVector& TmpSourceVector = vSources[l];
	  RealVector& TmpDestinationVector = vDestinations[l];
	  for (int i = firstComponent; i < LastComponent; ++i)
	    TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	}
      if (this->OneBodyTermFlag == true)
	{
	  for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	    {
	      m1 = this->OneBodyMValues[j];
	      m2 = this->OneBodyNValues[j];
	      TmpInteraction = this->OneBodyInteractionFactors[j];
	      for (int i = firstComponent; i < LastComponent; ++i)
		{
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			{
			  vDestinations[l][Index] += Coefficient * vSources[l][i];
			  vDestinations[l][i] += Coefficient * vSources[l][Index];
			}
		    }
		}
	    }
	}
      delete TmpParticles;
    }
  else
    {
      if (this->FastMultiplicationStep == 1)
	{
	  double* Coefficient2 = new double [nbrVectors];
	  double* TmpSum = new double [nbrVectors];
	  int* TmpIndexArray;
	  double* TmpCoefficientArray; 
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
		  TmpSum[l] = 0.0;
		  Coefficient2[l] = vSources[l][k];
		}
	      for (j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		      TmpSum[l] += Coefficient * vSources[l][Pos];
		    }
		}
	      for (int l = 0; l < nbrVectors; ++l)
		{
		  vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * Coefficient2[l];
		}
	      ++k;
	    }
	  delete[] Coefficient2;
	  delete[] TmpSum;
	}
      else
	{
	  if (this->DiskStorageFlag == false)
	    {
	      this->HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	  else
	    {
	      this->HermitianLowLevelMultipleAddMultiplyDiskStorage(vSources, vDestinations, nbrVectors, firstComponent, nbrComponent);
	    }
	}
   }
  if (this->L2Operator != 0)
    for (int l = 0; l < nbrVectors; ++l)
      this->L2Operator->HermitianLowLevelAddMultiply(vSources[l], vDestinations[l], firstComponent, nbrComponent);
  return vDestinations;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using partial fast multiply option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::HermitianLowLevelMultipleAddMultiplyPartialFastMultiply(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
											   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int Dim = this->Particles->GetHilbertSpaceDimension();
  double Coefficient;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int* TmpIndexArray;
  double* TmpCoefficientArray; 
  double* Coefficient2 = new double [nbrVectors];
  double* TmpSum = new double [nbrVectors];
  int j;
  int TmpNbrInteraction;
  firstComponent -= this->PrecalculationShift;
  LastComponent -= this->PrecalculationShift;
  int Pos2;
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
	  TmpSum[k] = 0.0;
	}
      for (j = 0; j < TmpNbrInteraction; ++j)
	{
	  Pos2 = TmpIndexArray[j];
	  Coefficient = TmpCoefficientArray[j];
	  for (int k = 0; k < nbrVectors; ++k)
	    {
	      vDestinations[k][Pos2] += Coefficient  * Coefficient2[k];
	      TmpSum[k] += Coefficient * vSources[k][Pos2];
	    }
	}
      for (int k = 0; k < nbrVectors; ++k)
	vDestinations[k][l] += TmpSum[k] + this->HamiltonianShift * Coefficient2[k];
      l += this->FastMultiplicationStep;
      ++Pos;
    }
  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double TmpInteraction;
  firstComponent += this->PrecalculationShift;
  LastComponent += this->PrecalculationShift;
  for (int k = 0; k < this->FastMultiplicationStep; ++k)
    if (PosMod != k)
      {	
	if (this->NbrM12Indices == 0)
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      TmpInteraction = this->InteractionFactors[j];
	      m4 = m1 + m2 - m3;
	      for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		{
		  Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < Dim)
		    {
		      Coefficient *= TmpInteraction;
		      for (int l = 0; l < nbrVectors; ++l)
			{
			  vDestinations[l][Index] += Coefficient * vSources[l][i];
			  vDestinations[l][i] += Coefficient * vSources[l][Index];
			}
		    }
		}
	    }
	else
	  {
	    double Coefficient2;
	    int SumIndices;
	    int TmpNbrM3Values;
	    int* TmpM3Values;
	    double* TmpCoefficients = new double[nbrVectors];
	    int ReducedNbrInteractionFactors;
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      {
		for (int l = 0; l < nbrVectors; ++l)
		  TmpSum[l] = 0.0;
		ReducedNbrInteractionFactors = 0;
		for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		  {
		    Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		    if (Coefficient != 0.0)
		      {
			SumIndices = this->M1Value[m1] + this->M2Value[m1];
			TmpNbrM3Values = this->NbrM3Values[m1];
			TmpM3Values = this->M3Values[m1];
			for (int l = 0; l < nbrVectors; ++l)
			  TmpCoefficients[l] = Coefficient * vSources[l][i];
			for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			  {
			    Index = TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			    if (Index < Dim)
			      {
				Coefficient2 *= this->InteractionFactors[ReducedNbrInteractionFactors];
				for (int l = 0; l < nbrVectors; ++l)
				  {
				    vDestinations[l][Index] += TmpCoefficients[l] * Coefficient2;
				    TmpSum[l] += Coefficient2 * Coefficient * vSources[l][Index];
				  }
			      }
			    ++ReducedNbrInteractionFactors;
			  }
		      }
		    else
		      ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		  }
		for (int l = 0; l < nbrVectors; ++l)
		  vDestinations[l][i] += TmpSum[l];
	      }
	    delete[] TmpCoefficients;
	  }
	for (int l = 0; l < nbrVectors; ++l)
	  {
	    RealVector& TmpSourceVector = vSources[l];
	    RealVector& TmpDestinationVector = vDestinations[l];
	    for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
	      TmpDestinationVector[i] += this->HamiltonianShift * TmpSourceVector[i];
	  }
	if (this->OneBodyTermFlag == true)
	  {
	    for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
	      {
		m1 = this->OneBodyMValues[j];
		m2 = this->OneBodyNValues[j];
		TmpInteraction = this->OneBodyInteractionFactors[j];
		for (int i = firstComponent + k; i < LastComponent; i += this->FastMultiplicationStep)
		  {
		    Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		    if (Index < Dim)
		      {
			Coefficient *= TmpInteraction;
			for (int l = 0; l < nbrVectors; ++l)
			  {
			    vDestinations[l][Index] += Coefficient * vSources[l][i];
			    vDestinations[l][i] += Coefficient * vSources[l][Index];
			  }
		      }
		  }
	      }
	  }
      }
  delete[] Coefficient2;
  delete[] TmpSum;
  delete TmpParticles;
  return vDestinations;
}

// multiply a et of vectors by the current hamiltonian for a given range of indices 
// and add result to another et of vectors, low level function (no architecture optimization)
// using disk storage option
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

RealVector* AbstractQHEOnSphereHamiltonian::HermitianLowLevelMultipleAddMultiplyDiskStorage(RealVector* vSources, RealVector* vDestinations, int nbrVectors, 
										   int firstComponent, int nbrComponent)
{
  double Coefficient;
  int* BufferIndexArray = new int [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* BufferCoefficientArray  = new double [this->BufferSize * this->MaxNbrInteractionPerComponent];
  double* Coefficient2 = new double [nbrVectors];
  double* TmpSum = new double [nbrVectors];
  int TmpNbrIteration = nbrComponent / this->BufferSize;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int TmpNbrInteraction;
  int k = firstComponent;
  int EffectiveHilbertSpaceDimension;
  int Pos;
  firstComponent -= this->PrecalculationShift;
  
  ifstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::in);
  File.read ((char*) &EffectiveHilbertSpaceDimension, sizeof(int));
  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    FileJump += (long) this->NbrInteractionPerComponent[i];
  FileJump *= sizeof(int);
  long FileOffset = 0;
  for (int i = this->DiskStorageStart; i < firstComponent; ++i)
    FileOffset += this->NbrInteractionPerComponent[i];
  File.seekg (((FileOffset + EffectiveHilbertSpaceDimension + 1) * sizeof(int)), ios::cur);
  FileJump += (sizeof(double) - sizeof(int)) * FileOffset;
  
  for (int i = 0; i < TmpNbrIteration; ++i)
    {
      int TmpPos = firstComponent;
      long ReadBlockSize = 0;
      for (int j = 0; j < this->BufferSize; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int m = 0; m < this->BufferSize; ++m)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      Coefficient2[l] = vSources[l][k];
	      TmpSum[l] = 0.0;
	    }
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		      TmpSum[l] += Coefficient * vSources[l][Pos];
		    }
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * Coefficient2[l];
	  ++k;
	  ++firstComponent;
	}
    }
  
  if ((TmpNbrIteration * this->BufferSize) != nbrComponent)
    {
      int TmpPos = firstComponent;
      int Lim =  nbrComponent % this->BufferSize;
      long ReadBlockSize = 0;
      for (int j = 0; j < Lim; ++j)
	{
	  ReadBlockSize += this->NbrInteractionPerComponent[TmpPos];
	  ++TmpPos;
	}		  
      File.read((char*) BufferIndexArray, sizeof(int) * ReadBlockSize);
      FileJump -= sizeof(int) * ReadBlockSize;
      File.seekg (FileJump, ios::cur);
      File.read((char*) BufferCoefficientArray, sizeof(double) * ReadBlockSize);		      
      FileJump += sizeof(double) * ReadBlockSize;
      File.seekg (-FileJump, ios::cur);
      
      TmpIndexArray = BufferIndexArray;
      TmpCoefficientArray = BufferCoefficientArray;
      for (int m = 0; m < Lim; ++m)
	{
	  TmpNbrInteraction = this->NbrInteractionPerComponent[firstComponent];
	  for (int l = 0; l < nbrVectors; ++l)
	    {
	      Coefficient2[l] = vSources[l][k];
	      TmpSum[l] = 0.0;
	    }
	  if (TmpNbrInteraction > 0)
	    {
	      for (int j = 0; j < TmpNbrInteraction; ++j)
		{
		  Pos = TmpIndexArray[j];
		  Coefficient = TmpCoefficientArray[j];
		  for (int l = 0; l < nbrVectors; ++l)
		    {
		      vDestinations[l][Pos] +=  Coefficient * Coefficient2[l];
		      TmpSum[l] += Coefficient * vSources[l][Pos];
		    }
		}
	      TmpIndexArray += TmpNbrInteraction;
	      TmpCoefficientArray += TmpNbrInteraction;
	    }
	  for (int l = 0; l < nbrVectors; ++l)
	    vDestinations[l][k] += TmpSum[l] + this->HamiltonianShift * Coefficient2[l];
	  ++k;
	  ++firstComponent;
	}
    }
  
  File.close();
  delete[] BufferIndexArray;
  delete[] BufferCoefficientArray;
  delete[] Coefficient2;
  delete[] TmpSum;
  return vDestinations;
}


 
// return a list of left interaction operators
//
// return value = list of left interaction operators

List<Matrix*> AbstractQHEOnSphereHamiltonian::LeftInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// return a list of right interaction operators
//
// return value = list of right interaction operators

List<Matrix*> AbstractQHEOnSphereHamiltonian::RightInteractionOperators()
{
  List<Matrix*> TmpList;
  return TmpList;
}

// get the preferred distribution over parallel execution in N tasks for parallel Hamiltonian-Vector multiplication
//
// nbrThreads = number of threads requested
// segmentIndices = array returning the reference to an array of the first index of each of the segments
// return value = true if no error occured

bool AbstractQHEOnSphereHamiltonian::GetLoadBalancing(int nbrTasks, long* &segmentIndices)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;

  if ((this->NbrInteractionPerComponent != 0) && (this->FastMultiplicationStep != 0))
    {
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;

      if ((this->LoadBalancingArray == 0) || (this->NbrBalancedTasks != nbrTasks))
	{
	  if (LoadBalancingArray!=0)
	    delete [] LoadBalancingArray;
	  long *SegmentSize = new long[nbrTasks];
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  this->NbrBalancedTasks = nbrTasks;
	  long TmpNbrElement = 0;
	  for (int i=0; i < ReducedSpaceDimension; ++i)
	    TmpNbrElement += this->NbrInteractionPerComponent[i];
	  long TmpNbrPerSegment = TmpNbrElement / nbrTasks;
	  TmpNbrElement = 0;
	  int Pos=0;
	  this->LoadBalancingArray[0] = MinIndex;
	  for (int i = 0; i < ReducedSpaceDimension; ++i)
	    {
	      TmpNbrElement += this->NbrInteractionPerComponent[i];
	      if (TmpNbrElement > TmpNbrPerSegment)
		{
		  SegmentSize[Pos] = TmpNbrElement;
		  this->LoadBalancingArray[Pos + 1]= MinIndex + (i * this->FastMultiplicationStep);
		  TmpNbrElement = 0;
		  ++Pos;
		}
	    }
	  while (Pos < (nbrTasks - 1))
	    {
	      LoadBalancingArray[Pos + 1] = MaxIndex + 1;
	      SegmentSize[Pos] = 0;
	      ++Pos;
	    }
	  LoadBalancingArray[nbrTasks] = MaxIndex + 1;
	  SegmentSize[nbrTasks - 1] = TmpNbrElement;
	  
	  cout << "LoadBalancingArray=[ ("<< LoadBalancingArray[1] - LoadBalancingArray[0] <<", "<<SegmentSize[0]<<")";
	  for (int i = 1; i < nbrTasks; ++i)
	    cout <<" ("<< LoadBalancingArray[i+1] - LoadBalancingArray[i] << ", " << SegmentSize[i] << ")";
	  cout << "]"<< endl;
	  delete [] SegmentSize;
	}
    }
  else
    {
      if ((LoadBalancingArray == 0) || (NbrBalancedTasks != nbrTasks))
	{
	  if (LoadBalancingArray!=0)
	    delete[] LoadBalancingArray;
	  this->LoadBalancingArray = new long[nbrTasks+1];
	  
	  int Step = EffectiveHilbertSpaceDimension / nbrTasks;
	  this->LoadBalancingArray[0] = MinIndex;
	  for (int i = 0; i < nbrTasks; ++i)
	    LoadBalancingArray[i]= MinIndex + (i * Step);
	  LoadBalancingArray[nbrTasks] = MaxIndex + 1;
	  NbrBalancedTasks = nbrTasks;
	}
    }
  segmentIndices = LoadBalancingArray;
  return true;
}


// test the amount of memory needed for fast multiplication algorithm
//
// allowedMemory = amount of memory that cam be allocated for fast multiplication
// return value = amount of memory needed

long AbstractQHEOnSphereHamiltonian::FastMultiplicationMemory(long allowedMemory)
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  
  this->NbrInteractionPerComponent = new int [EffectiveHilbertSpaceDimension];
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    this->NbrInteractionPerComponent[i] = 0;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;

  QHEParticlePrecalculationOperation Operation(this);
  Operation.ApplyOperation(this->Architecture);

  // reset load-balancing
  if (this->LoadBalancingArray!=0)
    delete [] this->LoadBalancingArray;
  this->LoadBalancingArray=0;
  this->NbrBalancedTasks=0;

  long Memory = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    Memory += this->NbrInteractionPerComponent[i];  

  cout << "nbr interaction = " << Memory << endl;
  long TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension;
  if ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
    {
      this->FastMultiplicationStep = 1;
      int ReducedSpaceDimension  = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
      while ((TmpMemory < 0) || ((TmpMemory / ((int) (sizeof (int) + sizeof(double)))) < Memory))
	{
	  ++this->FastMultiplicationStep;
	  ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
	  if (EffectiveHilbertSpaceDimension != (ReducedSpaceDimension * this->FastMultiplicationStep))
	    ++ReducedSpaceDimension;
	  TmpMemory = allowedMemory - (sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension;
	  Memory = 0;
	  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
	    Memory += this->NbrInteractionPerComponent[i];
	}
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * ReducedSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
      long ResidualMemory = allowedMemory - Memory;
      if (ResidualMemory > 0)
	{
	  if (this->DiskStorageFlag == false)
	    {
	      int TotalReducedSpaceDimension = ReducedSpaceDimension;
	      int* TmpNbrInteractionPerComponent = new int [TotalReducedSpaceDimension];
	      int i = 0;
	      int Pos = 0;
	      for (; i < ReducedSpaceDimension; ++i)
		{
		  TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[Pos];
		  Pos += this->FastMultiplicationStep;
		}
	      delete[] this->NbrInteractionPerComponent;
	      this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	    }
	}
      else
	if (this->DiskStorageFlag == false)
	  {
	    int* TmpNbrInteractionPerComponent = new int [ReducedSpaceDimension];
	    for (int i = 0; i < ReducedSpaceDimension; ++i)
	      TmpNbrInteractionPerComponent[i] = this->NbrInteractionPerComponent[i * this->FastMultiplicationStep];
	    delete[] this->NbrInteractionPerComponent;
	    this->NbrInteractionPerComponent = TmpNbrInteractionPerComponent;
	  }
    }
  else
    {
      Memory = ((sizeof (int*) + sizeof (int) + sizeof(double*)) * EffectiveHilbertSpaceDimension) + (Memory * (sizeof (int) + sizeof(double)));
      this->FastMultiplicationStep = 1;
    }

  cout << "reduction factor=" << this->FastMultiplicationStep << endl;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  return Memory;
}

// test the amount of memory needed for fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = number of components that has to be precalcualted
// return value = number of non-zero matrix element

long AbstractQHEOnSphereHamiltonian::PartialFastMultiplicationMemory(int firstComponent, int nbrComponent)
{
  int Index;
  double Coefficient;
  long Memory = 0;
  int m1;
  int m2;
  int m3;
  int m4;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();
  int LastComponent = nbrComponent + firstComponent;
  if (this->NbrM12Indices == 0)
    {
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  for (int j = 0; j < this->NbrInteractionFactors; ++j) 
	    {
	      m1 = this->M1Value[j];
	      m2 = this->M2Value[j];
	      m3 = this->M3Value[j];
	      m4 = m1 + m2 - m3;
	      Index = TmpParticles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
	      if (Index < this->Particles->GetHilbertSpaceDimension())
		{
		  ++Memory;
		  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		}
	    }    
	  if (this->OneBodyTermFlag == true)
	    {
	      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
		{
		  m1 = this->OneBodyMValues[j];
		  m2 = this->OneBodyNValues[j];
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }
	}
    }
  else
    {
      int SumIndices;
      int TmpNbrM3Values;
      int* TmpM3Values;
      for (int i = firstComponent; i < LastComponent; ++i)
	{
	  for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
	    {
	      Coefficient = TmpParticles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
	      if (Coefficient != 0.0)
		{
		  SumIndices = this->M1Value[m1] + this->M2Value[m1];
		  TmpM3Values = this->M3Values[m1];
		  TmpNbrM3Values = this->NbrM3Values[m1];
		  for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
		    {
		      if (TmpParticles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient) < this->Particles->GetHilbertSpaceDimension())
			{
			  ++Memory;
			  ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
			}
		    }    
		}
	    }
	  if (this->OneBodyTermFlag == true)
	    {
	      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
		{
		  m1 = this->OneBodyMValues[j];
		  m2 = this->OneBodyNValues[j];
		  Index = TmpParticles->AdA(i, m1, m2, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      ++Memory;
		      ++this->NbrInteractionPerComponent[i - this->PrecalculationShift];
		    }
		}
	    }
	}
    }
  delete TmpParticles;

  return Memory;
}

// enable fast multiplication algorithm
//

void AbstractQHEOnSphereHamiltonian::EnableFastMultiplication()
{
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  int ReducedSpaceDimension = EffectiveHilbertSpaceDimension / this->FastMultiplicationStep;
  if ((ReducedSpaceDimension * this->FastMultiplicationStep) != EffectiveHilbertSpaceDimension)
    ++ReducedSpaceDimension;
  this->InteractionPerComponentIndex = new int* [ReducedSpaceDimension];
  this->InteractionPerComponentCoefficient = new double* [ReducedSpaceDimension];

  // allocate all memory at the outset:
  long TotalPos = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; i += this->FastMultiplicationStep)
    {
      this->InteractionPerComponentIndex[TotalPos] = new int [this->NbrInteractionPerComponent[TotalPos]];
      this->InteractionPerComponentCoefficient[TotalPos] = new double [this->NbrInteractionPerComponent[TotalPos]];
      ++TotalPos;
    }

  QHEParticlePrecalculationOperation Operation(this, false);
  Operation.ApplyOperation(this->Architecture);
  
  this->FastMultiplicationFlag = true;
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// enable fast multiplication algorithm (partial evaluation)
//
// firstComponent = index of the first component that has to be precalcualted
// nbrComponent  = index of the last component that has to be precalcualted

void AbstractQHEOnSphereHamiltonian::PartialEnableFastMultiplication(int firstComponent, int nbrComponent)
{  
  int LastComponent = nbrComponent + firstComponent;
  ParticleOnSphere* TmpParticles = (ParticleOnSphere*) this->Particles->Clone();

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
      this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[Pos], 
							 this->InteractionPerComponentCoefficient[Pos], Pos);
    }

//   long TotalPos = ((firstComponent - this->PrecalculationShift - 1) / this->FastMultiplicationStep) + 1;
//   int InitialPos = ((firstComponent - 1) / this->FastMultiplicationStep) + 1;
//   InitialPos *= this->FastMultiplicationStep;
//   cout << "PartialEnableFastMultiplication: TotalPos="<<TotalPos<<", InitialPos="<<InitialPos<<endl;
//   for (int i = InitialPos; i < LastComponent; i += this->FastMultiplicationStep)
//     {
//       this->EvaluateMNTwoBodyFastMultiplicationComponent(TmpParticles, i, this->InteractionPerComponentIndex[TotalPos], 
// 							 this->InteractionPerComponentCoefficient[TotalPos], TotalPos);
//     }
  
  delete TmpParticles;
}

// enable fast multiplication algorithm using on disk cache 
//
// fileName = prefix of the name of the file where temporary matrix elements will be stored

void AbstractQHEOnSphereHamiltonian::EnableFastMultiplicationWithDiskStorage(char* fileName)
{
  if (this->FastMultiplicationStep == 1)
    {
      this->DiskStorageFlag = false;
      this->DiskStorageFileName = 0;
      this->EnableFastMultiplication();
      return;
    }
  this->DiskStorageFlag = true;
  this->DiskStorageFileName = new char [strlen(fileName) + 8];
  sprintf (this->DiskStorageFileName, "%s.ham", fileName);

  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  int EffectiveHilbertSpaceDimension = ((int) (MaxIndex - MinIndex)) + 1;
  this->DiskStorageStart = (int) MinIndex;
  int DiskStorageEnd = 1 + (int) MaxIndex;

  int Index;
  int m1;
  int m2;
  int m3;
  int m4;
  double Coefficient;
  int* TmpIndexArray;
  double* TmpCoefficientArray;
  int Pos;
  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  gettimeofday (&(TotalStartingTime2), 0);
  cout << "start" << endl;
  this->InteractionPerComponentIndex = 0;
  this->InteractionPerComponentCoefficient = 0;
  this->MaxNbrInteractionPerComponent = 0;

  int TotalPos = 0;
  ofstream File;
  File.open(this->DiskStorageFileName, ios::binary | ios::out);
 
  File.write((char*) &(EffectiveHilbertSpaceDimension), sizeof(int));
  File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
  File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * EffectiveHilbertSpaceDimension);

  long FileJump = 0;
  for (int i = 0; i < EffectiveHilbertSpaceDimension; ++i)
    {
      FileJump += (long) this->NbrInteractionPerComponent[i];
      if (this->MaxNbrInteractionPerComponent < this->NbrInteractionPerComponent[i])
	this->MaxNbrInteractionPerComponent = this->NbrInteractionPerComponent[i];
    }
  FileJump *= sizeof(int);

  TmpIndexArray = new int [this->MaxNbrInteractionPerComponent];
  TmpCoefficientArray = new double [this->MaxNbrInteractionPerComponent];      
  double Coefficient2;
  int SumIndices;
  int TmpNbrM3Values;
  int* TmpM3Values;
  int ReducedNbrInteractionFactors;

  for (int i = this->DiskStorageStart; i < DiskStorageEnd; ++i)
    {
      if (this->NbrInteractionPerComponent[TotalPos] > 0)
	{
	  Pos = 0;
	  if (this->NbrM12Indices == 0)
	    {
	      for (int j = 0; j < this->NbrInteractionFactors; ++j) 
		{
		  m1 = this->M1Value[j];
		  m2 = this->M2Value[j];
		  m3 = this->M3Value[j];
		  m4 = m1 + m2 - m3;
		  Index = this->Particles->AdAdAA(i, m1, m2, m3, m4, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[j];
		      ++Pos;
		    }
		}
	    }
	  else
	    {
	      ReducedNbrInteractionFactors = 0;
	      for (m1 = 0; m1 < this->NbrM12Indices; ++m1)
		{
		  Coefficient = this->Particles->AA(i, this->M1Value[m1], this->M2Value[m1]);	  
		  if (Coefficient != 0.0)
		    {
		      SumIndices = this->M1Value[m1] + this->M2Value[m1];
		      TmpM3Values = this->M3Values[m1];
		      TmpNbrM3Values = this->NbrM3Values[m1];
		      for (m3 = 0; m3 < TmpNbrM3Values; ++m3)
			{
			  Index = this->Particles->AdAd(TmpM3Values[m3], SumIndices - TmpM3Values[m3], Coefficient2);
			  if (Index < this->Particles->GetHilbertSpaceDimension())
			    {
			      TmpIndexArray[Pos] = Index;
			      TmpCoefficientArray[Pos] = Coefficient * this->InteractionFactors[ReducedNbrInteractionFactors] * Coefficient2;
			      ++Pos;
			    }
			  ++ReducedNbrInteractionFactors;
			}    
		    }
		  else
		    ReducedNbrInteractionFactors += this->NbrM3Values[m1];
		}	      
	    }
	  if (this->OneBodyTermFlag == true)
	    {
	      for (int j = 0; j < this->NbrOneBodyInteractionFactors; ++j)
		{
		  m1 = this->OneBodyMValues[j];
		  m2 = this->OneBodyNValues[j];
		  Index = this->Particles->AdA(i, m1, m2, Coefficient);
		  if (Index < this->Particles->GetHilbertSpaceDimension())
		    {
		      TmpIndexArray[Pos] = Index;
		      TmpCoefficientArray[Pos] = Coefficient * this->OneBodyInteractionFactors[j];
		      ++Pos;
		    }
		}
	    }
	  File.write((char*) TmpIndexArray, sizeof(int) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump -= sizeof(int) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(FileJump, ios::cur);
	  File.write((char*) TmpCoefficientArray, sizeof(double) * this->NbrInteractionPerComponent[TotalPos]);
	  FileJump += sizeof(double) * this->NbrInteractionPerComponent[TotalPos];
	  File.seekp(-FileJump, ios::cur);	  
	}
      ++TotalPos;
    }
  delete[] TmpIndexArray;
  delete[] TmpCoefficientArray;
  File.close();

  this->FastMultiplicationFlag = true;
  this->BufferSize = this->Memory / ((this->MaxNbrInteractionPerComponent * (sizeof(int) + sizeof(double))) + sizeof(int*) + sizeof(double*));

  gettimeofday (&(TotalEndingTime2), 0);
  cout << "------------------------------------------------------------------" << endl << endl;;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
}

// save precalculations in a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be stored
// return value = true if no error occurs

bool AbstractQHEOnSphereHamiltonian::SavePrecalculation (char* fileName)
{
  if (this->FastMultiplicationFlag)
    {
      ofstream File;
      File.open(fileName, ios::binary | ios::out);
      int Tmp = this->Particles->GetHilbertSpaceDimension();
      File.write((char*) &(Tmp), sizeof(int));
      File.write((char*) &(this->FastMultiplicationStep), sizeof(int));
      Tmp /= this->FastMultiplicationStep;
      if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
	++Tmp;
      File.write((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
	}
      for (int i = 0; i < Tmp; ++i)
	{
	  File.write((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
	}
      File.close();
      return true;
    }
  else
    {
      return false;
    }
}

// load precalculations from a file
// 
// fileName = pointer to a string containg the name of the file where precalculations have to be read
// return value = true if no error occurs

bool AbstractQHEOnSphereHamiltonian::LoadPrecalculation (char* fileName)
{
  ifstream File;
  File.open(fileName, ios::binary | ios::in);
  int Tmp;
  File.read((char*) &(Tmp), sizeof(int));
  if (Tmp != this->Particles->GetHilbertSpaceDimension())
    {
      File.close();
      return false;
    }
  File.read((char*) &(this->FastMultiplicationStep), sizeof(int));
  Tmp /= this->FastMultiplicationStep;
  if ((Tmp * this->FastMultiplicationStep) != this->Particles->GetHilbertSpaceDimension())
    ++Tmp;
  this->NbrInteractionPerComponent = new int [Tmp];
  File.read((char*) this->NbrInteractionPerComponent, sizeof(int) * Tmp);
  this->InteractionPerComponentIndex = new int* [Tmp];
  this->InteractionPerComponentCoefficient = new double* [Tmp];
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentIndex[i] = new int [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentIndex[i]), sizeof(int) * this->NbrInteractionPerComponent[i]);	  
    }
  for (int i = 0; i < Tmp; ++i)
    {
      this->InteractionPerComponentCoefficient[i] = new double [this->NbrInteractionPerComponent[i]];
      File.read((char*) (this->InteractionPerComponentCoefficient[i]), sizeof(double) * this->NbrInteractionPerComponent[i]);	  
    }
  File.close();
  this->FastMultiplicationFlag = true;
  return true;
}

// evaluate all interaction factors
//   

void AbstractQHEOnSphereHamiltonian::EvaluateInteractionFactors()
{
  this->InteractionFactors = 0;
  this->NbrInteractionFactors = 0;
  this->M1Value = 0;
  this->M2Value = 0;
  this->M3Value = 0;
  this->NbrM12Indices = 0;
  this->NbrM3Values = 0;
  this->M3Values = 0;
}

