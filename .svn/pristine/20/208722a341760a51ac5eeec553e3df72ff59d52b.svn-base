////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//             class density-density operator for particle with spin          //
//                                                                            //
//                        last modification : 10/12/2002                      //
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
#include "Operator/ParticleOnSphereWithSpinAllSzDensityDensityOddChannel.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationMomentumIndex1 = momentum index of the leftmost creation operator (from 0 to 2S)
// creationMomentumIndex2 = momentum index of the leftmost creation operator (from 0 to 2S)
// annihilationMomentumIndex1 = momentum index of the leftmost annihilation operator (from 0 to 2S)
// annihilationMomentumIndex2 = momentum index of the rightmost annihilation operator(from 0 to 2S)
ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::ParticleOnSphereWithSpinAllSzDensityDensityOddChannel(FermionOnSphereWithSpinAllSz* particle, int creationMomentumIndex1, int creationMomentumIndex2, int annihilationMomentumIndex1, int annihilationMomentumIndex2)
{
  this->Particle= particle;
  this->CreationMomentumIndex1 = creationMomentumIndex1;
  this->CreationMomentumIndex2 = creationMomentumIndex2;
  this->AnnihilationMomentumIndex1 = annihilationMomentumIndex1;
  this->AnnihilationMomentumIndex2 = annihilationMomentumIndex2;
}

// destructor
//

ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::~ParticleOnSphereWithSpinAllSzDensityDensityOddChannel()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (FermionOnSphereWithSpinAllSz*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int Index1, Index2;
  double Element = 0.0;
  int m = this->AnnihilationMomentumIndex1;
  int Q = this->AnnihilationMomentumIndex2;
  if ( (m != this->CreationMomentumIndex1) || (Q != this->CreationMomentumIndex2) )
   {
    cout << "error: difference in indices " << endl;
   }

  for (int i = 0; i < Dim; ++i)
   {

//****************** UP UP *************************
     //UpUpUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //UpUpUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += (-0.25) * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //UpUpDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //UpUpDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

//****************** UP DOWN *************************

     //UpDownUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += (-0.25) * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //UpDownUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //UpDownDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += (-0.25) * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //UpDownDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += (-0.25) * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

//****************** DOWN UP *************************

     //DownUpUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //DownUpUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	   Element += (-0.25) * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //DownUpDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //DownUpDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

//****************** DOWN DOWN *************************

     //DownDownUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //DownDownUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += (-0.25) * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //DownDownDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

     //DownDownDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	   Element += 0.25 * V1[Index] * V2[i] * Coefficient * Coefficient2;      
        }

   }
 
  return Complex(Element);
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  return Complex();
}
   
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSphereWithSpinAllSzDensityDensityOddChannel::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Coefficient2 = 0.0;
  int Sign;
  int m = this->AnnihilationMomentumIndex1;
  int Q = this->AnnihilationMomentumIndex2;
  if ( (m != this->CreationMomentumIndex1) || (Q != this->CreationMomentumIndex2) )
    {
      cout << "error: difference in indices " << endl;
    }
  
  for (int i = 0; i < Dim; ++i)
   {

//****************** UP UP *************************
     //UpUpUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //UpUpUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += -0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //UpUpDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //UpUpDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdu(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

//****************** UP DOWN *************************

     //UpDownUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += -0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //UpDownUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //UpDownDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += -0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //UpDownDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += -0.25 * vSource[i] * Coefficient * Coefficient2;
        }

//****************** DOWN UP *************************

     //DownUpUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //DownUpUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += -0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //DownUpDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //DownUpDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AduAdd(Q, m, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

//****************** DOWN DOWN *************************

     //DownDownUpUp
     Coefficient = this->Particle->AuAu(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //DownDownUpDown
     Coefficient = this->Particle->AuAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += -0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //DownDownDownUp
     Coefficient = this->Particle->AuAd(i, Q, m);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

     //DownDownDownDown
     Coefficient = this->Particle->AdAd(i, m, Q);
     if (Coefficient != 0.0)
       {
	 int Index = this->Particle->AddAdd(m, Q, Coefficient2);
	 if (Index != Dim)
	  vDestination[Index] += 0.25 * vSource[i] * Coefficient * Coefficient2;
        }

   }
 
  return vDestination;
}

