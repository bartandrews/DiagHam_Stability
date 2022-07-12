////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of Moore Read state wave function on sphere             //
//                                                                            //
//                        last modification : 19/09/2004                      //
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
#include "Tools/FQHEWaveFunction/NASSOnSphereWaveFunction.h"
#include "MathTools/BinomialCoefficients.h"
#include "MathTools/FactorialCoefficient.h"
#include "MathTools/SimplePermutations.h"
#include "Vector/RealVector.h"

#include "GeneralTools/OrderedList.h"

#include <iostream>
#include <cstdlib>
#include <cmath>


using std::cout;
using std::endl;


// constructor
//
// nbrParticlesPerCluster = number of particles per cluster (=N/k)
// fermionicStatistics = flag indicating whether the pfaffian should be multiplied by a squared Jastrow Factor
NASSOnSphereWaveFunction::NASSOnSphereWaveFunction(int nbrParticlesPerCluster, bool fermionicStatistics)
{
  this->NbrParticles = 4*nbrParticlesPerCluster;
  this->FermionicStatistics = fermionicStatistics;
  this->SpinorUCoordinates = new Complex[this->NbrParticles];
  this->SpinorVCoordinates = new Complex[this->NbrParticles];
  this->JastrowFactorElements = new Complex*[this->NbrParticles];
  this->JastrowFactorSquares = new Complex*[this->NbrParticles];  
  for (int i=0; i<NbrParticles; ++i)
    {
      this->JastrowFactorElements[i] = new Complex[this->NbrParticles];
      this->JastrowFactorSquares[i] = new Complex[this->NbrParticles];
    }
  this->EvaluatePermutations();
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

NASSOnSphereWaveFunction::NASSOnSphereWaveFunction(const NASSOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->FermionicStatistics = function.FermionicStatistics;
  this->SpinorUCoordinates = function.SpinorUCoordinates;
  this->SpinorVCoordinates = function.SpinorVCoordinates;
  this->Permutations = function.Permutations;
  this->NbrPermutations = function.NbrPermutations;
  this->JastrowFactorElements = function.JastrowFactorElements;
  this->JastrowFactorSquares = function.JastrowFactorSquares;
  this->Flag = function.Flag;
}

// destructor
//

NASSOnSphereWaveFunction::~NASSOnSphereWaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      for (unsigned i = 0; i < this->NbrPermutations; ++i)
	delete[] this->Permutations[i];
      delete[] this->Permutations;
      delete[] this->SpinorUCoordinates;
      delete[] this->SpinorVCoordinates;
      for (int i=0; i<NbrParticles; ++i)
	{
	  delete [] this->JastrowFactorSquares[i];
	  delete [] this->JastrowFactorElements[i];
	}
      delete [] this->JastrowFactorSquares;
      delete [] this->JastrowFactorElements;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* NASSOnSphereWaveFunction::Clone ()
{
  return new NASSOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex NASSOnSphereWaveFunction::operator ()(RealVector& x)
{  
  // CalculateSpinors
  double s,c;
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
      this->SpinorUCoordinates[i].Im = this->SpinorUCoordinates[i].Re;
      this->SpinorUCoordinates[i].Re *= (c=cos(0.5 * x[1 + (i << 1)]));
      this->SpinorUCoordinates[i].Im *= -(s=sin(0.5 * x[1 + (i << 1)]));
      this->SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
      this->SpinorVCoordinates[i].Im = this->SpinorVCoordinates[i].Re;
      this->SpinorVCoordinates[i].Re *= c;
      this->SpinorVCoordinates[i].Im *= s;
      //cout << "U["<<i<<"]="<<SpinorUCoordinates[i]<<", "<< "V["<<i<<"]="<<SpinorVCoordinates[i]<<endl;
    }

  return this->ComplexEvaluations();
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex NASSOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{  
  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
    }

  return this->ComplexEvaluations();
}  


// evaluate permutations required for the Moore-Read state evaluation
// using: two symmetric blocs, only permutations changing particles
//        between blocks are required

void NASSOnSphereWaveFunction::EvaluatePermutations()
{  
  SimplePermutations Generator(this->NbrParticles>>1);
  this->NbrPermutations = Generator.GetNbrPermutations();
  this->Permutations = Generator.CheckOutPermutations();
//   for (unsigned i=0; i<NbrPermutations; ++i)
//     {
//       cout << "Permutation["<<i<<"]= ["<<Permutations[i][0];
//       for (int j=1; j<NbrParticles>>1; ++j) cout << " " << Permutations[i][j];
//       cout << "]"<<endl;
//     }
  return;
}

// perform complex part of calculations
// uses internal spinor coordinates as input
//
Complex NASSOnSphereWaveFunction::ComplexEvaluations()
{
  Complex J=1.0, Jastrow=1.0, Tmp;
  if (FermionicStatistics)    
    for (int i = 0; i < this->NbrParticles; ++i)
      for (int j = 0; j < i; ++j)
	{
	  Tmp = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	  JastrowFactorElements[i][j] = Tmp;
	  JastrowFactorElements[j][i] = -Tmp;
	  J*=Tmp;
	}
  else
    for (int i = 0; i < this->NbrParticles; ++i)
      for (int j = 0; j < i; ++j)
	{
	  Tmp = ((this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]));
	  JastrowFactorElements[i][j] = Tmp;
	  JastrowFactorElements[j][i] = -Tmp;
	}
  Jastrow=J;

  unsigned *PUp;
  unsigned *PDown;
  int NFourth=NbrParticles>>2;
  int NHalf = NbrParticles>>1;
  Complex J2,Value(0.0,0.0);
  // cout << "Evaluating function"<<endl;
  for (unsigned i=0; i<NbrPermutations; ++i)
    {
      PUp=this->Permutations[i];
      J=1.0;
      for (int k=1; k<NFourth; ++k)
	for (int j=0; j<k; ++j)
	  {
	    J*=JastrowFactorElements[PUp[k]][PUp[j]];
	    J*=JastrowFactorElements[PUp[k+NFourth]][PUp[j+NFourth]];
	  }
      J*=J;

      for (unsigned i=0; i<NbrPermutations; ++i)
	{
	  PDown=this->Permutations[i];
	  J2=1.0;
	  for (int k=1; k<NFourth; ++k)
	    for (int j=0; j<k; ++j)
	      {
		J2*=JastrowFactorElements[NHalf+PDown[k]][NHalf+PDown[j]];
		J2*=JastrowFactorElements[NHalf+PDown[k+NFourth]][NHalf+PDown[j+NFourth]];
	      }
	  J2*=J2;
	  for (int k=0; k<NFourth; ++k)
	    for (int j=0; j<NFourth; ++j)
	      {
		J2*=JastrowFactorElements[PUp[k]][NHalf+PDown[j]];
		J2*=JastrowFactorElements[PUp[k+NFourth]][NHalf+PDown[j+NFourth]];
	      }
	  Value += J*J2;
	}
    }
  if (this->FermionicStatistics)
    Value *= Jastrow;
  //  cout << "Value ="<<Value<<endl;
  
  return Value;
}
