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
#include "Tools/FQHEWaveFunction/ExplicitMooreReadOnSphereWaveFunction.h"
#include "MathTools/BinomialCoefficients.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;

// constructor
//
// nbrParticlesPerCluster = number of particles per cluster
// nbrClusters = number of clusters
// fermionicStatistics = flag indicating whether the pfaffian should be multiplied by a squared Jastrow Factor
ExplicitMooreReadOnSphereWaveFunction::ExplicitMooreReadOnSphereWaveFunction(int nbrParticlesPerCluster, int nbrClusters, bool fermionicStatistics)
{
  this->NbrParticles = nbrClusters*nbrParticlesPerCluster;  
  this->ClusterSize = nbrParticlesPerCluster;
  this->NbrClusters = nbrClusters;
  this->FermionicStatistics = fermionicStatistics;
  this->Block = new ReadRezayiOnSphereBlock(nbrParticlesPerCluster, nbrClusters);
  this->Symmetrizer = new SymmetrizedComplexFunction(Block, NbrParticles, 2);
  this->SpinorUCoordinates = new Complex[this->NbrParticles];
  this->SpinorVCoordinates = new Complex[this->NbrParticles];
  this->RealCoordinates.Resize(2*NbrParticles);
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

ExplicitMooreReadOnSphereWaveFunction::ExplicitMooreReadOnSphereWaveFunction(const ExplicitMooreReadOnSphereWaveFunction& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ClusterSize = function.ClusterSize;
  this->NbrClusters = function.NbrClusters;
  this->FermionicStatistics = function.FermionicStatistics;
  this->Block = function.Block;
  this->Symmetrizer = function.Symmetrizer;
  this->SpinorUCoordinates = function.SpinorUCoordinates;
  this->SpinorVCoordinates = function.SpinorVCoordinates;
  this->RealCoordinates.Resize(2*NbrParticles);
  this->Flag = function.Flag;
}

// destructor
//

ExplicitMooreReadOnSphereWaveFunction::~ExplicitMooreReadOnSphereWaveFunction()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete[] this->SpinorUCoordinates;
      delete[] this->SpinorVCoordinates;
      delete Block;
      delete Symmetrizer;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* ExplicitMooreReadOnSphereWaveFunction::Clone ()
{
  return new ExplicitMooreReadOnSphereWaveFunction(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex ExplicitMooreReadOnSphereWaveFunction::operator ()(RealVector& x)
{  
  // for (int i = 0; i < this->NbrParticles; ++i)
  //     {
  //       SpinorUCoordinates[i].Re = cos(0.5 * x[i << 1]);
  //       SpinorUCoordinates[i].Im = SpinorUCoordinates[i].Re;
  //       SpinorUCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
  //       SpinorUCoordinates[i].Im *= sin(0.5 * x[1 + (i << 1)]);
  //       SpinorVCoordinates[i].Re = sin(0.5 * x[i << 1]);
  //       SpinorVCoordinates[i].Im = SpinorVCoordinates[i].Re;
  //       SpinorVCoordinates[i].Re *= cos(0.5 * x[1 + (i << 1)]);
  //       SpinorVCoordinates[i].Im *= -sin(0.5 * x[1 + (i << 1)]);
  //     }  
  
  Complex Value = (*Symmetrizer)(x);
  cout << "Explicit: Symmetric part:" << Value << endl;
  if (this->FermionicStatistics)
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
      for (int i = 0; i < this->NbrParticles; ++i)
	for (int j = 0; j < i; ++j)
	  Value *=  (this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]);
    }

  return Value;
}

// evaluate function at a given point
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)
Complex ExplicitMooreReadOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  // Import from spinors
  for (int i = 0; i < this->NbrParticles; ++i)
    {
      this->SpinorUCoordinates[i].Re = uv.Re(2*i);
      this->SpinorUCoordinates[i].Im = uv.Im(2*i);
      this->SpinorVCoordinates[i].Re = uv.Re(2*i+1);
      this->SpinorVCoordinates[i].Im = uv.Im(2*i+1);
      this->RealCoordinates[i<<1] = (2.0*acos(Norm(SpinorUCoordinates[i])));
      this->RealCoordinates[(i<<1)+1] = (Arg(SpinorVCoordinates[i])-Arg(SpinorUCoordinates[i]));
    }

  Complex Value = (*Symmetrizer)(RealCoordinates);

  if (this->FermionicStatistics)
    {      
      for (int i = 0; i < this->NbrParticles; ++i)
	for (int j = 0; j < i; ++j)
	  Value *=  (this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]);
    }
  
  return Value;
}  
