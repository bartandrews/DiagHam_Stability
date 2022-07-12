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
#include "Tools/FQHEWaveFunction/ReadRezayiOnSphereBlock.h"
#include "MathTools/BinomialCoefficients.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <math.h>


using std::cout;
using std::endl;


// constructor
//
// nbrParticlesPerCluster = number of particles per cluster (=N/2)
// nbrCluster = number of clusters
ReadRezayiOnSphereBlock::ReadRezayiOnSphereBlock(int nbrParticlesPerCluster, int nbrCluster)
{
  this->NbrParticles = nbrCluster*nbrParticlesPerCluster;
  this->ClusterSize = nbrParticlesPerCluster;
  this->NbrClusters = nbrCluster;
  this->SpinorUCoordinates = new Complex[this->NbrParticles];
  this->SpinorVCoordinates = new Complex[this->NbrParticles];
  this->Flag.Initialize();
}

// copy constructor
//
// function = reference on the wave function to copy

ReadRezayiOnSphereBlock::ReadRezayiOnSphereBlock(const ReadRezayiOnSphereBlock& function)
{
  this->NbrParticles = function.NbrParticles;
  this->ClusterSize = function.ClusterSize;
  this->NbrClusters = function.NbrClusters;
  this->SpinorUCoordinates = function.SpinorUCoordinates;
  this->SpinorVCoordinates = function.SpinorVCoordinates;
  this->Flag = function.Flag;
}

// destructor
//

ReadRezayiOnSphereBlock::~ReadRezayiOnSphereBlock()
{
  if ((this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete[] this->SpinorUCoordinates;
      delete[] this->SpinorVCoordinates;
    }
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* ReadRezayiOnSphereBlock::Clone ()
{
  return new ReadRezayiOnSphereBlock(*this);
}

// evaluate function at a given point
//
// x = point where the function has to be evaluated
// return value = function value at x  

Complex ReadRezayiOnSphereBlock::operator ()(RealVector& x)
{  
//   for (int i = 0; i < this->NbrParticles; ++i)
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

  Complex Value=1.0, Tmp;

  for (int c=0; c<NbrClusters; ++c)
    {
      for (int i=1+c*ClusterSize; i<(c+1)*ClusterSize; ++i)
	for (int j=c*ClusterSize; j<i; ++j)
	  {
	    Tmp=(this->SpinorUCoordinates[i] * this->SpinorVCoordinates[j]) - (this->SpinorUCoordinates[j] * this->SpinorVCoordinates[i]);
	    Value*=Tmp;
	  }
    }
  
  return Value*Value;  
}
