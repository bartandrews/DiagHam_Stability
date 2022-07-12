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
#include "Operator/ParticleOnSquareLatticeWithSpinBandDensityOperator.h"
#include "Output/MathematicaOutput.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

using std::cout;
using std::endl;

// constructor from default datas
//
// particle = hilbert space associated to the particles
// nbrOrbitalsPerBand = number of orbitals per band = nbr of unit cells
// bandIndex = selection of band to be considered: 0 / 1 
ParticleOnSquareLatticeWithSpinBandDensityOperator::ParticleOnSquareLatticeWithSpinBandDensityOperator(ParticleOnSphereWithSpin* particle, int nbrOrbitalsPerBand, int bandIndex)
{
  this->Particle = particle;
  this->NbrOrbitalsPerBand = nbrOrbitalsPerBand;
  this->BandIndex = bandIndex&1;
}

// destructor
//

ParticleOnSquareLatticeWithSpinBandDensityOperator::~ParticleOnSquareLatticeWithSpinBandDensityOperator()
{
}
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnSquareLatticeWithSpinBandDensityOperator::Clone ()
{
  return 0;
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnSquareLatticeWithSpinBandDensityOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnSphereWithSpin*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnSquareLatticeWithSpinBandDensityOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnSquareLatticeWithSpinBandDensityOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}
  
// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnSquareLatticeWithSpinBandDensityOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  if (this->BandIndex==0)
    {
      for (int i = (int) firstComponent; i < Dim; ++i)
	{
	  for (int OrbIndex=0; OrbIndex<NbrOrbitalsPerBand; ++OrbIndex)
	    {
	      int Index = this->Particle->AddAd(i, OrbIndex, OrbIndex, Coefficient);
	      if (Index != FullDim)
		Element += V1[Index] * V2[i] * Coefficient;
	    }
	}
    }
  else
    {
      for (int i = (int) firstComponent; i < Dim; ++i)
	{
	  for (int OrbIndex=0; OrbIndex<NbrOrbitalsPerBand; ++OrbIndex)
	    {
	      int Index = this->Particle->AduAu(i, OrbIndex, OrbIndex, Coefficient);
	      if (Index != FullDim)
		Element += V1[Index] * V2[i] * Coefficient;
	    }
	}
    }
  return Complex(Element);
}
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

RealVector& ParticleOnSquareLatticeWithSpinBandDensityOperator::LowLevelAddMultiply(RealVector& vSource, RealVector& vDestination, 
									 int firstComponent, int nbrComponent)
{
  int Last = firstComponent + nbrComponent;;
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  if (this->BandIndex==0)
    {
      for (int i = firstComponent; i < Last; ++i)
	  {
	    for (int OrbIndex=0; OrbIndex<NbrOrbitalsPerBand; ++OrbIndex)
	      {
		int Index = this->Particle->AddAd(i, OrbIndex, OrbIndex, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += vSource[i] * Coefficient;      
	      }
	  }
    }
  else
    {
      for (int i = firstComponent; i < Last; ++i)
	  {
	    for (int OrbIndex=0; OrbIndex<NbrOrbitalsPerBand; ++OrbIndex)
	      {
		int Index = this->Particle->AduAu(i, OrbIndex, OrbIndex, Coefficient);
		if (Index != Dim)
		  vDestination[Index] += vSource[i] * Coefficient;      
	      }
	  }
    }
  return vDestination;
}

