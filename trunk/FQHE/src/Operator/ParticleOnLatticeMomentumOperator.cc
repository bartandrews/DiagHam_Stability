////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of particle on lattice 1-body operator              //
//                                                                            //
//                        last modification : 09/04/2008                      //
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
#include "Operator/ParticleOnLatticeMomentumOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"
#include "GeneralTools/List.h"
#include "GeneralTools/ListIterator.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/PeriodicMomentumQuantumNumber.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
using std::cout;
using std::endl;
using std::sin;
using std::cos;

  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// lx = length in x-direction
// ly = length in y-direction
// subl = number of sublattices
// momentumX = X-component of the momentum
// momentumY = Y-component of the momentum
ParticleOnLatticeMomentumOperator::ParticleOnLatticeMomentumOperator(ParticleOnLattice* particle, int lx, int ly, int subl,
								     int momentumX, int momentumY, double offsetX, double offsetY)
{
  this->Particle = particle;
  this->Lx = lx;
  this->Ly = ly;
  this->NbrSubLattices = subl;
  if (this->Particle->GetHilbertSpaceAdditionalSymmetry()==ParticleOnLattice::NoSymmetry)
    {
      this->NbrTerms = Lx*Ly;
    }
  else if (this->Particle->GetHilbertSpaceAdditionalSymmetry()==ParticleOnLattice::XTranslations)
    {
      this->NbrTerms = Ly;
      List<AbstractQuantumNumber*> Numbers=this->Particle->GetQuantumNumbers();
      ListIterator<AbstractQuantumNumber*> Iterator(Numbers);
      AbstractQuantumNumber** Number;
      for (Number = Iterator();
	   (Number!=NULL) && ((*Number)->GetQuantumNumberType()!=AbstractQuantumNumber::PeriodicMomentum); Number = Iterator())
	{}
      if (Number == NULL)
	{
	  cout << "Error: Space has no periodic momentum"<<endl;
	  exit(1);
	}
      int TmpKx = ((PeriodicMomentumQuantumNumber*)Number)->GetMomentum();
      int TmpP = ((PeriodicMomentumQuantumNumber*)Number)->GetPeriod();
      if (momentumX%TmpP != TmpKx)
	{
	  cout << "Error: Momentum of space is not commensurate with requested momentum in ParticleOnLatticeMomentumOperator"<<endl;
	  exit(1);
	}
    }
  else if (this->Particle->GetHilbertSpaceAdditionalSymmetry()==ParticleOnLattice::YTranslations)
    {
      this->NbrTerms = Lx;
      List<AbstractQuantumNumber*> Numbers=this->Particle->GetQuantumNumbers ();
      ListIterator<AbstractQuantumNumber*> Iterator(Numbers);
      AbstractQuantumNumber** Number;
      for (Number = Iterator();
	   (Number!=NULL) && ((*Number)->GetQuantumNumberType()!=AbstractQuantumNumber::PeriodicMomentum); Number = Iterator())
	{}
      if (Number == NULL)
	{
	  cout << "Error: Space has no periodic momentum"<<endl;
	  exit(1);
	}
      int TmpKy = ((PeriodicMomentumQuantumNumber*)Number)->GetMomentum();
      int TmpP = ((PeriodicMomentumQuantumNumber*)Number)->GetPeriod();
      if (momentumY%TmpP != TmpKy)
	{
	  cout << "Error: Momentum of space is not commensurate with requested momentum in ParticleOnLatticeMomentumOperator"<<endl;
	  exit(1);
	}
    }
  this->NbrTerms *= this->NbrSubLattices;
  this->NbrTerms *= this->NbrTerms;
  this->CreationIndices = new int[NbrTerms];
  this->AnnihilationIndices = new int[NbrTerms];
  this->PhaseTableRe = new double[NbrTerms];
  this->PhaseTableIm = new double[NbrTerms];
  this->Flag.Initialize();
  this->SetMomentum(momentumX, momentumY);
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeMomentumOperator::ParticleOnLatticeMomentumOperator(const ParticleOnLatticeMomentumOperator& oper)
{
  this->Particle = oper.Particle;
  this->MomentumX = oper.MomentumX;
  this->MomentumY = oper.MomentumY;
  this->Lx = oper.Lx;
  this->Ly = oper.Ly;
  this->NbrSubLattices = oper.NbrSubLattices;
  this->NbrTerms = oper.NbrTerms;
  this->CreationIndices = oper.CreationIndices;
  this->AnnihilationIndices = oper.AnnihilationIndices;
  this->PhaseTableRe = oper.PhaseTableRe;
  this->PhaseTableIm = oper.PhaseTableIm;
  this->Flag = oper.Flag;
}

// destructor
//

ParticleOnLatticeMomentumOperator::~ParticleOnLatticeMomentumOperator()
{
  if (((this->NbrTerms) != 0) && (this->Flag.Used() == true) && (this->Flag.Shared() == false))
    {
      delete [] this->PhaseTableRe;
      delete [] this->PhaseTableIm;
      delete [] this->CreationIndices;
      delete [] this->AnnihilationIndices;
    }
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeMomentumOperator::Clone ()
{
  return new ParticleOnLatticeMomentumOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeMomentumOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnLattice*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeMomentumOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeMomentumOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}


// change values of momentum represented
// momentumX = X-component of the momentum
// momentumY = Y-component of the momentum
// offsetX = absolute offset of momentum values along x-axis
// offsetY = absolute offset of momentum values along y-axis
void ParticleOnLatticeMomentumOperator::SetMomentum (int momentumX, int momentumY, double offsetX, double offsetY)
{
  if (this->Particle->GetHilbertSpaceAdditionalSymmetry()==ParticleOnLattice::NoSymmetry)
    {
      while (momentumX<0) momentumX+=Lx;
      while (momentumY<0) momentumY+=Ly;
      this->MomentumX = momentumX%Lx;
      this->MomentumY = momentumY%Ly;
      double Kx=this->MomentumX*2.0*M_PI/Lx+offsetX;
      double Ky=this->MomentumY*2.0*M_PI/Ly+offsetY;
      int CreationIndex, AnnihilationIndex;
      Complex Tmp;
      int Pos=0;
      for (int xi=0; xi<Lx; ++xi)
	for (int yi=0; yi<Ly; ++yi)
	  for (int subi=0; subi<NbrSubLattices; ++subi)
	    {
	      AnnihilationIndex = Particle->EncodeQuantumNumber(xi, yi, subi, Tmp);
	      for (int xf=0; xf<Lx; ++xf)
		for (int yf=0; yf<Ly; ++yf)
		  for (int subf=0; subf<NbrSubLattices; ++subf)
		    {
		      CreationIndex = Particle->EncodeQuantumNumber(xf, yf, subf, Tmp);
		      CreationIndices[Pos] = CreationIndex;
		      AnnihilationIndices[Pos] = AnnihilationIndex;
		      PhaseTableRe[Pos]=cos(Kx*(xi-xf)+Ky*(yi-yf));
		      PhaseTableIm[Pos]=sin(Kx*(xi-xf)+Ky*(yi-yf));
		      ++Pos;
		    }
	    }
    }
  else if (this->Particle->GetHilbertSpaceAdditionalSymmetry()==ParticleOnLattice::XTranslations)
    {
      cout << "XTranslations need to be implemented in ParticleOnLatticeMomentumOperator::SetMomentum"<<endl;
      exit(1);
    }
  else if (this->Particle->GetHilbertSpaceAdditionalSymmetry()==ParticleOnLattice::YTranslations)
    {
      cout << "YTranslations need to be implemented in ParticleOnLatticeMomentumOperator::SetMomentum"<<endl;
      exit(1);
    }
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnLatticeMomentumOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      for (int j=0; j<NbrTerms; ++j)
	{
	  Index = TmpParticle->AdA(i, this->CreationIndices[j], this->AnnihilationIndices[j], Coefficient);
	  if ((Index<FullDim)&&(Coefficient!=0.0))
	    {
	      Coefficient *= V1[Index] * V2[i];
	      Element.Re += PhaseTableRe[j] * Coefficient;
	      Element.Im += PhaseTableIm[j] * Coefficient;
	    }
	}
    }
  delete TmpParticle;
  return Element;
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element
  
Complex ParticleOnLatticeMomentumOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->Particle->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  double Coefficient = 0.0;
  Complex Element;
  int Index;
  for (int i = (int) firstComponent; i < Dim; ++i)
    {
      for (int j=0; j<NbrTerms; ++j)
	{
	  Index = TmpParticle->AdA(i, this->CreationIndices[j], this->AnnihilationIndices[j], Coefficient);
	  if ((Index<FullDim)&&(Coefficient!=0.0))
	    Element += (Conj(V1[Index]) * Complex(PhaseTableRe[j], PhaseTableIm[j]) * V2[i] * Coefficient);
	}
    }
  delete TmpParticle;
  return Element; 
}   
  
// multiply a vector by the current operator for a given range of indices 
// and store result in another vector
//
// vSource = vector to be multiplied
// vDestination = vector where result has to be stored
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& ParticleOnLatticeMomentumOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->Particle->Clone();
  int Last = firstComponent + nbrComponent;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      for (int j=0; j<NbrTerms; ++j)
	{
	  Index = TmpParticle->AdA(i, this->CreationIndices[j], this->AnnihilationIndices[j], Coefficient);
	  if ((Index<Dim)&&(Coefficient!=0.0))
	    {
	      vDestination[Index].Re += (vSource[i].Re * PhaseTableRe[j] - vSource[i].Im * PhaseTableIm[j]) * Coefficient;
	      vDestination[Index].Im += (vSource[i].Im * PhaseTableRe[j] + vSource[i].Re * PhaseTableIm[j]) * Coefficient;
	    }
	}
    }
  delete TmpParticle;
  return vDestination;
}

