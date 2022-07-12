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
#include "Operator/ParticleOnLatticeVacancyOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"


#include <iostream>
using std::cout;
using std::endl;
  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// creationIndex = index of the creation operator
// annihilationIndex = index of the annihilation operator
//
ParticleOnLatticeVacancyOperator::ParticleOnLatticeVacancyOperator(ParticleOnLattice* source, ParticleOnLattice* target, LatticePhases *lattice,
								   int kx, int ky, int sublattice, bool particle)
{
  this->SourceSpace = source;
  this->TargetSpace = target;
  this->Lattice = lattice;
  this->NbrCells = Lattice->GetNbrCells();
  this->QuantumNbrs = new int[NbrCells];
  this->Phases = new Complex[NbrCells];
  this->SetMomenta (kx, ky, sublattice);
  this->ParticleFlag = particle;
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeVacancyOperator::ParticleOnLatticeVacancyOperator(const ParticleOnLatticeVacancyOperator& oper)
{
  this->SourceSpace = oper.SourceSpace;
  this->TargetSpace = oper.TargetSpace;
  this->Kx = oper.Kx;
  this->Ky = oper.Ky;
  this->SubLattice = oper.SubLattice;
  this->Lattice = oper.Lattice;
  this->NbrCells = Lattice->GetNbrCells();
  this->QuantumNbrs = new int[NbrCells];
  this->Phases = new Complex[NbrCells];
  this->SetMomenta (Kx, Ky, SubLattice);
  this->ParticleFlag = oper.ParticleFlag;
}

// destructor
//

ParticleOnLatticeVacancyOperator::~ParticleOnLatticeVacancyOperator()
{
  if (this->NbrCells>0)
    {
      delete [] this->Phases;
      delete [] this->QuantumNbrs;
    }
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeVacancyOperator::Clone ()
{
  return new ParticleOnLatticeVacancyOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeVacancyOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->SourceSpace = (ParticleOnLattice*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeVacancyOperator::GetHilbertSpace ()
{
  return this->SourceSpace;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeVacancyOperator::GetHilbertSpaceDimension ()
{
  return this->SourceSpace->GetHilbertSpaceDimension();
}


// change indices of creation / annihilation operators
// creationIndex = index of the creation operator
// annihilationIndex = index of the annihilation operator
void ParticleOnLatticeVacancyOperator::SetMomenta (int kx, int ky, int sublattice)
{
  this->Kx=kx;
  this->Ky=ky;
  this->SubLattice = sublattice;
  int Lx=Lattice->GetLatticeLength(0);
  int Ly=Lattice->GetLatticeLength(1);
  Complex TmpC;
  int Pos=0;
  //cout << "Setting momentum for kx="<<kx<<", ky="<<ky<<endl;
  //cout.flush();
  for (int x=0; x<Lx; ++x)
    for (int y=0; y<Ly; ++y)
      {
	this->QuantumNbrs[Pos] = SourceSpace->EncodeQuantumNumber(x, y, this->SubLattice, TmpC);
	this->Phases[Pos] = Polar(1.0, (((2.0*M_PI)/Lx) * kx * x) + (((2.0*M_PI)/Ly) * ky * y));
	++Pos;
      }
//   for (int i=0; i<NbrCells; ++i)
//     cout << "Phases["<<i<<"]="<<this->Phases[i]<<endl;
}


// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element

Complex ParticleOnLatticeVacancyOperator::PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->TargetSpace->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->SourceSpace->Clone();
  double Coefficient = 0.0;
  Complex Element = 0.0;
  int Index;
  if (ParticleFlag)
    for (int i = (int) firstComponent; i < Dim; ++i)
      {
	for (int q=0; q<NbrCells; ++q)
	  {
	    Index = TmpParticle->Ad(i, QuantumNbrs[q], Coefficient);
	    if ((Index<FullDim)&&(Coefficient!=0.0))
	      Element += V1[Index] * V2[i] * Coefficient * Phases[q];
	  }
      }
  else
    for (int i = (int) firstComponent; i < Dim; ++i)
      {
	for (int q=0; q<NbrCells; ++q)
	  {
	    Index = TmpParticle->A(i, QuantumNbrs[q], Coefficient);
	    if ((Index<FullDim)&&(Coefficient!=0.0))
	      Element += V1[Index] * V2[i] * Coefficient * Phases[q];
	  }
      }
  delete TmpParticle;
  return Complex(Element);
}

// evaluate part of the matrix element, within a given of indices
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = corresponding matrix element
  
Complex ParticleOnLatticeVacancyOperator::PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent)
{
  int Dim = (int) (firstComponent + nbrComponent);
  int FullDim = this->TargetSpace->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->SourceSpace->Clone();
  double Coefficient = 0.0;
  Complex Element;
  int Index;
  if (ParticleFlag)
    for (int i = (int) firstComponent; i < Dim; ++i)
      {
	for (int q=0; q<NbrCells; ++q)
	  {
	    Index = TmpParticle->Ad(i, QuantumNbrs[q], Coefficient);
	    if ((Index<FullDim)&&(Coefficient!=0.0))
	      Element += (Conj(V1[Index]) * V2[i] * Coefficient * Phases[q]);
	  }
      }
  else
    for (int i = (int) firstComponent; i < Dim; ++i)
      {
	for (int q=0; q<NbrCells; ++q)
	  {
	    Index = TmpParticle->A(i, QuantumNbrs[q], Coefficient);
	    if ((Index<FullDim)&&(Coefficient!=0.0))
	      Element += (Conj(V1[Index]) * V2[i] * Coefficient * Phases[q]);
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

ComplexVector& ParticleOnLatticeVacancyOperator::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int FullDim = this->TargetSpace->GetHilbertSpaceDimension();
  ParticleOnLattice* TmpParticle = (ParticleOnLattice*) this->SourceSpace->Clone();
  int Last = firstComponent + nbrComponent;;
  int Index;
  double Coefficient = 0.0;

  if (ParticleFlag)
    for (int i = (int) firstComponent; i < Last; ++i)
      {
	for (int q=0; q<NbrCells; ++q)
	  {
	    Index = TmpParticle->Ad(i, QuantumNbrs[q], Coefficient);
	    if ((Index<FullDim)&&(Coefficient!=0.0))
	      vDestination[Index] += vSource[i] * Coefficient * Phases[q];
	  }
      }
  else
    for (int i = (int) firstComponent; i < Last; ++i)
      {
	for (int q=0; q<NbrCells; ++q)
	  {
	    Index = TmpParticle->A(i, QuantumNbrs[q], Coefficient);
	    if ((Index<FullDim)&&(Coefficient!=0.0))
	      vDestination[Index] += vSource[i] * Coefficient * Phases[q];
	  }
      }
  delete TmpParticle;
  return vDestination;
}


