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
#include "Operator/ParticleOnLatticeFluxOperator.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "MathTools/Complex.h"

#include <iostream>
using std::cout;
using std::endl;
  
// constructor from default datas
//
// particle = hilbert space associated to the particles
// cellNumber = cell number where current operator is centered
ParticleOnLatticeFluxOperator::ParticleOnLatticeFluxOperator(ParticleOnLattice* particle, int cellNumber)
{
  this->Particle = (ParticleOnLattice*) (particle->Clone());
  this->CreationIndices = new int[8];
  this->AnnihilationIndices = new int[8];
  this->Coefficients = new Complex[8];
  this->CellNumber = cellNumber;
  int posx, posy, sublattice;
  this->Particle->DecodeQuantumNumber(cellNumber, posx, posy, sublattice);
  this->SetCellPosition(posx, posy, sublattice);  
}

// constructor from default datas
//
// particle = hilbert space associated to the particles
// posx, posy = coordinates of cell where current operator is centered
ParticleOnLatticeFluxOperator::ParticleOnLatticeFluxOperator(ParticleOnLattice* particle, int posx, int posy, int sublattice)
{
  this->Particle = (ParticleOnLattice*) (particle->Clone());
  this->CreationIndices = new int[8];
  this->AnnihilationIndices = new int[8];
  this->Coefficients = new Complex[8];
  this->SetCellPosition(posx, posy, sublattice);  
}

// copy constructor
//
// oper = reference on the operator to copy
 
ParticleOnLatticeFluxOperator::ParticleOnLatticeFluxOperator(const ParticleOnLatticeFluxOperator& oper)
{
  this->Particle = (ParticleOnLattice*) (oper.Particle->Clone());
  this->CreationIndices = oper.CreationIndices;
  this->AnnihilationIndices = oper.AnnihilationIndices;
  this->Flag = oper.Flag;
}

// destructor
//

ParticleOnLatticeFluxOperator::~ParticleOnLatticeFluxOperator()
{
  if ((this->Flag.Shared() == false) && (this->Flag.Used() == true))
    {
      delete [] CreationIndices;
      delete [] AnnihilationIndices;
    }
  delete this->Particle;
}
  
  
// clone operator without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractOperator* ParticleOnLatticeFluxOperator::Clone ()
{
  return new ParticleOnLatticeFluxOperator(*this);
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void ParticleOnLatticeFluxOperator::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  this->Particle = (ParticleOnLattice*) hilbertSpace;
}

// get Hilbert space on which operator acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* ParticleOnLatticeFluxOperator::GetHilbertSpace ()
{
  return this->Particle;
}

// return dimension of Hilbert space where operator acts
//
// return value = corresponding matrix elementdimension

int ParticleOnLatticeFluxOperator::GetHilbertSpaceDimension ()
{
  return this->Particle->GetHilbertSpaceDimension();
}


// change position of operator
// cellIndex = cell number where current operator is centered
void ParticleOnLatticeFluxOperator::SetCellIndex (int cellIndex)
{
  this->CellNumber = cellIndex;
  int posx, posy, sublattice;
  this->Particle->DecodeQuantumNumber(CellNumber, posx, posy, sublattice);
  this->SetCellPosition(posx, posy, sublattice);  
}

// change position of operator
// posx, posy = coordinates of cell where current operator is centered
void ParticleOnLatticeFluxOperator::SetCellPosition (int posx, int posy, int sublattice)
{
  Complex Phase;
  this->CellNumber=this->Particle->EncodeQuantumNumber(posx, posy, sublattice, Phase);
  if (this->Particle->GetLandauGaugeAxis()!='y')
    {
      cout << "Flux operators only defined for gauge along y-axis, sorry!"<<endl;
      exit(1);
    }
  Complex Phase1=Polar(1.0,2.0*M_PI*Particle->GetNbrFluxQuanta()*(double)posx);
  Complex Phase2=Polar(1.0,2.0*M_PI*Particle->GetNbrFluxQuanta()*(double)(posx+1));
  Complex TranslationPhase, TranslationPhase2;

  CreationIndices[0]=CellNumber;
  AnnihilationIndices[0]=Particle->EncodeQuantumNumber(posx+1, posy, 0, TranslationPhase);
  Coefficients[0]=-TranslationPhase;
  CreationIndices[1]=AnnihilationIndices[0];
  AnnihilationIndices[1]=CreationIndices[0];
  Coefficients[1]=Conj(TranslationPhase);

  CreationIndices[2]=AnnihilationIndices[0];
  AnnihilationIndices[2]=Particle->EncodeQuantumNumber(posx+1, posy+1, 0, TranslationPhase);
  Coefficients[2]=-Phase2*Conj(TranslationPhase);
  CreationIndices[3]=AnnihilationIndices[2];
  AnnihilationIndices[3]=CreationIndices[2];
  Coefficients[3]=-Conj(Coefficients[2]);

  CreationIndices[4]=AnnihilationIndices[2];
  AnnihilationIndices[4]=Particle->EncodeQuantumNumber(posx, posy+1, 0, TranslationPhase2);
  Coefficients[4]=-TranslationPhase*Conj(TranslationPhase2);
  CreationIndices[5]=AnnihilationIndices[4];
  AnnihilationIndices[5]=CreationIndices[4];
  Coefficients[5]=-Conj(Coefficients[4]);

  CreationIndices[6]=AnnihilationIndices[4];
  AnnihilationIndices[6]=CellNumber;
  Coefficients[6]=-Conj(Phase1)*TranslationPhase2;
  CreationIndices[7]=AnnihilationIndices[6];
  AnnihilationIndices[7]=CreationIndices[6];
  Coefficients[7]=-Conj(Coefficients[6]);
  
}


  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex ParticleOnLatticeFluxOperator::MatrixElement (RealVector& V1, RealVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  double Element = 0.0;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      for (int k=0; k<8; ++k)
	{
	  Index = this->Particle->AdA(i, this->CreationIndices[k], this->AnnihilationIndices[k], Coefficient);
	  if ((Index<Dim)&&(Coefficient!=0.0))
	    Element += V1[Index] * V2[i] * Coefficient * Coefficients[k].Re;
	}
    }
  return Complex(Element);
}


Complex ParticleOnLatticeFluxOperator::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  double Coefficient = 0.0;
  Complex Element;
  int Index;
  for (int i = 0; i < Dim; ++i)
    {
      for (int k=0; k<8; ++k)
	{
	  Index = this->Particle->AdA(i, this->CreationIndices[k], this->AnnihilationIndices[k], Coefficient);
	  if ((Index<Dim)&&(Coefficient!=0.0))
	    Element += (Conj(V1[Index]) * V2[i] * Coefficient*Coefficients[k]);
	}
    }
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

ComplexVector& ParticleOnLatticeFluxOperator::LowLevelMultiply(ComplexVector& vSource, ComplexVector& vDestination, int firstComponent, int nbrComponent)
{
  int Dim = this->Particle->GetHilbertSpaceDimension();
  int Last = firstComponent + nbrComponent;;
  int Index;
  double Coefficient = 0.0;
  for (int i = firstComponent; i < Last; ++i)
    {
      for (int k=0; k<8; ++k)
	{
	  Index = this->Particle->AdA(i, this->CreationIndices[k], this->AnnihilationIndices[k], Coefficient);
	  if ((Index<Dim)&&(Coefficient!=0.0))
	    {
	      vDestination[Index].Re = vSource[i].Re * Coefficient * Coefficients[k].Re - vSource[i].Im * Coefficient * Coefficients[k].Im;
	      vDestination[Index].Im = vSource[i].Im * Coefficient * Coefficients[k].Re + vSource[i].Re * Coefficient * Coefficients[k].Im;
	    }
	}
    }
  return vDestination;
}


