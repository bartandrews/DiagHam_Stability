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

#ifndef PARTICLEONLATTICEVACANCYOPERATOR_H
#define PARTICLEONLATTICEVACANCYOPERATOR_H

#include "config.h"
#include "GeneralTools/GarbageFlag.h"
#include "Operator/AbstractOperator.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Tools/FQHESpectrum/LatticePhases.h"
#include "MathTools/Complex.h"

class ParticleOnLatticeVacancyOperator : public AbstractOperator
{

 protected:

  // hilbert space associated to the source and target spaces
  ParticleOnLattice* SourceSpace;
  ParticleOnLattice* TargetSpace;
  
  // indices of the creation operator
  int Kx;
  // index of the annihilation operator
  int Ky;

  // sublattice-index shared by all operators
  int SubLattice;


  // flag indicating whether a particle or hole is considered
  bool ParticleFlag;

  // reference to generic lattice
  LatticePhases *Lattice;

  // number of cells / individual operators
  int NbrCells;
  
  // phase factors for current wavevector
  Complex *Phases;
  // corresponding quantum numbers
  int *QuantumNbrs;
 public:
  
  // constructor from default datas
  //
  // particle = hilbert space associated to the particles
  // creationIndex = index of the creation operator
  // annihilationIndex = index of the annihilation operator
  ParticleOnLatticeVacancyOperator(ParticleOnLattice* source, ParticleOnLattice* target, LatticePhases *lattice,
				   int kx=0, int ky=0, int sublattice=0, bool particle=false);

  // copy constructor
  //
  // oper = reference on the operator to copy
  ParticleOnLatticeVacancyOperator(const ParticleOnLatticeVacancyOperator& oper);

  // destructor
  //
  ~ParticleOnLatticeVacancyOperator();
  
  // clone operator without duplicating datas
  //
  // return value = pointer to cloned hamiltonian
  AbstractOperator* Clone ();

  // set Hilbert space
  //
  // hilbertSpace = pointer to source Hilbert space to use
  void SetHilbertSpace (AbstractHilbertSpace* hilbertSpace);

  // get Hilbert space on which operator acts
  //
  // return value = pointer to used Hilbert space
  AbstractHilbertSpace* GetHilbertSpace ();

  // return dimension of Hilbert space where operator acts
  //
  // return value = corresponding matrix elementdimension
  int GetHilbertSpaceDimension ();

  // change indices of creation / annihilation operators
  // creationIndex = index of the creation operator
  // annihilationIndex = index of the annihilation operator
  void SetMomenta (int kx, int ky, int sublattice); 
  
  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (RealVector& V1, RealVector& V2, long firstComponent, long nbrComponent);

  // evaluate part of the matrix element, within a given of indices
  //
  // V1 = vector to left multiply with current matrix
  // V2 = vector to right multiply with current matrix
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = corresponding matrix element
  Complex PartialMatrixElement (ComplexVector& V1, ComplexVector& V2, long firstComponent, long nbrComponent);
   
  // multiply a vector by the current operator for a given range of indices 
  // and store result in another vector
  //
  // vSource = vector to be multiplied
  // vDestination = vector where result has to be stored
  // firstComponent = index of the first component to evaluate
  // nbrComponent = number of components to evaluate
  // return value = reference on vector where result has been stored
  ComplexVector& LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
			       int firstComponent, int nbrComponent);
  
};

#endif 
