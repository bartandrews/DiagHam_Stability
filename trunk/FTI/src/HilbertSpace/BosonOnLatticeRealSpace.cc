////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                          class of bosons on lattice                        //
//                               in real space                                //
//                       class author: Antoine Sterdyniak                     //
//                                                                            //
//                        last modification : 09/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"

#include <math.h>
#include <cstdlib>
#include <fstream>

using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ifstream;
using std::ios;


// default constructor
// 

BosonOnLatticeRealSpace::BosonOnLatticeRealSpace ()
{
  this->NbrBosons = 0;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrSite = 0;
  this->LzMax = this->NbrSite;
  this->NbrLzValue = this->LzMax + 1;
  this->FermionBasis = 0;
  this->LargeHilbertSpaceDimension = 0l;
  this->HilbertSpaceDimension = 0;
  this->LargeHilbertSpaceDimension = 0;
}

// basic constructor
// 
// nbrFermions = number of fermions
// nbrSite = total number of sites 
// memory = amount of memory granted for precalculations

BosonOnLatticeRealSpace::BosonOnLatticeRealSpace (int nbrBosons, int nbrSite, unsigned long memory)
{  
  this->NbrBosons = nbrBosons;
  this->IncNbrBosons = this->NbrBosons + 1;
  this->TotalLz = 0;
  this->NbrSite = nbrSite;
  this->LzMax = this->NbrSite - 1;
  this->NbrLzValue = this->LzMax + 1;
  this->FermionBasis = new FermionOnLatticeRealSpace (this->NbrBosons, this->NbrSite+this->NbrBosons - 1);
  this->LargeHilbertSpaceDimension = this->EvaluateHilbertSpaceDimension(this->NbrBosons);
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  
  if (this->FermionBasis->GetHilbertSpaceDimension() != this->LargeHilbertSpaceDimension)
    {
      cout << "error while generating the Hilbert space, " << this->FermionBasis->GetHilbertSpaceDimension() << " generated states, should be " << this->LargeHilbertSpaceDimension << endl;
    }
  cout << "Hilbert space dimension = " << this->LargeHilbertSpaceDimension << endl;
  if (this->LargeHilbertSpaceDimension >= (1l << 30))
    this->HilbertSpaceDimension = 0;
  else
    this->HilbertSpaceDimension = (int) this->LargeHilbertSpaceDimension;
  if ( this->LargeHilbertSpaceDimension > 0l)
    {
      this->Flag.Initialize();
      this->TargetSpace = this;
    }
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeRealSpace::BosonOnLatticeRealSpace(const BosonOnLatticeRealSpace& bosons)
{
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->NbrSite = bosons.NbrSite;
  this->LzMax = bosons.LzMax;
  this->NbrLzValue = bosons.NbrLzValue;
  this->FermionBasis = (FermionOnLatticeRealSpace *) bosons.FermionBasis->Clone();
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
}

// destructor
//

BosonOnLatticeRealSpace::~BosonOnLatticeRealSpace ()
{
  if(this->FermionBasis != 0 )
    delete this->FermionBasis;
  this->FermionBasis = 0;
}

// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeRealSpace & BosonOnLatticeRealSpace::operator = (const BosonOnLatticeRealSpace & bosons)
{
  if (bosons.TargetSpace != &bosons)
    this->TargetSpace = bosons.TargetSpace;
  else
    this->TargetSpace = this;
  this->HilbertSpaceDimension = bosons.HilbertSpaceDimension;
  this->LargeHilbertSpaceDimension = bosons.LargeHilbertSpaceDimension;
  this->Flag = bosons.Flag;
  this->NbrBosons = bosons.NbrBosons;
  this->IncNbrBosons = bosons.IncNbrBosons;
  this->TotalLz = bosons.TotalLz;
  this->LzMax = bosons.LzMax;
  this->NbrSite = bosons.NbrSite;
  this->NbrLzValue = bosons.NbrLzValue;
  this->FermionBasis = (FermionOnLatticeRealSpace *) bosons.FermionBasis->Clone();
  this->TemporaryState = new unsigned long [this->NbrLzValue];
  this->ProdATemporaryState = new unsigned long [this->NbrLzValue];
  return *this;
}

// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeRealSpace::Clone()
{
  return new BosonOnLatticeRealSpace(*this);
}

// evaluate Hilbert space dimension
//
// nbrBosons = number of bosons
// return value = Hilbert space dimension
 long BosonOnLatticeRealSpace::EvaluateHilbertSpaceDimension(int nbrBosons)
{
  BinomialCoefficients binomials(this->NbrSite);
  long dimension = binomials(this->NbrSite+this->NbrBosons - 1, this->NbrBosons);
  return dimension;
}

// evaluate a density matrix of a subsystem of the whole system described by a given ground state, using particle partition. The density matrix is only evaluated in a given momentum sector.
// 
// nbrParticleSector = number of particles that belong to the subsytem 
// kxSector = kx sector in which the density matrix has to be evaluated 
// kySector = kx sector in which the density matrix has to be evaluated 
// groundState = reference on the total system ground state
// architecture = pointer to the architecture to use parallelized algorithm 
// return value = density matrix of the subsytem (return a wero dimension matrix if the density matrix is equal to zero)

HermitianMatrix BosonOnLatticeRealSpace::EvaluatePartialDensityMatrixParticlePartition (int nbrParticleSector, ComplexVector& groundState, AbstractArchitecture* architecture)
{
  if (nbrParticleSector == 0)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  if (nbrParticleSector == this->NbrBosons)
    {
      HermitianMatrix TmpDensityMatrix(1, true);
      TmpDensityMatrix(0, 0) = 1.0;
      return TmpDensityMatrix;
    }
  int ComplementaryNbrParticles = this->NbrBosons - nbrParticleSector;
  BosonOnLatticeRealSpace SubsytemSpace (nbrParticleSector, this->NbrSite);
  HermitianMatrix TmpDensityMatrix (SubsytemSpace.GetHilbertSpaceDimension(), true);
  BosonOnLatticeRealSpace ComplementarySpace (ComplementaryNbrParticles, this->NbrSite);
  cout << "subsystem Hilbert space dimension = " << SubsytemSpace.HilbertSpaceDimension << endl;

  FQHESphereParticleEntanglementSpectrumOperation Operation(this, &SubsytemSpace, &ComplementarySpace, groundState, TmpDensityMatrix);
  Operation.ApplyOperation(architecture);
  if (Operation.GetNbrNonZeroMatrixElements() > 0)	
    return TmpDensityMatrix;
  else
    {
      HermitianMatrix TmpDensityMatrixZero;
      return TmpDensityMatrixZero;
    }
}

