////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on lattice			              //
//       in real space with translation invariance in two directions          //
//                                                                            //
//                        class author: Antoine Sterdyniak                    //
//                                                                            //
//                        last modification : 12/09/2014                      //
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
#include "HilbertSpace/BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "QuantumNumber/SzQuantumNumber.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Vector/RealVector.h"
#include "FunctionBasis/AbstractFunctionBasis.h"
#include "MathTools/BinomialCoefficients.h"
#include "GeneralTools/UnsignedIntegerTools.h"
#include "MathTools/FactorialCoefficient.h"
#include "GeneralTools/Endian.h"
#include "GeneralTools/ArrayTools.h"
#include "Architecture/ArchitectureOperation/FQHESphereParticleEntanglementSpectrumOperation.h"
#include "GeneralTools/StringTools.h"

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

BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation::BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation ()
{
}

// basic constructor
// 
// nbrBosons = number of bosons
// nbrSite = number of sites
// xMomentum = momentum sector in the x direction
// maxYMomentum = maximum momentum in the x direction
// yMomentum = momentum sector in the y direction
// maxYMomentum = maximum momentum in the y direction 
// memory = amount of memory granted for precalculations

BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation::BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation(int nbrBosons,  int lx, int ly, int xMomentum, int  maxXMomentum,
														   int yMomentum, int maxYMomentum, unsigned long memory):
  BosonOnLatticeRealSpaceAnd2DTranslation (nbrBosons,lx*ly,xMomentum,maxXMomentum,yMomentum,maxYMomentum,memory), Lx(lx),Ly(ly)
{  
}

// copy constructor (without duplicating datas)
//
// bosons = reference on the hilbert space to copy to copy

BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation::BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation(const BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation & bosons) :  BosonOnLatticeRealSpaceAnd2DTranslation(bosons)
{
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
}

// destructor
//

BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation::~BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation ()
{
}


// assignement (without duplicating datas)
//
// fermions = reference on the hilbert space to copy to copy
// return value = reference on current hilbert space

BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation & BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation::operator = (const BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation& bosons)
{
  BosonOnLatticeRealSpaceAnd2DTranslation::operator = (bosons);
  this->Lx = bosons.Lx;
  this->Ly = bosons.Ly;
  return (*this);
}


// clone Hilbert space (without duplicating datas)
//
// return value = pointer to cloned Hilbert space

AbstractHilbertSpace* BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation::Clone()
{
  return new BosonOnLatticeRealSpaceOneOrbitalPerSiteAnd2DTranslation(*this);
}


