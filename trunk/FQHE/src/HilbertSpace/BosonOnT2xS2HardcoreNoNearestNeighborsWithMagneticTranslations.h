////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            Diagham  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2014 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of bosons on the 4D manifold T2 x S2,              //
// forbidding mutliple orbital occupations and nearest orbital occupations    //
//                       and using magnetic translations                      //
//                                                                            //
//                        last modification : 06/03/2017                      //
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


#ifndef BOSONONT2XS2HARDCORENONEARESTNEIGHBORSWITHMAGNETICTRANSLATIONS_H
#define BOSONONT2XS2HARDCORENONEARESTNEIGHBORSWITHMAGNETICTRANSLATIONS_H


#include "config.h"
#include "HilbertSpace/FermionOnTorusWithMagneticTranslations.h"

#include <iostream>



class BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations : public FermionOnTorusWithMagneticTranslations
{
  
 protected:

  // number of flux quanta piercing the torus
  int NbrFluxQuantumTorus;
  // number of flux quanta piercing the sphere
  int NbrFluxQuantumSphere;
  // number of orbitals for the sphere
  int NbrOrbitalSphere;
  // total momentum along the z axis 
  int TotalLz;
  // projection of the total angular momentum along the z axis, using the disk convention
  int ShiftedTotalLz;
  
 public:

  // default constructor
  // 
  BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations ();

  // basic constructor
  // 
  // nbrBosons = number of bosons
  // nbrFluxQuantumTorus = number of flux quanta piercing the torus
  // kxMomentum = momentum in the x direction for the torus
  // kyMomentum = momentum in the y direction for the torus
  // nbrFluxQuantumSphere = number of flux quanta piercing the sphere
  // totalLz = projection of the total angular momentum along the z axis for the sphere
  // memory = amount of memory granted for precalculations
  BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations (int nbrBosons, int nbrFluxQuantumTorus, int kxMomentum, int kyMomentum,
								  int nbrFluxQuantumSphere, int totalLz, unsigned long memory = 10000000);

  // copy constructor (without duplicating data)
  //
  // space = reference on the hilbert space to copy
  BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations(const BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations& space);

  // destructor
  //
  ~BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations & operator = (const BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations & bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // convert a state such that its components are now expressed in the unnormalized basis
  //
  // state = reference to the state to convert
  // reference = set which component has to be normalized to 1
  // symmetryFactor = if true also remove the symmetry factors
  // return value = converted state
  //  virtual RealVector& ConvertToUnnormalizedMonomial(RealVector& state, long reference = 0, bool symmetryFactor = true);

  // print a given State using the monomial notation
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  //  virtual ostream& PrintStateMonomial (ostream& Str, long state);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual bool HasPauliExclusions(int index, int pauliK, int pauliR);
 
  // get the occupation of a given orbital in a state, not assuming the coordinates are valid
  //
  // index = state index
  // xPosition = orbital x coordinates
  // yPosition = orbital y coordinates
  // return value = orbital occupation
  virtual unsigned long GetSafeOccupation(int index, int xPosition, int yPosition);

  // get the occupation of a given orbital in a state, not assuming the coordinates are valid, assume periodic bounadary condtions along x
  //
  // index = state index
  // xPosition = orbital x coordinates
  // yPosition = orbital y coordinates
  // return value = orbital occupation
  virtual unsigned long GetSafeOccupationWithPBC(int index, int xPosition, int yPosition);

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

 protected:

  // find canonical form of a state description and
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslations = reference on the number of translationsto obtain the canonical form of the resulting state
  // return value = canonical form of a state description
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  // apply a single translation for a state description
  //
  // stateDescription = reference on the state description  
  virtual void ApplySingleTranslation(unsigned long& stateDescription);

  // evaluate Hilbert space dimension
  //
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy);

  // generate all states corresponding to the constraints
  //
  // return value = Hilbert space dimension
  virtual long GenerateStates();

  // generate all states corresponding to the constraints
  // 
  // stateDescription = array that gives each state description
  // nbrBosons = number of bosons
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentFermionicPosition = current fermionic position within the state description
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrBosons, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int currentFermionicPosition, long pos);

  // request whether state with given index satisfies a general Pauli exclusion principle
  // index = state index
  // pauliK = number of particles allowed in consecutive orbitals
  // pauliR = number of consecutive orbitals
  virtual unsigned long FindCanonical(unsigned long state, int xPosition, int yPosition);

  // get a linearized index from the two momenta
  //
  // ky = momentum along the y direction for the torus
  // lz = z projection of the angular momentum for the sphere
  // return value = linearized index 
  virtual int GetLinearizedIndex(int ky, int lz);

  // get the two momenta associated to a given linearized index
  //
  // index = linearized index 
  // ky = reference on the momentum along the y direction for the torus
  // lz = reference on the z projection of the angular momentum for the sphere
  virtual void GetLinearizedIndex(int index, int& ky, int& lz);

};

// get a linearized index from the two momenta
//
// ky = momentum along the y direction for the torus
// lz = z projection of the angular momentum for the sphere
// return value = linearized index 

inline int BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::GetLinearizedIndex(int ky, int lz)
{
  return ((ky * this->NbrOrbitalSphere) + lz);
}

// get the two momenta associated to a given linearized index
//
// index = linearized index 
// ky = reference on the momentum along the y direction for the torus
// lz = reference on the z projection of the angular momentum for the sphere

inline void BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::GetLinearizedIndex(int index, int& ky, int& lz)
{
  lz = index % this->NbrOrbitalSphere;
  ky = index / this->NbrOrbitalSphere;
}

// get the occupation of a given orbital in a state, not assuming the coordinates are valid
//
// index = state index
// xPosition = orbital x coordinates
// yPosition = orbital y coordinates
// return value = orbital occupation

inline unsigned long BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::GetSafeOccupation(int index, int xPosition, int yPosition)
{
  if ((xPosition < 0) || (yPosition < 0) || (xPosition >= this->NbrFluxQuantumTorus) || (yPosition >= this->NbrOrbitalSphere))
    {
      return 0x0ul;
    }
  else
    {
      return ((this->StateDescription[index] >> ((xPosition * this->NbrOrbitalSphere) + yPosition)) & 0x1ul);
    }
}
// get the occupation of a given orbital in a state, not assuming the coordinates are valid
//
// index = state index
// xPosition = orbital x coordinates
// yPosition = orbital y coordinates
// return value = orbital occupation

inline unsigned long BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::GetSafeOccupationWithPBC(int index, int xPosition, int yPosition)
{
  if ((yPosition < 0) || (yPosition >= this->NbrOrbitalSphere))
    {
      return 0x0ul;      
    }
  while (xPosition < 0)
    {
      xPosition += this->NbrFluxQuantumTorus;
    }
  while (xPosition >= this->NbrFluxQuantumTorus)
    {
      xPosition -= this->NbrFluxQuantumTorus;
    }
  return ((this->StateDescription[index] >> ((xPosition * this->NbrOrbitalSphere) + yPosition)) & 0x1ul);
}

// find canonical form of a state description and
//
// stateDescription = unsigned integer describing the state
// nbrTranslations = reference on the number of translationsto obtain the canonical form of the resulting state
// return value = canonical form of a state description

inline unsigned long BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations)
{
  unsigned long CanonicalState = stateDescription;
  unsigned long TmpStateDescription;  
  nbrTranslations = 0;
  TmpStateDescription = stateDescription;
  for (int n = 1; n < this->MaxMomentum; ++n)
    {
      this->ApplySingleTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	}
    }
  return CanonicalState;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleTranslation(TmpStateDescription);      
    }
  return XSize;
}

// apply a single translation for a state description
//
// stateDescription = reference on the state description

inline void BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::ApplySingleTranslation(unsigned long& stateDescription)
{
  stateDescription = (stateDescription >> this->StateShift) | ((stateDescription & this->MomentumMask) << this->ComplementaryStateShift);
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool BosonOnT2xS2HardcoreNoNearestNeighborsWithMagneticTranslations::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleTranslation(TmpStateDescription);   
  while (stateDescription != TmpStateDescription)
    { 
      this->ApplySingleTranslation(TmpStateDescription);   
      ++XSize;
    }
  if (((this->XMomentum * XSize) % this->MomentumModulo) != 0)
    return false;
  return true;
}

#endif

