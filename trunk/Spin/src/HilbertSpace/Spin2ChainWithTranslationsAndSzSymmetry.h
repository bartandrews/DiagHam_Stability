////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 2 chain with translations                  //
//                           and the Sz<->-Sz symmetry                        //
//                                                                            //
//                        last modification : 12/11/2016                      //
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


#ifndef SPIN2CHAINWITHTRANSLATIONSANDSZSYMMETRY_H
#define SPIN2CHAINWITHTRANSLATIONSANDSZSYMMETRY_H


#include "config.h"
#include "HilbertSpace/Spin2ChainWithTranslations.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;


//  precalculation table used to apply the Sz<->-Sz symmetry
//
static  unsigned long Spin2ChainWithTranslationsSzInversionTable[] = {
0x124ul, 0x123ul, 0x122ul, 0x121ul, 0x120ul, 0x120ul, 0x120ul, 0x120ul, 0x11cul, 0x11bul, 0x11aul, 0x119ul, 0x118ul, 0x118ul, 0x118ul, 0x118ul, 
0x114ul, 0x113ul, 0x112ul, 0x111ul, 0x110ul, 0x110ul, 0x110ul, 0x110ul, 0x10cul, 0x10bul, 0x10aul, 0x109ul, 0x108ul, 0x108ul, 0x108ul, 0x108ul, 
0x104ul, 0x103ul, 0x102ul, 0x101ul, 0x100ul, 0x100ul, 0x100ul, 0x100ul, 0x104ul, 0x103ul, 0x102ul, 0x101ul, 0x100ul, 0x100ul, 0x100ul, 0x100ul, 
0x104ul, 0x103ul, 0x102ul, 0x101ul, 0x100ul, 0x100ul, 0x100ul, 0x100ul, 0x104ul, 0x103ul, 0x102ul, 0x101ul, 0x100ul, 0x100ul, 0x100ul, 0x100ul, 
0x0e4ul, 0x0e3ul, 0x0e2ul, 0x0e1ul, 0x0e0ul, 0x0e0ul, 0x0e0ul, 0x0e0ul, 0x0dcul, 0x0dbul, 0x0daul, 0x0d9ul, 0x0d8ul, 0x0d8ul, 0x0d8ul, 0x0d8ul, 
0x0d4ul, 0x0d3ul, 0x0d2ul, 0x0d1ul, 0x0d0ul, 0x0d0ul, 0x0d0ul, 0x0d0ul, 0x0ccul, 0x0cbul, 0x0caul, 0x0c9ul, 0x0c8ul, 0x0c8ul, 0x0c8ul, 0x0c8ul, 
0x0c4ul, 0x0c3ul, 0x0c2ul, 0x0c1ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 0x0c4ul, 0x0c3ul, 0x0c2ul, 0x0c1ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 
0x0c4ul, 0x0c3ul, 0x0c2ul, 0x0c1ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 0x0c4ul, 0x0c3ul, 0x0c2ul, 0x0c1ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 0x0c0ul, 
0x0a4ul, 0x0a3ul, 0x0a2ul, 0x0a1ul, 0x0a0ul, 0x0a0ul, 0x0a0ul, 0x0a0ul, 0x09cul, 0x09bul, 0x09aul, 0x099ul, 0x098ul, 0x098ul, 0x098ul, 0x098ul, 
0x094ul, 0x093ul, 0x092ul, 0x091ul, 0x090ul, 0x090ul, 0x090ul, 0x090ul, 0x08cul, 0x08bul, 0x08aul, 0x089ul, 0x088ul, 0x088ul, 0x088ul, 0x088ul, 
0x084ul, 0x083ul, 0x082ul, 0x081ul, 0x080ul, 0x080ul, 0x080ul, 0x080ul, 0x084ul, 0x083ul, 0x082ul, 0x081ul, 0x080ul, 0x080ul, 0x080ul, 0x080ul, 
0x084ul, 0x083ul, 0x082ul, 0x081ul, 0x080ul, 0x080ul, 0x080ul, 0x080ul, 0x084ul, 0x083ul, 0x082ul, 0x081ul, 0x080ul, 0x080ul, 0x080ul, 0x080ul, 
0x064ul, 0x063ul, 0x062ul, 0x061ul, 0x060ul, 0x060ul, 0x060ul, 0x060ul, 0x05cul, 0x05bul, 0x05aul, 0x059ul, 0x058ul, 0x058ul, 0x058ul, 0x058ul, 
0x054ul, 0x053ul, 0x052ul, 0x051ul, 0x050ul, 0x050ul, 0x050ul, 0x050ul, 0x04cul, 0x04bul, 0x04aul, 0x049ul, 0x048ul, 0x048ul, 0x048ul, 0x048ul, 
0x044ul, 0x043ul, 0x042ul, 0x041ul, 0x040ul, 0x040ul, 0x040ul, 0x040ul, 0x044ul, 0x043ul, 0x042ul, 0x041ul, 0x040ul, 0x040ul, 0x040ul, 0x040ul, 
0x044ul, 0x043ul, 0x042ul, 0x041ul, 0x040ul, 0x040ul, 0x040ul, 0x040ul, 0x044ul, 0x043ul, 0x042ul, 0x041ul, 0x040ul, 0x040ul, 0x040ul, 0x040ul, 
0x024ul, 0x023ul, 0x022ul, 0x021ul, 0x020ul, 0x020ul, 0x020ul, 0x020ul, 0x01cul, 0x01bul, 0x01aul, 0x019ul, 0x018ul, 0x018ul, 0x018ul, 0x018ul, 
0x014ul, 0x013ul, 0x012ul, 0x011ul, 0x010ul, 0x010ul, 0x010ul, 0x010ul, 0x00cul, 0x00bul, 0x00aul, 0x009ul, 0x008ul, 0x008ul, 0x008ul, 0x008ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 
0x024ul, 0x023ul, 0x022ul, 0x021ul, 0x020ul, 0x020ul, 0x020ul, 0x020ul, 0x01cul, 0x01bul, 0x01aul, 0x019ul, 0x018ul, 0x018ul, 0x018ul, 0x018ul, 
0x014ul, 0x013ul, 0x012ul, 0x011ul, 0x010ul, 0x010ul, 0x010ul, 0x010ul, 0x00cul, 0x00bul, 0x00aul, 0x009ul, 0x008ul, 0x008ul, 0x008ul, 0x008ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 
0x024ul, 0x023ul, 0x022ul, 0x021ul, 0x020ul, 0x020ul, 0x020ul, 0x020ul, 0x01cul, 0x01bul, 0x01aul, 0x019ul, 0x018ul, 0x018ul, 0x018ul, 0x018ul, 
0x014ul, 0x013ul, 0x012ul, 0x011ul, 0x010ul, 0x010ul, 0x010ul, 0x010ul, 0x00cul, 0x00bul, 0x00aul, 0x009ul, 0x008ul, 0x008ul, 0x008ul, 0x008ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 
0x024ul, 0x023ul, 0x022ul, 0x021ul, 0x020ul, 0x020ul, 0x020ul, 0x020ul, 0x01cul, 0x01bul, 0x01aul, 0x019ul, 0x018ul, 0x018ul, 0x018ul, 0x018ul, 
0x014ul, 0x013ul, 0x012ul, 0x011ul, 0x010ul, 0x010ul, 0x010ul, 0x010ul, 0x00cul, 0x00bul, 0x00aul, 0x009ul, 0x008ul, 0x008ul, 0x008ul, 0x008ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 
0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul, 0x004ul, 0x003ul, 0x002ul, 0x001ul, 0x000ul, 0x000ul, 0x000ul, 0x000ul};


class Spin2ChainWithTranslationsAndSzSymmetry : public Spin2ChainWithTranslations
{

 protected:

  // sign of the Sz<->-Sz symmetry sector
  double SzSymmetrySector;

  // mask needed for the Sz<->-Sz symmetry application
  unsigned long SzSymmetryMask;

 public:

  // default constructor
  //
  Spin2ChainWithTranslationsAndSzSymmetry ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslationsAndSzSymmetry (int chainLength, int momentum, int szSymmetrySector, unsigned long memory = 10000000);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momentum = total momentum of each state
  // sz = twice the value of total Sz component (should be equal to zero)
  // szSymmetrySector = Sz<->-Sz symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslationsAndSzSymmetry (int chainLength, int momentum, int szSymmetrySector, int sz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin2ChainWithTranslationsAndSzSymmetry (const Spin2ChainWithTranslationsAndSzSymmetry& chain);

  // destructor
  //
  ~Spin2ChainWithTranslationsAndSzSymmetry ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin2ChainWithTranslationsAndSzSymmetry& operator = (const Spin2ChainWithTranslationsAndSzSymmetry& chain);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  virtual AbstractHilbertSpace* Clone();

  // evaluate entanglement matrix of a subsystem of the whole system described by a given ground state. The entanglement matrix density matrix is only evaluated in a given Sz sector.
  // 
  // nbrSites = number of sites that are part of the A subsytem 
  // szSector = Sz sector in which the density matrix has to be evaluated 
  // groundState = reference on the total system ground state
  // architecture = pointer to the architecture to use parallelized algorithm 
  // return value = entanglement matrix of the subsytem (return a zero dimension matrix if the entanglement matrix is equal to zero)
  virtual ComplexMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, ComplexVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // generate all states with constraints 
  //
  // return value = number of generated states
  virtual long GenerateStates();

  // factorized code that is used to symmetrize the result of any operator action
  //
  // state = reference on the state that has been produced with the operator action
  // nbrStateInOrbit = original number of states in the orbit before the operator action
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // return value = index of the destination state  
  virtual int SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, int& nbrTranslations);

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
  // return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint
  virtual unsigned long FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign);

  //  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // return value = true if the state satisfies the momentum constraint
  virtual bool TestMomentumConstraint(unsigned long stateDescription);

  // find the size of the orbit for a given state
  //
  // return value = orbit size
  virtual int FindOrbitSize(unsigned long stateDescription);

  // compute the rescaling factors
  //
  virtual void ComputeRescalingFactors();

  // apply the Sz<->-Sz symmetry to a state description
  //
  // stateDescription = reference on the state on which the Sz<->-Sz symmetry has to be applied
  virtual void ApplySzSymmetry(unsigned long& stateDescription);

};

// factorized code that is used to symmetrize the result of any operator action
//
// state = reference on the state that has been produced with the operator action
// nbrStateInOrbit = original number of states in the orbit before the operator action
// coefficient = reference on the double where the multiplicative factor has to be stored
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// return value = index of the destination state  

inline int Spin2ChainWithTranslationsAndSzSymmetry::SymmetrizeResult(unsigned long& state, int nbrStateInOrbit, double& coefficient, 
								     int& nbrTranslations)
{
  double TmpSign;
  state = this->FindCanonicalForm(state, nbrTranslations, TmpSign);
  int TmpMaxMomentum = 3 * this->ChainLength;
  while (((state >> TmpMaxMomentum) == 0x0ul) && (TmpMaxMomentum > 0))
    --TmpMaxMomentum;
  int TmpIndex = this->FindStateIndex(state, TmpMaxMomentum);
  if (TmpIndex < this->HilbertSpaceDimension)
    {
      coefficient *= TmpSign;
      coefficient *= this->RescalingFactors[nbrStateInOrbit][this->NbrStateInOrbit[TmpIndex]];
      nbrTranslations = (this->MaxXMomentum - nbrTranslations) % this->MaxXMomentum;
     }
  return TmpIndex;
}

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// szSymmetrySign = reference on the additional sign coming from the Sz<->-Sz symmetry
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long Spin2ChainWithTranslationsAndSzSymmetry::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign)
{
  unsigned long CanonicalState = stateDescription;
  nbrTranslations = 0;
  szSymmetrySign = 1.0;
  unsigned long TmpStateDescription = stateDescription;  
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);      
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	}
    }
  TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslations = 0;	      
      szSymmetrySign = this->SzSymmetrySector;      
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);    
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	  szSymmetrySign = this->SzSymmetrySector;      
	}
    }
  return CanonicalState;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int Spin2ChainWithTranslationsAndSzSymmetry::FindOrbitSize(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    return XSize;  
  int XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  
  if (XSize2 != XSize)
    {
      return XSize;
    }
  return (2 * XSize);
}

//  test if the state and its translated version can be used to create a state corresponding to the momentum constraint
//
// stateDescription = unsigned integer describing the state
// return value = true if the state satisfies the momentum constraint

inline bool Spin2ChainWithTranslationsAndSzSymmetry::TestMomentumConstraint(unsigned long stateDescription)
{
  unsigned long TmpStateDescription = stateDescription;
  int XSize = 1;
  this->ApplySingleXTranslation(TmpStateDescription);   
  while (stateDescription != TmpStateDescription)
    {
      ++XSize;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }
  if (((this->Momentum * XSize) % this->MaxXMomentum) != 0)
    return false;
  TmpStateDescription = stateDescription;
  this->ApplySzSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      if (this->SzSymmetrySector < 0.0)
	return false;
      else
	return true;
    }
  int XSize2 = 1;
  this->ApplySingleXTranslation(TmpStateDescription);      
  while ((stateDescription != TmpStateDescription) && (XSize2 < XSize))
    {
      ++XSize2;
      this->ApplySingleXTranslation(TmpStateDescription);      
    }  

  if (XSize2 < XSize)
    {
      if (this->SzSymmetrySector < 0.0)
	{
	  if ((((this->Momentum * XSize2 * 2) + this->MaxXMomentum) % (2 * this->MaxXMomentum)) != 0)
	    return false;
	  else
	    return true;
	}
      if ((((this->Momentum * XSize2)) % this->MaxXMomentum) != 0)
	return false;
      else
	return true;  
    }
  return true;
}

// apply the Sz<->-Sz symmetry to a state description
//
// stateDescription = reference on the state on which the Sz<->-Sz symmetry has to be applied

inline void Spin2ChainWithTranslationsAndSzSymmetry::ApplySzSymmetry(unsigned long& stateDescription)
{
#ifdef __64_BITS__
  unsigned long TmpState = Spin2ChainWithTranslationsSzInversionTable[stateDescription & 0x1fful];
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 9) & 0x1fful] << 9;
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 18) & 0x1fful] << 18;
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 27) & 0x1fful] << 27;
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 36) & 0x1fful] << 36;
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 45) & 0x1fful] << 45;
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 54) & 0x1fful] << 54;
#else
  unsigned long TmpState = Spin2ChainWithTranslationsSzInversionTable[stateDescription & 0x1fful];
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 9) & 0x1fful] << 9;
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[(stateDescription >> 18) & 0x1fful] << 18;
  TmpState |= Spin2ChainWithTranslationsSzInversionTable[stateDescription >> 27] << 27;
#endif	
  TmpState &= this->SzSymmetryMask;
  stateDescription = TmpState;
}

#endif


