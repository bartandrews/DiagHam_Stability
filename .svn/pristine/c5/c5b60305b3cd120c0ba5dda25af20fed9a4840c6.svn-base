////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                   class of spin 2 chain with translations                  //
//                          and the inversion symmetry                        //
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


#ifndef SPIN2CHAINWITHTRANSLATIONSANDINVERSIONSYMMETRY_H
#define SPIN2CHAINWITHTRANSLATIONSANDINVERSIONSYMMETRY_H


#include "config.h"
#include "HilbertSpace/Spin2ChainWithTranslationsAndSzSymmetry.h"

#include <iostream>


using std::ostream;
using std::cout;
using std::endl;
using std::hex;
using std::dec;


//  precalculation table used to apply the inversion symmetry
//
static  unsigned long Spin2ChainWithTranslationsInversionTable[] = {
0x000ul, 0x040ul, 0x080ul, 0x0c0ul, 0x100ul, 0x140ul, 0x180ul, 0x1c0ul, 0x008ul, 0x048ul, 0x088ul, 0x0c8ul, 0x108ul, 0x148ul, 0x188ul, 0x1c8ul, 
0x010ul, 0x050ul, 0x090ul, 0x0d0ul, 0x110ul, 0x150ul, 0x190ul, 0x1d0ul, 0x018ul, 0x058ul, 0x098ul, 0x0d8ul, 0x118ul, 0x158ul, 0x198ul, 0x1d8ul, 
0x020ul, 0x060ul, 0x0a0ul, 0x0e0ul, 0x120ul, 0x160ul, 0x1a0ul, 0x1e0ul, 0x028ul, 0x068ul, 0x0a8ul, 0x0e8ul, 0x128ul, 0x168ul, 0x1a8ul, 0x1e8ul, 
0x030ul, 0x070ul, 0x0b0ul, 0x0f0ul, 0x130ul, 0x170ul, 0x1b0ul, 0x1f0ul, 0x038ul, 0x078ul, 0x0b8ul, 0x0f8ul, 0x138ul, 0x178ul, 0x1b8ul, 0x1f8ul, 
0x001ul, 0x041ul, 0x081ul, 0x0c1ul, 0x101ul, 0x141ul, 0x181ul, 0x1c1ul, 0x009ul, 0x049ul, 0x089ul, 0x0c9ul, 0x109ul, 0x149ul, 0x189ul, 0x1c9ul, 
0x011ul, 0x051ul, 0x091ul, 0x0d1ul, 0x111ul, 0x151ul, 0x191ul, 0x1d1ul, 0x019ul, 0x059ul, 0x099ul, 0x0d9ul, 0x119ul, 0x159ul, 0x199ul, 0x1d9ul, 
0x021ul, 0x061ul, 0x0a1ul, 0x0e1ul, 0x121ul, 0x161ul, 0x1a1ul, 0x1e1ul, 0x029ul, 0x069ul, 0x0a9ul, 0x0e9ul, 0x129ul, 0x169ul, 0x1a9ul, 0x1e9ul, 
0x031ul, 0x071ul, 0x0b1ul, 0x0f1ul, 0x131ul, 0x171ul, 0x1b1ul, 0x1f1ul, 0x039ul, 0x079ul, 0x0b9ul, 0x0f9ul, 0x139ul, 0x179ul, 0x1b9ul, 0x1f9ul, 
0x002ul, 0x042ul, 0x082ul, 0x0c2ul, 0x102ul, 0x142ul, 0x182ul, 0x1c2ul, 0x00aul, 0x04aul, 0x08aul, 0x0caul, 0x10aul, 0x14aul, 0x18aul, 0x1caul, 
0x012ul, 0x052ul, 0x092ul, 0x0d2ul, 0x112ul, 0x152ul, 0x192ul, 0x1d2ul, 0x01aul, 0x05aul, 0x09aul, 0x0daul, 0x11aul, 0x15aul, 0x19aul, 0x1daul, 
0x022ul, 0x062ul, 0x0a2ul, 0x0e2ul, 0x122ul, 0x162ul, 0x1a2ul, 0x1e2ul, 0x02aul, 0x06aul, 0x0aaul, 0x0eaul, 0x12aul, 0x16aul, 0x1aaul, 0x1eaul, 
0x032ul, 0x072ul, 0x0b2ul, 0x0f2ul, 0x132ul, 0x172ul, 0x1b2ul, 0x1f2ul, 0x03aul, 0x07aul, 0x0baul, 0x0faul, 0x13aul, 0x17aul, 0x1baul, 0x1faul, 
0x003ul, 0x043ul, 0x083ul, 0x0c3ul, 0x103ul, 0x143ul, 0x183ul, 0x1c3ul, 0x00bul, 0x04bul, 0x08bul, 0x0cbul, 0x10bul, 0x14bul, 0x18bul, 0x1cbul, 
0x013ul, 0x053ul, 0x093ul, 0x0d3ul, 0x113ul, 0x153ul, 0x193ul, 0x1d3ul, 0x01bul, 0x05bul, 0x09bul, 0x0dbul, 0x11bul, 0x15bul, 0x19bul, 0x1dbul, 
0x023ul, 0x063ul, 0x0a3ul, 0x0e3ul, 0x123ul, 0x163ul, 0x1a3ul, 0x1e3ul, 0x02bul, 0x06bul, 0x0abul, 0x0ebul, 0x12bul, 0x16bul, 0x1abul, 0x1ebul, 
0x033ul, 0x073ul, 0x0b3ul, 0x0f3ul, 0x133ul, 0x173ul, 0x1b3ul, 0x1f3ul, 0x03bul, 0x07bul, 0x0bbul, 0x0fbul, 0x13bul, 0x17bul, 0x1bbul, 0x1fbul, 
0x004ul, 0x044ul, 0x084ul, 0x0c4ul, 0x104ul, 0x144ul, 0x184ul, 0x1c4ul, 0x00cul, 0x04cul, 0x08cul, 0x0ccul, 0x10cul, 0x14cul, 0x18cul, 0x1ccul, 
0x014ul, 0x054ul, 0x094ul, 0x0d4ul, 0x114ul, 0x154ul, 0x194ul, 0x1d4ul, 0x01cul, 0x05cul, 0x09cul, 0x0dcul, 0x11cul, 0x15cul, 0x19cul, 0x1dcul, 
0x024ul, 0x064ul, 0x0a4ul, 0x0e4ul, 0x124ul, 0x164ul, 0x1a4ul, 0x1e4ul, 0x02cul, 0x06cul, 0x0acul, 0x0ecul, 0x12cul, 0x16cul, 0x1acul, 0x1ecul, 
0x034ul, 0x074ul, 0x0b4ul, 0x0f4ul, 0x134ul, 0x174ul, 0x1b4ul, 0x1f4ul, 0x03cul, 0x07cul, 0x0bcul, 0x0fcul, 0x13cul, 0x17cul, 0x1bcul, 0x1fcul, 
0x005ul, 0x045ul, 0x085ul, 0x0c5ul, 0x105ul, 0x145ul, 0x185ul, 0x1c5ul, 0x00dul, 0x04dul, 0x08dul, 0x0cdul, 0x10dul, 0x14dul, 0x18dul, 0x1cdul, 
0x015ul, 0x055ul, 0x095ul, 0x0d5ul, 0x115ul, 0x155ul, 0x195ul, 0x1d5ul, 0x01dul, 0x05dul, 0x09dul, 0x0ddul, 0x11dul, 0x15dul, 0x19dul, 0x1ddul, 
0x025ul, 0x065ul, 0x0a5ul, 0x0e5ul, 0x125ul, 0x165ul, 0x1a5ul, 0x1e5ul, 0x02dul, 0x06dul, 0x0adul, 0x0edul, 0x12dul, 0x16dul, 0x1adul, 0x1edul, 
0x035ul, 0x075ul, 0x0b5ul, 0x0f5ul, 0x135ul, 0x175ul, 0x1b5ul, 0x1f5ul, 0x03dul, 0x07dul, 0x0bdul, 0x0fdul, 0x13dul, 0x17dul, 0x1bdul, 0x1fdul, 
0x006ul, 0x046ul, 0x086ul, 0x0c6ul, 0x106ul, 0x146ul, 0x186ul, 0x1c6ul, 0x00eul, 0x04eul, 0x08eul, 0x0ceul, 0x10eul, 0x14eul, 0x18eul, 0x1ceul, 
0x016ul, 0x056ul, 0x096ul, 0x0d6ul, 0x116ul, 0x156ul, 0x196ul, 0x1d6ul, 0x01eul, 0x05eul, 0x09eul, 0x0deul, 0x11eul, 0x15eul, 0x19eul, 0x1deul, 
0x026ul, 0x066ul, 0x0a6ul, 0x0e6ul, 0x126ul, 0x166ul, 0x1a6ul, 0x1e6ul, 0x02eul, 0x06eul, 0x0aeul, 0x0eeul, 0x12eul, 0x16eul, 0x1aeul, 0x1eeul, 
0x036ul, 0x076ul, 0x0b6ul, 0x0f6ul, 0x136ul, 0x176ul, 0x1b6ul, 0x1f6ul, 0x03eul, 0x07eul, 0x0beul, 0x0feul, 0x13eul, 0x17eul, 0x1beul, 0x1feul, 
0x007ul, 0x047ul, 0x087ul, 0x0c7ul, 0x107ul, 0x147ul, 0x187ul, 0x1c7ul, 0x00ful, 0x04ful, 0x08ful, 0x0cful, 0x10ful, 0x14ful, 0x18ful, 0x1cful, 
0x017ul, 0x057ul, 0x097ul, 0x0d7ul, 0x117ul, 0x157ul, 0x197ul, 0x1d7ul, 0x01ful, 0x05ful, 0x09ful, 0x0dful, 0x11ful, 0x15ful, 0x19ful, 0x1dful, 
0x027ul, 0x067ul, 0x0a7ul, 0x0e7ul, 0x127ul, 0x167ul, 0x1a7ul, 0x1e7ul, 0x02ful, 0x06ful, 0x0aful, 0x0eful, 0x12ful, 0x16ful, 0x1aful, 0x1eful, 
0x037ul, 0x077ul, 0x0b7ul, 0x0f7ul, 0x137ul, 0x177ul, 0x1b7ul, 0x1f7ul, 0x03ful, 0x07ful, 0x0bful, 0x0fful, 0x13ful, 0x17ful, 0x1bful, 0x1fful};


class Spin2ChainWithTranslationsAndInversionSymmetry : public Spin2ChainWithTranslationsAndSzSymmetry
{

 protected:

  // sign of the inversion symmetry sector
  double InversionSector;

  // bit shift to apply before performing the inversion symmetry
  int InversionShift;
  // bit shift to apply after performing the inversion symmetry
  int InversionUnshift;

 public:

  // default constructor
  //
  Spin2ChainWithTranslationsAndInversionSymmetry ();

  // constructor for complete Hilbert space with no restriction on total spin projection Sz
  //
  // chainLength = number of spin
  // momemtum = total momentum of each state
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslationsAndInversionSymmetry (int chainLength, int momentum, int inversionSector, unsigned long memory = 10000000);

  // constructor for complete Hilbert space corresponding to a given total spin projection Sz
  //
  // chainLength = number of spin 1
  // momentum = total momentum of each state
  // sz = twice the value of total Sz component (should be equal to zero)
  // inversionSector = inversion symmetry sector (can be either +1 or -1)
  // memory = amount of memory granted for precalculations
  Spin2ChainWithTranslationsAndInversionSymmetry (int chainLength, int momentum, int inversionSector, int sz, unsigned long memory = 10000000);

  // copy constructor (without duplicating datas)
  //
  // chain = reference on chain to copy
  Spin2ChainWithTranslationsAndInversionSymmetry (const Spin2ChainWithTranslationsAndInversionSymmetry& chain);

  // destructor
  //
  ~Spin2ChainWithTranslationsAndInversionSymmetry ();

  // assignement (without duplicating datas)
  //
  // chain = reference on chain to copy
  // return value = reference on current chain
  Spin2ChainWithTranslationsAndInversionSymmetry& operator = (const Spin2ChainWithTranslationsAndInversionSymmetry& chain);

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
  virtual RealMatrix EvaluatePartialEntanglementMatrix (int nbrSites, int szSector, RealVector& groundState, AbstractArchitecture* architecture = 0);

 protected:

  // find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
  //
  // stateDescription = unsigned integer describing the state
  // nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
  // szSymmetrySign = reference on the additional sign coming from the inversion symmetry
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

  // apply the inversion symmetry to a state description
  //
  // stateDescription = reference on the state on which the inversion symmetry has to be applied
  virtual void ApplyInversionSymmetry(unsigned long& stateDescription);

};

// find canonical form of a state description and if test if the state and its translated version can be used to create a state corresponding to themomentum constraint
//
// stateDescription = unsigned integer describing the state
// nbrTranslations = reference on the number of translations to obtain the canonical form of the resulting state
// szSymmetrySign = reference on the additional sign coming from the inversion symmetry
// return value = canonical form of a state description and -1 in nbrTranslationX if the state does not fit the momentum constraint

inline unsigned long Spin2ChainWithTranslationsAndInversionSymmetry::FindCanonicalForm(unsigned long stateDescription, int& nbrTranslations, double& szSymmetrySign)
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
//  cout << "inversion " << hex << TmpStateDescription << " ";
  this->ApplyInversionSymmetry(TmpStateDescription);
//  cout << TmpStateDescription << dec << endl;
  if (TmpStateDescription < CanonicalState)
    {
      CanonicalState = TmpStateDescription;
      nbrTranslations = 0;	      
      szSymmetrySign = this->InversionSector;      
    }
  for (int n = 1; n < this->MaxXMomentum; ++n)
    {
      this->ApplySingleXTranslation(TmpStateDescription);    
      if (TmpStateDescription < CanonicalState)
	{
	  CanonicalState = TmpStateDescription;
	  nbrTranslations = n;	      
	  szSymmetrySign = this->InversionSector;      
	}
    }
  return CanonicalState;
}

// find the size of the orbit for a given state
//
// return value = orbit size

inline int Spin2ChainWithTranslationsAndInversionSymmetry::FindOrbitSize(unsigned long stateDescription)
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
  this->ApplyInversionSymmetry(TmpStateDescription);
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

inline bool Spin2ChainWithTranslationsAndInversionSymmetry::TestMomentumConstraint(unsigned long stateDescription)
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
  this->ApplyInversionSymmetry(TmpStateDescription);
  if (stateDescription == TmpStateDescription)
    {
      if (this->InversionSector < 0.0)
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
      if (this->InversionSector < 0.0)
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

// apply the inversion symmetry to a state description
//
// stateDescription = reference on the state on which the inversion symmetry has to be applied

inline void Spin2ChainWithTranslationsAndInversionSymmetry::ApplyInversionSymmetry(unsigned long& stateDescription)
{
#ifdef __64_BITS__
  unsigned long InitialState = stateDescription & 0xfffffffffffffff8ul;  
#else
  unsigned long InitialState = stateDescription & 0xfffffff8ul;  
#endif	
  InitialState <<= this->InversionShift;
#ifdef __64_BITS__
  unsigned long TmpState = Spin2ChainWithTranslationsInversionTable[InitialState & 0x1fful] << 54;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 9) & 0x1fful] << 45;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 18) & 0x1fful] << 36;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 27) & 0x1fful] << 27;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 36) & 0x1fful] << 18;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 45) & 0x1fful] << 9;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 54) & 0x1fful]; 
#else
  unsigned long TmpState = Spin2ChainWithTranslationsInversionTable[InitialState & 0x1fful] << 27;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 9) & 0x1fful] << 18;
  TmpState |= Spin2ChainWithTranslationsInversionTable[(InitialState >> 18) & 0x1fful] << 9;
  TmpState |= Spin2ChainWithTranslationsInversionTable[InitialState >> 27];
#endif	
  TmpState >>= this->InversionUnshift;
  TmpState |= stateDescription & 0x7ul;
  stateDescription = TmpState;
}

#endif


