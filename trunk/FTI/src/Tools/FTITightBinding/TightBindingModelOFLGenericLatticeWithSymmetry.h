////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2014 Nicolas Regnault                  //
//                                                                            //
//                          class author: Gunnar Möller                       //
//                                                                            //
//      class of Hamiltonian of bosons (currently) in a generic optical       //
//                    flux lattice with extended Brillouin zone               //
//                                                                            //
//                        last modification : 01/07/2014                      //
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


#ifndef TIGHTBINDINGMODELOFLGENERICLATTICEWITHSYMMETRY_H
#define TIGHTBINDINGMODELOFLGENERICLATTICEWITHSYMMETRY_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "MathTools/Complex.h"

#include "Options/Options.h"

#include <iostream>
using std::cout;
using std::endl;

// flag to switch additional debug output
#define DEBUG_OUTPUT

class TightBindingModelOFLGenericLatticeWithSymmetry : public Abstract2DTightBindingModel
{
  friend class ParticleOnLatticeOFLGenericLatticeWithSymmetryTwoBandHamiltonian;

 public:
  enum CutOffModes
  {
    Periodic,
    Square,
    Circular   
  } Mode;

 protected:

  // number of k-points in the x-/y-direction
  int NbrPoints1;
  int NbrPoints2;
  
  // maximum offset in reciprocal space terms of lattice vectors
  int NMax1;
  int NMax2;
  double CutOffMomentum;
  double LatticeDepth;
  double DeltaK1;
  double DeltaK2;

  RealVector LatticeVector1;
  RealVector LatticeVector2;

  // number of states in reciprocal unit cell
  int NbrSubLattices;
  RealMatrix SubLatticeVectors;
  // symmetry extended number of sublattice vectors
  int ExtNbrSubLattices;
  // type of flavour associated with each of the sublattice vectors
  int* SubLatticeFlavours;
  // number of distinct flavours
  int NbrSubLatticeFlavours;
  // number of sublattices per flavour
  int NbrSubLatticesPerFlavour;

  // SymmetryMultipliers: enlargement of the unit cell along LatticeVectors 1,2
  int SymmetryMultiplier1;
  int SymmetryMultiplier2;

  // pattern for evolution of flavour indices along k1, k2
  int FlavourOffset1;
  int FlavourOffset2;

  // external parameters in tight binding model
  int NbrExtParameters;
  double *ExtParameters;

  // flag indicating whether Time-reversal pairs shall be generated
  bool TRSymmetrize;
  
  // list of the jump terms in the tight binding model
  int NbrJumpTerms;

  int *DeltaG1;
  int *DeltaG2;

  int *InitialSublattice;
  int *FinalSublattice;

  Complex *JumpAmplitudes;

  // internal number of bands / sublattices in reciprocal space
  int FullNbrBands;
  
  // name of the lattice
  char *Descriptor;

 public:

  // default constructor
  //
  // nbrPointsX = number of k-points in the 1-direction
  // nbrPointsY = number of k-points in the 2-direction
  // gamma1 = boundary condition twisting angle along x
  // gamma2 = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // cutOffMode = 0: Circular, 1: Square, 2: Periodic
  // cutOffMomentum = maximum value of momenta considered (circular: norm; Square: NMax, Periodic: NMax)
  // nMax1 = maximum unit cell index in symmetry enlarged k1 direction
  // nMax2 = maximum unit cell index in symmetry enlarged k2 direction
  // latticeDepth = parameter for the depth of the optical lattice in recoil energies
  // nbrBandsToKeep = number of bands that should be calculated (total number is unbounded above)
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelOFLGenericLatticeWithSymmetry(int nbrPoints1, int nbrPoints2, double gamma1, double gamma2, AbstractArchitecture* architecture, CutOffModes cutOffMode, double cutOffMomentum, int nMax1, int nMax2, double latticeDepth, int nbrBandsToKeep=-1, double symmetryThreshold=1.0e-6, bool storeOneBodyMatrices=true);

  // destructor
  //
  ~TightBindingModelOFLGenericLatticeWithSymmetry();
  
  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // get the lattice descriptor
  char* GetDescriptor(){return this->Descriptor;}

  // get the symmetry multipliers:
  int GetSymmetryMultiplier1(){return SymmetryMultiplier1;}
  int GetSymmetryMultiplier2(){return SymmetryMultiplier2;}

  // get the number of reciprocal lattice points
  int GetNbrReciprocalVectors1(){return 2*NMax1+1;}
  int GetNbrReciprocalVectors2(){return 2*NMax2+1;}

  // get the number of spin flavours
  int GetNbrSubLatticeFlavours(){return NbrSubLatticeFlavours;}

  // get the number of sublattices per flavour
  int GetNbrSubLatticesPerFlavour(){return NbrSubLatticesPerFlavour;}

  // get the number of sublattices
  int GetNbrSublattices(){return NbrSubLattices;}

  // get the number of symmetry extended sublattices
  int GetExtNbrSublattices(){return ExtNbrSubLattices;}

  // add an option group containing all options related to the LatticeGeometry options
  //
  // manager = pointer to the option manager
  static void AddOptionGroup(OptionManager* manager);

  // pointer to the option manager
  static OptionManager* Options;

  
 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex = 0l, long nbrStates= 0l);
  
  // get linearized indices for the full reciprocal simulation cell up to CutOffMomentum
  //
  // momentum1 = linear index for offset in reciprocal lattice vector G1 (0 .. 2*NMax1)
  // momentum2 = linear index for offset in reciprocal lattice vector G2 (0 .. 2*NMax1)
  // crystalMomentum1 = crystal momentum offset in units of G1
  // crystalMomentum2 = crystal momentum offset in units of G2
  // subLattice = joint reciprocal sublattice / spin index
  // norm = norm of corresponding reciprocal lattice vector
  // inBounds = flag indicating if the site falls within the cutoff momentum
  // index for the linearized momentum, or -1 (if outside simulation cell)
  int LinearizedReciprocalSpaceIndex(int momentum1, int momentum2, double crystalMomentum1, double crystalMomentum2, int subLattice, double &norm, bool &inBounds);

 public:
  // purely periodic version of the previous function
  int PeriodicLinearizedReciprocalSpaceIndex(int momentum1, int momentum2, int extSubLattice);

  // get linearized momentum indices for k-points
  //
  int GetLinearizedMomentumIndex(int k1, int k2);

  // get the sublattice from the flavour and flavour sublattice
  // flavour = spin flavour
  // flavourSublattice = sublattice index per flavour
  int TotalSublatticeIndex(int flavour, int flavourSublattice);

  // periodize index within fundamental interval along direction d
  // coordinate = number to periodize
  // symmetryFactor = length of periodicity
  // shift = translation of coordinate necessary to end up in unit cell, in units of lattice vectors
  int Periodize(int coordinate, int symmetryFactor, int &shift);

 protected:


  // inline, see below
  int LinearizedSublatticeIndex(int Sublattice, int n1, int n2);

  // Linearized indices for final bands after symmetry expansion
  // n = fundamental band index
  // s1 = symmetry related copy along k1
  // s2 = symmetry related copy along k2
  int LinearizedBandIndex(int n, int s1, int s2);

  // inline, see below
  void ElementarySubLatticeIndices(int linearized_index, int &Sublattice, int &n1, int &n2);



  // express a hopping term for the effective extended BZ in terms of a translation by n1 v1 + n2 v2 with respect to the fundamental BZ
  //
  void ExpressInExtendedCell(int dG1, int dG2, int sub1, int sub2, int n1, int n2, int &dG1Eff, int &dG2Eff, int& sub1Eff, int& sub2Eff);

  

  // remap a vector for a translation of a given momentum
  //
  // vectorIn = input vector
  // vectorOut = input vector
  // dG1 = amount of translation along G1
  // dG2 = amount of translation along G2
  void TranslateVector(ComplexVector &vectorIn, ComplexVector &vectorOut, int dG1, int dG2);


};

// get linearized indices for the full reciprocal simulation cell up to CutOffMomentum
//
// momentum1 = linear index for offset in symmetry enhanced reciprocal lattice vector ~G1 (0 .. 2*NMax1)
// momentum2 = linear index for offset in symmetry enhanced reciprocal lattice vector ~G2 (0 .. 2*NMax1)
// crystalMomentum1 = crystal momentum offset within the Brillouin zone in units of G1
// crystalMomentum2 = crystal momentum offset within the Brillouin zone in units of G2
// extSubLattice = joint reciprocal sublattice / spin (or symmetry related points) index
// norm = norm of corresponding reciprocal lattice vector
// inBounds = flag indicating if the site falls within the cutoff momentum
// index for the linearized momentum, or -1 (if outside simulation cell)
inline int TightBindingModelOFLGenericLatticeWithSymmetry::LinearizedReciprocalSpaceIndex(int momentum1, int momentum2, double crystalMomentum1, double crystalMomentum2, int extSubLattice, double &norm, bool &inBounds)
{
  int SignedMomentum1, SignedMomentum2, Shift;
  // calculate absolute value of lattice vector
  RealVector TmpVector (2,true);
  SignedMomentum1 = momentum1 - NMax1;
  SignedMomentum2 = momentum2 - NMax2;
  int ElementarySublattice, S1, S2;
  this->ElementarySubLatticeIndices(extSubLattice, ElementarySublattice, S1, S2);
  TmpVector.AddLinearCombination((double)(SignedMomentum1*SymmetryMultiplier1+S1)+crystalMomentum1,this->LatticeVector1);
  TmpVector.AddLinearCombination((double)(SignedMomentum2*SymmetryMultiplier2+S2)+crystalMomentum2,this->LatticeVector2);
  TmpVector+=this->SubLatticeVectors[ElementarySublattice];
  norm = TmpVector.Norm();
  //cout << "this->LatticeVector1="<<endl<<this->LatticeVector1<<"this->LatticeVector2="<<endl<<this->LatticeVector2<<"TmpVector="<<endl<<TmpVector;

  switch (this->Mode)
    {
    case Circular:
      
      if ((abs(SignedMomentum1)>NMax1)||(abs(SignedMomentum2)>NMax2))
        {
          inBounds = false;
          return -1;
        }
      inBounds = (norm<this->CutOffMomentum);
      return ( ExtNbrSubLattices*(momentum1 * (2*NMax2+1) + momentum2 ) + extSubLattice );
      
    case Square:
      SignedMomentum1 = momentum1 - NMax1;
      SignedMomentum2 = momentum2 - NMax2;
      if ((abs(SignedMomentum1)>NMax1)||(abs(SignedMomentum2)>NMax2))
        {
          inBounds = false;
          return -1;
        }  
      inBounds = true;
      return ( ExtNbrSubLattices*(momentum1 * (2*NMax2+1) + momentum2 ) + extSubLattice );

    case Periodic:
      momentum1 = this->Periodize(momentum1, 2*NMax1+1, Shift);
      momentum2 = this->Periodize(momentum2, 2*NMax2+1, Shift);
      inBounds = true;
      return ( ExtNbrSubLattices*(momentum1 * (2*NMax2+1) + momentum2 ) + extSubLattice );

    default:
      std::cout << "Error: unknown mode of periodicity in TightBindingModelOFLGenericLatticeWithSymmetry"<<std::endl;
      exit(1);
      break;

    }
    
}


// get linearized indices for the full reciprocal simulation cell up to CutOffMomentum
//
// momentum1 = linear index for offset in symmetry enhanced reciprocal lattice vector ~G1 (0 .. 2*NMax1)
// momentum2 = linear index for offset in symmetry enhanced reciprocal lattice vector ~G2 (0 .. 2*NMax1)
// crystalMomentum1 = crystal momentum offset within the Brillouin zone in units of G1
// crystalMomentum2 = crystal momentum offset within the Brillouin zone in units of G2
// extSubLattice = joint reciprocal sublattice / spin (or symmetry related points) index
// norm = norm of corresponding reciprocal lattice vector
// inBounds = flag indicating if the site falls within the cutoff momentum
// index for the linearized momentum, or -1 (if outside simulation cell)
inline int TightBindingModelOFLGenericLatticeWithSymmetry::PeriodicLinearizedReciprocalSpaceIndex(int momentum1, int momentum2, int extSubLattice)
{
  int Shift;
  momentum1 = this->Periodize(momentum1, 2*NMax1+1, Shift);
  momentum2 = this->Periodize(momentum2, 2*NMax2+1, Shift);

  return ( ExtNbrSubLattices*(momentum1 * (2*NMax2+1) + momentum2 ) + extSubLattice );    
}


// get the sublattice from the flavour and flavour sublattice
// flavour = spin flavour
// flavourSublattice = sublattice index per flavour
inline int TightBindingModelOFLGenericLatticeWithSymmetry::TotalSublatticeIndex(int flavour, int flavourSublattice)
{
  int sub=0;
  for (int i=0; i<ExtNbrSubLattices; ++i)
    if (SubLatticeFlavours[i]==flavour)
      {
	if (sub == flavourSublattice)
	  {
	    //cout << "flavour="<<flavour<<", flavourSublattice="<<flavourSublattice<<", Total Sublattice ="<< i<<endl;
	    return i;
	  }
	else
	  ++sub;
      }
  std::cout << "Error: flavourSublattice could not be decoded"<<endl;
  exit(1);
  return -1;
}


//Indices for symmetry related translations of the fundamental reciprocal lattice vectors in BT
// Sublattice = fundamental sublattice
// n1 = translation along LatticeVector1
// n2 = translation along LatticeVector2
inline int TightBindingModelOFLGenericLatticeWithSymmetry::LinearizedSublatticeIndex(int Sublattice, int n1, int n2)
{
  return ( (n2*this->SymmetryMultiplier1 + n1)*this->NbrSubLattices + Sublattice);
}

// Linearized indices for final bands after symmetry expansion
// n = fundamental band index
// s1 = symmetry related copy along k1
// s2 = symmetry related copy along k2
inline int TightBindingModelOFLGenericLatticeWithSymmetry::LinearizedBandIndex(int n, int s1, int s2)
{
  return ( (n*this->SymmetryMultiplier1 + s1)*this->SymmetryMultiplier2 + s2);
}

// ElementarySubLatticeIndices : 
/// Param[in] linearized_index combined sublattice index
/// Param[out] Sublattice = fundamental sublattice
/// Param[out] n1 = translation along LatticeVector1
/// Param[out] n2 = translation along LatticeVector2
inline void TightBindingModelOFLGenericLatticeWithSymmetry::ElementarySubLatticeIndices(int linearized_index, int &Sublattice, int &n1, int &n2)
{
  Sublattice = linearized_index % NbrSubLattices;
  linearized_index/=NbrSubLattices;
  n1 = linearized_index % this->SymmetryMultiplier1;
  n2 = linearized_index / this->SymmetryMultiplier1;
}

// periodize index within fundamental interval along direction d
// coordinate = number to periodize
// symmetryFactor = length of periodicity
// shift = translation of coordinate necessary to end up in unit cell, in units of lattice vectors<
inline int TightBindingModelOFLGenericLatticeWithSymmetry::Periodize(int coordinate, int symmetryFactor, int &shift)
{
  int result;
  shift=0;
  //std::cout << "Raw value: "<<coordinate<<", symm "<<symmetryFactor;
  while (coordinate<0)
    {
      coordinate += symmetryFactor;
      shift += symmetryFactor;
    }
  shift += (result=(coordinate%symmetryFactor)) - coordinate;
  //std::cout << ", shift="<<shift<<", result="<<result<<std::endl;
  return result;
}


// express a hopping term for the effective extended BZ in terms of a translation by n1 v1 + n2 v2 with respect to the fundamental BZ
//
// dG1, dG2 = offset in terms of elementary reciprocal lattice vectors
// sub1 = initial sublattice index of the elementary BZ
// sub2 = final sublattice index in terms of the elementary BZ
// n1, n2 = offset along k1, k2 inside the extended BZ
// dG1Eff, dG2Eff [output] offset in terms of symmetry extended lattice vectors
// sub1Eff = initial sublattice index of the extended BZ
// sub2Eff = final sublattice index in terms of the extended BZ
inline void TightBindingModelOFLGenericLatticeWithSymmetry::ExpressInExtendedCell(int dG1, int dG2, int sub1, int sub2, int n1, int n2, int &dG1Eff, int &dG2Eff, int& sub1Eff, int& sub2Eff)
{
#ifdef DEBUG_OUTPUT
  std::cout << "Begin ExpressInCell (dG1="<<dG1<<", dG2="<<dG2<<", sub1="<<sub1<<", sub2="<<sub2<<", n1="<<n1<<", n2="<<n2<<")"<<std::endl;
#endif
  int n1eff = this->Periodize(dG1+n1, SymmetryMultiplier1, dG1Eff);
  int n2eff = this->Periodize(dG2+n2, SymmetryMultiplier2, dG2Eff);
  sub1Eff = this->LinearizedSublatticeIndex(sub1, n1, n2);
  sub2Eff = this->LinearizedSublatticeIndex(sub2, n1eff, n2eff);
  dG1Eff = - dG1Eff/SymmetryMultiplier1;
  dG2Eff = - dG2Eff/SymmetryMultiplier2;
#ifdef DEBUG_OUTPUT
  std::cout << "intermediate: n1eff="<<n1eff<<", n2eff="<<n2eff<<std::endl;
  std::cout << "checking: dG1="<<dG1<<", dG2="<<dG2<<", sub1="<<sub1<<", sub2="<<sub2<<", dG1Eff="<<dG1Eff<<", dG2Eff="<<dG2Eff<<", sub1Eff="<<sub1Eff<<", sub2Eff="<<sub2Eff<<std::endl;
#endif
}



// linearized momentum index in the enlarged Brillouin zone
// k1 = 0 ... NbrPoint1 * SymmetryMultiplier1
// k2 = 0 ... NbrPoint2 * SymmetryMultiplier2
// return = linearized index
inline int TightBindingModelOFLGenericLatticeWithSymmetry::GetLinearizedMomentumIndex(int k1, int k2)
{
  return (k1*this->NbrPoints2*this->SymmetryMultiplier2+k2);
}
#endif
