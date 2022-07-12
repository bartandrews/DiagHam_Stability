////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the Checkerboard lattice       //
//                                                                            //
//                        last modification : 08/05/2012                      //
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


#ifndef TIGHTBINDINGMODELOFLGENERICLATTICE_H
#define TIGHTBINDINGMODELOFLGENERICLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"
#include "MathTools/Complex.h"

#include "Options/Options.h"


class TightBindingModelOFLGenericLattice : public Abstract2DTightBindingModel
{

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
  // cutOffMomentum = maximum (absolute value) of momenta considered
  // latticeDepth = parameter for the depth of the optical lattice in recoil energies
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelOFLGenericLattice(int nbrPoints1, int nbrPoints2, double gamma1, double gamma2, AbstractArchitecture* architecture, double cutOffMomentum, double latticeDepth, bool storeOneBodyMatrices=true);

  // destructor
  //
  ~TightBindingModelOFLGenericLattice();
  
  // write the energy spectrum in an ASCII file
  //
  // fileName = name of the ASCII file 
  // return value = true if no error occured
  virtual bool WriteAsciiSpectrum(char* fileName);

  // get the lattice descriptor
  char* GetDescriptor(){return this->Descriptor;}
  
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
  
  // get linearized momentum indices for k-points
  //
  int GetLinearizedMomentumIndex(int k1, int k2);

 private:
  RealVector TmpVector;

};

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
inline int TightBindingModelOFLGenericLattice::LinearizedReciprocalSpaceIndex(int momentum1, int momentum2, double crystalMomentum1, double crystalMomentum2, int subLattice, double &norm, bool &inBounds)
{
  this->TmpVector.ClearVector();
  
  int SignedMomentum1 = momentum1 - NMax1;
  int SignedMomentum2 = momentum2 - NMax2;
  TmpVector.AddLinearCombination((double)SignedMomentum1+crystalMomentum1,this->LatticeVector1);
  TmpVector.AddLinearCombination((double)SignedMomentum2+crystalMomentum2,this->LatticeVector2);
  TmpVector+=this->SubLatticeVectors[subLattice];

  if ((abs(SignedMomentum1)>NMax1)||(abs(SignedMomentum2)>NMax2))
    {
      inBounds = false;
      return -1;
    }
  norm = TmpVector.Norm();
  inBounds = (norm<this->CutOffMomentum);
    
  return ( NbrSubLattices*(momentum1 * (2*NMax2+1) + momentum2 ) + subLattice );
}


inline int TightBindingModelOFLGenericLattice::GetLinearizedMomentumIndex(int k1, int k2)
{
  return (k1*this->NbrPoints2+k2);
}
#endif
