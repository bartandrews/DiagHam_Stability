////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                    Copyright (C) 2001-2011 Nicolas Regnault                //
//                                                                            //
//                                                                            //
//                   class of fermions on square lattice with spin            //
//                                  in momentum space                         //
//                                                                            //
//                        last modification : 16/02/2011                      //
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


#ifndef FERMIONONSQUARELATTICEWITHSPINMOMENTUMSPACETRUNCATED_H
#define FERMIONONSQUARELATTICEWITHSPINMOMENTUMSPACETRUNCATED_H

#include "config.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"

#ifdef HAVE_FTI
#include "Tools/FTITightBinding/AbstractTightBindingModel.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#endif

#include <iostream>



class FermionOnSquareLatticeWithSpinMomentumSpaceTruncated : public FermionOnSphereWithSpin
{

 protected:

  // number of sites in the x direction
  int NbrSiteX;
  // number of sites in the y direction
  int NbrSiteY;
  // momentum along the x direction
  int KxMomentum;
  // momentum along the y direction
  int KyMomentum;

  // flag to indicate that the Hilbert space should preserve Sz
  bool SzFlag;

  // pointer to the energies of the corresponding tight-binding model
  double *TightBindingEnergies;

  // Cut-Off in Energy (absolute units)
  double EnergyCutoff;

 public:

#ifdef HAVE_FTI
  // constructors relying on FTI technology...
  // options as below, except:
  //
  // tightBindingModel = pointer to the relevant tight-binding model
  // cutoff = energy cut-off in configuration space, with respect to the groundstate energy
  //
  FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, Abstract2DTightBindingModel *tightBindingModel, double cutoff, unsigned long memory = 10000000);

  FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, Abstract2DTightBindingModel *tightBindingModel, double cutoff, unsigned long memory = 10000000);
  
#endif

  
  // basic constructor
  // 
  // nbrFermions = number of fermions
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // tightBindingEnergies = pointer to the energies of a tight binding model of the corresponding band structure
  // cutoff = energy cut-off in configuration space, in absolute units
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, double *tightBindingEnergies, double cutoff, unsigned long memory = 10000000);

  // basic constructor when Sz is preserved
  // 
  // nbrFermions = number of fermions
  // nbrSpinUp = number of particles with spin up
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // kxMomentum = momentum along the x direction
  // kyMomentum = momentum along the y direction
  // tightBindingEnergies = pointer to the energies of a tight binding model of the corresponding band structure
  // cutoff = energy cut-off in configuration space, in absolute units
  // memory = amount of memory granted for precalculations
  FermionOnSquareLatticeWithSpinMomentumSpaceTruncated (int nbrFermions, int nbrSpinUp, int nbrSiteX, int nbrSiteY, int kxMomentum, int kyMomentum, double *tightBindingEnergies, double cutoff, unsigned long memory = 10000000);


  // copy constructor (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  FermionOnSquareLatticeWithSpinMomentumSpaceTruncated(const FermionOnSquareLatticeWithSpinMomentumSpaceTruncated& fermions);

  // destructor
  //
  ~FermionOnSquareLatticeWithSpinMomentumSpaceTruncated ();

  // assignement (without duplicating datas)
  //
  // fermions = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  FermionOnSquareLatticeWithSpinMomentumSpaceTruncated& operator = (const FermionOnSquareLatticeWithSpinMomentumSpaceTruncated& fermions);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

 protected:

  // evaluate Hilbert space dimension
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentEnergy = current energy of particles in configuration
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, double currentEnergy);

  // evaluate Hilbert space dimension with a fixed number of fermions with spin up
  //
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // nbrSpinUp = number of fermions with spin up
  // currentEnergy = current energy of particles in configuration
  // return value = Hilbert space dimension
  virtual long EvaluateHilbertSpaceDimension(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp, double currentEnergy);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // currentEnergy = current energy of particles in configuration
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, double currentEnergy, long pos);

  // generate all states corresponding to the constraints
  // 
  // nbrFermions = number of fermions
  // currentKx = current momentum along x for a single particle
  // currentKy = current momentum along y for a single particle
  // currentTotalKx = current total momentum along x
  // currentTotalKy = current total momentum along y
  // nbrSpinUp = number of fermions with spin up
  // currentEnergy = current energy of particles in configuration
  // pos = position in StateDescription array where to store states
  // return value = position from which new states have to be stored
  virtual long GenerateStates(int nbrFermions, int currentKx, int currentKy, int currentTotalKx, int currentTotalKy, int nbrSpinUp, double currentEnergy, long pos);

  // find state index (needs to be overridden to correctly evaluate operators leading to a truncated state)
  //
  // stateDescription = unsigned integer describing the state
  // lzmax = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  
  int FindStateIndex(unsigned long stateDescription, int lzmax);


 private:
  // find energy of a given state:
  // kx = Kx momentum
  // ky = Ky momentum
  // band = band index (down=0, up=1)
  inline double GetEnergy(int kx, int ky, int band);
    

};

inline double FermionOnSquareLatticeWithSpinMomentumSpaceTruncated::GetEnergy(int kx, int ky, int band)
{
  return this->TightBindingEnergies[(((kx * this->NbrSiteY) + ky) << 1) + band];
}



#endif


