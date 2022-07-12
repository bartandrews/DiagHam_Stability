////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//  class of bosons on linearly indexed system for a system size such that    //
//               NbrStates  < 63 or 31 (64 bits or 32bits systems)            //
//                                                                            //
//                  states indexed by y-momentum and x-position               //
//                    (Landau gauge, A(r) = \alpha x \vec e_y)                //
//                                                                            //
//                        last modification : 29/09/2008                      //
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


#ifndef HARDCOREBOSONONLATTICEKY_H
#define HARDCOREBOSONONLATTICEKY_H


#include "config.h"
#include "HilbertSpace/ParticleOnLattice.h"
#include "Matrix/RealSymmetricMatrix.h"

#include <iostream>


using std::cout;
using std::endl;
using std::dec;
using std::hex;


class HardCoreBosonOnLatticeKy : public ParticleOnLattice
{

 protected:

  // array describing each state
  unsigned long* StateDescription;
  // array giving the highest bit reached for a bososn in a given state
  int* StateHighestBit;

  // maximum shift used for searching a position in the look-up table
  int MaximumLookUpShift;
  // memory used for the look-up table in a given HighestBit sector
  unsigned long LookUpTableMemorySize;
  // shift used in each lzmax sector
  int* LookUpTableShift;
  // look-up table with two entries : the first one used lzmax value of the state an the second 
  int** LookUpTable;

  // number of bosons
  int NbrBosons;
  // total number of states
  int NbrStates;
  // length in x-direction
  int Lx;
  // length in y-direction
  int Ly;
  // many-body momentum k_y constraint
  int Ky;
  // number of sites per translation-cell
  int TranslationCell;
  // maximum many-body momentum (Ly / Kmod)
  int Kmax;
  // number of sublattices
  int NbrSublattices;
  // number of flux quanta piercing the simulation cell
  int NbrFluxQuanta;  
  // flux density (flux quanta per unit cell)
  double FluxDensity;

  // phases occurred when translating the system by a full system length (at y=1, or x=1)
  // in the x-direction
  Complex LxTranslationPhase;
  // and in the y-direction
  Complex LyTranslationPhase;  

  // pointer to an integer which indicate which coordinates are kept for the next time step iteration
  int* KeptCoordinates;
  // minors of permanents used for the time coherent wave function evaluation
  Complex** Minors;
  

  // temporary state used when applying ProdA operator
  unsigned long ProdATemporaryState;
  // Lz maximum value associated to temporary state used when applying ProdA operator
  int ProdAHighestBit;

  // pointer to the target space when an index is require after applying basic operation
  HardCoreBosonOnLatticeKy* TargetSpace;
  
 public:

  // default constructor
  //
  HardCoreBosonOnLatticeKy ();

  // basic constructor -> yields a square lattice in Landau gauge
  // 
  // nbrBosons = number of bosons
  // lx = length of simulation cell in x-direction
  // ly = length of simulation cell in y-direction
  // ky = many-body momentum in y-direction
  // nbrFluxQuanta = number of flux quanta piercing the simulation cell
  // memory = memory that can be allocated for precalculations
  HardCoreBosonOnLatticeKy (int nbrBosons, int lx, int ly, int ky, int nbrFluxQuanta, unsigned long memory = 10000000);


  // copy constructor (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  HardCoreBosonOnLatticeKy(const HardCoreBosonOnLatticeKy& bosons);

  // destructor
  //
  virtual ~HardCoreBosonOnLatticeKy ();

  // assignement (without duplicating datas)
  //
  // bosons = reference on the hilbert space to copy to copy
  // return value = reference on current hilbert space
  HardCoreBosonOnLatticeKy& operator = (const HardCoreBosonOnLatticeKy& bosons);

  // clone Hilbert space (without duplicating datas)
  //
  // return value = pointer to cloned Hilbert space
  AbstractHilbertSpace* Clone();

  // get the particle statistic 
  //
  // return value = particle statistic
  virtual int GetParticleStatistic();

  // get the quantization axis 
  //
  // return value = particle statistic
  virtual char GetLandauGaugeAxis();

  // get information about any additional symmetry of the Hilbert space
  //
  // return value = symmetry id  
  int GetHilbertSpaceAdditionalSymmetry();

  // return a list of all possible quantum numbers 
  //
  // return value = pointer to corresponding quantum number
  virtual List<AbstractQuantumNumber*> GetQuantumNumbers ();

  // return quantum number associated to a given state
  //
  // index = index of the state
  // return value = pointer to corresponding quantum number
  virtual AbstractQuantumNumber* GetQuantumNumber (int index);

  // extract subspace with a fixed quantum number
  //
  // q = quantum number value
  // converter = reference on subspace-space converter to use
  // return value = pointer to the new subspace
  virtual AbstractHilbertSpace* ExtractSubspace (AbstractQuantumNumber& q, 
						 SubspaceSpaceConverter& converter);

  // in presence of translation symmetries, it is NOT possible to change the flux through the simulation cell
  // this function is provided for compatibility with the interface, only
  // nbrFluxQuanta = number of quanta of flux piercing the simulation cell
  virtual void SetNbrFluxQuanta(int nbrFluxQuanta);

  // get maximum possible momentum for this geometry
  // return = maximum value of Ky
  virtual int GetMaximumKy();
  
  // apply creation operator to a word, using the conventions
  // for state-coding and quantum numbers of this space
  // state = word to be acted upon
  // q = quantum number of boson to be added
  virtual unsigned long Ad (unsigned long state, int q, double &coefficient);

  // apply a^+_q1 a^+_q2 a_r1 a_r2 operator to a given state (with q1+q2=r1+r2)
  //
  // index = index of the state on which the operator has to be applied
  // q1 = first index for creation operator
  // q2 = second index for creation operator
  // r1 = first index for annihilation operator
  // r2 = second index for annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAdAA (int index, int q1, int q2, int r1, int r2, double& coefficient);

  // apply a_r1 a_r2 operator to a given state. Warning, the resulting state may not belong to the current Hilbert subspace. It will be keep in cache until next AdAd call
  //
  // index = index of the state on which the operator has to be applied
  // r1 = first index for annihilation operator
  // r2 = second index for annihilation operator
  // return value =  multiplicative factor 
  virtual double AA (int index, int r1, int r2);

  // apply a^+_q1 a^+_q2 operator to the state produced using AA method (without destroying it)
  //
  // q1 = first index for creation operator
  // q2 = second index for creation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  virtual int AdAd (int q1, int q2, double& coefficient);

  // apply a^+_m a_m operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // m = index of the creation and annihilation operator
  // return value = coefficient obtained when applying a^+_m a_m
  virtual double AdA (int index, int m);

  // apply a^+_q a_r operator to a given state 
  //
  // index = index of the state on which the operator has to be applied
  // q = index of the creation operator
  // r = index of the annihilation operator
  // coefficient = reference on the double where the multiplicative factor has to be stored
  // return value = index of the destination state 
  int AdA (int index, int q, int r, double& coefficient);

  // apply \sum q U_q a^+_q a_q ( a^+_q a_q - 1 )
  // index = index of the state on which the operator has to be applied
  // nbrInteraction = number of q-values in sum, if equals NbrStates, ordered sequence 0,...,NbrStates-1 assumed
  // qValues = array of quantum numbers where an interaction is present
  // interactionPerQ = coefficient U_q of the interaction
  //
  virtual double AdAdAADiagonal(int index, int nbrInteraction, double *interactionPerQ, int *qValues);

  // code set of quantum numbers posx, posy into a single integer
  // posx = position along x-direction
  // ky = momentum in y-direction
  // sublattice = sublattice index
  // translationPhase = returns phase occurred from translating the
  //                    site to the fundamental region [0,Lx-1] x [0,Ly-1]
  virtual int EncodeQuantumNumber(int posx, int ky, int sublattice, Complex &translationPhase); 

  // decode a single encoded quantum number q to the set of quantum numbers posx, posy, sublattice
  // posx = position along x-direction
  // ky = momentum in y-direction
  // sublattice = sublattice index
  virtual void DecodeQuantumNumber(int q, int &posx, int &ky, int &sublattice);  

  // ky = true momentum in y-direction
  // fluxSubLattice = 'sublattice' index remaining after translation symmetry
  virtual int EncodeCompositeMomentum(int ky, int fluxSubLattice);

  // ky = true momentum in y-direction
  // fluxSubLattice = 'sublattice' index remaining after translation symmetry
  virtual void DecodeCompositeMomentum(int q, int &ky, int &fluxSubLattice);

  // extract the momentum ky from a quantum number q
  // return: momentum ky (in range 0...Kmax-1)
  int DecodeKy(int q);
  
  // obtain a list of quantum numbers in state
  // index = index of many-body state to be considered
  // quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
  // normalization = indicating the multiplicity of the state for bosonic spaces
  virtual void ListQuantumNumbers(int index, int *quantumNumbers, double &normalization);

  // obtain a list of quantum numbers in state
  // index = index of many-body state to be considered
  // quantumNumbers = integer array of length NbrParticles, to be written with quantum numbers of individual particles
  virtual void ListQuantumNumbers(int index, int *quantumNumbers);

  // translate a state by a multiple of the lattice vectors
  // shiftX = length of translation in x-direction
  // shiftY = length of translation in y-direction
  // translationPhase = returns phase inccurred by translation
  // return value = index of translated state
  virtual int TranslateState(int index, int shiftX, int shiftY, Complex &translationPhase);

  // find whether there is a translation vector from state i to state f
  // i = index of initial state
  // f = index of final state
  // shiftX = length of translation in x-direction
  // shiftY = length of translation in y-direction
  // return value = final state can be reached by translation
  virtual bool IsTranslation(int i, int f, int &shiftX, int &shiftY);


  // print a given State
  //
  // Str = reference on current output stream 
  // state = ID of the state to print
  // return value = reference on current output stream 
  virtual ostream& PrintState (ostream& Str, int state);

  // carefully test whether state is in Hilbert-space and find corresponding state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum nonzero bit reached by a particle in the state (can be given negative, if not known)
  // return value = corresponding index, or dimension of space, if not found
  virtual int CarefulFindStateIndex(unsigned long stateDescription, int highestBit);

 protected:

  // find state index
  //
  // stateDescription = unsigned integer describing the state
  // highestBit = maximum Lz value reached by a fermion in the state
  // return value = corresponding index
  int FindStateIndex(unsigned long stateDescription, int highestBit);

  // generate all states corresponding to the constraints
  // 
  // nbrBosons = number of bosons
  // maxQ = maximum value for the quantum number of a boson in the state
  // currentMaxQ = current max value for the quantum number of bosons that are still to be placed
  // pos = position in StateDescription array where to store states
  // currentMomentum = current value of the momentum
  // return value = position from which new states have to be stored
  int GenerateStates(int nbrBosons, int maxQ, int currentMaxQ, int pos, int currentMomentum, int debugLevel);


  // recursively evaluate Hilbert space dimension 
  //
  // nbrBosons = number of bosons
  // maxQ = maximum value for the quantum number of a boson in the state
  // currentMaxQ = current max value for the quantum number of bosons that are still to be placed
  // currentMomentum = current value of the momentum
  // return value = Hilbert space dimension
  int EvaluateHilbertSpaceDimension(int nbrBosons, int maxQ, int currentMaxQ, int currentMomentum);

  // generate look-up table associated to current Hilbert space
  // 
  // memory = memory size that can be allocated for the look-up table  
  void GenerateLookUpTable(unsigned long memory);


};


// get the particle statistic 
//
// return value = particle statistic

inline int HardCoreBosonOnLatticeKy::GetParticleStatistic()
{
  return ParticleOnLattice::BosonicStatistic;
}

// get the quantization axis 
//
// return value = particle statistic
inline char HardCoreBosonOnLatticeKy::GetLandauGaugeAxis()
{
  return 'y';
}

#endif
