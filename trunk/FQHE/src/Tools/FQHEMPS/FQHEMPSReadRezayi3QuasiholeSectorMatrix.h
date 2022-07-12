////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//              class of MPS matrix for the Read-Rezayi k=3 state             //
//                           in its quasihole sector                          //
//                                                                            //
//                        last modification : 22/02/2013                      //
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


#ifndef FQHEMPSREADREZAYI3QUASIHOLESECTORMATRIX_H
#define FQHEMPSREADREZAYI3QUASIHOLESECTORMATRIX_H


#include "config.h"
#include "Tools/FQHEMPS/FQHEMPSReadRezayi3Matrix.h"


class FQHEMPSReadRezayi3QuasiholeSectorMatrix : public FQHEMPSReadRezayi3Matrix
{

 protected:


 public:
  
  // constructor 
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital)
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  // architecture = architecture to use for precalculation
  FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, int nbrBMatrices = 2, bool bosonicVersion = false, bool useRational = true, 
					  bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
					  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0,
					  AbstractArchitecture* architecture = 0);

  // constructor 
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // nbrBMatrices = number of B matrices to compute (max occupation per orbital)
  // cftDirectory = path to the directory where all the pure CFT matrices are stored
  // bosonicVersion = use a version of the code that is compatible with bosonic wave functions
  // useRational = use arbitrary precision numbers for all the CFT calculations
  // trimChargeIndices = trim the charge indices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  // architecture = architecture to use for precalculation
  FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, int nbrBMatrices, char* cftDirectory, bool bosonicVersion = false, bool useRational = true, 
					  bool trimChargeIndices = false, bool cylinderFlag = false, double kappa = 1.0, 
					  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0, 
					  AbstractArchitecture* architecture = 0);

  // constructor from stored B matrices
  //
  // laughlinIndex = power of the Laughlin part minus 1 (i.e.  laughlinIndex=1 for the fermionic RR state)  
  // pLevel = |P| level truncation
  // fileName = name of the file that contains the B matrices
  // cylinderFlag = true if B_0 has to be normalized on the cylinder geometry
  // kappa = cylinder aspect ratio
  // torusFlag = true the torus geometry should be used instead of a genus-0 surface
  // nbrFluxQuanta = number of flux quanta piercing the torus
  // aspectRatio = aspect ratio of the torus(norm of tau)
  // angle = angle between the two vectors (i.e. 1 and tau) that span the torus (in pi unit)
  // fluxInsertion = flux insertion along the tau direction
  FQHEMPSReadRezayi3QuasiholeSectorMatrix(int laughlinIndex, int pLevel, char* fileName, bool cylinderFlag = false, double kappa = 1.0, 
					  bool torusFlag = false, int nbrFluxQuanta = 0, double aspectRatio = 1.0, double angle = 0.5, double fluxInsertion = 0.0);

  // destructor
  //
  ~FQHEMPSReadRezayi3QuasiholeSectorMatrix();
  
  // get the name describing the B matrices 
  // 
  // return value = name 
  virtual char* GetName ();

  // create the B matrices for the laughlin state
  //
  // cftDirectory = an optional path to the directory where all the CFT matrices are stored
  // architecture = architecture to use for precalculation
  virtual void CreateBMatrices (char* cftDirectory = 0, AbstractArchitecture* architecture = 0);

  // get the Q sector shift for a given CFT sector compared to the x=0 CFT sector
  //
  // cftSector = index of the CFT sector
  // return value = Q sector shift
  virtual int GetQValueCFTSectorShift(int cftSector);

  // compute the charge index range at a given truncation level
  // 
  // pLevel = tuncation level
  // cftSector = CFT sector
  // minQ = reference on the lowest charge index
  // maxQ = reference on the lowest charge index
  virtual void ComputeChargeIndexRange(int pLevel, int cftSector, int& minQ, int& maxQ);

 protected:



};


#endif
