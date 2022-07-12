////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//            class of tight binding model for the three-dimensional          //
//           second order TI with chiral hinge states, corresponding          //
//              to dipole pumping of the topological quadrupole               //
//            Schindler et al. Science advances 4, eaat0346 (2018).           //
//                                                                            //
//                        last modification : 18/03/2020                      //
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


#ifndef TIGHTBINDINGMODEL3DFCHARGEPUMPEDCHIRALHINGECUBICLATTICE_H
#define TIGHTBINDINGMODEL3DFCHARGEPUMPEDCHIRALHINGECUBICLATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"


class TightBindingModel3DChargePumpedChiralHingeCubicLattice : public Abstract3DTightBindingModel
{

 protected:


  // neareast neighbor hopping amplitude within the unit cell
  double NNHoppingIntraUnitCell;
  // neareast neighbor hopping amplitude between unit cells in the XY plane
  double NNHoppingInterUnitCellXY;
  // neareast neighbor hopping amplitude along the Z direction
  double NNHoppingAlongZ;
  // next neareast neighbor hopping amplitude along the Z direction
  double NextNNHoppingAlongZ;
  
  // sublattice chemical potential on even sites
  double MuS;

public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // nbrSiteZ = number of sites in the z direction
  // nNHoppingIntraUnitCell = neareast neighbor hopping amplitude within the unit cell
  // nNHoppingInterUnitCellXY = neareast neighbor hopping amplitude between unit cells in the XY plane
  // nNHoppingAlongZ = neareast neighbor hopping amplitude along the Z direction
  // nextNNHoppingAlongZ = next neareast neighbor hopping amplitude along the Z direction
  // mus = sublattice chemical potential on even sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // gammaZ = boundary condition twisting angle along z
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModel3DChargePumpedChiralHingeCubicLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
							 double nNHoppingIntraUnitCell, double nNHoppingInterUnitCellXY,
							 double nNHoppingAlongZ, double nextNNHoppingAlongZ, double mus,
							 double gammaX, double gammaY,  double gammaZ, 
							 AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);
  
  // destructor
  //
  ~TightBindingModel3DChargePumpedChiralHingeCubicLattice();

 protected :

  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

  // get the accumulated phase due to the twisted boundray conditions along a straight line connecting two lattice points from A to B
  //
  // unitCellCoordinateA = array containing the unit cell coordinates for point A
  // embeddingA  = array containing the site coordinates within the unit cell for point A
  // unitCellCoordinateB = array containing the unit cell coordinates for point A
  // embeddingB  = array containing the site coordinates within the unit cell for point A
  // return value = accumulated phase 
  virtual Complex GetTwistedBoundaryConditionPhase(int* unitCellCoordinateA, double* embeddingA, int* unitCellCoordinateB, double* embeddingB);

  // add the contribution from one oriented link (from A to B) to the connected orbital arrays
  //
  // hoppingAmplitude = reference on the hopping amplitude
  // aSite = intra unit cell index of the A site
  // bSite = intra unit cell index of the A site
  // unitCellCoordinateA = array containing the unit cell coordinates for point A
  // unitCellCoordinateB = array containing the unit cell coordinates for point B
  // embeddings = array containing all the intra unit cell site embeddings
  // tmpIndex = array of indices for the connected orbital arrays
  virtual void AddOrbitalConnection (Complex hoppingAmplitude, int aSite, int bSite, int* unitCellCoordinateA, int* unitCellCoordinateB, double* embeddings, int* tmpIndex);

};


// get the accumulated phase due to the twisted boundray conditions along a straight line connecting two lattice points from A to B
//
// unitCellCoordinateA = array containing the unit cell coordinates for point A
// embeddingA = array containing the site coordinates within the unit cell for point A
// unitCellCoordinateB = array containing the unit cell coordinates for point B
// embeddingB = array containing the site coordinates within the unit cell for point B
// return value = accumulated phase 

inline Complex TightBindingModel3DChargePumpedChiralHingeCubicLattice::GetTwistedBoundaryConditionPhase(int* unitCellCoordinateA, double* embeddingA, int* unitCellCoordinateB, double* embeddingB)
{
  return Phase(this->KxFactor * this->GammaX * (((double) (unitCellCoordinateB[0] - unitCellCoordinateA[0])) + embeddingB[0] - embeddingA[0])
	       + this->KyFactor * this->GammaY * (((double) (unitCellCoordinateB[1] - unitCellCoordinateA[1])) + embeddingB[1] - embeddingA[1])
	       + this->KzFactor * this->GammaZ * (((double) (unitCellCoordinateB[2] - unitCellCoordinateA[2])) + embeddingB[2] - embeddingA[2]));
}

// add the contribution from one oriented link (from A to B) to the connected orbital arrays
//
// hoppingAmplitude = reference on the hopping amplitude
// aSite = intra unit cell index of the A site
// bSite = intra unit cell index of the A site
// unitCellCoordinateA = array containing the unit cell coordinates for point A
// unitCellCoordinateB = array containing the unit cell coordinates for point B
// embeddings = array containing all the intra unit cell site embeddings
// tmpIndex = array of indices for the connected orbital arrays

inline void TightBindingModel3DChargePumpedChiralHingeCubicLattice::AddOrbitalConnection (Complex hoppingAmplitude, int aSite, int bSite,
											  int* unitCellCoordinateA, int* unitCellCoordinateB, double* embeddings, int* tmpIndex)
{
  this->ConnectedOrbitalIndices[aSite][tmpIndex[aSite]] = bSite;
  this->ConnectedOrbitalSpatialIndices[aSite][tmpIndex[aSite] * 3] = unitCellCoordinateB[0];
  this->ConnectedOrbitalSpatialIndices[aSite][(tmpIndex[aSite] * 3) + 1] = unitCellCoordinateB[1];
  this->ConnectedOrbitalSpatialIndices[aSite][(tmpIndex[aSite] * 3) + 2] = unitCellCoordinateB[2];
  this->ConnectedOrbitalHoppingAmplitudes[aSite][tmpIndex[aSite]] = hoppingAmplitude * this->GetTwistedBoundaryConditionPhase(unitCellCoordinateA, &(embeddings[aSite * 3]),
																	  &(this->ConnectedOrbitalSpatialIndices[aSite][tmpIndex[aSite] * 3]),
																	  &(embeddings[bSite * 3]));
  ++tmpIndex[aSite];
}

#endif
