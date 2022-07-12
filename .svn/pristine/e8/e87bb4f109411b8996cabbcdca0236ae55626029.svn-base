////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//       class of tight binding model for the Dice lattice at flux 1/2        //
//                                                                            //
//                        last modification : 07/12/2017                      //
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


#ifndef TIGHTBINDINGMODELDICELATTICE_H
#define TIGHTBINDINGMODELDICELATTICE_H


#include "config.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"


class TightBindingModelDiceLattice : public Abstract2DTightBindingModel
{

 protected:

  // hopping amplitude between neareast neighbor sites
  double NNHopping;
  // hopping amplitude between next neareast neighbor sites
  double NextNNHopping;
  // spin orbit coupling to neareast neighbor sites
  double NNSpinOrbit;
  // spin orbit coupling to next neareast neighbor sites
  double NextNNSpinOrbit;
  
  // four times the sublattice staggered chemical potential 
  double MuS;
  
  // use the Bloch form instead of the the traditional form
  bool BlochFormFlag;

 public:

  // default constructor
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t1 = real part of the hopping amplitude between neareast neighbor sites
  // t2 = real part of the hopping amplitude between next neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A1 sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelDiceLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double lambda1, double lambda2, double mus, 
				 double gammaX, double gammaY, 
				 AbstractArchitecture* architecture, bool storeOneBodyMatrices = true, bool blochFormFlag = false);
  
  // constructor for a tilted lattice
  //
  // nbrSiteX = number of sites in the x direction
  // nbrSiteY = number of sites in the y direction
  // t1 = real part of the hopping amplitude between neareast neighbor sites
  // t2 = real part of the hopping amplitude between next neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
  // lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
  // mus = sublattice chemical potential on A1 sites
  // gammaX = boundary condition twisting angle along x
  // gammaY = boundary condition twisting angle along y
  // architecture = pointer to the architecture
  // storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
  TightBindingModelDiceLattice(int nbrSiteX, int nbrSiteY, int nx1, int ny1, int nx2, int ny2, int offset, double t1, double t2, double lambda1, double lambda2, double mus, 
				 double gammaX, double gammaY, 
				 AbstractArchitecture* architecture, bool storeOneBodyMatrices = true);


  // destructor
  //
  ~TightBindingModelDiceLattice();

  // compute the Bloch hamiltonian at a point of the Brillouin zone
  //
  // kx = momentum along the x axis
  // ky = momentum along the x axis
  // return value = Bloch hamiltonian
  virtual HermitianMatrix ComputeBlochHamiltonian(double kx, double ky);

  // get the high symmetry points 
  //
  // pointNames = name of each high symmetry point
  // pointCoordinates = coordinates in the first Brillouin zone of the high symmetry points
  // return value = number of high symmetry points
  virtual int GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates);

  // compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
  //
  // kx1 = momentum of the first point along the x axis
  // ky1 = momentum of the first point along the y axis
  // kx2 = reference on the momentum of the second point along the x axis
  // ky2 = reference on the momentum of the second point along the y axis
  // return value = distance between the two points
  virtual double GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2);

 protected :

  // core part that compute the band structure
  //
  // minStateIndex = minimum index of the state to compute
  // nbrStates = number of states to compute
  virtual void CoreComputeBandStructure(long minStateIndex, long nbrStates);
  
  // find the orbitals connected to those located at the origin unit cell
  // 
  virtual void FindConnectedOrbitals();

};


#endif
