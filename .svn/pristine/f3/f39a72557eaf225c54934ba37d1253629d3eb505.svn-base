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


#include "config.h"
#include "Tools/FTITightBinding/TightBindingModel3DChargePumpedChiralHingeCubicLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;



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

TightBindingModel3DChargePumpedChiralHingeCubicLattice::TightBindingModel3DChargePumpedChiralHingeCubicLattice(int nbrSiteX, int nbrSiteY, int nbrSiteZ,
													       double nNHoppingIntraUnitCell, double nNHoppingInterUnitCellXY,
													       double nNHoppingAlongZ, double nextNNHoppingAlongZ, double mus, 
													       double gammaX, double gammaY, double gammaZ,
													       AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->NbrSiteZ = nbrSiteZ;
  this->NbrSiteYZ = this->NbrSiteY *  this->NbrSiteZ;
  this->Nx1 = this->NbrSiteX;
  this->Ny1 = 0;
  this->Nz1 = 0;
  this->Nx2 = 0;
  this->Ny2 = this->NbrSiteY;
  this->Nz2 = 0;
  this->Nx3 = 0;
  this->Ny3 = 0;
  this->Nz3 = this->NbrSiteZ;
  this->OffsetReal = 0;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->KzFactor = 2.0 * M_PI / ((double) this->NbrSiteZ);
  this->NNHoppingIntraUnitCell = nNHoppingIntraUnitCell;
  this->NNHoppingInterUnitCellXY = nNHoppingInterUnitCellXY;
  this->NNHoppingAlongZ = nNHoppingAlongZ;
  this->NextNNHoppingAlongZ = nextNNHoppingAlongZ;
  this->MuS = mus;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->GammaZ = gammaZ;
  this->NbrBands = 4;
  this->NbrBandsAs2DModel = this->NbrBands * this->NbrSiteZ;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ;
  this->Architecture = architecture;
  
  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
  else
    {
      this->OneBodyBasis = 0;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }

  this->FindConnectedOrbitals();
  this->ComputeAllProjectedMomenta();
  this->ComputeBandStructure();
}



// destructor
//

TightBindingModel3DChargePumpedChiralHingeCubicLattice::~TightBindingModel3DChargePumpedChiralHingeCubicLattice()
{
}

// find the orbitals connected to those located at the origin unit cell
// 
  
void TightBindingModel3DChargePumpedChiralHingeCubicLattice::FindConnectedOrbitals()
{
  if (this->NbrConnectedOrbitals != 0)
    {
      return;
    }
      
  double* Embedding = new double [12];
  // site 1 within the unit cell
  Embedding[0] = 0.0;
  Embedding[1] = 0.0;
  Embedding[2] = 0.0;
  // site 2 within the unit cell
  Embedding[3] = 0.5;
  Embedding[4] = 0.0;
  Embedding[5] = 0.0;
  // site 3 within the unit cell
  Embedding[6] = 0.5;
  Embedding[7] = 0.5;
  Embedding[8] = 0.0;
  // site 4 within the unit cell
  Embedding[9] = 0.0;
  Embedding[10] = 0.5;
  Embedding[11] = 0.0;
  
  
  int ASite = 0;
  int BSite = 0;
  int TmpUnitCellA[3];
  this->GetRealSpaceIndex(0, 0, 0, TmpUnitCellA);
  int TmpUnitCellB[3];

  this->NbrConnectedOrbitals = new int [this->NbrBands];
  this->ConnectedOrbitalIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalSpatialIndices = new int* [this->NbrBands];
  this->ConnectedOrbitalHoppingAmplitudes = new Complex* [this->NbrBands];
  this->NbrConnectedOrbitals[0] = 11; 
  this->NbrConnectedOrbitals[1] = 11;
  this->NbrConnectedOrbitals[2] = 11;
  this->NbrConnectedOrbitals[3] = 11;
  
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->ConnectedOrbitalIndices[i] = new int[this->NbrConnectedOrbitals[i]];
      this->ConnectedOrbitalSpatialIndices[i] = new int[3 * this->NbrConnectedOrbitals[i]];
      this->ConnectedOrbitalHoppingAmplitudes[i] = new Complex[this->NbrConnectedOrbitals[i]];
    }
  
  this->NbrConnectedOrbitals[0] = 0;
  this->NbrConnectedOrbitals[1] = 0;
  this->NbrConnectedOrbitals[2] = 0;
  this->NbrConnectedOrbitals[3] = 0;
  
  // staggered chemical potential
  this->GetRealSpaceIndex(0, 0, 0, TmpUnitCellB);
  ASite = 0;
  BSite = 0;
  this->ConnectedOrbitalIndices[ASite][this->NbrConnectedOrbitals[ASite]] = BSite;
  this->ConnectedOrbitalSpatialIndices[ASite][this->NbrConnectedOrbitals[ASite] * 3] = TmpUnitCellB[0];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) + 1] = TmpUnitCellB[1];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) + 2] = TmpUnitCellB[2];
  this->ConnectedOrbitalHoppingAmplitudes[ASite][this->NbrConnectedOrbitals[ASite]] = this->MuS;
  ++this->NbrConnectedOrbitals[ASite];
  this->GetRealSpaceIndex(0, 0, 0, TmpUnitCellB);
  ASite = 1;
  BSite = 1;
  this->ConnectedOrbitalIndices[ASite][this->NbrConnectedOrbitals[ASite]] = BSite;
  this->ConnectedOrbitalSpatialIndices[ASite][this->NbrConnectedOrbitals[ASite] * 3] = TmpUnitCellB[0];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) +  1] = TmpUnitCellB[1];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) + 2] = TmpUnitCellB[2];
  this->ConnectedOrbitalHoppingAmplitudes[ASite][this->NbrConnectedOrbitals[ASite]] = -this->MuS;
  ++this->NbrConnectedOrbitals[ASite];
  this->GetRealSpaceIndex(0, 0, 0, TmpUnitCellB);
  ASite = 2;
  BSite = 2;
  this->ConnectedOrbitalIndices[ASite][this->NbrConnectedOrbitals[ASite]] = BSite;
  this->ConnectedOrbitalSpatialIndices[ASite][this->NbrConnectedOrbitals[ASite] * 3] = TmpUnitCellB[0];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) +  1] = TmpUnitCellB[1];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) + 2] = TmpUnitCellB[2];
  this->ConnectedOrbitalHoppingAmplitudes[ASite][this->NbrConnectedOrbitals[ASite]] = this->MuS;
  ++this->NbrConnectedOrbitals[ASite];
  this->GetRealSpaceIndex(0, 0, 0, TmpUnitCellB);
  ASite = 3;
  BSite = 3;
  this->ConnectedOrbitalIndices[ASite][this->NbrConnectedOrbitals[ASite]] = BSite;
  this->ConnectedOrbitalSpatialIndices[ASite][this->NbrConnectedOrbitals[ASite] * 3] = TmpUnitCellB[0];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) +  1] = TmpUnitCellB[1];
  this->ConnectedOrbitalSpatialIndices[ASite][(this->NbrConnectedOrbitals[ASite] * 3) + 2] = TmpUnitCellB[2];
  this->ConnectedOrbitalHoppingAmplitudes[ASite][this->NbrConnectedOrbitals[ASite]] = -this->MuS;
  ++this->NbrConnectedOrbitals[ASite];
  
  // intra-unit cell nearest neighbor hopping amplitude
  this->GetRealSpaceIndex(0, 0, 0, TmpUnitCellB);
  this->AddOrbitalConnection(this->NNHoppingIntraUnitCell, 0, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(-this->NNHoppingIntraUnitCell, 0, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(this->NNHoppingIntraUnitCell, 1, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(this->NNHoppingIntraUnitCell, 1, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(this->NNHoppingIntraUnitCell, 2, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(this->NNHoppingIntraUnitCell, 2, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(this->NNHoppingIntraUnitCell, 3, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(-this->NNHoppingIntraUnitCell, 3, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  
  // inter-unit cell nearest neighbor hopping amplitude in the XY plane
  this->GetRealSpaceIndex(-1, 0, 0, TmpUnitCellB);
  this->AddOrbitalConnection(this->NNHoppingInterUnitCellXY, 0, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, -1, 0, TmpUnitCellB);
  this->AddOrbitalConnection(-this->NNHoppingInterUnitCellXY, 0, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(1, 0, 0, TmpUnitCellB);
  this->AddOrbitalConnection(this->NNHoppingInterUnitCellXY, 1, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, -1, 0, TmpUnitCellB);
  this->AddOrbitalConnection(this->NNHoppingInterUnitCellXY, 1, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 1, 0, TmpUnitCellB);
  this->AddOrbitalConnection(this->NNHoppingInterUnitCellXY, 2, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(1, 0, 0, TmpUnitCellB);
  this->AddOrbitalConnection(this->NNHoppingInterUnitCellXY, 2, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(-1, 0, 0, TmpUnitCellB);
  this->AddOrbitalConnection(this->NNHoppingInterUnitCellXY, 3, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 1, 0, TmpUnitCellB);
  this->AddOrbitalConnection(-this->NNHoppingInterUnitCellXY, 3, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  
  // inter-unit cell nearest neighbor hopping amplitude along Z
  Complex TmpFactor = -I() * (0.5 * this->NNHoppingAlongZ);
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(TmpFactor, 0, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(Conj(TmpFactor), 0, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(Conj(TmpFactor), 1, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(TmpFactor, 1, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(TmpFactor, 2, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(Conj(TmpFactor), 2, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(Conj(TmpFactor), 3, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(TmpFactor, 3, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  
  // inter-unit cell next nearest neighbor hopping amplitude along Z
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 0, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(0.5 * this->NextNNHoppingAlongZ, 0, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 0, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(0.5 * this->NextNNHoppingAlongZ, 0, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 1, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 1, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 1, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 1, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 2, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 2, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 2, 1, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 2, 3, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, 1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 3, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(0.5 * this->NextNNHoppingAlongZ, 3, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->GetRealSpaceIndex(0, 0, -1, TmpUnitCellB);
  this->AddOrbitalConnection(-0.5 * this->NextNNHoppingAlongZ, 3, 2, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);
  this->AddOrbitalConnection(0.5 * this->NextNNHoppingAlongZ, 3, 0, TmpUnitCellA, TmpUnitCellB, Embedding, this->NbrConnectedOrbitals);   


}

