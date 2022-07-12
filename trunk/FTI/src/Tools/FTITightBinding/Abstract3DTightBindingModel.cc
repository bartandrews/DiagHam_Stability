////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                  class of abstract 3D tight binding model                  //
//                                                                            //
//                        last modification : 01/05/2012                      //
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
#include "Tools/FTITightBinding/Abstract3DTightBindingModel.h"
#include "GeneralTools/Endian.h"

#include <fstream>


using std::ofstream;
using std::endl;


// default constructor
//

Abstract3DTightBindingModel::Abstract3DTightBindingModel()
{
  this->EmbeddingZ = RealVector();
  this->Nz1 = 0;
  this->Nz2 = 0;
  this->Nx3 = 0;
  this->Ny3 = 0;
  this->Nz3 = 1;
  this->LatticeVector1 = RealVector(3, true);
  this->LatticeVector2 = RealVector(3, true);
  this->LatticeVector3 = RealVector(3, true);
  this->LatticeVector1[0] = 1.0;
  this->LatticeVector2[1] = 1.0;
  this->LatticeVector3[2] = 1.0;
  this->NbrBandsAs2DModel = 1;
}

// destructor
//

Abstract3DTightBindingModel::~Abstract3DTightBindingModel()
{
}

// get the size (length / area / volume ... ) of the unit cell
//
// return value =  size
double Abstract3DTightBindingModel::GetUnitCellSize()
{
  cout << "Unit cell size (volume) not yet defined in 3D"<<endl;
  return 0.0;
}


// write an header that describes the tight binding model
// 
// output = reference on the output stream
// return value  = reference on the output stream

ofstream& Abstract3DTightBindingModel::WriteHeader(ofstream& output)
{
  int Dimension = 3;
  int HeaderSize = ((2 * Dimension) * sizeof(double)) + ((Dimension + 1) * sizeof(int));
  WriteLittleEndian(output, HeaderSize);
  WriteLittleEndian(output, Dimension);
  WriteLittleEndian(output, this->NbrSiteX);
  WriteLittleEndian(output, this->KxFactor);
  WriteLittleEndian(output, this->GammaX);
  if (this->EmbeddingX.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingX[i]);
  }
  WriteLittleEndian(output, this->NbrSiteY);
  WriteLittleEndian(output, this->KyFactor);
  WriteLittleEndian(output, this->GammaY);
  if (this->EmbeddingY.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingY[i]);
  }
  WriteLittleEndian(output, this->NbrSiteZ);
  WriteLittleEndian(output, this->KzFactor);
  WriteLittleEndian(output, this->GammaZ);
  if (this->EmbeddingZ.GetVectorDimension() != this->NbrBands)
  {
      double Tmp = 0.0;
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, Tmp);
  }
  else
  {
      for (int i = 0; i < this->NbrBands; ++i)
          WriteLittleEndian(output, this->EmbeddingZ[i]);
  }
  return output; 
}

// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool Abstract3DTightBindingModel::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  File << "# kx    ky    kz";
  for (int i = 0; i < this->NbrBands; ++i)
    File <<  "    E_" << i;
  File << endl;
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    {
      for (int ky = 0; ky < this->NbrSiteY; ++ky)
	{
	  for (int kz = 0; kz < this->NbrSiteZ; ++kz)
	    {
	      int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(kx, ky, kz);
	      File << kx << " " << ky << " " << kz; 
	      for (int i = 0; i < this->NbrBands; ++i)
		File << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex];
	      File << endl;
	    }
	}
    }
  File.close();
  return true;
}

//computes all the values of the projected momentum and stores them in a double array
//

void Abstract3DTightBindingModel::ComputeAllProjectedMomenta()
{
  this->ProjectedMomenta = new double* [this->NbrStatePerBand];
  for (int i = 0; i < this->NbrStatePerBand; ++i)
    this->ProjectedMomenta[i] = new double [3];
  double projectedMomentum1;
  double projectedMomentum2;
  double projectedMomentum3;
  if (this->NbrConnectedOrbitals != 0)
    {
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
		{
		  double kx_trans = kx + this->Offset * ky;
		  double ky_trans = ky;
		  double kz_trans = kz;
		  projectedMomentum1 = 2.0 * M_PI * (((double) kx_trans * (double) this->Ny2) - ((double) ky_trans * (double) this->Ny1)) / ((double) (this->NbrSiteX * this->NbrSiteY));
		  projectedMomentum2 = 2.0 * M_PI * (((double) kx_trans * (double) (-this->Nx2)) + ((double) ky_trans * (double)this->Nx1)) / ((double) (this->NbrSiteX * this->NbrSiteY));
		  projectedMomentum3 = 2.0 * M_PI * ((double) kz_trans) / ((double) this->NbrSiteZ);
		  this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky, kz)][0] = projectedMomentum1;
		  this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky, kz)][1] = projectedMomentum2;
		  this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky, kz)][2] = projectedMomentum3;
		}
	    }
	}
    }
  else
    {
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
		{
		  double kx_trans = kx + this->Offset * ky + this->GammaX;
		  double ky_trans = ky + this->GammaY;
		  double kz_trans = kz + this->GammaZ;
		  projectedMomentum1 = 2.0 * M_PI * ((double) kx_trans * (double) this->Ny2 - (double) ky_trans * (double) this->Ny1) / ((double) (this->NbrSiteX * this->NbrSiteY));
		  projectedMomentum2 = 2.0 * M_PI * ((double) kx_trans * (double) (-this->Nx2) + (double) ky_trans * (double)this->Nx1) / ((double) (this->NbrSiteX * this->NbrSiteY));
		  projectedMomentum3 = 2.0 * M_PI * ((double) kz_trans) / ((double) this->NbrSiteZ);
		  this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky, kz)][0] = projectedMomentum1;
		  this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky, kz)][1] = projectedMomentum2;
		  this->ProjectedMomenta[this->GetLinearizedMomentumIndex(kx, ky, kz)][2] = projectedMomentum3;
		}
	    }
	}
    }
}

// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix Abstract3DTightBindingModel::BuildTightBindingHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  HermitianMatrix TmpHamiltonian(this->NbrBands * this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ, true);
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      for (int j = 0; j < this->NbrSiteY; ++j)
	{
	  for (int m = 0; m < this->NbrSiteZ; ++m)
	    {
	      for (int k = 0; k < this->NbrBands; ++k)
		{
		  int Index2 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, m, k);
		  for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
		    {
		      int Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(spatialIndices[k][l * 3] + i, spatialIndices[k][(l * 3) + 1] + j,
										     spatialIndices[k][(l * 3) + 2] + m, orbitalIndices[k][l]);
		      if(Index1 >= Index2)
			{
			  TmpHamiltonian.AddToMatrixElement(Index1, Index2, hoppingAmplitudes[k][l]);
			}
		    }
		}
	    }
	}      
    }
  return TmpHamiltonian;
}

// build the tight binding hamiltonian in real space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions but without assuming its hermiticiy
//
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

ComplexMatrix Abstract3DTightBindingModel::BuildTightBindingNonHermitianHamiltonianRealSpace(int* nbrConnectedOrbitals, int** orbitalIndices, int** spatialIndices, Complex** hoppingAmplitudes)
{
  ComplexMatrix TmpHamiltonian(this->NbrBands * this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ,
			       this->NbrBands * this->NbrSiteX * this->NbrSiteY * this->NbrSiteZ, true);
  for (int i = 0; i < this->NbrSiteX; ++i)
    {
      for (int j = 0; j < this->NbrSiteY; ++j)
	{
	  for (int m = 0; m < this->NbrSiteZ; ++m)
	    {
	      for (int k = 0; k < this->NbrBands; ++k)
		{
		  int Index2 = this->GetRealSpaceTightBindingLinearizedIndexSafe(i, j, m, k);
		  for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
		    {
		      int Index1 = this->GetRealSpaceTightBindingLinearizedIndexSafe(spatialIndices[k][l * 3] + i, spatialIndices[k][(l * 3) + 1] + j,
										     spatialIndices[k][(l * 3) + 2] + m, orbitalIndices[k][l]);
		      TmpHamiltonian.AddToMatrixElement(Index1, Index2, hoppingAmplitudes[k][l]);
		    }
		}
	    }
	}      
    }
  return TmpHamiltonian;
}

// build the tight binding hamiltonian in recirpocal space from the hopping parameters of the unit cell located at the origin, assuming periodic boundary conditions 
//
// kx = momentum along the x direction (in 2pi /N_x unit) for which the hamiltonian in recirpocal space has to be computed
// ky = momentum along the y direction (in 2pi /N_y unit) for which the hamiltonian in recirpocal space has to be computed
// ky = momentum along the z direction (in 2pi /N_z unit) for which the hamiltonian in recirpocal space has to be computed
// nbrConnectedOrbitals = array that gives the number of connected orbitals for each orbital within the unit cell located at the origin
// orbitalIndices = array that gives the orbital indices of the connected orbitals
// spatialIndices = array that gives the coordinates of the connected orbitals (each coordinate being a consecutive series of d integers where d is the space dimension)
// hoppingAmplitudes = array that gives the hopping amplitudes for each pair of connected orbitals
// return value = tight binding hamiltonian in real space 

HermitianMatrix Abstract3DTightBindingModel::BuildTightBindingHamiltonianReciprocalSpace(int kx, int ky, int kz, int* nbrConnectedOrbitals, int** orbitalIndices, 
											 int** spatialIndices, Complex** hoppingAmplitudes)
{
  HermitianMatrix TmpHamiltonian(this->NbrBands, true);
  double TmpKx = this->GetProjectedMomentum(kx, ky, kz, 0);
  double TmpKy = this->GetProjectedMomentum(kx, ky, kz, 1);
  double TmpKz = this->GetProjectedMomentum(kx, ky, kz, 2);
  int p;
  int q;
  int r;
  
  for (int k = 0; k < this->NbrBands; ++k)
    {
      for (int l = 0; l < nbrConnectedOrbitals[k]; ++l)
	{
	  this->GetRealSpaceIndex(spatialIndices[k][l * 3], spatialIndices[k][(l * 3) + 1], spatialIndices[k][(l * 3) + 2], p, q, r);
	  double TmpPhase = ((TmpKx * ((double) p)) + (TmpKy * ((double) q)) + (TmpKz * ((double) r)));
	  if (k >= orbitalIndices[k][l])
	    {
	      TmpHamiltonian.AddToMatrixElement(k, orbitalIndices[k][l], Conj(hoppingAmplitudes[k][l]) * Phase(TmpPhase));
	    }
	}
    }
 return TmpHamiltonian; 
}

// compute the band structure at a single point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the y axis
// kz = momentum along the z axis
// energies = array where the energies will be stored

void Abstract3DTightBindingModel::ComputeBandStructureSinglePoint(double kx, double ky, double kz, double* energies)
{
  HermitianMatrix TmpOneBodyHamiltonian = this->ComputeBlochHamiltonian(kx, ky, kz);
  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
  for (int i = 0; i < this->NbrBands; ++i)
    energies[i] = TmpDiag(i, i);
}

// compute the Bloch hamiltonian at a point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the y axis
// kz = momentum along the z axis
// return value = Bloch hamiltonian

HermitianMatrix Abstract3DTightBindingModel::ComputeBlochHamiltonian(double kx, double ky, double kz)
{
  HermitianMatrix TmpOneBodyHamiltonian;
  return TmpOneBodyHamiltonian;
}


// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void Abstract3DTightBindingModel::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  this->FindConnectedOrbitals();
  if (this->NbrConnectedOrbitals != 0)
    {
      if (nbrStates == 0l)
	nbrStates = this->NbrStatePerBand;
      long MaxStateIndex = minStateIndex + nbrStates;
      double KX;
      double KY;
      for (int kx = 0; kx < this->NbrSiteX; ++kx)
	{
	  for (int ky = 0; ky < this->NbrSiteY; ++ky)
	    {
	      for (int kz = 0; kz < this->NbrSiteZ; ++kz)
		{
		  int Index = this->GetLinearizedMomentumIndex(kx, ky, kz);
		  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
		    {
		      HermitianMatrix TmpOneBodyHamiltonian = this->BuildTightBindingHamiltonianReciprocalSpace(kx, ky, kz, this->NbrConnectedOrbitals, this->ConnectedOrbitalIndices,
														this->ConnectedOrbitalSpatialIndices, this->ConnectedOrbitalHoppingAmplitudes);
		      if (this->OneBodyBasis != 0)
			{
			  ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
			  TmpMatrix.SetToIdentity();
			  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
			  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
			  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
			  this->OneBodyBasis[Index] = TmpMatrix;
			  for (int i = 0; i < this->NbrBands; ++i)
			    {
			      this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
			    }
			}
		      else
			{
			  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
			  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
			  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
			  for (int i = 0; i < this->NbrBands; ++i)
			    {
			      this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
			    }
			}
		    }
		}
	    }
	}
    }
}

