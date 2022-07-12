////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                         Class author : Cecile Repellin                     //
//                                                                            //
//         class of tight binding model for the Kagome lattice                //
//                     Time Reversal Invariant Model                          //
//                   last modification : 16/04/2013                           //
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
#include "Tools/FTITightBinding/TightBindingModelTimeReversalKagomeLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

using std::cout;
using std::endl;
// default constructor
//
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// t1 = real part of the hopping amplitude between neareast neighbor sites
// t2 = real part of the hopping amplitude between next neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between neareast neighbor sites
// lambda1 = imaginary part of the hopping amplitude between next neareast neighbor sites
// mixingTerm12 = mixing term coupling the two copies of the kagome lattice (sites 1 and 2)
// mixingTerm13 = mixing term coupling the two copies of the kagome lattice (sites 1 and 3)
// mixingTerm23 = mixing term coupling the two copies of the kagome lattice (sites 2 and 3)
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelTimeReversalKagomeLattice::TightBindingModelTimeReversalKagomeLattice(int nbrSiteX, int nbrSiteY, double t1, double t2, double lambda1, double lambda2, double mixingTerm12, double mixingTerm13, double mixingTerm23, double gammaX, double gammaY, AbstractArchitecture* architecture, bool timeReversalFlag, bool storeOneBodyMatrices)
{
   
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->NNHopping = t1;
  this->NextNNHopping = t2;
  this->NNSpinOrbit = lambda1;
  this->NextNNSpinOrbit = lambda2;

  this->MixingTerm12 = mixingTerm12;
  this->MixingTerm13 = mixingTerm13;
  this->MixingTerm23 = mixingTerm23;
 
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 6;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
  this->Architecture = architecture;
  this->TimeReversal = timeReversalFlag;

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
  this->ComputeBandStructure();
}

// destructor
//

TightBindingModelTimeReversalKagomeLattice::~TightBindingModelTimeReversalKagomeLattice()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelTimeReversalKagomeLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
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
	  int Index = (kx * this->NbrSiteY) + ky;
	  
	  HermitianMatrix TmpOneBodyHamiltonian = this->ComputeBlochHamiltonian(this->KxFactor* ((double) kx), this->KyFactor * ((double) ky));
	  	  
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
		this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
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
		this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
	    }
	}
    }
}



// compute the Bloch hamiltonian at a point of the Brillouin zone
//
// kx = momentum along the x axis
// ky = momentum along the x axis
// return value = Bloch hamiltonian

HermitianMatrix TightBindingModelTimeReversalKagomeLattice::ComputeBlochHamiltonian(double kx, double ky)
{
  double KX = kx + (this->KxFactor * this->GammaX);
  double KY = ky + (this->KyFactor * this->GammaY);
  double InvKX = -kx + (this->KxFactor * this->GammaX);
  double InvKY = -ky + (this->KyFactor * this->GammaY);
  double Sign = -1.0;

  HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
  Complex HAB (-this->NNHopping, -this->NNSpinOrbit);
  HAB *= 1 + Phase(KX);
  Complex HAC(-this->NNHopping, this->NNSpinOrbit);
  HAC *= 1 + Phase(KY);
  Complex HBC(- this->NNHopping, - this->NNSpinOrbit);
  HBC *= 1 + Phase(KY - KX);
  
  
  Complex InvHAB = Complex(- this->NNHopping, - this->NNSpinOrbit);
  InvHAB *= 1 + Phase(InvKX);
  Complex InvHAC = Complex(- this->NNHopping,  this->NNSpinOrbit);
  InvHAC *= 1 + Phase(InvKY);
  Complex InvHBC = Complex(- this->NNHopping, - this->NNSpinOrbit);
  InvHBC *= 1 + Phase(InvKY - InvKX);
  
  TmpOneBodyHamiltonian.SetMatrixElement(0, 1, HAB);
  TmpOneBodyHamiltonian.SetMatrixElement(0, 2, HAC);
  TmpOneBodyHamiltonian.SetMatrixElement(1, 2, HBC);
  if (this->TimeReversal == true)
    {
      TmpOneBodyHamiltonian.SetMatrixElement(3, 4, Conj(InvHAB));
      TmpOneBodyHamiltonian.SetMatrixElement(3, 5, Conj(InvHAC));
      TmpOneBodyHamiltonian.SetMatrixElement(4, 5, Conj(InvHBC));
    }
  else
    {
      TmpOneBodyHamiltonian.SetMatrixElement(3, 4, HAB);
      TmpOneBodyHamiltonian.SetMatrixElement(3, 5, HAC);
      TmpOneBodyHamiltonian.SetMatrixElement(4, 5, HBC);
    }
  
  TmpOneBodyHamiltonian.SetMatrixElement(0, 4, this->MixingTerm12);
  TmpOneBodyHamiltonian.SetMatrixElement(0, 5, this->MixingTerm13);
  TmpOneBodyHamiltonian.SetMatrixElement(1, 3, Sign * this->MixingTerm12);
  TmpOneBodyHamiltonian.SetMatrixElement(1, 5, this->MixingTerm23);
  TmpOneBodyHamiltonian.SetMatrixElement(2, 3, Sign * this->MixingTerm13);
  TmpOneBodyHamiltonian.SetMatrixElement(2, 4, Sign * this->MixingTerm23);
  return TmpOneBodyHamiltonian;
}

// get the high symmetry points 
//
// pointNames = name of each high symmetry point
// pointCoordinates = coordinates in the first Brillouin zone of the high symmetry points
// return value = number of high symmetry points

int TightBindingModelTimeReversalKagomeLattice::GetHighSymmetryPoints(char**& pointNames, double**& pointCoordinates)
{
  int NbrHighSymmetryPoints = 3;
  pointNames = new char*[NbrHighSymmetryPoints];
  pointCoordinates = new double*[NbrHighSymmetryPoints];

  pointNames[0] = new char[16];
  sprintf (pointNames[0], "Gamma");
  pointCoordinates[0] = new double[2];
  pointCoordinates[0][0] = 0.0; 
  pointCoordinates[0][1] = 0.0; 

  pointNames[1] = new char[16];
  sprintf (pointNames[1], "K");
  pointCoordinates[1] = new double[2];
  pointCoordinates[1][0] = 4.0 * M_PI / 3.0; 
  pointCoordinates[1][1] = 2.0 * M_PI / 3.0; 

  pointNames[2] = new char[16];
  sprintf (pointNames[2], "M");
  pointCoordinates[2] = new double[2];
  pointCoordinates[2][0] = M_PI; 
  pointCoordinates[2][1] = 0.0; 

  return NbrHighSymmetryPoints;
}

// compute the distance between two points in the first Brillouin zone, changing the coordinates the second one by a reciprocal lattice vector if needed
//
// kx1 = momentum of the first point along the x axis
// ky1 = momentum of the first point along the y axis
// kx2 = reference on the momentum of the second point along the x axis
// ky2 = reference on the momentum of the second point along the y axis
// return value = distance between the two points

double TightBindingModelTimeReversalKagomeLattice::GetDistanceReciprocalSpace(double kx1, double ky1, double& kx2, double& ky2)
{
  double AngleFactor = 2.0 * cos(2.0 * M_PI / 3.0);
  double DiffKx = kx1 - kx2;
  double DiffKy = ky1 - ky2;
  double MinDistance = sqrt ((DiffKx * DiffKx) + (DiffKy * DiffKy) + (AngleFactor * DiffKx * DiffKy));
  double MinKx2 = kx2;
  double MinKy2 = ky2;
  for (int i = -1; i <= 1; ++i)
    {
      double TmpKx2  = kx2 + (2.0 * ((double) i) * M_PI);
      for (int j = -1; j <= 1; ++j)
	{
	  double TmpKy2  = ky2 + (2.0 * ((double) j) * M_PI);	  
	  double DiffKx = kx1 - TmpKx2;
	  double DiffKy = ky1 - TmpKy2;
	  double TmpDistance =  sqrt ((DiffKx * DiffKx) + (DiffKy * DiffKy) + (AngleFactor * DiffKx * DiffKy));
	  if (TmpDistance < MinDistance)
	    {
	      MinDistance = TmpDistance;
	      MinKx2 = TmpKx2;
	      MinKy2 = TmpKy2;
	    }
	}
    }
  kx2 = MinKx2;
  ky2 = MinKy2;
  return MinDistance;
}
