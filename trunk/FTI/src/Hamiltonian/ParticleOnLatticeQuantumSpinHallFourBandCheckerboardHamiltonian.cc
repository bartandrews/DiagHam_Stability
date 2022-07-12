////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//               class of quatum spin Hall restricted to four bands           //
//                         using the checkerboard model                       //
//                                                                            //
//                        last modification : 25/12/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include <iostream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ostream;


// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = strength of the repulsive two body neareast neighbor interaction with identical spin
// vPotential = strength of the repulsive on site two body interaction with opposite spin
// wPotential = strength of the repulsive two body neareast neighbor interaction between opposite spins
// t1 = hoping amplitude between neareast neighbor sites
// t2 = hoping amplitude between next neareast neighbor sites
// t2p = hoping amplitude between second next neareast neighbor sites
// mixingTermNorm = norm of the mixing term coupling the two copies of the checkerboard lattice
// mixingTermArgv = argument of the mixing term coupling the two copies of the checkerboard lattice
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = maximum amount of memory that can be allocated for fast multiplication (negative if there is no limit)

ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian(ParticleOnSphereWithSU4Spin* particles, int nbrParticles, int nbrSiteX, 
															       int nbrSiteY, double uPotential, double vPotential, double wPotential, double t1, double t2, double t2p, double mixingTermNorm, double mixingTermArg, double gammaX, double gammaY, bool flatBandFlag, AbstractArchitecture* architecture, long memory)
{
  this->Particles = particles;
  this->NbrParticles = nbrParticles;
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->LzMax = nbrSiteX * nbrSiteY - 1;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->HamiltonianShift = 0.0;
  this->NNHoping = t1;
  this->NextNNHoping = t2;
  this->SecondNextNNHoping = t2p;
  this->MixingTerm = mixingTermNorm * Phase(mixingTermArg);
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->FlatBand = flatBandFlag;
  this->UPotential = uPotential;
  this->VPotential = vPotential;
  this->WPotential = wPotential;

  this->Architecture = architecture;
  this->Memory = memory;
  this->FastMultiplicationFlag = false;
  this->HermitianSymmetryFlag = false;
  long MinIndex;
  long MaxIndex;
  this->Architecture->GetTypicalRange(MinIndex, MaxIndex);
  this->PrecalculationShift = (int) MinIndex;  
  this->EvaluateInteractionFactors();
  if (memory > 0)
    {
      long TmpMemory = this->FastMultiplicationMemory(memory);
      if (TmpMemory < 1024)
	cout  << "fast = " <<  TmpMemory << "b ";
      else
	if (TmpMemory < (1 << 20))
	  cout  << "fast = " << (TmpMemory >> 10) << "kb ";
	else
	  if (TmpMemory < (1 << 30))
	    cout  << "fast = " << (TmpMemory >> 20) << "Mb ";
	  else
	    {
	      cout  << "fast = " << (TmpMemory >> 30) << ".";
	      TmpMemory -= ((TmpMemory >> 30) << 30);
	      TmpMemory *= 100l;
	      TmpMemory >>= 30;
	      if (TmpMemory < 10l)
		cout << "0";
	      cout  << TmpMemory << " Gb ";
	    }
      this->EnableFastMultiplication();
    }
}

// destructor
//

ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::~ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian()
{
}
  
// compute the matrix element for the two body interaction between two sites A and B  belonging to the same layer
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ComputeTwoBodyMatrixElementAUpBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) - (this->KyFactor * ((double) (ky4 - ky2)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) + (this->KyFactor * ((double) (ky4 - ky2)))))));
}

// compute the matrix element for the two body interaction between two sites A and B with down spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ComputeTwoBodyMatrixElementADownBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) - (this->KyFactor * ((double) (ky4 - ky2)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) + (this->KyFactor * ((double) (ky4 - ky2)))))));
}
  
// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin down
// ky1 = momentum along y for the creation operator on A site with spin down
// kx2 = momentum along x for the creation operator on B site with spin up
// ky2 = momentum along y for the creation operator on B site with spin up
// kx3 = momentum along x for the annihilation operator on A site with spin down
// ky3 = momentum along y for the annihilation operator on A site with spin down
// kx4 = momentum along x for the annihilation operator on B site with spin up
// ky4 = momentum along y for the annihilation operator on B site with spin up
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ComputeTwoBodyMatrixElementADownBUp(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx3 - kx1))) - (this->KyFactor * ((double) (ky3 - ky1)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx3 - kx1))) + (this->KyFactor * ((double) (ky3 - ky1)))))));
}

// compute the matrix element for the two body interaction between two sites A and B with opposite spins
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ComputeTwoBodyMatrixElementAUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return (2.0 * (cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) - (this->KyFactor * ((double) (ky4 - ky2)))))
		 + cos (0.5 * ((this->KxFactor * ((double) (kx4 - kx2))) + (this->KyFactor * ((double) (ky4 - ky2)))))));
}

// compute the matrix element for the two body interaction between two sites A with opposite spins 
//
// kx1 = momentum along x for the creation operator on A site with spin up
// ky1 = momentum along y for the creation operator on A site with spin up
// kx2 = momentum along x for the creation operator on A site with spin down
// ky2 = momentum along y for the creation operator on A site with spin down
// kx3 = momentum along x for the annihilation operator on A site with spin up
// ky3 = momentum along y for the annihilation operator on A site with spin up
// kx4 = momentum along x for the annihilation operator on A site with spin down
// ky4 = momentum along y for the annihilation operator on A site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ComputeTwoBodyMatrixElementAUpADown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  Complex Tmp = 1.0;
  return Tmp;
}

// compute the matrix element for the two body interaction between two sites B with opposite spins 
//
// kx1 = momentum along x for the creation operator on B site with spin up
// ky1 = momentum along y for the creation operator on B site with spin up
// kx2 = momentum along x for the creation operator on B site with spin down
// ky2 = momentum along y for the creation operator on B site with spin down
// kx3 = momentum along x for the annihilation operator on B site with spin up
// ky3 = momentum along y for the annihilation operator on B site with spin up
// kx4 = momentum along x for the annihilation operator on B site with spin down
// ky4 = momentum along y for the annihilation operator on B site with spin down
// return value = corresponding matrix element

Complex ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ComputeTwoBodyMatrixElementBUpBDown(int kx1, int ky1, int kx2, int ky2, int kx3, int ky3, int kx4, int ky4)
{
  return  Phase (0.5 * ((this->KxFactor * ((double) (kx4 + kx3 - kx2 - kx1))) +
			(this->KyFactor * ((double) (ky4 + ky3 - ky2 - ky1)))));
}


// compute the one body hamiltonians related to the band stucture contribution
//
// oneBodyHamiltonians = array of one body hamiltonians

void ParticleOnLatticeQuantumSpinHallFourBandCheckerboardHamiltonian::ComputeOneBodyHamiltonian(HermitianMatrix* oneBodyHamiltonians)
{
  for (int kx = 0; kx < this->NbrSiteX; ++kx)
    for (int ky = 0; ky < this->NbrSiteY; ++ky)
      {
	HermitianMatrix TmpOneBodyHamiltonian(4, true);
	int Index = ((kx * this->NbrSiteY) + ky);
	double TmpKx = (((double) kx) + this->GammaX) * this->KxFactor;
	double TmpKy = (((double) ky) + this->GammaY) * this->KyFactor;
	Complex B1 = 4.0 * this->NNHoping * Complex (cos (0.5 * TmpKx) * cos (0.5 * TmpKy) * cos(M_PI * 0.25), 
						     sin (0.5 * TmpKx) * sin (0.5 * TmpKy) * sin(M_PI * 0.25));
	double d1 = 4.0 * this->SecondNextNNHoping * cos (TmpKx) * cos (TmpKy);
	double d3 = 2.0 * this->NextNNHoping * (cos (TmpKx) - cos (TmpKy));
	TmpOneBodyHamiltonian.SetMatrixElement(0, 0, d1 + d3);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 1, B1);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 1, d1 - d3);
	TmpKx = (((double) -kx) + this->GammaX) * this->KxFactor;
	TmpKy = (((double) -ky) + this->GammaY) * this->KyFactor;
	B1 = 4.0 * this->NNHoping * Complex (cos (0.5 * TmpKx) * cos (0.5 * TmpKy) * cos(M_PI * 0.25), 
						     sin (0.5 * TmpKx) * sin (0.5 * TmpKy) * sin(M_PI * 0.25));
	d1 = 4.0 * this->SecondNextNNHoping * cos (TmpKx) * cos (TmpKy);
	d3 = 2.0 * this->NextNNHoping * (cos (TmpKx) - cos (TmpKy));
	TmpOneBodyHamiltonian.SetMatrixElement(2, 2, d1 + d3);
	TmpOneBodyHamiltonian.SetMatrixElement(2, 3, Conj(B1));
	TmpOneBodyHamiltonian.SetMatrixElement(3, 3, d1 - d3);
	TmpOneBodyHamiltonian.SetMatrixElement(0, 3, - I() * this->MixingTerm);
	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, I() * this->MixingTerm);
	//	TmpOneBodyHamiltonian.SetMatrixElement(1, 2, - I() * this->MixingTerm);
	//	TmpOneBodyHamiltonian.SetMatrixElement(0, 3, I() * this->MixingTerm);
	
	if (this->FlatBand == false)
	  {
	    oneBodyHamiltonians[Index] = TmpOneBodyHamiltonian;
	  }
	else
	  {
	    ComplexMatrix OneBodyBasis(4, 4);
	    OneBodyBasis.SetToIdentity();
	    oneBodyHamiltonians[Index].Copy(TmpOneBodyHamiltonian);
	    RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
	    TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, OneBodyBasis);
#else
	    TmpOneBodyHamiltonian.Diagonalize(TmpDiag, OneBodyBasis);
#endif   
	    for (int i = 0; i < 4; ++i)
	      {
		cout << TmpDiag(i, i) << " ";
	      }
	    cout << endl;
	    TmpOneBodyHamiltonian.ClearMatrix();
	    double Tmp = 1.0;
	    TmpOneBodyHamiltonian.SetMatrixElement(2, 2, Tmp);
	    TmpOneBodyHamiltonian.SetMatrixElement(3, 3, Tmp);
	    oneBodyHamiltonians[Index] = TmpOneBodyHamiltonian.InvConjugate(OneBodyBasis);
	  }
      }
}
