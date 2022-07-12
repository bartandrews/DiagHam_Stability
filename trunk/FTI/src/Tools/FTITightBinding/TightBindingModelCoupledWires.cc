////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                        Copyright (C) 2012-2014 Bin Xu                      //
//                                                                            //
//                                                                            //
//            class of tight binding model for the coupled wired model        //
//                         described in  arXiv:1403.1791                      //
//                                                                            //
//                        last modification : 22/07/2014                      //
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

#include "Tools/FTITightBinding/TightBindingModelCoupledWires.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"


// default constructor
//
// nbrSiteX = number of sites along the wire
// nbrSiteY = number of unit cell of wires, i.e., number of wires / 4
// tx = hopping amplitude along a wire
// phiX = complex phase angle of the hopping amplitude along a wire
// m = mass of electron
// t1 = electron-electron and hole-hole hopping amplitudes
// t2 = electron-hole hopping amplitudes
// architecture = pointer to the architecture
// gammaX = boundary condition twisting angle along x
// gammaY = boundary condition twisting angle along y
// storeOneBodyMatrices = flag flag to indicate if the
//   one body transformation matrices have to be computed and stored

TightBindingModelCoupledWires::TightBindingModelCoupledWires
                              (int nbrSiteX, int nbrSiteY, double tx,
                              double phiX, double t1, double t2, double m,
                              double gammaX, double gammaY,
                              AbstractArchitecture* architecture,
                              bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSiteX;
  this->NbrSiteY = nbrSiteY;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->KyFactor = 2.0 * M_PI / ((double) this->NbrSiteY);
  this->Hopping1 = t1;
  this->Hopping2 = t2;
  this->HoppingX = tx;
  this->PhiX = phiX;
  this->mass = m;
  this->GammaX = gammaX;
  this->GammaY = gammaY;
  this->NbrBands = 4;
  this->NbrStatePerBand = this->NbrSiteX * this->NbrSiteY;
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
  this->ComputeBandStructure();
}

// destructor
//

TightBindingModelCoupledWires::~TightBindingModelCoupledWires()
{
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelCoupledWires::CoreComputeBandStructure(long minStateIndex,
                                                             long nbrStates)
{
  if (nbrStates == 0)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  
  for (int kx = 0; kx < this->NbrSiteX; kx++) 
    {
      for (int ky = 0; ky < this->NbrSiteY; ky++) 
	{
	  int Index = this->GetLinearizedMomentumIndex(kx, ky);
	  if ((Index >= minStateIndex) && (Index <= MaxStateIndex)) 
	    {
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      
	      //Compute a few matrix elements, ref. arXiv: 1403.1791 eq. (8)
	      double txm = 2.0 * this->HoppingX * cos((2.0 * M_PI * double(kx)/((double) this->NbrSiteX) + this->PhiX/2.0 ));
	      double txp = 2.0 * this->HoppingX * cos((2.0 * M_PI * double(kx)/((double) this->NbrSiteX) - this->PhiX/2.0 ));
	      Complex Hopping2C = this->Hopping2 * Complex(cos(4.0 * 2.0 * M_PI * double(ky)/((double)this->NbrSiteY)),
							   sin(4.0 * 2.0 * M_PI * double(ky)/((double)this->NbrSiteY)));
	      Complex Hopping2Ccc = this->Hopping2 * Complex(cos(4.0 * 2.0 * M_PI * double(ky)/((double)this->NbrSiteY)),
							   -sin(4.0 * 2.0 * M_PI * double(ky)/((double)this->NbrSiteY)));
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 0, this->mass - txm);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 1, this->mass - txp);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 2, -this->mass + txp);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 3, -this->mass + txm);//a few of them
        
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 0, -this->Hopping1);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 1, -this->Hopping2);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 2, -this->Hopping1);//a few of them
	      
              TmpOneBodyHamiltonian.SetMatrixElement(0, 1, -this->Hopping1);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(1, 2, -this->Hopping2);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(2, 3, -this->Hopping1);//a few of them
	      
	      TmpOneBodyHamiltonian.SetMatrixElement(0, 3, -Hopping2C);//a few of them
	      TmpOneBodyHamiltonian.SetMatrixElement(3, 0, -Hopping2Ccc);//a few of them
	      
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
		  for (int i = 0; i < this->NbrBands; i++) 
		    {
		      this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
                    }
                }
            }
        }
    }
}
