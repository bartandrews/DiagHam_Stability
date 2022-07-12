////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//           class of tight binding model for the simple square lattice       //
//                                                                            //
//                        last modification : 01/01/2016                      //
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
#include "Tools/FTITightBinding/TightBindingModelSSH.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include <iostream>

using std::cout;
using std::endl;
using std::ostream;



// default constructor
//
// nbrSite = number of sites in the x direction
// t1 = hoping amplitude between neareast neighbor sites
// epsilon = square confinement amplitude 
// gammaX = boundary condition twisting angle along x
// architecture = pointer to the architecture
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored

TightBindingModelSSH::TightBindingModelSSH (int nbrSite, double delta, bool cylinder, AbstractArchitecture* architecture, bool storeOneBodyMatrices)
{
  this->NbrSiteX = nbrSite;
  this->KxFactor = 2.0 * M_PI / ((double) this->NbrSiteX);
  this->Delta = delta;
  this->CylinderFlag = cylinder;
  this->NbrBands = 2;
  this->NbrStatePerBand = this->NbrSiteX;
  this->Architecture = architecture;
}


// destructor
//

TightBindingModelSSH::~TightBindingModelSSH ()
{
}

// find the orbitals connected to those located at the origin unit cell
// 
  
HermitianMatrix TightBindingModelSSH::GetRealSpaceTightBindingHamiltonian()
{
  HermitianMatrix RealSpaceTightBindingHamiltonian (2*this->NbrSiteX,true);

  for(int i = 0; i <  this->NbrSiteX ; i++)
    {
      RealSpaceTightBindingHamiltonian.SetMatrixElement(2*i,2*i+1,(1.0-this->Delta));
    }
  
  for(int i = 0; i <  this->NbrSiteX-1 ; i++)
    {
      RealSpaceTightBindingHamiltonian.SetMatrixElement(2*i+1,2*(i+1),(1.0+this->Delta));
    }
  
  if(this->CylinderFlag == false) 
    {
      RealSpaceTightBindingHamiltonian.SetMatrixElement(2*this->NbrSiteX-1,0,(1.0+this->Delta));
    }

  return  RealSpaceTightBindingHamiltonian;
}

 
