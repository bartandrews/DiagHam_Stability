////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                 class of quantum spin Hall restricted to four bands        //
//                                                                            //
//                        last modification : 26/08/2011                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian.h"
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


// default constructor
//

ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian()
{
}

// destructor
//

ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian::~ParticleOnLatticeQuantumSpinHallFullFourBandHamiltonian()
{
  if (this->InteractionFactorsSigma != 0)
    {
       for (int sigma1 = 0; sigma1 < 4; ++sigma1)
	{
	  for (int sigma2 = sigma1; sigma2 < 4; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < 4; ++sigma3)
		{
		  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
		    {
		      delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma3][i];
		    }
		  delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma3] ;
		  for (int sigma4 = sigma3 + 1; sigma4 < 4; ++sigma4)
		    {			  
		      for (int i = 0; i < this->NbrInterSectorSums; ++i)
			{
			  delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma4][i];
			}
		      delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma4] ;
		    }
		}
	    }
	}
    }
}



