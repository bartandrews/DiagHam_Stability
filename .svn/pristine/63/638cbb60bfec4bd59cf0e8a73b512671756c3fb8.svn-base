////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                                                                            //
//                 class of quantum spin Hall restricted to two bands         //
//                  with a fully SU(2) symmetry breaking interaction          //
//                            and real matrix elements                        //
//                                                                            //
//                        last modification : 04/05/2020                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullTwoBandRealHamiltonian.h"
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

ParticleOnLatticeQuantumSpinHallFullTwoBandRealHamiltonian::ParticleOnLatticeQuantumSpinHallFullTwoBandRealHamiltonian()
{
  this->OneBodyInteractionFactorsSigma = 0;
  this->NbrInternalIndices = 2;
}

// destructor
//

ParticleOnLatticeQuantumSpinHallFullTwoBandRealHamiltonian::~ParticleOnLatticeQuantumSpinHallFullTwoBandRealHamiltonian()
{
  if (this->InteractionFactorsSigma != 0)
    {
       for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)
	    {
	      for (int sigma3 = 0; sigma3 < this->NbrInternalIndices; ++sigma3)
		{
		  if (this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma3] != 0)
		    {
		      for (int i = 0; i < this->NbrIntraSectorSums; ++i)
			{
			  delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma3][i];
			}
		      delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma3] ;
		    }
		  for (int sigma4 = sigma3 + 1; sigma4 < this->NbrInternalIndices; ++sigma4)
		    {
		      if (this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma4] != 0)
			{
			  for (int i = 0; i < this->NbrInterSectorSums; ++i)
			    {
			      delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma4][i];
			    }
			  delete[] this->InteractionFactorsSigma[sigma1][sigma2][sigma3][sigma4];
			}
		    }
		}
	    }
	}
    }
  if (this->OneBodyInteractionFactorsSigma != 0)
    {
       for (int sigma1 = 0; sigma1 < this->NbrInternalIndices; ++sigma1)
	{
	  for (int sigma2 = sigma1; sigma2 < this->NbrInternalIndices; ++sigma2)
	    {
	      if (this->OneBodyInteractionFactorsSigma[sigma1][sigma2] != 0)
		{
		  delete[] this->OneBodyInteractionFactorsSigma[sigma1][sigma2];
		}	      
	    }
	  delete[] this->OneBodyInteractionFactorsSigma[sigma1];
	}
       delete[] this->OneBodyInteractionFactorsSigma;
    }
}



