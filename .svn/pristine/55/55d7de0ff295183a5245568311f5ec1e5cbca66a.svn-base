////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                        class author: Nicolas Regnault                      //
//                                                                            //
//                class of quantum spin Hall restricted to three bands        //
//                  with a fully SU(3) symmetry breaking interaction          //
//                                                                            //
//                        last modification : 30/06/2012                      //
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
#include "Hamiltonian/ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian.h"
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

ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian()
{
}

// destructor
//

ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian::~ParticleOnLatticeQuantumSpinHallFullThreeBandHamiltonian()
{
  if (this->InteractionFactors1122 != 0)
    {
      if (this->NbrIntraSectorSums > 0)
	{
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      delete[] this->InteractionFactors1122[i];
	      delete[] this->InteractionFactors1133[i];
	      delete[] this->InteractionFactors2211[i];
	      delete[] this->InteractionFactors2233[i];
	      delete[] this->InteractionFactors3311[i];
	      delete[] this->InteractionFactors3322[i];
	    }
	  delete[] this->InteractionFactors1122;
	  delete[] this->InteractionFactors1133;
	  delete[] this->InteractionFactors2211;
	  delete[] this->InteractionFactors2233;
	  delete[] this->InteractionFactors3311;
	  delete[] this->InteractionFactors3322;
	}
    }
  if (this->InteractionFactors1213 != 0)
    {
      if (this->NbrInterSectorSums > 0)
	{
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      delete[] InteractionFactors1213[i];
	      delete[] InteractionFactors1223[i];
	      delete[] InteractionFactors1312[i];
	      delete[] InteractionFactors1323[i];
	      delete[] InteractionFactors2312[i];
	      delete[] InteractionFactors2313[i];
	    }
	  delete[] this->InteractionFactors1213;
	  delete[] this->InteractionFactors1223;
	  delete[] this->InteractionFactors1312;
	  delete[] this->InteractionFactors1323;
	  delete[] this->InteractionFactors2312;
	  delete[] this->InteractionFactors2313;
	}
    }
  if (this->InteractionFactors1112 != 0)
    {
      if (this->NbrIntraSectorSums > 0)
	{
	  for (int i = 0; i < this->NbrInterSectorSums; ++i)
	    {
	      delete[] InteractionFactors1112[i];
	      delete[] InteractionFactors1113[i];
	      delete[] InteractionFactors1123[i];
	      delete[] InteractionFactors2212[i];
	      delete[] InteractionFactors2213[i];
	      delete[] InteractionFactors2223[i];
	      delete[] InteractionFactors3312[i];
	      delete[] InteractionFactors3313[i];
	      delete[] InteractionFactors3323[i];
	    }
	  delete[] this->InteractionFactors1112;
	  delete[] this->InteractionFactors1113;
	  delete[] this->InteractionFactors1123;
	  delete[] this->InteractionFactors2212;
	  delete[] this->InteractionFactors2213;
	  delete[] this->InteractionFactors2223;
	  delete[] this->InteractionFactors3312;
	  delete[] this->InteractionFactors3313;
	  delete[] this->InteractionFactors3323;
	}
    }
  if (this->InteractionFactors1211 != 0)
    {
      if (this->NbrInterSectorSums > 0)
	{
	  for (int i = 0; i < this->NbrIntraSectorSums; ++i)
	    {
	      delete[] InteractionFactors1211[i];
	      delete[] InteractionFactors1311[i];
	      delete[] InteractionFactors2311[i];
	      delete[] InteractionFactors1222[i];
	      delete[] InteractionFactors1322[i];
	      delete[] InteractionFactors2322[i];
	      delete[] InteractionFactors1233[i];
	      delete[] InteractionFactors1333[i];
	      delete[] InteractionFactors2333[i];
	    }
	  delete[] this->InteractionFactors1211;
	  delete[] this->InteractionFactors1311;
	  delete[] this->InteractionFactors2311;
	  delete[] this->InteractionFactors1222;
	  delete[] this->InteractionFactors1322;
	  delete[] this->InteractionFactors2322;
	  delete[] this->InteractionFactors1233;
	  delete[] this->InteractionFactors1333;
	  delete[] this->InteractionFactors2333;
	}
    }
}

