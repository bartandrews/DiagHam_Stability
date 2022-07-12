////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//         class of SU(4) generalized Halperine wave function on sphere       //
//                                                                            //
//                        last modification : 12/11/2006                      //
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
#include "Tools/FQHEWaveFunction/SU4HalperinOnSphereWaveFunction.h"
#include "GeneralTools/ConfigurationParser.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <cstring>

using std::cout;
using std::endl;


// constructor
//
// nbrSpinUpIsospinPlusParticles = number of particles with spin up and isopsin plus
// nbrSpinUpIsospinMinusParticles = number of particles with spin up and isopsin minus
// nbrSpinDownIsospinPlusParticles = number of particles with spin down and isopsin plus
// nbrSpinDownIsospinMinusParticles = number of particles with spin down and isopsin minus
// mupIndex = power of the laughlin-like part for spin up - isopsin plus
// mumIndex = power of the laughlin-like part for spin up - isopsin minus
// mdpIndex = power of the laughlin-like part for spin down - isopsin plus
// mdmIndex = power of the laughlin-like part for spin down - isopsin minus
// nIntraIsospinIndex = power of the intra-isospin part (i.e (z_u{pm} -z_d{pm}))
// nIntraSpinIndex = power of the intra-spin part (i.e (z_{ud}p -z_{ud}m))
// nCrossSpinIsospinIndex = power of the cross spin-isospin part (i.e (z_up -z_dm) and (z_um -z_dp))

SU4HalperinOnSphereWaveFunction::SU4HalperinOnSphereWaveFunction(int nbrSpinUpIsospinPlusParticles, int nbrSpinUpIsospinMinusParticles, 
								 int nbrSpinDownIsospinPlusParticles, int nbrSpinDownIsospinMinusParticles,
								 int mupIndex, int mumIndex, int mdpIndex, int mdmIndex,
								 int nIntraIsospinIndex, int nIntraSpinIndex, int nCrossSpinIsospinIndex)
{
  this->NbrSpinUpIsospinPlusParticles = nbrSpinUpIsospinPlusParticles;
  this->NbrSpinUpIsospinMinusParticles = nbrSpinUpIsospinMinusParticles;
  this->NbrSpinDownIsospinPlusParticles = nbrSpinDownIsospinPlusParticles;
  this->NbrSpinDownIsospinMinusParticles = nbrSpinDownIsospinMinusParticles;
  this->MupIndex = mupIndex;
  this->MumIndex = mumIndex;
  this->MdpIndex = mdpIndex;
  this->MdmIndex = mdmIndex;
  this->NIntraIsospinIndex = nIntraIsospinIndex;
  this->NIntraSpinIndex = nIntraSpinIndex;
  this->NCrossSpinIsospinIndex = nCrossSpinIsospinIndex;
  this->TotalNbrParticles = (this->NbrSpinUpIsospinPlusParticles + this->NbrSpinUpIsospinMinusParticles + 
			     this->NbrSpinDownIsospinPlusParticles + this->NbrSpinDownIsospinMinusParticles);
  this->NbrSpinUpParticles = this->NbrSpinUpIsospinPlusParticles +  this->NbrSpinUpIsospinMinusParticles;
  this->PartialNbrParticles = this->TotalNbrParticles - this->NbrSpinDownIsospinMinusParticles;
  cout << this->NbrSpinUpIsospinPlusParticles << " "  << this->NbrSpinUpIsospinMinusParticles << " "
       <<  this->NbrSpinDownIsospinPlusParticles << " " << this->NbrSpinDownIsospinMinusParticles << endl;
  cout <<  this->MupIndex << " " << this->MumIndex << " " << this->MdpIndex << " " << this->MdmIndex << " " << this->NIntraIsospinIndex 
       << " " <<  this->NIntraSpinIndex<< " " << this->NCrossSpinIsospinIndex << endl;
}

// copy constructor
//
// function = reference on the wave function to copy

SU4HalperinOnSphereWaveFunction::SU4HalperinOnSphereWaveFunction(const SU4HalperinOnSphereWaveFunction& function)
{
  this->NbrSpinUpIsospinPlusParticles = function.NbrSpinUpIsospinPlusParticles;
  this->NbrSpinUpIsospinMinusParticles = function.NbrSpinUpIsospinMinusParticles;
  this->NbrSpinDownIsospinPlusParticles = function.NbrSpinDownIsospinPlusParticles;
  this->NbrSpinDownIsospinMinusParticles = function.NbrSpinDownIsospinMinusParticles;
  this->MupIndex = function.MupIndex;
  this->MumIndex = function.MumIndex;
  this->MdpIndex = function.MdpIndex;
  this->MdmIndex = function.MdmIndex;
  this->NIntraIsospinIndex = function.NIntraIsospinIndex;
  this->NIntraSpinIndex = function.NIntraSpinIndex;
  this->NCrossSpinIsospinIndex = function.NCrossSpinIsospinIndex;
  this->TotalNbrParticles = function.TotalNbrParticles;
  this->NbrSpinUpParticles = this->NbrSpinUpIsospinPlusParticles +  this->NbrSpinUpIsospinMinusParticles;
  this->PartialNbrParticles = this->TotalNbrParticles - this->NbrSpinDownIsospinMinusParticles;
}

// constructor from configuration file
//
// configuration = reference on the configuration file parser
// errorFlag = reference on the error flag that is set to false if an error occured while retriving the configuration
// nbrParticles = reference on the total number of particles (computed from the configuration file datas)
// lzMax = twice the maximum angular momentum for a particle (computed from the configuration file datas)

SU4HalperinOnSphereWaveFunction::SU4HalperinOnSphereWaveFunction(ConfigurationParser& configuration, bool& errorFlag, int& nbrParticles, int& lzMax)
{
  errorFlag = true;
  if ((configuration["WaveFunction"] == 0) || (strcmp ("SU4Halperin", configuration["WaveFunction"]) != 0))
    {
      errorFlag = false;
      return;
    }
  if ((configuration.GetAsSingleInteger("NbrSpinUpIsospinPlusParticles", this->NbrSpinUpIsospinPlusParticles) == false) ||
      (configuration.GetAsSingleInteger("NbrSpinUpIsospinMinusParticles", this->NbrSpinUpIsospinMinusParticles) == false) || 
      (configuration.GetAsSingleInteger("NbrSpinDownIsospinPlusParticles", this->NbrSpinDownIsospinPlusParticles) == false) ||
      (configuration.GetAsSingleInteger("NbrSpinDownIsospinMinusParticles", this->NbrSpinDownIsospinMinusParticles) == false))
    {
      errorFlag = false;
      return;
    }
  cout << this->NbrSpinUpIsospinPlusParticles << " "  << this->NbrSpinUpIsospinMinusParticles << " "
       <<  this->NbrSpinDownIsospinPlusParticles << " " << this->NbrSpinDownIsospinMinusParticles << endl;
  this->TotalNbrParticles = (this->NbrSpinUpIsospinPlusParticles + this->NbrSpinUpIsospinMinusParticles + 
			     this->NbrSpinDownIsospinPlusParticles + this->NbrSpinDownIsospinMinusParticles);
  this->NbrSpinUpParticles = this->NbrSpinUpIsospinPlusParticles +  this->NbrSpinUpIsospinMinusParticles;
  this->PartialNbrParticles = this->TotalNbrParticles - this->NbrSpinDownIsospinMinusParticles;
  nbrParticles = this->TotalNbrParticles;
  if ((this->TotalNbrParticles <= 0) || (this->NbrSpinUpIsospinPlusParticles < 0) || (this->NbrSpinUpIsospinMinusParticles < 0) ||
      (this->NbrSpinDownIsospinPlusParticles < 0) || (this->NbrSpinDownIsospinMinusParticles < 0))
    {
      errorFlag = false;
      return;
    }
  if ((configuration.GetAsSingleInteger("MUpPlus", this->MupIndex) == false) ||
      (configuration.GetAsSingleInteger("MUpMinus", this->MumIndex) == false) ||
      (configuration.GetAsSingleInteger("MDownPlus", this->MdpIndex) == false) ||
      (configuration.GetAsSingleInteger("MDownMinus", this->MdmIndex) == false) ||
      (configuration.GetAsSingleInteger("NIntraIsospin", this->NIntraIsospinIndex) == false) ||
      (configuration.GetAsSingleInteger("NIntraSpin", this->NIntraSpinIndex) == false) ||
      (configuration.GetAsSingleInteger("NCrossSpinIsospin", this->NCrossSpinIsospinIndex) == false) ||
      (this->MupIndex < 0) || (this->MdpIndex < 0) || (this->MdmIndex < 0) || (this->MumIndex < 0) || 
      (this->NIntraIsospinIndex < 0) || (this->NIntraSpinIndex < 0) || (this->NCrossSpinIsospinIndex < 0))
    {
      errorFlag = false;
      return;
    }
  cout <<  this->MupIndex << " " << this->MumIndex << " " << this->MdpIndex << " " << this->MdmIndex << " " << this->NIntraIsospinIndex 
       << " " <<  this->NIntraSpinIndex<< " " << this->NCrossSpinIsospinIndex << endl;
  lzMax = ((this->MupIndex * this->NbrSpinUpIsospinPlusParticles) + (this->NIntraIsospinIndex * this->NbrSpinDownIsospinPlusParticles) 
	   + (this->NIntraSpinIndex * this->NbrSpinUpIsospinMinusParticles) + (this->NCrossSpinIsospinIndex * this->NbrSpinDownIsospinMinusParticles));


  cout <<lzMax  << " " << ((this->NIntraIsospinIndex * this->NbrSpinUpIsospinPlusParticles) + (this->MumIndex * this->NbrSpinDownIsospinPlusParticles) 
			   + (this->NCrossSpinIsospinIndex * this->NbrSpinUpIsospinMinusParticles) + (this->NIntraSpinIndex * this->NbrSpinDownIsospinMinusParticles))
       << " " << ((this->NIntraSpinIndex * this->NbrSpinUpIsospinPlusParticles) + (this->NCrossSpinIsospinIndex * this->NbrSpinDownIsospinPlusParticles) 
		 + (this->MumIndex * this->NbrSpinUpIsospinMinusParticles) + (this->NIntraIsospinIndex * this->NbrSpinDownIsospinMinusParticles))
       << " " << ((this->NCrossSpinIsospinIndex * this->NbrSpinUpIsospinPlusParticles) + (this->NIntraSpinIndex * this->NbrSpinDownIsospinPlusParticles) 
		 + (this->NIntraIsospinIndex * this->NbrSpinUpIsospinMinusParticles) + (this->MdmIndex * this->NbrSpinDownIsospinMinusParticles)) << " " << endl;

  if ((lzMax != ((this->NIntraIsospinIndex * this->NbrSpinUpIsospinPlusParticles) + (this->MumIndex * this->NbrSpinDownIsospinPlusParticles) 
		 + (this->NCrossSpinIsospinIndex * this->NbrSpinUpIsospinMinusParticles) + (this->NIntraSpinIndex * this->NbrSpinDownIsospinMinusParticles))) ||
      (lzMax != ((this->NIntraSpinIndex * this->NbrSpinUpIsospinPlusParticles) + (this->NCrossSpinIsospinIndex * this->NbrSpinDownIsospinPlusParticles) 
		 + (this->MumIndex * this->NbrSpinUpIsospinMinusParticles) + (this->NIntraIsospinIndex * this->NbrSpinDownIsospinMinusParticles))) ||
      (lzMax != ((this->NCrossSpinIsospinIndex * this->NbrSpinUpIsospinPlusParticles) + (this->NIntraSpinIndex * this->NbrSpinDownIsospinPlusParticles) 
		 + (this->NIntraIsospinIndex * this->NbrSpinUpIsospinMinusParticles) + (this->MdmIndex * this->NbrSpinDownIsospinMinusParticles))))
    {
      errorFlag = false;
      return;
    }
}

// destructor
//

SU4HalperinOnSphereWaveFunction::~SU4HalperinOnSphereWaveFunction()
{
}

// clone function 
//
// return value = clone of the function 

Abstract1DComplexFunction* SU4HalperinOnSphereWaveFunction::Clone ()
{
  return new SU4HalperinOnSphereWaveFunction(*this);
}

// evaluate function at a given point (the first 2*nbrSpinUpIsospinPlusParticles coordinates correspond to the position of the spin up - isposin plus particles, 
//                                     the following 2*nbrSpinDownIsospinPlusParticles coordinates correspond to the position of spin down - isposin plus particles,
//                                     the other coordinates obey to the same scheme for isospin minus)
//
//
// uv = ensemble of spinor variables on sphere describing point
//      where function has to be evaluated
//      ordering: u[i] = uv [2*i], v[i] = uv [2*i+1]
// return value = function value at (uv)

Complex SU4HalperinOnSphereWaveFunction::CalculateFromSpinorVariables(ComplexVector& uv)
{
  Complex Tmp;
  Complex TmpU;
  Complex TmpV;
  double Factor = M_PI * 0.5;
  Complex TotalWaveFunction(1.0);
  if (this->MupIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpIsospinPlusParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrSpinUpIsospinPlusParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->MupIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->MumIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = this->NbrSpinUpIsospinPlusParticles; i < this->NbrSpinUpParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->NbrSpinUpParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->MumIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->MdpIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = this->NbrSpinUpParticles; i < this->PartialNbrParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->PartialNbrParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->MdpIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->MdmIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = this->PartialNbrParticles; i < this->TotalNbrParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = i + 1; j < this->TotalNbrParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->MdmIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }


  if (this->NIntraIsospinIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpIsospinPlusParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = this->NbrSpinUpParticles; j < this->PartialNbrParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = this->NbrSpinUpIsospinPlusParticles; i < this->NbrSpinUpParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = this->PartialNbrParticles; j < this->TotalNbrParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->NIntraIsospinIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->NIntraSpinIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpIsospinPlusParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = this->NbrSpinUpIsospinPlusParticles; j < this->NbrSpinUpParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = this->NbrSpinUpParticles; i < this->PartialNbrParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = this->PartialNbrParticles; j < this->TotalNbrParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->NIntraSpinIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }
  if (this->NCrossSpinIsospinIndex > 0)
    {
      Complex WaveFunction(1.0);
      for (int i = 0; i < this->NbrSpinUpIsospinPlusParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = this->PartialNbrParticles; j < this->TotalNbrParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = this->NbrSpinUpIsospinPlusParticles; i < this->NbrSpinUpParticles; ++i)
	{
	  TmpU = uv[2 * i];
	  TmpV = uv[2 * i + 1];
	  for (int j = this->NbrSpinUpParticles; j < this->PartialNbrParticles; ++j)
	    WaveFunction *=  Factor * ((TmpU * uv[2 * j + 1]) - (TmpV * uv[2 * j]));
	}
      for (int i = 0; i < this->NCrossSpinIsospinIndex; ++i)
	TotalWaveFunction *= WaveFunction;
    }

  return TotalWaveFunction;
}
