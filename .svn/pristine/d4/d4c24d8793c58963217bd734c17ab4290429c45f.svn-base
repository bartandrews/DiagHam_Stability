////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                              class of XYZ chain                            //
//                                                                            //
//                        last modification : 16/12/2013                      //
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


#include "Hamiltonian/SpinChainLongRangeXYZHamiltonian.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/Complex.h"
#include "Output/MathematicaOutput.h"

#include <iostream>


using std::cout;
using std::endl;
using std::ostream;


// default constructor
//

SpinChainLongRangeXYZHamiltonian::SpinChainLongRangeXYZHamiltonian()
{
}

// constructor from default data
//
// chain = reference on Hilbert space of the associated system
// nbrSpin = number of spin
// jxFactor = coupling along the x direction
// jyFactor = coupling along the y direction
// jzFactor = coupling along the z direction
// powerLawXX/YY/ZZ = power law decay for XX, YY, ZZ couplings
// hFactor = Zeeman term 

SpinChainLongRangeXYZHamiltonian::SpinChainLongRangeXYZHamiltonian(Spin1_2ChainWithTranslations* chain, int nbrSpin, double jxFactor, double jyFactor, double jzFactor, double powerLawXX, double powerLawYY, double powerLawZZ, double hFactor)
{
  this->Chain = chain;
  this->NbrSpin = nbrSpin;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->JxFactor = jxFactor;
  this->JyFactor = jyFactor;
  this->JzFactor = jzFactor;
  this->PowerLawXX = powerLawXX;
  this->PowerLawYY = powerLawYY;
  this->PowerLawZZ = powerLawZZ;
  this->FFactors = new double[this->NbrSpin];
  for (int i = 0; i < this->NbrSpin; ++i)
    this->FFactors[i] = hFactor;
  this->EvaluateDiagonalMatrixElements();
  this->EvaluateCosinusTable();
}

// destructor
//

SpinChainLongRangeXYZHamiltonian::~SpinChainLongRangeXYZHamiltonian() 
{
}

// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void SpinChainLongRangeXYZHamiltonian::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
  delete[] this->SzSzContributions;
  this->Chain = (Spin1_2ChainWithTranslations*) hilbertSpace;
  this->SzSzContributions = new double [this->Chain->GetHilbertSpaceDimension()];
  this->EvaluateDiagonalMatrixElements();
}

// get Hilbert space on which Hamiltonian acts
//
// return value = pointer to used Hilbert space

AbstractHilbertSpace* SpinChainLongRangeXYZHamiltonian::GetHilbertSpace ()
{
  return this->Chain;
}

// return dimension of Hilbert space where Hamiltonian acts
//
// return value = corresponding matrix elementdimension

int SpinChainLongRangeXYZHamiltonian::GetHilbertSpaceDimension ()
{
  return this->Chain->GetHilbertSpaceDimension();
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void SpinChainLongRangeXYZHamiltonian::ShiftHamiltonian (double shift)
{
  for (int i = 0; i < this->Chain->GetHilbertSpaceDimension(); i ++)
    this->SzSzContributions[i] += shift;
}


// evaluate all cosinus/sinus that are needed when computing matrix elements
//

void SpinChainLongRangeXYZHamiltonian::EvaluateCosinusTable()
{
  this->CosinusTable = new double [this->NbrSpin];
  this->SinusTable = new double [this->NbrSpin];
  this->ExponentialTable = new Complex [this->NbrSpin];
  double Coef = 2.0 * M_PI / ((double) this->NbrSpin) * ((double) this->Chain->GetMomentum());
  for (int i = 0; i < this->NbrSpin ; ++i)
    {
      this->CosinusTable[i] = cos(Coef * ((double) i));
      this->SinusTable[i] = sin(Coef * ((double) i));
      this->ExponentialTable[i].Re = this->CosinusTable[i];
      this->ExponentialTable[i].Im = this->SinusTable[i];
    }
}

  
// multiply a vector by the current hamiltonian for a given range of indices 
// and add result to another vector, low level function (no architecture optimization)
//
// vSource = vector to be multiplied
// vDestination = vector at which result has to be added
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = reference on vector where result has been stored

ComplexVector& SpinChainLongRangeXYZHamiltonian::LowLevelAddMultiply(ComplexVector& vSource, ComplexVector& vDestination, 
							  int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int NbrTranslation;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;
  Complex Magnitude;
  double Jxx, Jyy;
  for (int i = firstComponent; i < LastComponent; ++i)
    {
      vDestination[i] += this->SzSzContributions[i] * vSource[i];

      // J part of Hamiltonian      
      for (int j1 = 0; j1 <= MaxPos; ++j1)
      	for (int j2 = 0; j2 < j1; ++j2)
	      {
           	pos = this->Chain->SmiSpj(j1, j2, i, coef, NbrTranslation);
	  		if (pos < dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = vSource[i] * (Jxx + Jyy);

	      			vDestination[pos] += Magnitude * coef * this->ExponentialTable[NbrTranslation];
	    		}	
	  		pos = this->Chain->SmiSpj(j2, j1, i, coef, NbrTranslation);
	  		if (pos < dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = vSource[i] * (Jxx + Jyy);

	      			vDestination[pos] += Magnitude * coef * this->ExponentialTable[NbrTranslation];
	    		}
	  		pos = this->Chain->SpiSpj(j1, j2, i, coef, NbrTranslation);
	  		if (pos < dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = vSource[i] * (Jxx - Jyy);

	      			vDestination[pos] += Magnitude * coef * this->ExponentialTable[NbrTranslation];
	    		}
	  		pos = this->Chain->SmiSmj(j1, j2, i, coef, NbrTranslation);
	  		if (pos < dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = vSource[i] * (Jxx - Jyy);

	      			vDestination[pos] += Magnitude * coef * this->ExponentialTable[NbrTranslation];
	    		}
		  }
      
    }
  return vDestination;
}

// multiply a set of vectors by the current hamiltonian for a given range of indices 
// and add result to another set of vectors, low level function (no architecture optimization)
//
// vSources = array of vectors to be multiplied
// vDestinations = array of vectors at which result has to be added
// nbrVectors = number of vectors that have to be evaluated together
// firstComponent = index of the first component to evaluate
// nbrComponent = number of components to evaluate
// return value = pointer to the array of vectors where result has been stored

ComplexVector* SpinChainLongRangeXYZHamiltonian::LowLevelMultipleAddMultiply(ComplexVector* vSources, ComplexVector* vDestinations, int nbrVectors, 
									   int firstComponent, int nbrComponent)
{
  int LastComponent = firstComponent + nbrComponent;
  int dim = this->Chain->GetHilbertSpaceDimension();
  double coef;
  double coef2;
  int NbrTranslation;
  Complex Magnitude;
  double Jxx, Jyy;
  int pos;
  int pos2;
  int MaxPos = this->NbrSpin - 1;

  for (int k = 0; k < nbrVectors; ++k)
    {
      ComplexVector& TmpSource = vSources[k];
      ComplexVector& TmpDestination = vDestinations[k];
      for (int i = firstComponent; i < LastComponent; ++i)
		{
	  		TmpDestination[i] += this->SzSzContributions[i] * TmpSource[i];
		}
    }

  for (int i = firstComponent; i < LastComponent; ++i)
    {
      // J part of Hamiltonian      
      for (int j1 = 0; j1 <= MaxPos; ++j1)
      	for (int j2 = 0; j2 < j1; ++j2)
		{
		  	pos = this->Chain->SmiSpj(j1, j2, i, coef, NbrTranslation);
	  		if (pos != dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = (Jxx + Jyy);

	      			for (int k = 0; k < nbrVectors; ++k)
		     			vDestinations[k][pos] += Magnitude * coef * vSources[k][i] * this->ExponentialTable[NbrTranslation];
	    		}		
	  		pos = this->Chain->SmiSpj(j2, j1, i, coef, NbrTranslation);
	  		if (pos != dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = (Jxx + Jyy);

	      			for (int k = 0; k < nbrVectors; ++k)
		    			vDestinations[k][pos] += Magnitude * coef * vSources[k][i] * this->ExponentialTable[NbrTranslation];
	    		}
	  		pos = this->Chain->SpiSpj(j1, j2, i, coef, NbrTranslation);
	  		if (pos != dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = (Jxx - Jyy);

	      			for (int k = 0; k < nbrVectors; ++k)
		    			vDestinations[k][pos] += Magnitude * coef * vSources[k][i] * this->ExponentialTable[NbrTranslation];
	    		}
	  		pos = this->Chain->SmiSmj(j1, j2, i, coef, NbrTranslation);
	  		if (pos != dim)
	    		{
	    			Jxx = this->JxFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawXX);
	    			Jyy = this->JyFactor/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawYY);
  					Magnitude = (Jxx - Jyy);

	      			for (int k = 0; k < nbrVectors; ++k)
		     			vDestinations[k][pos] += Magnitude * coef * vSources[k][i] * this->ExponentialTable[NbrTranslation];
	    		}
		}

    }
  return vDestinations;
}

// evaluate diagonal matrix elements
// 

void SpinChainLongRangeXYZHamiltonian::EvaluateDiagonalMatrixElements()
{
  int dim = this->Chain->GetHilbertSpaceDimension();
  double Magnitude;

  for (int i = 0; i < dim; i++)
    {
      // SzSz part
      double Tmp = 0.0;
      double Tmp2 = 0.0;
      
      for (int j = 0; j < this->NbrSpin; ++j)
	   {
	     Tmp2 += this->FFactors[j] * this->Chain->Szi(j, i);
	   }
	   
      for (int j1 = 0; j1 < this->NbrSpin; ++j1)
        for (int j2 = 0; j2 < j1; ++j2)
          {
            Magnitude = 1.0/pow((this->NbrSpin/M_PI) * sin((j1-j2) * M_PI/this->NbrSpin), this->PowerLawZZ);
	        Tmp += Magnitude * this->JzFactor * 4.0 * this->Chain->SziSzj(j1, j2, i);
	      }

      this->SzSzContributions[i] = Tmp + Tmp2;
    }
}

