////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2005 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//        class of hamiltonian associated quantum dots in 3 dimensions        //
//                                                                            //
//                      last modification : 10/13/2005                        //
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
#include "Hamiltonian/QuantumWellHamiltonianInMagneticField.h"
 
#include <math.h>
#include <iostream>
#include <stdlib.h>


using std::cout;
using std::endl;


#define HBARE_M0 0.115767635
#define HBAR 1.05457168e-34
#define HBAR_E 6.582119138e-16
#define ECHARGE 1.60217653e-19
#define M0 9.1093826e-31


// constructor from default data
//
// xSize = system dimension in the x direction (in Angstrom unit)
// ySize = system dimension in the y direction (in Angstrom unit)
// zSize = system dimension in the z direction (in Angstrom unit)
// mass = effective mass in the x direction (in electron mass unit)
// bField = B field value (in Tesla)
// zEnergy1 = z confinement in the first subband
// zEnergy2 = z confinement in the second subband
// landauIndex1 = Landau index of the first subband
// landauIndex2 = Landau index of the second subband
// mailleParameter =
// bandOffset = conduction band offset between GaAs and InAs
// inDopage = In/Ga dopage ratio (=x with Ga_(1-x) In_x As)
// potentialDescription = name of the file that contains the potential description (null if the potential has to be evaluated)

QuantumWellHamiltonianInMagneticField::QuantumWellHamiltonianInMagneticField(double xSize, double ySize, double zSize, double mass, double bField, double zEnergy1, double zEnergy2,
									     int landauIndex1, int landauIndex2, double mailleParameter, double bandOffset, double inDopage, 
									     char* potentialDescription)
{
  this->XSize = xSize;
  this->YSize = ySize;
  this->ZSize = zSize;
  this->Mass = mass;
  this->BField = bField;
  this->ZEnergy1 = zEnergy1;
  this->ZEnergy2 = zEnergy2;
  this->LandauIndex1 = landauIndex1;
  this->LandauIndex2 = landauIndex2;
  this->MailleParameter = mailleParameter;
  this->BandOffset = bandOffset;
  this->InDopage = inDopage;

  this->MagneticLength = 1.0e10 * sqrt(HBAR_E / this->BField);
  cout << this->MagneticLength << endl;
  this->CyclotronEnergy = HBARE_M0 * this->BField / this->Mass;
  cout << this->CyclotronEnergy << endl;
  this->LandauDegeneracy = (int) ((this->XSize * this->YSize) / (2.0 * M_PI * this->MagneticLength * this->MagneticLength));
  cout << this->LandauDegeneracy << endl;
  this->NbrCells = (int) ((4.0 * this->XSize * this->YSize * this->ZSize) / (this->MailleParameter * this->MailleParameter * this->MailleParameter));
  cout << this->NbrCells << endl;
  this->GaXPosition = new double [this->NbrCells];
  this->GaYPosition = new double [this->NbrCells];
  this->GaZPosition = new double [this->NbrCells];
  this->InXPosition = new double [this->NbrCells];
  this->InYPosition = new double [this->NbrCells];
  this->InZPosition = new double [this->NbrCells];
  HermitianMatrix TmpHamiltonian(this->LandauDegeneracy * 2, true);
  this->Hamiltonian = TmpHamiltonian;
  this->NbrXCells = (int)  (2.0 * this->XSize / this->MailleParameter);
  this->NbrYCells = (int)  (2.0 * this->YSize / this->MailleParameter);
  this->NbrZCells = (int)  (this->ZSize / this->MailleParameter);
  this->Potential = new BinaryThreeDConstantCellPotential(NbrXCells, NbrYCells, NbrZCells);
  if (potentialDescription == 0)
    {
      int Threshold = (int) (this->InDopage * ((double) RAND_MAX));
      double GaCoefficient = this->MailleParameter * this->MailleParameter * this->MailleParameter * 0.25 * (2.0 / (this->YSize * this->ZSize)) * this->BandOffset;
      double InCoefficient = GaCoefficient * (this->InDopage - 1.0);
      GaCoefficient *= this->InDopage;
      for (int k = 0; k < this->NbrZCells; ++k)
	for (int j = 0; j < this->NbrYCells; ++j)
	  for (int i = 0; i < this->NbrXCells; ++i)
	    {
	      if (rand() < Threshold)
		this->Potential->SetPotential(i, j, k, InCoefficient);
	      else
		this->Potential->SetPotential(i, j, k, GaCoefficient);
	    }
    }
  else
    {
      this->Potential = new BinaryThreeDConstantCellPotential(NbrXCells, NbrYCells, NbrZCells);
      this->Potential->LoadBinaryPotential(potentialDescription);
    }
  this->EvaluateInteractionFactors();
}

// copy constructor (without duplicating datas)
//
// hamiltonian = reference on hamiltonian to copy

QuantumWellHamiltonianInMagneticField::QuantumWellHamiltonianInMagneticField(const QuantumWellHamiltonianInMagneticField& hamiltonian)
{
}

// destructor
//

QuantumWellHamiltonianInMagneticField::~QuantumWellHamiltonianInMagneticField()
{
  delete this->Potential;
}

// clone hamiltonian without duplicating datas
//
// return value = pointer to cloned hamiltonian

AbstractHamiltonian* QuantumWellHamiltonianInMagneticField::Clone ()
{
  return new QuantumWellHamiltonianInMagneticField(*this);
}


// set Hilbert space
//
// hilbertSpace = pointer to Hilbert space to use

void QuantumWellHamiltonianInMagneticField::SetHilbertSpace (AbstractHilbertSpace* hilbertSpace)
{
}

// shift Hamiltonian from a given energy
//
// shift = shift value

void QuantumWellHamiltonianInMagneticField::ShiftHamiltonian (double shift)
{
}


// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumWellHamiltonianInMagneticField::MatrixElement (RealVector& V1, RealVector& V2)
{
  Complex Tmp;
  return Tmp;
}
  
// evaluate matrix element
//
// V1 = vector to left multiply with current matrix
// V2 = vector to right multiply with current matrix
// return value = corresponding matrix element

Complex QuantumWellHamiltonianInMagneticField::MatrixElement (ComplexVector& V1, ComplexVector& V2)
{
  Complex Tmp;
  return Tmp;
}

// store Hamiltonian into an hermitian matrix
//
// M = reference on matrix where Hamiltonian has to be stored
// return value = reference on  corresponding hermitian matrix

HermitianMatrix& QuantumWellHamiltonianInMagneticField::GetHamiltonian (HermitianMatrix& M)
{
  M = this->Hamiltonian;
  return M;
}
  
// evaluate all interaction factors
//   

void QuantumWellHamiltonianInMagneticField::EvaluateInteractionFactors()
{
  int Lim = 2 * this->LandauDegeneracy;
  double DiagonalTerm1 = this->ZEnergy1 + ((0.5 + ((double) this->LandauIndex1)) * this->CyclotronEnergy);
  double DiagonalTerm2 = this->ZEnergy2 + ((0.5 + ((double) this->LandauIndex2)) * this->CyclotronEnergy);
  
  for (int i = 0; i < Lim; i += 2)
    {      
      this->Hamiltonian.SetMatrixElement(i, i, DiagonalTerm1);
      this->Hamiltonian.SetMatrixElement(i + 1, i + 1, DiagonalTerm2);
    }

  double XInc = this->MailleParameter * 0.5 / this->MagneticLength;
  double YInc = this->MailleParameter * 0.5 / this->MagneticLength;
  double ZInc = this->MailleParameter;
  double X = 0.5 * XInc;
  double Y = 0.5 * YInc;
  double Z = 0.5 * ZInc;
  double Coefficient;
  double KCoeffcient = this->MagneticLength * 2.0 * M_PI / this->YSize;
  double LandauPrefactor = pow(M_PI * this->MagneticLength * this->MagneticLength, -0.25);  
  double* TmpSin1 = new double [this->NbrZCells];
  double* TmpSin2 = new double [this->NbrZCells];
  double* TmpLandau1 = new double [this->LandauDegeneracy];
  double* TmpLandau2 = new double [this->LandauDegeneracy];
  double Coefficient11;
  double Coefficient12;
  double Coefficient22;
  Complex Tmp11;
  Complex TmpB11;
  Complex TmpB12;
  Complex TmpB21;
  Complex TmpB22;  
  for (int k = 0; k < this->NbrZCells; ++k)
    {
      TmpSin1[k] = sin (2.0 * M_PI * Z / this->ZSize);
      TmpSin2[k] = sin (M_PI * Z / this->ZSize);
      Z += ZInc;
    }
  for (int i = 0; i < this->NbrXCells; ++i)
    {
      for (int m = 0; m < this->LandauDegeneracy; ++m)
	{
	  double ShiftXM = X - (KCoeffcient * m);
	  TmpLandau1[m] = LandauPrefactor * exp (-0.5 * (ShiftXM * ShiftXM));
	  TmpLandau2[m] = sqrt(0.5) * ((2.0 * ShiftXM * ShiftXM) - 1.0) * TmpLandau1[m];
	}
      for (int m = 0; m < this->LandauDegeneracy; ++m)
	{
	  for (int n = m + 1; n < this->LandauDegeneracy; ++n)
	    {	
	      Y = 0.5 * YInc;
	      TmpB11 = 0.0;
	      TmpB12 = 0.0;
	      TmpB22 = 0.0;
	      for (int j = 0; j < this->NbrYCells; ++j)
		{		      
		  Tmp11.Re = cos(Y* KCoeffcient *((double) (n - m)));
		  Tmp11.Im = -sin (Y* KCoeffcient *((double) (n - m)));
		  Coefficient11 = 0.0;
		  Coefficient12 = 0.0;
		  Coefficient22 = 0.0;
		  for (int k = 0; k < this->NbrZCells; ++k)
		    {
		      Coefficient = this->Potential->GetPotential(i, j, k);
		      Coefficient11 += Coefficient * TmpSin1[k] * TmpSin1[k];
		      Coefficient12 += Coefficient * TmpSin1[k] * TmpSin2[k];
		      Coefficient22 += Coefficient * TmpSin2[k] * TmpSin2[k];
		    }
		  TmpB11.Re += Tmp11.Re * Coefficient11;
		  TmpB11.Im += Tmp11.Im * Coefficient11;
		  TmpB12.Re += Tmp11.Re * Coefficient12;
		  TmpB12.Im += Tmp11.Im * Coefficient12;
		  TmpB22.Re += Tmp11.Re * Coefficient22;		  
		  TmpB22.Im += Tmp11.Im * Coefficient22;		  
		  Y += YInc;
		}
	      TmpB21 = TmpB12;
	      TmpB11 *= TmpLandau1[m] * TmpLandau1[n];
	      TmpB21 *= TmpLandau2[m] * TmpLandau1[n];
	      TmpB12 *= TmpLandau1[m] * TmpLandau2[n];
	      TmpB22 *= TmpLandau2[m] * TmpLandau2[n];
	      this->Hamiltonian.AddToMatrixElement(2 * m, 2 * n, TmpB11);
	      this->Hamiltonian.AddToMatrixElement(2 * m, 2 * n + 1, TmpB12);
	      this->Hamiltonian.AddToMatrixElement(2 * m + 1, 2 * n, TmpB21);
	      this->Hamiltonian.AddToMatrixElement(2 * m + 1, 2 * n + 1, TmpB22);
	    }
	  Coefficient11 = 0.0;
	  Coefficient12 = 0.0;
	  Coefficient22 = 0.0;
	  for (int k = 0; k < this->NbrZCells; ++k)
	    {
	      Coefficient = 0.0;
	      for (int j = 0; j < this->NbrYCells; ++j)
		Coefficient += this->Potential->GetPotential(i, j, k);
	      Coefficient11 += Coefficient * TmpSin1[k] * TmpSin1[k];
	      Coefficient12 += Coefficient * TmpSin1[k] * TmpSin2[k];
	      Coefficient22 += Coefficient * TmpSin2[k] * TmpSin2[k];
	    }
	  this->Hamiltonian.AddToMatrixElement(2 * m, 2 * m, Coefficient11  * TmpLandau1[m] * TmpLandau1[m]);
	  this->Hamiltonian.AddToMatrixElement(2 * m, 2 * m + 1, Coefficient12 * TmpLandau1[m] * TmpLandau2[m]);
	  this->Hamiltonian.AddToMatrixElement(2 * m + 1, 2 * m + 1, Coefficient22 * TmpLandau2[m] * TmpLandau2[m]);		  
	}
      X += XInc;
    }
  delete[] TmpSin1;
  delete[] TmpSin2;  
  delete[] TmpLandau1;
  delete[] TmpLandau2;
}

// save potential on disk
// 
// filename = name of the file (with path) where potential has to be saved
// return value = true if no error occured

bool QuantumWellHamiltonianInMagneticField::SavePotential(char* filename)
{
  this->Potential->SaveBinaryPotential(filename);
  return true;
}
