////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for absorption spectra of quantum well            //
//                                in magnetic field                           //
//                                                                            //
//                        last modification : 06/12/2005                      //
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

#include "Tools/Spectra/QuantumWellBFieldAbsorptionSpectra.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include <fstream>
#include <math.h>
#include <stdlib.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


// constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
//
// nbrFiles=  number of files
// nbrInitialStates = number of initial states per sample
// initialStateSpectrumFiles = array of names of the file containing the initial state spectrum
// initialStateEigenstateFiles = pointers to arrays that contains names of the eigenvectors associated to each spectrum (for a given spectrum, 
//                               eigenvectors have to be sorted in the same manner as the eigenvalues)
// nbrFinalStates = number of final states per sample
// finalStateSpectrumFiles = array of names of the file containing the  final state spectrum
// finalStateEigenstateFiles = pointers to arrays that contains names of the eigenvectors associated to each spectrum (for a given spectrum, 
//                             eigenvectors have to be sorted in the same manner as the eigenvalues)
// thetaPolarizationAngle = angle between the z axis and the polarization vector
// phiPolarizationAngle = angle between the x axis and the projection of the polarization vector onto the x-y plane
// zSize = system dimension in the z direction (in Angstrom unit)
// gamma = lorentzian broadening parameter
// beta = beta factor (inverse of the temperature in spectrum data energy unit) 
// eMin = photon minimum energy (must use same unit than the spectrum datas)
// eMax = photon maximum energy (must use same unit than the spectrum datas)
// deltaE = photon energy step (must use same unit than the spectrum datas)

QuantumWellBFieldAbsorptionSpectra::QuantumWellBFieldAbsorptionSpectra(int nbrFiles, int nbrInitialStates, char** initialStateSpectrumFiles, char*** initialStateEigenstateFiles, 	  
								       int nbrFinalStates, char** finalStateSpectrumFiles, char*** finalStateEigenstateFiles,
								       double thetaPolarizationAngle, double phiPolarizationAngle, double zSize, 
								       double gamma, double beta, double eMin, double eMax, double deltaE)
{
  this->Gamma = gamma;
  this->Beta = beta;
  this->ZSize = zSize;
  this->NbrInitialStates = nbrInitialStates;
  this->NbrFinalStates = nbrFinalStates;

  int N = (int) ((eMax - eMin) / deltaE);
  double* Energy = new double [N];
  double tmp1 = eMin; 
  for (int i = 0; i < N; ++i)
    {
      Energy[i] = tmp1;
      tmp1 += deltaE;
    }

  this->AxeX = new RealVector(Energy, N);
  this->AxeY = new RealVector(N, true);

  this->ComputeOscillatorStrengthMatrix(thetaPolarizationAngle, phiPolarizationAngle, zSize);

  for (int i = 0; i < nbrFiles; ++i)
    {
      cout << "evaluate contribution of " << initialStateSpectrumFiles[i] << " and " << finalStateSpectrumFiles[i] << endl;
      this->AddSample(nbrInitialStates, initialStateSpectrumFiles[i], initialStateEigenstateFiles[i],
		      nbrFinalStates, finalStateSpectrumFiles[i], finalStateEigenstateFiles[i]);      
    }
  double tmp3 = 1.0 / ((double) nbrFiles);
  for (int i = 0; i < N; ++i)
    (*(this->AxeY))[i] *= tmp3;

  this->PointNumber = N;
}

// add contribution of a given sample to the absorption spectrum
//
// nbrInitialStates = number of initial states
// initialStateSpectrumFileName = name of the file containing the initial state spectrum
// initialStateEigenstateFile = array that contains names of the eigenvectors associated to the spectrum (for a given spectrum, 
//                              eigenvectors have to be sorted in the same manner as the eigenvalues)
// nbrFinalStates = number of final states
// finalStateSpectrumFiles = name of the file containing the final state spectrum
// finalStateEigenstateFileName = array that contains names of the eigenvectors associated to the spectrum (for a given spectrum, 
//                                eigenvectors have to be sorted in the same manner as the eigenvalues)

void QuantumWellBFieldAbsorptionSpectra::AddSample (int nbrInitialStates, char* initialStateSpectrumFileName, char** initialEigenstateFileNames,
						    int nbrFinalStates, char* finalStateSpectrumFileName, char** finalEigenstateFileNames)
{
  int NbrFinalStates = 2 * nbrInitialStates;
  double* InitialEnergies = new double [nbrInitialStates];
  double* FinalEnergies = new double [NbrFinalStates];

  if (this->ReadSpectrum(initialStateSpectrumFileName, InitialEnergies, nbrInitialStates) == false)
    {
      delete[] InitialEnergies;
      delete[] FinalEnergies;
      return;
    }
  if (this->ReadSpectrum(finalStateSpectrumFileName, FinalEnergies, nbrFinalStates) == false)
    {
      delete[] InitialEnergies;
      delete[] FinalEnergies;
      return;
    }
      
  ComplexVector TmpInitialVector(nbrInitialStates);
  ComplexVector TmpInitialVector2(nbrInitialStates);
  ComplexVector TmpFinalVector(NbrFinalStates);
  double Tmp;
  double Tmp2;
  int TmpNbrValues = this->AxeX->GetVectorDimension();
  double Factor = this->Gamma / (2.0 * M_PI);
  double GammaSqr = 0.25 * this->Gamma * this->Gamma;
  double TmpEnergy;
  for (int i = 0; i < nbrFinalStates; ++i)
    {
      TmpFinalVector.ReadVector(finalEigenstateFileNames[i]);   
      TmpInitialVector2.Multiply(this->OscillatorStrength, TmpFinalVector);
      Tmp2 = Factor * exp (-this->Beta * FinalEnergies[i]);
      for (int j = 0; j < nbrInitialStates; ++j)
	if (FinalEnergies[i] > InitialEnergies[j])
	  {
	    TmpInitialVector.ReadVector(initialEigenstateFileNames[j]);
	    Tmp = SqrNorm(TmpInitialVector * TmpInitialVector2) * Tmp2;
	    for (int k = 0; k < TmpNbrValues; ++k)
	      {
		TmpEnergy = FinalEnergies[i] - InitialEnergies[j] - (*(this->AxeX))[k];
		(*(this->AxeY))[k] += Tmp / ((TmpEnergy *TmpEnergy) + GammaSqr);
	      }
	}
    }
  delete[] InitialEnergies;
  delete[] FinalEnergies;
}

// compute the oscillator strength matrix 
//
// thetaPolarizationAngle = angle between the z axis and the polarization vector
// phiPolarizationAngle = angle between the x axis and the projection of the polarization vector onto the x-y plane
// zSize = system dimension in the z direction (in Angstrom unit)

void QuantumWellBFieldAbsorptionSpectra::ComputeOscillatorStrengthMatrix(double thetaPolarizationAngle, double phiPolarizationAngle, double zSize)
{
  this->OscillatorStrength = ComplexMatrix (this->NbrInitialStates, this->NbrFinalStates, true);
  int HalfNbrFinalStates = this->NbrFinalStates / 2;

  double ZFactor = cos (thetaPolarizationAngle) * 2.0 * zSize / (M_PI * M_PI);
  double Tmp = 0.0;
  int Diff = 2 - 1;
  int Sum = 2 + 1;
  if ((Diff & 1) != 0)
    Tmp += 2.0 / ((double) (Diff * Diff));
  if ((Sum & 1) != 0)
    Tmp -= 2.0 / ((double) (Sum * Sum));
  ZFactor *= Tmp;
  cout << "ZFactor = " << ZFactor << endl;
  for (int m = 0; m < HalfNbrFinalStates; ++m)
    for (int n = 0; n < this->NbrInitialStates; ++n)
      {
	this->OscillatorStrength.SetMatrixElement(n, 2 * m, ZFactor);
      }
}


