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


#ifndef QUANTUMWELLBFIELDABSORPTIONSPECTRA_H
#define QUANTUMWELLBFIELDABSORPTIONSPECTRA_H

#include "config.h"

#include "Tools/Spectra/Spectra.h"
#include "Matrix/ComplexMatrix.h"


class QuantumWellBFieldAbsorptionSpectra : public Spectra
{

  protected:
  
  // lorentzian broadening parameter
  double Gamma;
  // beta factor (inverse of the temperature in spectrum data energy unit) 
  double Beta;
  // system dimension in the z direction (in Angstrom unit)
  double ZSize;

  // number of initial states per sample
  int NbrInitialStates;
  // number of final states per sample
  int NbrFinalStates;

  // representation of the oscillator strength operator (with a given polarizarion) onto the initial and final state basis
  ComplexMatrix OscillatorStrength;

 public:
  
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
  QuantumWellBFieldAbsorptionSpectra(int NbrFiles, int nbrInitialStates, char** initialStateSpectrumFiles, char*** initialStateEigenstateFiles, 	  
				     int nbrFinalStates, char** finalStateSpectrumFiles, char*** finalStateEigenstateFiles,
				     double thetaPolarizationAngle, double phiPolarizationAngle, double zSize, 
				     double gamma, double beta, double eMin, double eMax, double deltaE);
    
  
 protected:

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
  void AddSample (int nbrInitialStates, char* initialStateSpectrumFileName, char** initialEigenstateFileNames,
		  int nbrFinalStates, char* finalStateSpectrumFileName, char** finalEigenstateFileNames);


  // compute the oscillator strength matrix 
  //
  // thetaPolarizationAngle = angle between the z axis and the polarization vector
  // phiPolarizationAngle = angle between the x axis and the projection of the polarization vector onto the x-y plane
  // zSize = system dimension in the z direction (in Angstrom unit)
  void ComputeOscillatorStrengthMatrix(double thetaPolarizationAngle, double phiPolarizationAngle, double zSize);

};


#endif
