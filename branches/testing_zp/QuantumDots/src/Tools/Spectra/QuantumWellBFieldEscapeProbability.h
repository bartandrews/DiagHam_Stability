////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2004 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                    class for escape probability of quantum well            //
//                                in magnetic field                           //
//                                                                            //
//                        last modification : 29/05/2006                      //
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


#ifndef QUANTUMWELLBFIELDESCAPEPROBABILITY_H
#define QUANTUMWELLBFIELDESCAPEPROBABILITY_H

#include "config.h"

#include "Tools/Spectra/Spectra.h"
#include "Matrix/RealMatrix.h"


class QuantumWellBFieldEscapeProbability : public Spectra
{

  protected:
  
  // time step value (in hbar/E units)
  double TimeStep;
  // number of time steps
  int NbrTimeSteps;  

  // index of the initial state (-1 if probality has to evaluated for all possible states)
  int InitialStateIndex;

  // number of states per sample and per Landau level
  int NbrStates;
  // array that contains the probability as a function of time for each possible state in the lowest Landau level
  RealMatrix Probabilities;

  // flag to indicate if it has to compute the probability to stay in the same subband that the initial state instead of the probability to stay in the initial state
  bool SubbandSumFlag;
  // flag to indicate if the -log of the probability has to be plotted instead of the probability
  bool LogarithmicPlotFlag;

 public:
  
  // constructor from a set of energy files
  //
  // nbrFiles=  number of files
  // nbrStates = number of states per sample and per Landau level
  // stateSpectrumFiles = array of names of the file containing the state spectrum
  // stateEigenstateFiles = pointers to arrays that contains names of the eigenvectors associated to each spectrum (for a given spectrum, 
  //                        eigenvectors have to be sorted in the same manner as the eigenvalues)
  // timeStep = time step value (in hbar/E units)
  // nbrTimeSteps = number of time steps
  // initialStateIndex = index of the initial state (-1 if probality has to evaluated for all possible states)
  // subbandSum =  true if it has to compute the probability to stay in the same subband that the initial state instead of the probability to stay in the initial state
  // logarithmicPlot = true if the -log of the probability has to be plotted instead of the probability
  QuantumWellBFieldEscapeProbability(int nbrFiles, int nbrStates, char** stateSpectrumFiles, char*** stateEigenstateFiles, 	  
				     double timeStep, int nbrTimeSteps, int initialStateIndex, 
				     bool subbandSum = false, bool logarithmicPlot = false);
    
  
  // virtual method to write the spectrum in a file in ASCII mode
  //
  // fileName = name of the file where the spectrum will be stored
  // return = true if no error occurs
  bool WriteSpectra(char * fileName);

 protected:

  // add contribution of a given sample to the escape probability
  //
  // nbrStates = number of states per Landau level
  // stateSpectrumFileName = name of the file containing the state spectrum
  // eigenstateFile = array that contains names of the eigenvectors associated to the spectrum (for a given spectrum, 
  //                              eigenvectors have to be sorted in the same manner as the eigenvalues)
  void AddSample (int nbrStates, char* stateSpectrumFileName, char** eigenstateFileNames);


};


#endif
