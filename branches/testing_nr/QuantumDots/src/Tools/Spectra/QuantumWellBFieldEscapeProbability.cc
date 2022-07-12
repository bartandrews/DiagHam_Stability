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


#include "config.h"

#include "Tools/Spectra/QuantumWellBFieldEscapeProbability.h"
#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include <fstream>
#include <math.h>
#include <stdlib.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


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

QuantumWellBFieldEscapeProbability::QuantumWellBFieldEscapeProbability(int nbrFiles, int nbrStates, char** stateSpectrumFiles, char*** stateEigenstateFiles, 	  
								       double timeStep, int nbrTimeSteps, int initialStateIndex, bool subbandSum, bool logarithmicPlot)

{
  this->TimeStep = timeStep;
  this->NbrTimeSteps = nbrTimeSteps;
  this->NbrStates = nbrStates;
  this->InitialStateIndex = initialStateIndex;
  this->LogarithmicPlotFlag = logarithmicPlot;
  this->SubbandSumFlag = subbandSum;
  
  double* TimeValues = new double [this->NbrTimeSteps];
  double tmp1 = 0.0; 
  for (int i = 0; i < this->NbrTimeSteps; ++i)
    {
      TimeValues[i] = tmp1;
      tmp1 += this->TimeStep;
    }

  this->AxeX = new RealVector(TimeValues, this->NbrTimeSteps);
  this->AxeY = new RealVector(this->NbrTimeSteps, true);
  
  this->Probabilities = RealMatrix(this->NbrTimeSteps, this->NbrStates, true);

//  nbrFiles = 1;
  for (int i = 0; i < nbrFiles; ++i)
    {
      cout << "evaluate contribution of " << stateSpectrumFiles[i] << endl;
      this->AddSample(nbrStates, stateSpectrumFiles[i], stateEigenstateFiles[i]);      
    }
  double tmp3 = 1.0 / ((double) nbrFiles);
//   for (int i = 0; i < this->NbrTimeSteps; ++i)
//     for (int j = 0; j < this->NbrStates; ++j)
  this->Probabilities *= tmp3;

  this->PointNumber = this->NbrTimeSteps;
}

// add contribution of a given sample to the escape probability
//
// nbrStates = number of states per Landau level
// stateSpectrumFileName = name of the file containing the state spectrum
// eigenstateFile = array that contains names of the eigenvectors associated to the spectrum (for a given spectrum, 
//                              eigenvectors have to be sorted in the same manner as the eigenvalues)

void QuantumWellBFieldEscapeProbability::AddSample (int nbrStates, char* stateSpectrumFileName, char** eigenstateFileNames)
{
  int TotalHilbertSpaceSize = 2 * nbrStates;

  double* Energies = new double [TotalHilbertSpaceSize];
  if (this->ReadSpectrum(stateSpectrumFileName, Energies, TotalHilbertSpaceSize) == false)
    {
      delete[] Energies;
      return;
    }

  ComplexVector* TmpVectors = new ComplexVector [TotalHilbertSpaceSize];
  for (int i = 0; i < TotalHilbertSpaceSize; ++i)
    {
      TmpVectors[i].ReadVector(eigenstateFileNames[i]);   
    }
     
  int Shift = 0;
  HermitianMatrix ReducedHamiltonian (nbrStates, true);
  for (int n = 0; n < TotalHilbertSpaceSize; ++n)
    {
      ComplexVector& TmpVector2 = TmpVectors[n];
      for (int i1 = 0; i1 < nbrStates; ++i1)
 	for (int i2 = 0; i2 < nbrStates; ++i2)
 	  ReducedHamiltonian.AddToMatrixElement(i1, i2, Energies[n] * Conj(TmpVector2[(2 * i1) + Shift]) * TmpVector2[(2 * i2) + Shift]);
    }
  RealDiagonalMatrix DiagonalizedHamiltonian (ReducedHamiltonian.GetNbrRow());
  ComplexMatrix Eigenvectors(ReducedHamiltonian.GetNbrRow(), ReducedHamiltonian.GetNbrRow());
  ReducedHamiltonian.Diagonalize(DiagonalizedHamiltonian, Eigenvectors);


  int n = 0;
  int Lim = nbrStates;
  if ((this->InitialStateIndex >= 0) && (this->InitialStateIndex < nbrStates))
    {
      n = this->InitialStateIndex;
      Lim = this->InitialStateIndex + 1;
    }
  if (this->SubbandSumFlag == false)
    {
      double* TmpCoefficients = new double [TotalHilbertSpaceSize];
      for (; n < Lim; ++n)
	{
	  ComplexVector& TmpVector3 = Eigenvectors[n];
	  for (int i = 0; i < TotalHilbertSpaceSize; ++i)
	    {
	      Complex Tmp = 0.0;
	      ComplexVector& TmpVector2 = TmpVectors[i];
	      for (int j = 0; j < nbrStates; ++j)
		Tmp += Conj(TmpVector2[(2 * j) + Shift]) * TmpVector3[j];
	      TmpCoefficients[i] = SqrNorm(Tmp);
	    }
	  
	  double Time = 0.0;
	  for (int t = 0; t < this->NbrTimeSteps; ++t)
	    {
	      Complex Tmp = 0.0;
	      for (int i = 0; i < TotalHilbertSpaceSize; ++i)
		{
		  Tmp.Re += TmpCoefficients[i] * cos (Time * Energies[i]);
		  Tmp.Im += TmpCoefficients[i] * sin (Time * Energies[i]);
		}
	      this->Probabilities(t, n) += SqrNorm(Tmp);
	      Time += this->TimeStep;
	    }
	}
      delete[] TmpCoefficients;      
    }
  else
    {
      ComplexMatrix ProjectedScalarProducts(TotalHilbertSpaceSize, TotalHilbertSpaceSize, true);
      for (int i = 0; i < TotalHilbertSpaceSize; ++i)
	{
	  ComplexVector& TmpVector2 = TmpVectors[i];
	  Complex Tmp;
	  for (int j = i; j < TotalHilbertSpaceSize; ++j)
	    {
	      ProjectedScalarProducts.SetMatrixElement(i, j, TmpVector2.PartialScalarProduct(TmpVectors[j], Shift, nbrStates, 2));	  
	    }
	}
      Complex* TmpCoefficients = new Complex [TotalHilbertSpaceSize];
      for (; n < Lim; ++n)
	{
	  ComplexVector& TmpVector3 = Eigenvectors[n];
	  for (int i = 0; i < TotalHilbertSpaceSize; ++i)
	    {
	      Complex Tmp = 0.0;
	      ComplexVector& TmpVector2 = TmpVectors[i];
	      for (int j = 0; j < nbrStates; ++j)
		Tmp += Conj(TmpVector2[(2 * j) + Shift]) * TmpVector3[j];
	      TmpCoefficients[i] = Tmp;
	    }
	  
	  double Time = 0.0;
	  Complex Tmp2;
	  for (int t = 0; t < this->NbrTimeSteps; ++t)
	    {
	      double Tmp = 0.0;
	      for (int i = 0; i < TotalHilbertSpaceSize; ++i)
		{
		  Tmp += SqrNorm(TmpCoefficients[i]) * ProjectedScalarProducts[i][i].Re;
		  for (int j = i + 1; j < TotalHilbertSpaceSize; ++j)
		    {
		      Tmp2 = ProjectedScalarProducts[j][i] * Conj(TmpCoefficients[i]) * (TmpCoefficients[j]);		      
		      Tmp += 2.0 * cos (Time * (Energies[i] - Energies[j])) * Tmp2.Re;
		      Tmp -= 2.0 * sin (Time * (Energies[i] - Energies[j])) * Tmp2.Im;
		    }
		}
	      this->Probabilities(t, n) += Tmp;
	      Time += this->TimeStep;
	    }
	}
      delete[] TmpCoefficients;      
    }


  delete[] Energies;
  delete[] TmpVectors;
}

// virtual method to write the spectrum in a file in ASCII mode
//
// fileName = name of the file where the spectrum will be stored
// return = true if no error occurs

bool QuantumWellBFieldEscapeProbability::WriteSpectra(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::out);
  File.precision(14);
  for (int i = 0; i < this->PointNumber; ++i)
    {
      File << (*AxeX)[i] << '\t';
      if (this->LogarithmicPlotFlag == false)
	if (this->InitialStateIndex < 0)
	  for (int j = 0; j < NbrStates; ++j)
	    File << '\t' <<  this->Probabilities(i, j);
	else
	  File << '\t' <<  this->Probabilities(i, this->InitialStateIndex);	
      else
	if (this->InitialStateIndex < 0)
	  for (int j = 0; j < NbrStates; ++j)
	    File << '\t' <<  -log (this->Probabilities(i, j));
	else
	  File << '\t' <<  -log (this->Probabilities(i, this->InitialStateIndex));	
      File << endl;
    }
  File.close();
  return true;
}

