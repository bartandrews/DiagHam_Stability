////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                            base class for spectra                          //
//                                                                            //
//                        last modification : 09/15/2003                      //
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


#include "Tools/Spectra/OverlapSpectra.h"
#include "Vector/RealVector.h"

#include <iostream>
#include <fstream>
#include <strstream>

using std::ifstream;
using std::cout;
using std::endl;
using std::ios;


// default constructor
//

OverlapSpectra::OverlapSpectra()
{
  this->AxeX = 0;
  this->AxeY = 0;
  this->PointNumber = 0;
}

// constructor from a BINARY data file which contains the dimension at each new line
//
// ElectronStateFile, ElectronEnergyFile, ElectronNumber = state file, energy file and number of states for electrons
// HoleStateFile, HoleEnergyFile, HoleNumber = state file, energy file and number of states for holes 

OverlapSpectra::OverlapSpectra(char* ElectronStateFile, char* ElectronEnergyFile ,int ElectronNumber, char* HoleStateFile, char* HoleEnergyFile, int HoleNumber)
{
  int N = HoleNumber * ElectronNumber;
  this->PointNumber = N;
  double* eenergy = new double[ElectronNumber];
  double* henergy = new double[HoleNumber];
  ifstream EState(ElectronStateFile); ifstream HState(HoleStateFile);
  ifstream EEnergy(ElectronEnergyFile); ifstream HEnergy(HoleEnergyFile);

  RealVector* estate = new RealVector[ElectronNumber];
  RealVector* hstate = new RealVector[HoleNumber];
  for (int i = 0; i < ElectronNumber; ++i)
    {
      EState >> estate[i];
      EEnergy >> eenergy[i];
    }
  for (int i = 0; i < HoleNumber; ++i)
    {
      HState >> hstate[i];
      HEnergy >> henergy[i];
    }

  double * Energy = new double [N];
  double * Overlap = new double [N];

  int tmp = 0;

  for (int i = 0; i < ElectronNumber; ++i)
    for (int j = 0; j < HoleNumber; ++j)
      {
	tmp = i * HoleNumber + j;
	Energy[tmp] = eenergy[i] + henergy[j];
	Overlap[tmp] = (estate[i] * hstate[j]);
      }
  this->AxeX = new RealVector(Energy, N);
  this->AxeY = new RealVector(Overlap, N);

  delete[] estate; delete[] hstate;
  EState.close(); HState.close();
  EEnergy.close();HEnergy.close();
}

// constructor from ASCII data file
//
// ElectronStateFile, ElectronEnergyFile, ElectronNumber = state file, energy file and number of states for electrons
// HoleStateFile, HoleEnergyFile, HoleNumber = state file, energy file and number of states for holes   
// NumberState = number of states

OverlapSpectra::OverlapSpectra(char* ElectronStateFile, char* ElectronEnergyFile ,int ElectronNumber, char* HoleStateFile, char* HoleEnergyFile, int HoleNumber, int NumberState)
{
  int N = HoleNumber * ElectronNumber;
  PointNumber = N;
  double* eenergy = new double[ElectronNumber];
  double* henergy = new double[HoleNumber];

  RealVector* estate = new RealVector[ElectronNumber];
  RealVector* hstate = new RealVector[HoleNumber];
  ifstream EState(ElectronStateFile); ifstream HState(HoleStateFile);
  ifstream EEnergy(ElectronEnergyFile); ifstream HEnergy(HoleEnergyFile);

  for (int i = 0; i < ElectronNumber; ++i)
    {
      estate[i] = RealVector(NumberState);
      for (int j = 0; j < NumberState; ++j)
	EState >> estate[i][j];
      EEnergy >> eenergy[i];
    }
  for (int i = 0; i < HoleNumber; ++i)
    {
      hstate[i] = RealVector(NumberState);
      for (int j = 0; j < NumberState; ++j)
	HState >> hstate[i][j];
      HEnergy >> henergy[i];
    }

  double * Energy = new double [N];
  double * Overlap = new double [N];

  int tmp = 0;

  for (int i = 0; i < ElectronNumber; ++i)
    for (int j = 0; j < HoleNumber; ++j)
      {
	tmp = i * HoleNumber + j;
	Energy[tmp] = eenergy[i] + henergy[j];
	Overlap[tmp] = estate[i] * hstate[j];
      }
  AxeX = new RealVector(Energy, N);
  AxeY = new RealVector(Overlap, N);

  delete[] estate; delete[] hstate;
  EState.close(); HState.close();
  EEnergy.close();HEnergy.close();
}

// method to write in ASCII mode with 2 columns: recombination energy and square overlap
//
// fileName =  name of the file where the spectrum will be stored
// return = true if no error occurs

bool OverlapSpectra::WriteSquareOverlap(char* fileName)
{
  ofstream File;
  File.open(fileName, ios::binary | ios::out | ios::app);
  File.precision(14);
  for (int i = 0; i < this->PointNumber; ++i)
    {
      File << (*AxeX)[i] << '\t';
      File << ((*AxeY)[i]) * ((*AxeY)[i]) << '\n';
    }
  File.close();
  return true;
}

// evaluate the overlap between two functions in periodic basis
//
// ElectronStateFile & HoleStateFile: state files
// Dimension: dimension of the Hilbert space
// Real: reference to the real part of overlap
// Imaginary: reference to the imaginary part of overlap

void OverlapSpectra::GetOverlap(char* ElectronStateFile, char* HoleStateFile, int Dimension, double& Real, double& Imaginary)
{
  ifstream EState(ElectronStateFile); ifstream HState(HoleStateFile);
  double TmpReE = 0.0; double TmpReH = 0.0; double TmpImE = 0.0; double TmpImH = 0.0;
  Real = 0.0; Dimension = 0.0;
  for (int n = 0; n < Dimension; ++n)
    {
      EState >> TmpReE >> TmpImE; 
      HState >> TmpReH >> TmpImH; 
      Real += (TmpReE * TmpReH - TmpImE * TmpImH);
      Imaginary += (TmpImE * TmpReH + TmpReE * TmpImH);
    }
  EState.close(); HState.close(); 
}
