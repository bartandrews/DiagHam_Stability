////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                     Copyright (C) 2003 Duc-Phuong Nguyen                   //
//                                                                            //
//                  class for spectra of time resolved experiment             //
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



#include "Tools/Spectra/TimeResolvedPLSpectra.h"

#include "Vector/RealVector.h"
#include "Matrix/RealMatrix.h"

#include <iostream>
#include <fstream>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


// constructor from a set of squared overlap files
//
// FileNumber = number of squared overlap files
// FileName = name of square overlap files
// ElectronState, HoleState = number of states in each file
// Emin, Emax = range of energy observed
// T, NT = time and number of time intervals

TimeResolvedPLSpectra::TimeResolvedPLSpectra(int FileNumber, char** FileName, int* ElectronState, int* HoleState, double Emin, double Emax, double T, int NT):Spectra(NT)
{
  RealMatrix* Overlap1 = new RealMatrix[FileNumber];
  RealMatrix* Overlap2 = new RealMatrix[FileNumber];
  RealMatrix* Energy = new RealMatrix[FileNumber];
  RealVector* Electron = new RealVector[FileNumber];
  RealVector* Hole = new RealVector[FileNumber];
  RealVector* IE = new RealVector[FileNumber];
  RealVector* IH = new RealVector[FileNumber];
  int NE = 0; int NH = 0; double tmp = 0.0;
  for (int n = 0; n < FileNumber; ++n)
    {
      NE = ElectronState[n]; NH = HoleState[n];
      Overlap1[n] = RealMatrix(NE, NH);
      Overlap2[n] = RealMatrix(NH, NE);
      Energy[n] = RealMatrix(NE, NH);
      Electron[n] = RealVector(NE);
      Hole[n] = RealVector(NH);
      IE[n] = RealVector(NE);
      IH[n] = RealVector(NH);
      for (int i = 0; i < NE; ++i)
	{
	  Electron[n][i] = 1.0;
	  IE[n][i] = 1.0;
	}
      for (int j = 0; j < NH; ++j)
	{
	  Hole[n][j] = 1.0;
	  IH[n][j] = 1.0;
	}
      ifstream File;
      File.open(FileName[n], ios::in);
      for (int i = 0; i < NE; ++i)
	for (int j = 0; j < NH; ++j)
	  {
	    File >> tmp;	   
	    Energy[n](i, j) = tmp;	    
	    File >> tmp;	    
	    Overlap2[n](j, i) = tmp;
	    Overlap1[n](i, j) = tmp;
	  }      
      File.close();
    }
  int counter = 0;
  for (int n = 0; n < FileNumber; ++n)
    for (int i = 0; i < ElectronState[n]; ++i)
      for (int j = 0; j < HoleState[n]; ++j)
	if ((Emin < Energy[n](i,j)) && (Emax > Energy[n](i,j)))
	  ++counter;
  int* ListF = new int[counter];
  int* ListE = new int[counter];
  int* ListH = new int[counter];
  counter = 0; double deltaT = T/NT; double signal = 0.0;
  for (int n = 0; n < FileNumber; ++n)
    for (int i = 0; i < ElectronState[n]; ++i)
      for (int j = 0; j < HoleState[n]; ++j)
	if ((Emin < Energy[n](i,j)) && (Emax > Energy[n](i,j)))
	  {
	    ListF[counter] = n;
	    ListE[counter] = i;
	    ListH[counter] = j;
	    signal += Overlap1[n](i, j);
	    ++counter;
	  }
  (*AxeX)[0] = 0.0;
  (*AxeY)[0] = signal;

  int tF; int tE; int tH; RealVector tvE(NE); RealVector tvH(NH);
  double signal10 = signal/10.0;
  for (int nt = 1; nt < NT; ++nt)
    {
      signal = 0.0;
      for (int n = 0; n < FileNumber; ++n)
	{
	  for (int indice = 0; indice < NE; ++indice)
	    Electron[n][indice] *= (1.0 - ((tvE.Multiply(Overlap1[n], Hole[n])[indice]) * deltaT));
	  for (int indice = 0; indice < NH; ++indice)
	    Hole[n][indice] *= (1.0 - ((tvH.Multiply(Overlap2[n], Electron[n])[indice]) * deltaT));

	  // Electron[n] *= (IE[n] - ((tv.Multiply(Hole[n], Overlap1[n])) * deltaT));
	  // Hole[n] *= (IH[n] - ((tv.Multiply(Electron[n], Overlap2[n])) * deltaT));       
	}
      for (int i = 0; i < counter; ++i)
	{
	  tF = ListF[i]; tE = ListE[i]; tH = ListH[i];
	  signal += Electron[tF][tE] * Overlap1[tF](tE,tH) * Hole[tF][tH];
	}
      (*AxeX)[nt] = nt * deltaT;
      (*AxeY)[nt] = signal;    
      if ((signal < signal10) && (tmp > signal10))
	T10 = nt * deltaT;
      tmp = signal;
    }  
}
