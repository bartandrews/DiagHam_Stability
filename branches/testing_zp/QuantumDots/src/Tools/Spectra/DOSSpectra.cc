#include "Tools/Spectra/DOSSpectra.h"
#include "Vector/RealVector.h"

#include <fstream>
#include <math.h>
#include <stdlib.h>

using std::ifstream;
using std::ios;
using std::cout;
using std::endl;


// constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
//
// FileNumber: number of files, Files: name of files
// StateNumber: integer array containing number of states in each file
// Gamma: FWHM
// Emin, Emax, dE: two energy bounds and the step in energy

DOSSpectra::DOSSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE)
{
  int N = (int) ((Emax - Emin) / dE);
  double * Energy = new double [N];
  double * DOS = new double [N]; 
  double tmp1; 
  double tmp2 = 0.0; 
  double g = Gamma * Gamma / 4;

  for (int i = 0; i < N; ++i)
    {
      Energy[i] = Emin + dE * i;
      DOS[i] = 0.0;
    }

  for (int i = 0; i < FileNumber; ++i)
    {
      int n = StateNumber[i];
      double* tmp = new double [n];
      ifstream file;
      file.open(Files[i],ios::out);
      if (!file.is_open())
        {
	  cout << "Error in open the file: " << Files[i] << "Exit now" << endl;
	  exit(0);
	}
      for (int j = 0; j < n; ++j)
	{
	  file >> tmp[j];
	}
      file.close();

      for (int j = 0; j < N; ++j)
	{
	  tmp1 = Energy[j]; tmp2 = 0.0;
	  for (int k = 0; k < n; ++k)
	    tmp2 += 1.0 / ((tmp1 - tmp[k]) * (tmp1 - tmp[k]) + g);
	  DOS[j] += tmp2;
	}
      tmp = 0;
      delete[] tmp;
    }
  double tmp3 = Gamma/(2.0 * M_PI * ((double) FileNumber));
  for (int i = 0; i < N; ++i)
    DOS[i] *= tmp3;

  this->AxeX = new RealVector(Energy, N);
  this->AxeY = new RealVector(DOS, N);
  this->PointNumber = N;
}

// constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
//
// nbrFiles=  number of files
// spectrumFiles = array of names of the file containing the state spectrum
// nbrStates = number of initial states per sample
// gamma = lorentzian broadening parameter
// eMin = photon minimum energy (must use same unit than the spectrum datas)
// eMax = photon maximum energy (must use same unit than the spectrum datas)
// deltaE = photon energy step (must use same unit than the spectrum datas)

DOSSpectra::DOSSpectra(int nbrFiles, char** spectrumFiles, int nbrStates,  double gamma, double eMin, double eMax, double deltaE)
{
  int N = (int) ((eMax - eMin) / deltaE);
  double* Energy = new double [N];
  double tmp1 = eMin; 
  double tmp2;
  double g = gamma * gamma / 4.0;
  for (int i = 0; i < N; ++i)
    {
      Energy[i] = tmp1;
      tmp1 += deltaE;
    }

  this->AxeX = new RealVector(Energy, N);
  this->AxeY = new RealVector(N, true);

  double* TmpSpectrum = 0;
  if (nbrStates == 0)
    {
      TmpSpectrum = new double[nbrStates];
    }
  for (int i = 0; i < nbrFiles; ++i)
    {
      this->ReadSpectrum(spectrumFiles[i], TmpSpectrum, nbrStates);
      for (int j = 0; j < N; ++j)
	{
	  tmp1 = Energy[j];
	  tmp2 = 0.0;
	  for (int k = 0; k < nbrStates; ++k)
	    tmp2 += 1.0 / ((tmp1 - TmpSpectrum[k]) * (tmp1 - TmpSpectrum[k]) + g);
	  (*(this->AxeY))[j] += tmp2;
	}
    }
  delete[] TmpSpectrum;

  double tmp3 = gamma / (2.0 * M_PI * ((double) nbrFiles));
  for (int i = 0; i < N; ++i)
    (*(this->AxeY))[i] *= tmp3;

  this->PointNumber = N;
}

// read spectrum raw data from a file
// 
// filename = name of  the file that conatins the spectrum (with optional relative/absolute path)
// energies = reference ont the array where energy values will be stored
// nbrValues = reference on the number of energy values to retrieve from the file (0 it has to be automatically determined, and so is energy array allocation)
// return value = true if no error occured

bool DOSSpectra::ReadSpectrum(char* filename, double*& energies, int& nbrValues)
{
  if (nbrValues == 0)
    {
      ifstream File2;
      File2.open(filename, ios::out);
      if (!File2.is_open())
	{
	  cout << "error while opening file : " << filename << endl;
	  return false;
	}
      double Dummy;
      while (File2.tellg() >= 0)
	{
	  File2 >> Dummy;
	  ++nbrValues;
	}      
      --nbrValues;
      File2.close();
      energies = new double[nbrValues];
    }
  ifstream File;
  File.open(filename, ios::out);
  if (!File.is_open())
    {
      cout << "error while opening file : " << filename << endl;
      return false;
    }
  for (int j = 0; j < nbrValues; ++j)
    {
      if (File.tellg() < 0)
	{
	  cout << filename <<  " has to few eigenvalues" << endl;
	  File.close();
	  return false;
	}
      else
	File >> energies[j];
    }
  File.close();
  return true;  
}
