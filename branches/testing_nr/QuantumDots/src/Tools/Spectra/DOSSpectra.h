#ifndef DOSSPECTRA_H
#define DOSSPECTRA_H

#include "config.h"

#include "Tools/Spectra/Spectra.h"

class DOSSpectra : public Spectra
{

 public:

  // constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
  //
  // FileNumber: number of files, Files: name of files
  // StateNumber: integer array containing number of states in each file
  // Gamma: FWHM
  // Emin, Emax, dE: two energy bounds and the step in energy
  DOSSpectra(int FileNumber, char** Files, int * StateNumber, double Gamma, double Emin, double Emax, double dE);

  // constructor from a set of energy files. Each peak is assimilated to a Lorentzian function.
  //
  // nbrFiles=  number of files
  // spectrumFiles = array of names of the file containing the state spectrum
  // nbrStates = number of initial states per sample
  // gamma = lorentzian broadening parameter
  // eMin = photon minimum energy (must use same unit than the spectrum datas)
  // eMax = photon maximum energy (must use same unit than the spectrum datas)
  // deltaE = photon energy step (must use same unit than the spectrum datas)
  DOSSpectra(int nbrFiles, char** spectrumFiles, int nbrStates, double gamma, double eMin, double eMax, double deltaE);

 private:

  // read spectrum raw data from a file
  // 
  // filename = name of  the file that conatins the spectrum (with optional relative/absolute path)
  // energies = reference ont the array where energy values will be stored
  // nbrValues = reference on the number of energy values to retrieve from the file (0 it has to be automatically determined, and so is energy array allocation)
  // return value = true if no error occured
  bool ReadSpectrum(char* filename, double*& energies, int& nbrValues);

};

#endif
