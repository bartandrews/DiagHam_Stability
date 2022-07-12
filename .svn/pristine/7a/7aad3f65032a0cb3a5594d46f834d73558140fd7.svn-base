#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "HilbertSpace/PeriodicXYReflexionZPeriodicThreeDOneParticle.h"

#include "Tools/Spectra/Spectra.h"
#include "Tools/Spectra/XYReflexionSymmetricPeriodicSpectra.h"

#include <iostream>
#include <fstream>

using std::cout;
using std::ifstream;
using std::ofstream;
using std::ios;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);  
  OptionManager Manager ("DefaultXYReflexionSymmetryOscillatorForce" , "0.01");
  OptionGroup* HilbertSpaceGroup = new OptionGroup ("Hilbert space options");
  OptionGroup* FileGroup =  new OptionGroup ("File options");
  OptionGroup* SampleGroup =  new OptionGroup ("Sample options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += HilbertSpaceGroup;
  Manager += FileGroup;
  Manager += SampleGroup;
  Manager += MiscGroup;

  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statex", "number of states in x direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statey", "number of states in y direction of the full Hilbert space (no symmetry reduction)", 81);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "nbr-statez", "number of states in z direction", 21);
  (*HilbertSpaceGroup) += new SingleIntegerOption ('\n', "lowz", "lower impulsion in z direction", -10);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairX1", "pair function in X direction of the first particle", false);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairY1", "pair function in Y direciton of the first particle", false);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairX2", "pair function in X direction of the second particle", false);
  (*HilbertSpaceGroup) += new BooleanOption ('\n', "pairY2", "pair function in Y direciton of the second particle", false);

  (*FileGroup) += new SingleStringOption ('\n', "state", "name of the file containing the first state", "eigenvector.0");
  (*FileGroup) += new SingleStringOption ('\n', "energy", "name of the file energy of the first state", "eigenvalues");
  (*FileGroup) += new SingleIntegerOption ('\n', "last", "number of the last destination state", 1);

  (*SampleGroup) += new SingleDoubleOption ('X', "sizeX", "size of the sample (in Angstrom unit)", 500);
  (*SampleGroup) += new SingleDoubleOption ('Y', "sizeY", "size of the sample (in Angstrom unit)", 500);
  (*SampleGroup) += new SingleDoubleOption ('Z', "sizeZ", "size of the sample (in Angstrom unit)", 500);

  (*MiscGroup) += new BooleanOption ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsLaplacianDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrStateX = ((SingleIntegerOption*) Manager["nbr-statex"])->GetInteger();
  int NbrStateY = ((SingleIntegerOption*) Manager["nbr-statey"])->GetInteger();
  int NbrStateZ = ((SingleIntegerOption*) Manager["nbr-statez"])->GetInteger();
  int LowImpulsionZ = ((SingleIntegerOption*) Manager["lowz"])->GetInteger();
  bool PairX1 = ((BooleanOption*) Manager["pairX1"])->GetBoolean();
  bool PairY1 = ((BooleanOption*) Manager["pairY1"])->GetBoolean();
  bool PairX2 = ((BooleanOption*) Manager["pairX2"])->GetBoolean();
  bool PairY2 = ((BooleanOption*) Manager["pairY2"])->GetBoolean();

  char* State = ((SingleStringOption*) Manager["state"])->GetString();
  char* Energy = ((SingleStringOption*) Manager["energy"])->GetString();
  int Last = ((SingleIntegerOption*) Manager["last"])->GetInteger();
 
  double SizeX = ((SingleDoubleOption*) Manager["sizeX"])->GetDouble();
  double SizeY = ((SingleDoubleOption*) Manager["sizeY"])->GetDouble();
  double SizeZ = ((SingleDoubleOption*) Manager["sizeZ"])->GetDouble();

  // define Hilbert space   
  PeriodicXYReflexionZPeriodicThreeDOneParticle* Space1 = new PeriodicXYReflexionZPeriodicThreeDOneParticle (NbrStateX, PairX1, NbrStateY, PairY1, NbrStateZ, LowImpulsionZ);
  PeriodicXYReflexionZPeriodicThreeDOneParticle* Space2 = new PeriodicXYReflexionZPeriodicThreeDOneParticle (NbrStateX, PairX2, NbrStateY, PairY2, NbrStateZ, LowImpulsionZ);

  XYReflexionSymmetricPeriodicSpectra spectra(Space1, State);  

  double ReX, ImX, ReY, ImY, ReZ, ImZ;
  // polarization in the Z direction
  if ((PairX1 == PairX2) && (PairY1 == PairY2))
    {
      ofstream polarization ("Polarization.txt");  
      ofstream PX("PolarizationX.txt");
      ofstream PY("PolarizationY.txt");  
      ofstream PZ("PolarizationZ.txt");
      ifstream energy ("eigenvalues");

      char** Files = new char* [Last + 1];
      double tmpE;  
      double fundamental;
      energy >> fundamental;
      for (int i = 1; i <= Last; ++i)
	{
	  Files[i] = new char[80];
	  AddString(Files[i], "eigenvector.", i, "");
	  spectra.GetImpulsion(Space2, Files[i], SizeX, SizeY, SizeZ, ReX, ImX, ReY, ImY, ReZ, ImZ);
	  energy >> tmpE;
	  polarization << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << '\n';
	  PX << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << '\n';
	  PY << tmpE - fundamental << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << '\n';
	  PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << '\n';	  
	}
      PX.close(); PY.close(); PZ.close();
      energy.close(); polarization.close();
      delete[] Files;
    }
  // polarization in the plane
  else
    {
      ofstream polarization ("Polarization.txt");  
      ofstream PX("PolarizationX.txt");
      ofstream PY("PolarizationY.txt");  
      ofstream PZ("PolarizationZ.txt");
      ifstream energybis (Energy);
      double fundamental;
      energybis >> fundamental;
      energybis.close();
      ifstream energy ("eigenvalues");

      char** Files = new char* [Last + 1];
      double tmpE;  
      for (int i = 0; i <= Last; ++i)
	{
	  Files[i] = new char[80];
	  AddString(Files[i], "eigenvector.", i, "");
	  spectra.GetImpulsion(Space2, Files[i], SizeX, SizeY, SizeZ, ReX, ImX, ReY, ImY, ReZ, ImZ);
	  energy >> tmpE;
	  polarization << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << '\n';
	  PX << tmpE - fundamental << '\t' << ((ReX * ReX) + (ImX * ImX)) / (tmpE - fundamental) << '\n';
	  PY << tmpE - fundamental << '\t' << ((ReY * ReY) + (ImY * ImY)) / (tmpE - fundamental) << '\n';
	  PZ << tmpE - fundamental << '\t' << ((ReZ * ReZ) + (ImZ * ImZ)) / (tmpE - fundamental) << '\n';	  
	}
      PX.close(); PY.close(); PZ.close();
      energy.close(); polarization.close();
      delete[] Files;
    }

  return 1;
}
