#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereEntanglementSpectrumChargeCumulants" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "density-matrix", "file containing the reduced density matrix");  
  (*SystemGroup) += new SingleIntegerOption ('n', "nbr-particles", "number of particles in the A part", 0l);
  (*SystemGroup) += new SingleIntegerOption ('l', "nbr-orbitals", "number of orbitals in the A part", 0l);
  (*SystemGroup) += new SingleIntegerOption ('c', "nbr-cumulants", "number of cumulant that have to be computed", 2l);
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output name for the entanglement spectrum (default name replace density-matrix full.ent extension with la_x_na_y.charge)");
  (*SystemGroup) += new BooleanOption ('\n', "realspace-entanglement", "compute the charge cumulants from the real space entanglement spectrum");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereEntanglementSpectrumChargeCumulants -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrOrbitalsInPartition = Manager.GetInteger("nbr-orbitals");
  int NbrParticlesInPartition = Manager.GetInteger("nbr-particles");
  int NbrCumulant = Manager.GetInteger("nbr-cumulants");
  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int TotalSz = 0;
  bool Statistics = true;
  
  if (Manager.GetString("density-matrix") == 0)
    {
      cout << "a reduced density matrix has to be provided, see man page for option syntax or type FQHESphereEntanglementSpectrumChargeCumulants -h" << endl;
      return -1;
    }

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("density-matrix"), NbrParticles, NbrFluxQuanta, TotalLz, Statistics) == false)
    {
      cout << "can't retrieve system informations from the reduced density matrix file name" << endl;
      return -1;
    } 

  MultiColumnASCIIFile DensityMatrix;
  if (DensityMatrix.Parse(Manager.GetString("density-matrix")) == false)
    {
      DensityMatrix.DumpErrors(cout);
      return -1;
    }

  int MinLza = 1 << 30;
  int MaxLza = -MinLza; 
  int* LzaValueArray = 0;

  double* CentralMoments  = new double[NbrCumulant + 1];
  double* Moments = new double[NbrCumulant + 1];
  double* Cumulants = new double[NbrCumulant + 1];
  for (int i = 0; i <= NbrCumulant; ++i)
    {
      CentralMoments[i] = 0.0;
      Moments[i] = 0.0;
      Cumulants[i] = 0.0;
    }
  int* NaValues = 0;
  double* Coefficients = 0;
  long MinIndex = 0l;
  long MaxIndex = DensityMatrix.GetNbrLines();
  char* OutputFileName = Manager.GetString("output");
   
  if (Manager.GetBoolean("realspace-entanglement") == false)
    {
      if (DensityMatrix.GetNbrColumns() < 4)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      int* LaValues = DensityMatrix.GetAsIntegerArray(0);
      NaValues = DensityMatrix.GetAsIntegerArray(1);
      Coefficients = DensityMatrix.GetAsDoubleArray(3);
      while ((MinIndex < MaxIndex) && (LaValues[MinIndex] != NbrOrbitalsInPartition))
	++MinIndex;
      long TmpIndex = MinIndex;
      while ((TmpIndex < MaxIndex) && (LaValues[TmpIndex] == NbrOrbitalsInPartition))
	{
	  ++TmpIndex;
	}
      MaxIndex = TmpIndex;
      if (OutputFileName == 0)
	{
	  char* TmpExtension = new char[256];
	  sprintf(TmpExtension, "la_%d_na_%d.charge", NbrOrbitalsInPartition, NbrParticlesInPartition);
	  if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent", TmpExtension);
	    }
	  else
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent.bz2", TmpExtension);
	    }
	}
    }      
  else
    {
      if (DensityMatrix.GetNbrColumns() < 3)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      NaValues = DensityMatrix.GetAsIntegerArray(0);
      Coefficients = DensityMatrix.GetAsDoubleArray(2);
      if (OutputFileName == 0)
	{
	  char* TmpExtension = new char[256];
	  sprintf(TmpExtension, "rsescharge");
	  if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.rses", TmpExtension);
	    }
	  else
	    {
	      OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.rses.bz2", TmpExtension);
	    }
	}
    }
  
  long TmpIndex = MinIndex;
  while (TmpIndex < MaxIndex)
    {
      Moments[0] += Coefficients[TmpIndex];
      Moments[1] += ((double) NaValues[TmpIndex]) * Coefficients[TmpIndex];
      ++TmpIndex;
    }
  CentralMoments[0] = Moments[0];
  CentralMoments[1] = 0.0;
  Cumulants[0] = 0.0;
  Cumulants[1] = Moments[1];
  
  TmpIndex = MinIndex;
  while (TmpIndex < MaxIndex)
    {
      double Tmp1 = ((double) NaValues[TmpIndex]) - Moments[1];
      double Tmp2 = Tmp1 * Tmp1;
      double Tmp3 = ((double) NaValues[TmpIndex]);
      double Tmp4 = Tmp3 * Tmp3;
      for (int i = 2; i <= NbrCumulant; ++i)
	{
	  CentralMoments[i] += Tmp2 * Coefficients[TmpIndex];
	  Tmp2 *= Tmp1;
	  Moments[i] += Tmp4 * Coefficients[TmpIndex];
	  Tmp4 *= Tmp3;
	}
      ++TmpIndex;
    }

  BinomialCoefficients Binomials (NbrCumulant - 1);

  for (int i = 2; i <= NbrCumulant; ++i)
    {
      Cumulants[i] = CentralMoments[i];
      for (int j = 1; j < (i - 1); ++j)
	Cumulants[i] -= Binomials.GetNumericalCoefficient(i - 1, i - j - 1) * Cumulants[i - j] * CentralMoments[j];
    }


  ofstream File;
  File.open(OutputFileName, ios::out);
  File.precision(14);
  File << "# i <(N_A-<N_A>)^i> <N_A^i> <<N_A^i>>" << endl;
  for (int i = 0; i <= NbrCumulant; ++i)
    {
      File << i << " " << CentralMoments[i] << " " << Moments[i] << " " << Cumulants[i] << endl;
    }
  File.close();	      

  return 0;
}

