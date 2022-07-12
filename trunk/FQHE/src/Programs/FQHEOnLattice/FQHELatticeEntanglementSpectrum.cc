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

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"

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
  OptionManager Manager ("FQHELatticeEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "density-matrix", "file containing the reduced density matrix");  
  (*SystemGroup) += new SingleIntegerOption ('n', "nbr-particles", "number of particles in the A part", 0l);
  //  (*SystemGroup) += new SingleIntegerOption ('l', "nbr-orbitals", "number of orbitals in the A part", 0l);
  //  (*SystemGroup) += new SingleIntegerOption ('s', "sz-value", "twice the Sz value of A part (SU(2) mode only)", 0l);
  //  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin (override autodetection from the reduced density matrix file name)");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output name for the entanglement spectrum (default name replace density-matrix full.ent extension with la_x_na_y.entspec)");
  (*SystemGroup) += new SingleDoubleOption  ('e', "eigenvalue-error", "lowest acceptable reduced density matrix eigenvalue", 1e-14);   (*SystemGroup) += new SingleDoubleOption  ('\n', "xi-error", "minus log of the lowest acceptable reduced density matrix eigenvalue (o if error control relies on the eigenvalue-error option)", 0);  
  (*SystemGroup) += new BooleanOption ('\n', "show-minmaxkya", "show minimum an maximum Ky value that can be reached");
  (*SystemGroup) += new BooleanOption ('\n', "show-counting", "show degeneracy counting for each Ky value");
  (*SystemGroup) += new BooleanOption ('\n', "particle-entanglement", "compute particle entanglement spectrum");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky-periodic", "set the periodicity for for the ky momentum (0 if non-periodic )", 0);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  //  int NbrOrbitalsInPartition = Manager.GetInteger("nbr-orbitals");
  int NbrParticlesInPartition = Manager.GetInteger("nbr-particles");
  //  int TotalSzInPartition = Manager.GetInteger("sz-value");
  double Error = Manager.GetDouble("eigenvalue-error");
  if (Manager.GetDouble("xi-error") != 0.0)
    {
      Error = exp(-Manager.GetDouble("xi-error"));
    }
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  bool Statistics = true;
  int Modulo = Manager.GetInteger("ky-periodic");
  
  if (Manager.GetString("density-matrix") == 0)
    {
      cout << "a reduced density matrix has to be provided, see man page for option syntax or type FQHELatticeEntanglementSpectrum -h" << endl;
      return -1;
    }

//   bool SU2SpinFlag = false;
//   if (strcasestr(Manager.GetString("density-matrix"), "_su2_") != 0)
//     {
//       SU2SpinFlag = true;
//     }
//   else
//     {
//       SU2SpinFlag = Manager.GetBoolean("su2-spin");      
//     }
  
  double Mass = 0.0;
  if (FQHEOnSquareLatticeFindSystemInfoFromFileName(Manager.GetString("density-matrix"),
							  NbrParticles, NbrSiteX, NbrSiteY, Statistics) == false)
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

  int MinKa = 1 << 30;
  int MaxKa = -MinKa; 
  int* KaValueArray = 0;


  if (Manager.GetBoolean("particle-entanglement") == false)
    {
    }
  else
    {
      if (DensityMatrix.GetNbrColumns() < 4)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      
      int* NaValues = DensityMatrix.GetAsIntegerArray(0);
      int* KxValues = DensityMatrix.GetAsIntegerArray(1);
      int* KyValues = DensityMatrix.GetAsIntegerArray(2);
      long Index = 0l;
      long MaxIndex = DensityMatrix.GetNbrLines();
      while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
	++Index;

      if (Index < MaxIndex)
	{
	  double* Coefficients = DensityMatrix.GetAsDoubleArray(3);
	  char* OutputFileName = Manager.GetString("output");
	  if (OutputFileName == 0)
	    {
	      char* TmpExtension = new char[256];
	      sprintf(TmpExtension, "na_%d.parentspec", NbrParticlesInPartition);
	      if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
		{
		  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent", TmpExtension);
		}
	      else
		{
		  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent.bz2", TmpExtension);
		}
	    }
	  ofstream File;
	  File.open(OutputFileName, ios::out);
	  File.precision(14);
	  File << "# na kx ky linearized_k lambda -log(lambda)";
	  File << endl;
	  int TmpIndex = Index;
	  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
	    {
	      double Tmp = Coefficients[Index];
	      if (Tmp > Error)
		{
		  int TmpKxa = KxValues[Index];
		  int TmpKya = KyValues[Index];
		  int TmpLinearizedK = TmpKxa + (TmpKya * NbrSiteX);
		  File << NbrParticlesInPartition << " " << TmpKxa << " " << TmpKya << " " << TmpLinearizedK << " " << Tmp << " " << (-log(Tmp));
		  File << endl;
 		  if (TmpLinearizedK < MinKa)
 		    MinKa = TmpLinearizedK;
 		  if (TmpLinearizedK > MaxKa)
 		    MaxKa = TmpLinearizedK;
		}
	      ++Index;
	    }
          if (MaxKa>=MinKa)
          {
              KaValueArray = new int[(MaxKa - MinKa + 1)];
              for (int i = MinKa; i <= MaxKa; ++i)
                KaValueArray[(i - MinKa)] = 0; 
              Index = TmpIndex;
              while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
                {
                  if (Coefficients[Index] > Error)
                    {
                      KaValueArray[(KxValues[Index] + (KyValues[Index] * NbrSiteX) - MinKa)]++; 
                    }
                  ++Index;
                }
          }
	  File.close();	      
	}
      else
	{
	  cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
	  return -1;
	}
    }

//   if (Manager.GetBoolean("show-minmaxkya"))
//     {
//       cout << "min Kya = " << MinKa << endl;
//       cout << "max Kya = " << MaxKa << endl;
//     }
  
  if ((Manager.GetBoolean("show-counting")) && (KaValueArray != 0))
     {
       long TotalDegenracy = 0l;
       for (int i = MinKa; i <= MaxKa; ++i)
 	{
 	  TotalDegenracy += KaValueArray[(i - MinKa)]; 
 	}
       cout << "total degeneracy counting " << TotalDegenracy << endl;
       cout << "degeneracy counting : " << endl;
       for (int i = MinKa; i <= MaxKa; ++i)
 	{
 	  cout << i << " (kx=" << (i % NbrSiteX) << ", ky=" << (i / NbrSiteX) << ") = "<< KaValueArray[(i - MinKa)] << endl; 
 	}
     }

  return 0;
}

