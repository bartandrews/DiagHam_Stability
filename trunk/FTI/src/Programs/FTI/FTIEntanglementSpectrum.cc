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
#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"

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
  OptionManager Manager ("FQHETopInsulatorEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "density-matrix", "file containing the reduced density matrix");  
  (*SystemGroup) += new SingleIntegerOption ('n', "nbr-particles", "number of particles in the A part", 0l);
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output name for the entanglement spectrum (default name replace density-matrix full.parent extension with la_x_na_y.entspec)");
  (*SystemGroup) += new SingleDoubleOption  ('e', "eigenvalue-error", "lowest acceptable reduced density matrix eigenvalue", 1e-14);   (*SystemGroup) += new SingleDoubleOption  ('\n', "xi-error", "minus log of the lowest acceptable reduced density matrix eigenvalue (o if error control relies on the eigenvalue-error option)", 0);  
  (*SystemGroup) += new BooleanOption ('\n', "show-minmaxkya", "show minimum an maximum Ky value that can be reached");
  (*SystemGroup) += new BooleanOption ('\n', "show-counting", "show degeneracy counting for each Ky value");
  (*SystemGroup) += new BooleanOption ('\n', "show-countingSz", "show degeneracy counting for each sz value (in decoupled mode)");
  (*SystemGroup) += new BooleanOption ('\n', "particle-entanglement", "compute particle entanglement spectrum");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky-periodic", "set the periodicity for for the ky momentum (0 if non-periodic )", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "3d", "consider a 3d model instead of a 2d model");
  (*SystemGroup) += new BooleanOption ('\n', "decoupled", "assume that the FTI states are made of two decoupled FCI copies");
  (*SystemGroup) += new BooleanOption  ('\n', "Wannier", "use Wannier basis");
  (*SystemGroup) += new BooleanOption  ('\n', "Wannier-block-diagonal", "use 'block diagonal' Wannier basis ");
  (*SystemGroup) += new BooleanOption  ('\n', "real-space", "consider a model written in real space basis");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHETopInsulatorEntanglementSpectrum -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

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
  int NbrSiteZ = 0;
  bool Statistics = true;
  int Modulo = Manager.GetInteger("ky-periodic");
  bool Flag3d = Manager.GetBoolean("3d");
  bool FlagDecoupled = Manager.GetBoolean("decoupled");
  
  if (Manager.GetString("density-matrix") == 0)
    {
      cout << "a reduced density matrix has to be provided, see man page for option syntax or type FQHETopInsulatorEntanglementSpectrum -h" << endl;
      return -1;
    }

  if (Flag3d == true)
    {
      if (FQHEOnCubicLatticeFindSystemInfoFromFileName(Manager.GetString("density-matrix"),
						       NbrParticles, NbrSiteX, NbrSiteY, NbrSiteZ, Statistics) == false)
	{
	  cout << "can't retrieve system informations from the reduced density matrix file name" << endl;
	  return -1;
	} 
    }
  else
    {
      NbrSiteZ = 1;
      double Mass = 0.0;
      if (Manager.GetBoolean("real-space") == false)
      {
	if (FQHEOnSquareLatticeFindSystemInfoFromFileName(Manager.GetString("density-matrix"),
							NbrParticles, NbrSiteX, NbrSiteY, Statistics) == false)
	  {
	    cout << "can't retrieve system informations from the reduced density matrix file name" << endl;
	    return -1;
	  }
      }
      else
      {
	int NbrSites;
	bool GutzwillerFlag;
	if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("density-matrix"), NbrParticles, NbrSites, Statistics, GutzwillerFlag) == false)
	{
	  cout << "can't retrieve system informations from the reduced density matrix file name" << endl;
	  return -1;
	}
      }
    }

  MultiColumnASCIIFile DensityMatrix;
  if (DensityMatrix.Parse(Manager.GetString("density-matrix")) == false)
    {
      DensityMatrix.DumpErrors(cout);
      return -1;
    }

  int MinKa = 1 << 30;
  int MaxKa = -MinKa; 
  int MinSza = Manager.GetInteger("nbr-particles") % 2;
  int MaxSza = Manager.GetInteger("nbr-particles") / 2;
  int* KaValueArray = 0;
  int* SzaValueArray = 0;


  if (Manager.GetBoolean("particle-entanglement") == false)
    {
    }
  else
    {
      if (DensityMatrix.GetNbrColumns() < 4 )
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      
      int* NaValues = DensityMatrix.GetAsIntegerArray(0);
      int*  KxValues = DensityMatrix.GetAsIntegerArray(1);
      int*  KyValues = DensityMatrix.GetAsIntegerArray(2);
      int* KzValues = 0;
      if (Flag3d == true)
	KzValues = DensityMatrix.GetAsIntegerArray(3);
      int* SzValues = 0;
      if (FlagDecoupled == true)
	SzValues = DensityMatrix.GetAsIntegerArray(3);
      long Index = 0l;
      long MaxIndex = DensityMatrix.GetNbrLines();
      while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
	++Index;

      if (Index < MaxIndex)
	{
	  double* Coefficients = 0;
	  if ((Flag3d == true) || (FlagDecoupled == true))
	    Coefficients = DensityMatrix.GetAsDoubleArray(4);
	  else
	    Coefficients = DensityMatrix.GetAsDoubleArray(3);
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
	      if (OutputFileName==0)
		{
		  cout << "Problem generating output file name"<<endl;
		  return -1;
		}
	    }
	  ofstream File;
	  File.open(OutputFileName, ios::out);
	  File.precision(14);
	  if (Flag3d == true)
	    File << "# na kx ky kz linearized_k lambda -log(lambda)";
	  else
	    {
	      if (FlagDecoupled == false)
		{
		  File << "# na kx ky linearized_k lambda -log(lambda)";
		}
	      else
		{
		  File << "# na kx ky linearized_k sz lambda -log(lambda)";
		}
	    }
	  File << endl;
	  int TmpIndex = Index;
	  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
	    {
	      double Tmp = Coefficients[Index];
	      if (Tmp > Error)
		{
		  int TmpKxa = KxValues[Index];
		  int TmpKya = KyValues[Index];
		  int TmpLinearizedK = 0;
		  if (Flag3d == true)
		    {
		      int TmpKza = KzValues[Index];
		      TmpLinearizedK = (TmpKxa + ((TmpKya + (TmpKza * NbrSiteY)) * NbrSiteX));
		      File << NbrParticlesInPartition << " " << TmpKxa << " " << TmpKya << " " << TmpKza << " " << TmpLinearizedK << " " << Tmp << " " << (-log(Tmp));
		    }
		  else
		    {
		      TmpLinearizedK = (TmpKxa + (TmpKya * NbrSiteX));
		      if (FlagDecoupled == false)
			{
			  File << NbrParticlesInPartition << " " << TmpKxa << " " << TmpKya << " " << TmpLinearizedK << " " << Tmp << " " << (-log(Tmp));
			}
		      else
			{
			  File << NbrParticlesInPartition << " " << TmpKxa << " " << TmpKya << " " << TmpLinearizedK << " " << SzValues[Index] << " " << Tmp << " " << (-log(Tmp));
			}
		    }
		  File << endl;
 		  if (TmpLinearizedK < MinKa)
 		    MinKa = TmpLinearizedK;
 		  if (TmpLinearizedK > MaxKa)
 		    MaxKa = TmpLinearizedK;
		}
	      ++Index;
	    }
          if (MaxKa >= MinKa)
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
	    	    
	    if (FlagDecoupled == true)
	    {
              SzaValueArray = new int[(NbrParticlesInPartition + 1)];
              for (int i = 0; i <= NbrParticlesInPartition ; ++i)
                SzaValueArray[i] = 0; 
              Index = TmpIndex;
              while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
                {
                  if (Coefficients[Index] > Error)
                    {
                      SzaValueArray[(SzValues[Index] + NbrParticlesInPartition) / 2]++; 
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
       if (Flag3d == false)
	 {
	   for (int i = MinKa; i <= MaxKa; ++i)
	     {
	       cout << i << " (kx=" << (i % NbrSiteX) << ", ky=" << (i / NbrSiteX) << ") = "<< KaValueArray[(i - MinKa)] << endl; 
	     }
	 }
       else
	 {
	   for (int i = MinKa; i <= MaxKa; ++i)
	     {
	       int TmpKx = i % (NbrSiteX * NbrSiteY);
	       int TmpKz = i / (NbrSiteX * NbrSiteY);
	       int TmpKy = (i / NbrSiteX) % NbrSiteY;
	       cout << i << " (kx=" << TmpKx << ", ky=" << TmpKy << ", kz=" 
		    << TmpKz << ") = "<< KaValueArray[(i - MinKa)] << endl; 
	     }
	 }
     }

     
     
     if ((Manager.GetBoolean("show-countingSz")) && ( Manager.GetBoolean("decoupled") == true))
     {
       long TotalDegenracy = 0l;
       for (int i = MinKa; i <= MaxKa; ++i)
 	{
 	  TotalDegenracy += KaValueArray[(i - MinKa)]; 
 	}
       cout << "total degeneracy counting " << TotalDegenracy << endl;
       cout << "degeneracy counting : " << endl;
       for (int i = 0; i <= NbrParticlesInPartition; ++i)
	     {
	       cout << "Sza = " << (2*i - NbrParticlesInPartition) << " : " << SzaValueArray[i] << endl; 
	     }
	
     }
  return 0;
}

