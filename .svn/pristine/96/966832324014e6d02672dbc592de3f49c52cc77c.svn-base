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
  OptionManager Manager ("FQHESphereEntanglementSpectrum" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('\0', "density-matrix", "file containing the reduced density matrix");  
  (*SystemGroup) += new SingleIntegerOption ('n', "nbr-particles", "number of particles in the A part", 0l);
  (*SystemGroup) += new SingleIntegerOption ('l', "nbr-orbitals", "number of orbitals in the A part", 0l);
  (*SystemGroup) += new SingleIntegerOption ('s', "sz-value", "twice the Sz value of A part (SU(2) mode only)", 0l);
  (*SystemGroup) += new BooleanOption  ('\n', "su2-spin", "consider particles with SU(2) spin (override autodetection from the reduced density matrix file name)");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output name for the entanglement spectrum (default name replace density-matrix full.ent extension with la_x_na_y.entspec)");
  (*SystemGroup) += new SingleDoubleOption  ('e', "eigenvalue-error", "lowest acceptable reduced density matrix eignvalue", 1e-14);  
  (*SystemGroup) += new BooleanOption ('\n', "show-minmaxlza", "show minimum an maximum Lz value that can be reached");
  (*SystemGroup) += new BooleanOption ('\n', "show-counting", "show degeneracy counting for each Lz value");
  (*SystemGroup) += new BooleanOption ('\n', "particle-entanglement", "compute particle entanglement spectrum");
  (*SystemGroup) += new BooleanOption ('\n', "discrete-symmetries", "the reduced density matrix spectrum includes discrete symmetry information such as Sz<->-Sz");
  (*SystemGroup) += new BooleanOption ('\n', "ls-sorted", "for the particle entanglement spectrum with su2-spin on the sphere, sort the spectrum with respect to L and S");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereEntanglementSpectrum -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrOrbitalsInPartition = Manager.GetInteger("nbr-orbitals");
  int NbrParticlesInPartition = Manager.GetInteger("nbr-particles");
  int TotalSzInPartition = Manager.GetInteger("sz-value");
  double Error = Manager.GetDouble("eigenvalue-error");
  int NbrParticles = 0;
  int NbrFluxQuanta = 0;
  int TotalLz = 0;
  int TotalSz = 0;
  bool Statistics = true;
  
  if (Manager.GetString("density-matrix") == 0)
    {
      cout << "a reduced density matrix has to be provided, see man page for option syntax or type FQHESphereEntanglementSpectrum -h" << endl;
      return -1;
    }

  bool SU2SpinFlag = false;
  if (strcasestr(Manager.GetString("density-matrix"), "_su2_") != 0)
    {
      SU2SpinFlag = true;
    }
  else
    {
      SU2SpinFlag = Manager.GetBoolean("su2-spin");      
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


  if (Manager.GetBoolean("particle-entanglement") == false)
    {
      if (DensityMatrix.GetNbrColumns() < 4)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      
      int* LaValues = DensityMatrix.GetAsIntegerArray(0);
      int* NaValues = DensityMatrix.GetAsIntegerArray(1);
      int* LzValues = DensityMatrix.GetAsIntegerArray(2);
      long Index = 0l;
      long MaxIndex = DensityMatrix.GetNbrLines();
      while ((Index < MaxIndex) && (LaValues[Index] != NbrOrbitalsInPartition))
	++Index;
      
      if (Index < MaxIndex)
	{
	  if (SU2SpinFlag == false)
	    {
	      if (Index < MaxIndex)
		{
		  double* Coefficients = DensityMatrix.GetAsDoubleArray(3);
		  char* OutputFileName = Manager.GetString("output");
		  if (OutputFileName == 0)
		    {
		      char* TmpExtension = new char[256];
		      sprintf(TmpExtension, "la_%d_na_%d.entspec", NbrOrbitalsInPartition, NbrParticlesInPartition);
		      if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent.bz2", TmpExtension);
			}
		    }
		  ofstream File;
		  File.open(OutputFileName, ios::out);
		  File.precision(14);
		  File << "# la na lz shifted_lz lambda -log(lambda)" << endl;
		  if (NbrParticlesInPartition == 0)
		    {
		      int TmpIndex = Index;
		      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
			{
			  double Tmp = Coefficients[Index];
			  if (Tmp > Error)
			    {
			      int TmpLza = (- (LzValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
			      File << NbrOrbitalsInPartition << " " << NaValues[Index] << " " << LzValues[Index] << " " <<  (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp)) << endl;
			      if (TmpLza < MinLza)
				MinLza = TmpLza;
			      if (TmpLza > MaxLza)
				MaxLza = TmpLza;
			    }
			  ++Index;
			}
		      LzaValueArray = new int[((MaxLza - MinLza) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			LzaValueArray[(i - MinLza) >> 1] = 0; 
		      Index = TmpIndex;
		      while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition))
			{
			  if (Coefficients[Index] > Error)
			    {
			      int TmpLza = (- (LzValues[Index] + ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NaValues[Index])));
			      LzaValueArray[(TmpLza - MinLza) >> 1]++; 
			    }
			  ++Index;
			}
		    }
		  else
		    {
		      while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
			++Index;
		      if (Index < MaxIndex)
			{
			  int TmpIndex = Index;
			  int Shift = ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NbrParticlesInPartition);
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition))
			    {
			      double Tmp = Coefficients[Index];
			      if (Tmp > Error)
				{
				  int TmpLza = (-(LzValues[Index] + Shift));
				  File << NbrOrbitalsInPartition << " " << NbrParticlesInPartition << " " << LzValues[Index] << " " << (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp)) << endl;
				  if (TmpLza < MinLza)
				    MinLza = TmpLza;
				  if (TmpLza > MaxLza)
				    MaxLza = TmpLza;
				}
			      ++Index;
			    }	      
			  LzaValueArray = new int[((MaxLza - MinLza) >> 1) + 1];
			  for (int i = MinLza; i <= MaxLza; i += 2)
			    LzaValueArray[(i - MinLza) >> 1] = 0; 
			  Index = TmpIndex;
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition))
			    {
			      if (Coefficients[Index] > Error)
				{
				  int TmpLza = (-(LzValues[Index] + Shift));
				  LzaValueArray[(TmpLza - MinLza) >> 1]++; 
				}
			      ++Index;
			    }
			}
		      else
			{
			  cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
			  return -1;
			}
		    }
		  File.close();
		}
	    }
	  else
	    {
	      if (((TotalSzInPartition ^ NbrParticlesInPartition) & 1) == 0)
		{
		  int* SzaValues = DensityMatrix.GetAsIntegerArray(3);
		  double* Coefficients = DensityMatrix.GetAsDoubleArray(4);
		  char* OutputFileName = Manager.GetString("output");
		  if (OutputFileName == 0)
		    {
		      char* TmpExtension = new char[256];
		      sprintf(TmpExtension, "la_%d_na_%d_sza_%d.entspec", NbrOrbitalsInPartition, NbrParticlesInPartition, TotalSzInPartition);
		      if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.ent.bz2", TmpExtension);
			}
		    }
		  ofstream File;
		  File.open(OutputFileName, ios::out);
		  File.precision(14);
		  File << "# la sza na lz shifted_lz lambda -log(lambda)" << endl;
		  
		  while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
		    ++Index;
		  if (Index < MaxIndex)
		    {
		      while ((Index < MaxIndex) && (SzaValues[Index] != TotalSzInPartition))
			++Index;
		      if (Index < MaxIndex)
			{
			  double Shift = ((NbrOrbitalsInPartition - 1 - NbrFluxQuanta) * NbrParticlesInPartition);
			  while ((Index < MaxIndex) && (LaValues[Index] == NbrOrbitalsInPartition) && (NaValues[Index] == NbrParticlesInPartition) && (SzaValues[Index] == TotalSzInPartition))
			    {
			      double Tmp = Coefficients[Index];
			      if (Tmp > Error)
				{
				  int TmpLza = (-(LzValues[Index] + Shift));
				  File << LaValues[Index] << " " << SzaValues[Index] << " " << NaValues[Index] << " " << LzValues[Index] << " " <<  (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp)) << endl;
				  if (TmpLza < MinLza)
				    MinLza = TmpLza;
				  if (TmpLza > MaxLza)
				    MaxLza = TmpLza;
				}
			      ++Index;
			    }	      
			}
		      else
			{
			  cout << "error, no entanglement spectrum can be computed from current data (invalid Sz value)" << endl;
			  return -1;
			}
		    }
		  else
		    {
		      cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
		      return -1;
		    }
		  File.close();
		}
	      else
		{
		  cout << "error, no entanglement spectrum can be computed from current data (invalid Sz value)" << endl;
		  return -1;
		}
	    }
	}
      else
	{
	  cout << "error, no entanglement spectrum can be computed from current data (invalid number of orbitals)" << endl;
	  return -1;
	}
      
    }
  else
    {
      if (DensityMatrix.GetNbrColumns() < 3)
	{
	  cout << "wrong number of columns in " << Manager.GetString("density-matrix") << endl;
	  return -1;
	}
      
      int* NaValues = DensityMatrix.GetAsIntegerArray(0);

      if (SU2SpinFlag == false)
	{
	  int* LzValues = DensityMatrix.GetAsIntegerArray(1);
	  double* LValues = 0;
	  double* L2Values = 0;
	  if (DensityMatrix.GetNbrColumns() > 3)
	    {
	      L2Values = DensityMatrix.GetAsDoubleArray(3);
	      LValues = DensityMatrix.GetAsDoubleArray(4);
	    }
	  long Index = 0l;
	  long MaxIndex = DensityMatrix.GetNbrLines();
	  while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
	    ++Index;
	  
	  if (Index < MaxIndex)
	    {
	      double* Coefficients = DensityMatrix.GetAsDoubleArray(2);
	      char* OutputFileName = Manager.GetString("output");
	      if (OutputFileName == 0)
		{
		  char* TmpExtension = new char[256];
		  if (strcasestr(Manager.GetString("density-matrix"), "realent") == 0)
		    {
		      sprintf(TmpExtension, "na_%d.parentspec", NbrParticlesInPartition);
		    }
		  else
		    {
		      sprintf(TmpExtension, "na_%d.realentspec", NbrParticlesInPartition);
		    }
		  if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
		    {
		      if (strcasestr(Manager.GetString("density-matrix"), "realent") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.realent", TmpExtension);
			}
		    }
		  else
		    {
		      if (strcasestr(Manager.GetString("density-matrix"), "realent") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent.bz2", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.realent.bz2", TmpExtension);
			}
		    }
		}
	      ofstream File;
	      File.open(OutputFileName, ios::out);
	      File.precision(14);
	      File << "# na lz lambda -log(lambda)";
	      if (LValues != 0)
		{
		  File << " L^2 L";
		}
	      File << endl;
	      int TmpIndex = Index;
	      while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		{
		  double Tmp = Coefficients[Index];
		  if (Tmp > Error)
		    {
		      int TmpLza = LzValues[Index];
		      File << NbrParticlesInPartition << " " << (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp));
		      if (LValues != 0)
			{
			  File << " " << L2Values[Index] << " " << LValues[Index];
			}
		      File << endl;
		      if (TmpLza < MinLza)
			MinLza = TmpLza;
		      if (TmpLza > MaxLza)
			MaxLza = TmpLza;
		    }
		  ++Index;
		}
	      LzaValueArray = new int[((MaxLza - MinLza ) >> 1) + 1];
	      for (int i = MinLza; i <= MaxLza; i += 2)
		LzaValueArray[(i - MinLza) >> 1] = 0; 
	      Index = TmpIndex;
	      while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		{
		  if (Coefficients[Index] > Error)
		    {
		      LzaValueArray[(LzValues[Index] - MinLza) >> 1]++; 
		    }
		  ++Index;
		}
	      File.close();	      
	    }
	  else
	    {
	      cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
	      return -1;
	    }
	}
      else
	{
	  int* LzValues = 0;
	  if (Manager.GetBoolean("discrete-symmetries") == false)
	    {
	      LzValues = DensityMatrix.GetAsIntegerArray(4);
	    }
	  else
	    {
	      LzValues = DensityMatrix.GetAsIntegerArray(5);
	    }
	  int* SzSymmetryValues = 0;
	  if (Manager.GetBoolean("discrete-symmetries") == true)
	    {
	      SzSymmetryValues = DensityMatrix.GetAsIntegerArray(2);
	    }
	  int* SzValues = DensityMatrix.GetAsIntegerArray(1);
	  int** SzaLzaValueArray = 0;
	  long Index = 0l;
	  long MaxIndex = DensityMatrix.GetNbrLines();
	  while ((Index < MaxIndex) && (NaValues[Index] != NbrParticlesInPartition))
	    ++Index;
	  
	  if (Index < MaxIndex)
	    {
	      double* Coefficients = 0;
	      if (Manager.GetBoolean("discrete-symmetries") == false)
		{
		  Coefficients = DensityMatrix.GetAsDoubleArray(5);
		}
	      else
		{
		  Coefficients = DensityMatrix.GetAsDoubleArray(6);
		}
	      char* OutputFileName = Manager.GetString("output");
	      if (OutputFileName == 0)
		{
		  char* TmpExtension = new char[256];
		  if (strcasestr(Manager.GetString("density-matrix"), "realent") == 0)
		    {
		      sprintf(TmpExtension, "na_%d.parentspec", NbrParticlesInPartition);
		    }
		  else
		    {
		      sprintf(TmpExtension, "na_%d.realentspec", NbrParticlesInPartition);
		    }
		  if (strcasestr(Manager.GetString("density-matrix"), "bz2") == 0)
		    {
		      if (strcasestr(Manager.GetString("density-matrix"), "realent") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.realent", TmpExtension);
			}
		    }
		  else
		    {
		      if (strcasestr(Manager.GetString("density-matrix"), "realent") == 0)
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.parent.bz2", TmpExtension);
			}
		      else
			{
			  OutputFileName = ReplaceExtensionToFileName(Manager.GetString("density-matrix"), "full.realent.bz2", TmpExtension);
			}
		    }
		}
	      ofstream File;
	      File.open(OutputFileName, ios::out);
	      File.precision(14);
	      if (Manager.GetBoolean("ls-sorted") == false)
		{
		  if (Manager.GetBoolean("discrete-symmetries") == false)
		    {
		      File << "# na sz lz lambda -log(lambda)";
		    }
		  else
		    {
		      File << "# na sz sz<->-sz lz lambda -log(lambda)";
		    }
		  File << endl;
		  int TmpIndex = Index;
		  if (Manager.GetBoolean("discrete-symmetries") == false)
		    {
		      while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
			{
			  double Tmp = Coefficients[Index];
			  if (Tmp > Error)
			    {
			      int TmpLza = LzValues[Index];
			      File << NbrParticlesInPartition << " " << (0.5 * SzValues[Index]) << " " << (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp));
			      File << endl;
			      if (TmpLza < MinLza)
				MinLza = TmpLza;
			      if (TmpLza > MaxLza)
				MaxLza = TmpLza;
			    }
			  ++Index;
			}
		    }
		  else
		    {
		      while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
			{
			  double Tmp = Coefficients[Index];
			  if (Tmp > Error)
			    {
			      int TmpLza = LzValues[Index];
			      File << NbrParticlesInPartition << " " << (0.5 * SzValues[Index]) << " " << SzSymmetryValues[Index] << " " 
				   << (0.5 * TmpLza) << " " << Tmp << " " << (-log(Tmp));
			      File << endl;
			      if (TmpLza < MinLza)
				MinLza = TmpLza;
			      if (TmpLza > MaxLza)
				MaxLza = TmpLza;
			    }
			  ++Index;
			}
		    }
		  SzaLzaValueArray = new int* [NbrParticlesInPartition + 1];
		  for (int TmpSza = 0; TmpSza <= NbrParticlesInPartition; ++TmpSza)
		    {
		      SzaLzaValueArray[TmpSza] = new int[((MaxLza - MinLza ) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			SzaLzaValueArray[TmpSza][(i - MinLza) >> 1] = 0; 
		    }
		  Index = TmpIndex;
		  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		    {
		      if (Coefficients[Index] > Error)
			{
			  SzaLzaValueArray[(SzValues[Index] + NbrParticlesInPartition) >> 1][(LzValues[Index] - MinLza) >> 1]++; 
			}
		      ++Index;
		    }
		  if (Manager.GetBoolean("show-minmaxlza"))
		    {
		      cout << "min Lza = " << MinLza << endl;
		      cout << "max Lza = " << MaxLza << endl;
		    }
		  
		  if ((Manager.GetBoolean("show-counting")) && (SzaLzaValueArray != 0))
		    {
		      cout << "degeneracy counting (SzA Lza NbrStates) : " << endl;
		      for (int j = 0; j <= NbrParticlesInPartition; ++j)
			{
			  for (int i = MinLza; i <= MaxLza; i += 2)
			    {
			      cout << (j - (NbrParticlesInPartition * 0.5)) << " " << (0.5 * i) << " " << SzaLzaValueArray[j][(i - MinLza) >> 1] << endl; 
			    }
			}
		    }	      
		}
	      else
		{
		  File << "# na s l lambda -log(lambda)";
		  File << endl;
		  int MinSza = 1 << 30;
		  int MaxSza = -MinSza; 
		  int TmpIndex = Index;
		  int LargestSzSector = SzValues[Index];
		  int LargestLzSector = LzValues[Index];
		  while ((Index < MaxIndex) && (NaValues[Index] == NbrParticlesInPartition))
		    {
		      if (abs(LzValues[Index]) < abs(LargestLzSector))
			LargestLzSector = LzValues[Index];
		      if (abs(SzValues[Index]) < abs(LargestSzSector))
			LargestSzSector = SzValues[Index];
		      double Tmp = Coefficients[Index];
		      if (Tmp > Error)
			{
			  int TmpLza = LzValues[Index];
			  if (TmpLza < MinLza)
			    MinLza = TmpLza;
			  if (TmpLza > MaxLza)
			    MaxLza = TmpLza;
			  int TmpSza = SzValues[Index];
			  if (TmpSza < MinSza)
			    MinSza = TmpSza;
			  if (TmpSza > MaxSza)
			    MaxSza = TmpSza;
			}
		      ++Index;
		    }
		  MaxIndex = Index;
		  Index = TmpIndex;
		  int** NbrSzLzSortedSpectrum = new int* [((MaxSza - MinSza) >> 1) + 1];
		  for (int j = MinSza; j <= MaxSza; j += 2)
		    {
		      NbrSzLzSortedSpectrum[(j - MinSza) >> 1] = new int[((MaxLza - MinLza ) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			NbrSzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1] = 0; 
		    }
		  while (Index < MaxIndex)
		    {
		      if (Coefficients[Index] > Error)
			{
			  NbrSzLzSortedSpectrum[(SzValues[Index] - MinSza) >> 1][(LzValues[Index] - MinLza) >> 1]++; 		       
			}
		      ++Index;
		    }
		  double*** SzLzSortedSpectrum = new double** [((MaxSza - MinSza) >> 1) + 1];
		  for (int j = MinSza; j <= MaxSza; j += 2)
		    {
		      SzLzSortedSpectrum[(j - MinSza) >> 1] = new double*[((MaxLza - MinLza ) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			{
			  int& TmpIndex2 = NbrSzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1];
			  if (TmpIndex2 > 0)
			    {
			      SzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1] = new double [TmpIndex2];
			      TmpIndex2 = 0; 
			    }
			}
		    }
		  Index = TmpIndex;
		  while (Index < MaxIndex)
		    {
		      if (Coefficients[Index] > Error)
			{
			  int& TmpIndex2 = NbrSzLzSortedSpectrum[(SzValues[Index] - MinSza) >> 1][(LzValues[Index] - MinLza) >> 1];
			  SzLzSortedSpectrum[(SzValues[Index] - MinSza) >> 1][(LzValues[Index] - MinLza) >> 1][TmpIndex2] = Coefficients[Index];
			  ++TmpIndex2; 		       
			}
		      ++Index;
		    }

		  for (int j = MaxSza; j >= 0; j -= 2)
		    {
		      for (int i = MaxLza; i >= ((MaxLza & 1) + 1); i -= 2)
			{
			  double* Spectrum1  = SzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1];
			  int NbrValues1 = NbrSzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1];
			  for (int k = 0; k < NbrValues1; ++k)
			    {
			      double TmpValue = Spectrum1[k];
			      if (TmpValue >= 0.0)
				{
				  for (int i2 = i - 2; i2 >= (MaxLza & 1); i2 -= 2)
				    {
				      double* Spectrum2  = SzLzSortedSpectrum[(j - MinSza) >> 1][(i2 - MinLza) >> 1];
				      int NbrValues2 = NbrSzLzSortedSpectrum[(j - MinSza) >> 1][(i2 - MinLza) >> 1];
				      double MinError = 1.0;
				      int TmpPosition = -1;
				      for (int l = 0; l < NbrValues2; ++l)
					{
					  if ((Spectrum2[l] >= 0.0) && (fabs(Spectrum2[l] - TmpValue) < MinError))
					    {
					      MinError = fabs(Spectrum2[l] - TmpValue);
					      TmpPosition = l;
					    }
					}		
				      if (TmpPosition >= 0)	
					Spectrum2[TmpPosition] = -1.0;
				    }
				}
			    }
			}
		    }

		  for (int i = MaxLza; i >= 0; i -= 2)
		    {
		      for (int j = MaxSza; j >= ((MaxSza & 1) + 1); j -= 2)
			{
			  double* Spectrum1  = SzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1];
			  int NbrValues1 = NbrSzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1];
			  for (int k = 0; k < NbrValues1; ++k)
			    {
			      double TmpValue = Spectrum1[k];
			      if (TmpValue >= 0.0)
				{
				  for (int j2 = j - 2; j2 >= (MaxSza & 1); j2 -= 2)
				    {
				      double* Spectrum2  = SzLzSortedSpectrum[(j2 - MinSza) >> 1][(i - MinLza) >> 1];
				      int NbrValues2 = NbrSzLzSortedSpectrum[(j2 - MinSza) >> 1][(i - MinLza) >> 1];
				      double MinError = 1.0;
				      int TmpPosition = -1;
				      for (int l = 0; l < NbrValues2; ++l)
					{
					  if ((Spectrum2[l] >= 0.0) && (fabs(Spectrum2[l] - TmpValue) < MinError))
					    {
					      MinError = fabs(Spectrum2[l] - TmpValue);
					      TmpPosition = l;
					    }
					}		
				      if (TmpPosition >= 0)	
					Spectrum2[TmpPosition] = -1.0;
				    }
				}
			    }
			}
		    }

		  SzaLzaValueArray = new int* [NbrParticlesInPartition + 1];
		  for (int TmpSza = 0; TmpSza <= NbrParticlesInPartition; ++TmpSza)
		    {
		      SzaLzaValueArray[TmpSza] = new int[((MaxLza - MinLza ) >> 1) + 1];
		      for (int i = MinLza; i <= MaxLza; i += 2)
			SzaLzaValueArray[TmpSza][(i - MinLza) >> 1] = 0; 
		    }

		  for (int j = MaxSza; j >= 0; j -= 2)
		    {
		      for (int i = MaxLza; i >= (MaxLza & 1); i -= 2)
			{
			  double* Spectrum  = SzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1];
			  int NbrValues = NbrSzLzSortedSpectrum[(j - MinSza) >> 1][(i - MinLza) >> 1];
			  for (int k = 0; k < NbrValues; ++k)
			    {
			      double Tmp = Spectrum[k];
			      if (Tmp >= 0.0)
				{
				  SzaLzaValueArray[(j - MinSza) >> 1][(i - MinLza) >> 1]++;
				  File << NbrParticlesInPartition << " " << (0.5 * j) << " " << (0.5 * i) << " " << Tmp << " " << (-log(Tmp));
				  File << endl;
				}
			    }
			}
		    }
		
		  if (Manager.GetBoolean("show-minmaxlza"))
		    {
		      cout << "min Lza = " << MinLza << endl;
		      cout << "max Lza = " << MaxLza << endl;
		    }
		  
		  if ((Manager.GetBoolean("show-counting")) && (SzaLzaValueArray != 0))
		    {
		      cout << "degeneracy counting (Sa La NbrStates) : " << endl;
		      for (int j = MaxSza; j >= 0; j -= 2)
			{
			  for (int i = MaxLza; i >= 0; i -= 2)
			    {
			      if (SzaLzaValueArray[(j - MinSza) >> 1][(i - MinLza) >> 1] > 0)
				cout << (j * 0.5) << " " << (0.5 * i) << " " << SzaLzaValueArray[(j - MinSza) >> 1][(i - MinLza) >> 1] << endl; 
			    }
			}
		    }		  
		}
	      File.close();	      

	      return 0;	      
	    }
	  else
	    {
	      cout << "error, no entanglement spectrum can be computed from current data (invalid number of particles)" << endl;	      
	      return -1;
	    }
	}
    }

  if (Manager.GetBoolean("show-minmaxlza"))
    {
      cout << "min Lza = " << MinLza << endl;
      cout << "max Lza = " << MaxLza << endl;
    }
  
  if ((Manager.GetBoolean("show-counting")) && (LzaValueArray != 0))
    {
      cout << "degeneracy counting : " << endl;
      for (int i = MinLza; i <= MaxLza; i += 2)
	{
	  cout << (0.5 * i) << " " << LzaValueArray[(i - MinLza) >> 1] << endl; 
	}
    }

  return 0;
}

