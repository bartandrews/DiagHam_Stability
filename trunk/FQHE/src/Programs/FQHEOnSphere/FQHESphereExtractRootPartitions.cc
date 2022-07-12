#include "config.h"

#include "Vector/RealVector.h"
#include "Vector/LongRationalVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"
#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/ParticleOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"
#include "HilbertSpace/FermionOnSphereWithSU3Spin.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/BosonOnSphereWithSpin.h"
#include "HilbertSpace/BosonOn4DSphere.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


// resuffle an array of vectors such that the first vector is the one with the first non-zero component having the smallest index
//
// vectors = array of vectors
// nbrVectors = numerb of vectors
// componentError = error threshold to assume a component is zero
// return value = index of the first non-zero component for the first vector 
int ReshuffleVectors (RealVector* vectors, int nbrVectors, double componentError);

// resuffle an array of vectors such that the first vector is the one with the first non-zero component having the smallest index
//
// vectors = array of vectors
// nbrVectors = numerb of vectors
// return value = index of the first non-zero component for the first vector 
int ReshuffleVectors (LongRationalVector* vectors, int nbrVectors);

// get the quantum numbers Jz, Kz of a N-particle 4D state 
//
//totalJz = array that gives the sector Jz for all values of the linearized index
//totalKz = array that gives the sector Kz for all values of the linearized index
void GetSectorsJzKzFromLinearizedIndex(int nbrParticles, int nbrFluxQuanta, int* totalJz, int* totalKz);

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereExtractRootPartitions" , "0.01");
  OptionGroup* MainGroup = new OptionGroup ("main options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += MiscGroup;
  Manager += MainGroup;
 
  (*MainGroup) += new SingleStringOption  ('\0', "input-file", "name of the file which contains the spectrum");
  (*MainGroup) += new SingleStringOption  ('\n', "eigenstate-list", "use a list of eigenstates instead of extracting the information from the spectrum");
  (*MainGroup) += new SingleDoubleOption  ('\n', "energy-error", "energy below which a state is considered as a zero energy states", 1e-10);
  (*MainGroup) += new SingleDoubleOption  ('\n', "component-error", "value below which norm of a vector component is assumed to be zero", 1e-8);
  (*MainGroup) += new SingleIntegerOption  ('\n', "output-precision", "numerical display precision", 14, true, 2, true, 14);
  (*MainGroup) += new SingleStringOption  ('o', "output", "output name for the root partition list (default name uses input-file, replacing dat extension with root)");
  (*MainGroup) += new BooleanOption  ('\n', "save-states", "save the states corresponding to each root configuration");
  (*MainGroup) += new BooleanOption  ('\n', "rational" , "use rational numbers instead of double precision floating point numbers");
  (*MainGroup) += new BooleanOption  ('\n', "4-D" , "use a Hilbert space with SO(5) symmetry (only available in eigenstate-list mode)");
  (*MainGroup) += new BooleanOption  ('\n', "PES-eigenstate" , "use a list of reduced density matrix eigenstates instead of energy eigenstates");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereExtractRootPartitions -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("input-file") == 0) && (Manager.GetString("eigenstate-list") == 0))
    {
      cout << "no input spectrum nor eigenstate list" << endl << "see man page for option syntax or type FQHESphereExtractRootPartitions -h" << endl;
      return -1;
    }

  double EnergyError = Manager.GetDouble("energy-error");
  double VectorError = Manager.GetDouble("component-error");
  int NbrFluxQuanta = 0;
  int NbrParticles = 0;
  int NbrTotParticles = 0;
  bool FermionFlag = true;
  bool SU2SpinFlag = false;
  bool SU3SpinFlag = false;
  bool SU4SpinFlag = false;
  bool FourDFlag = Manager.GetBoolean("4-D");
  bool PESFlag = Manager.GetBoolean("PES-eigenstate");
  int TotalSz = 0;
  int TotalTz = 0;
  int TotalY = 0;
  int TotalIz = 0;
  int TotalPz = 0;
  int TotalJz = 0;
  int TotalKz = 0;
  int TotalMaxLz = 0;
  int MaxLzValue = 0;
  int MaxJzValue = 0;
  int MaxKzValue = 0;
  int MinLzValue = 0;
  int* LzDegeneracy = 0;
  int* totalJz = 0;
  int* totalKz = 0;
  long TotalNbrZeroEnergyStates = 0l;
  bool SaveStates = Manager.GetBoolean("save-states");

  if (Manager.GetString("input-file") != 0)
    {
      if (FQHEOnSphereFindSystemInfoFromFileName(Manager.GetString("input-file"), NbrParticles, NbrFluxQuanta, FermionFlag) == false)
	{
	  cout << "can't retrieve system information form file name " << Manager.GetString("input-file") << endl;
	  return -1;
	}

      FQHEOnSphereFindInternalSymmetryGroupFromFileName(Manager.GetString("input-file"), SU2SpinFlag, SU3SpinFlag, SU4SpinFlag);    
      FQHEOnSphereLzSortedSpectrum Spectrum (Manager.GetString("input-file"), EnergyError);
      if (Spectrum.IsSpectrumValid() == false)
	{
	  cout << "Spectrum " << Manager.GetString("input-file") << " is unreadable or is not valid" << endl;
	  return -1;           
	}
      TotalMaxLz = Spectrum.GetMaxLzValue();
      if ((TotalMaxLz & 1) != 0)
	MaxLzValue = 1;
      LzDegeneracy = new int[TotalMaxLz + 1];
      while ((MaxLzValue <= TotalMaxLz) && (fabs(Spectrum.GetEnergy(MaxLzValue, 0)) < EnergyError))
	{
	  LzDegeneracy[MaxLzValue / 2] = Spectrum.GetDegeneracy(MaxLzValue, 0);
	  if (MaxLzValue == 0)
	    TotalNbrZeroEnergyStates += LzDegeneracy[MaxLzValue / 2];
	  else
	    TotalNbrZeroEnergyStates += 2l * LzDegeneracy[MaxLzValue / 2];
	  MaxLzValue += 2;
	}
      MaxLzValue -= 2;
      MinLzValue = (MaxLzValue & 1);
    }
  else
    {
      MultiColumnASCIIFile EigenstateListFile;
      if (EigenstateListFile.Parse(Manager.GetString("eigenstate-list")) == false)
	{
	  EigenstateListFile.DumpErrors(cout);
	  return -1;
	}
      if (FourDFlag == false)
	{
	  if (FQHEOnSphereFindSystemInfoFromVectorFileName(EigenstateListFile(0,0), NbrParticles, NbrFluxQuanta, MaxLzValue, FermionFlag) == false)  
	    {
	      cout << "can't retrieve system information from file name " << EigenstateListFile(0,0) << endl;
	      return -1;
	    }
	  
	  FQHEOnSphereFindInternalSymmetryGroupFromFileName(EigenstateListFile(0,0), SU2SpinFlag, SU3SpinFlag, SU4SpinFlag);          
	  MinLzValue = MaxLzValue;
	  TotalMaxLz = (MaxLzValue - MinLzValue) >> 1;
	  LzDegeneracy = new int[TotalMaxLz + 1];
	  TotalNbrZeroEnergyStates = EigenstateListFile.GetNbrLines();
	  for (int i = 0 ; i < TotalMaxLz; ++i)
	    LzDegeneracy[i] = 0;
	  LzDegeneracy[TotalMaxLz] = TotalNbrZeroEnergyStates;
	}
      else
	{
	  if (PESFlag == false)
	    {
	      if (FQHEOn4DSphereFindSystemInfoFromVectorFileName(EigenstateListFile(0,0),NbrParticles, NbrFluxQuanta, MaxJzValue, MaxKzValue,FermionFlag) == false)
		{
		  cout << "can't retrieve system information from 4d file name " << EigenstateListFile(0,0) << endl;
		  return -1;
		}
	    }
	  else
	    {
	      if (FQHEOn4DSphereFindSystemInfoFromPESVectorFileName(EigenstateListFile(0,0),NbrTotParticles,NbrParticles, NbrFluxQuanta, MaxJzValue, MaxKzValue,FermionFlag) == false)
		{
		  cout << "can't retrieve system information from 4d file name " << EigenstateListFile(0,0) << endl;
		  return -1;
		}
	    }
	  int TotalNbrSectors = 2*NbrParticles*NbrFluxQuanta*(NbrParticles*NbrFluxQuanta + 1) + 1;
	  // 	cout << TotalNbrSectors << endl;
	  totalJz = new int[TotalNbrSectors];
	  totalKz = new int[TotalNbrSectors];
	  GetSectorsJzKzFromLinearizedIndex(NbrParticles, NbrFluxQuanta, totalJz, totalKz);
	  if (MaxJzValue >= 0) 
	    MaxLzValue = MaxKzValue + NbrParticles*NbrFluxQuanta + 2*NbrParticles*NbrFluxQuanta*MaxJzValue - MaxJzValue*(MaxJzValue - 1);
	  else
	    MaxLzValue = -(MaxKzValue + NbrParticles*NbrFluxQuanta + MaxJzValue + (2*NbrParticles*NbrFluxQuanta + MaxJzValue + 1)*(-MaxJzValue - 1)) - 1;
	  // 	MaxLzValue = MaxLzValue + NbrParticles*NbrFluxQuanta*NbrParticles*NbrFluxQuanta ;
	  MinLzValue = MaxLzValue;
	  TotalMaxLz = (MaxLzValue - MinLzValue) >> 1;
	  LzDegeneracy = new int[TotalMaxLz + 1];
	  TotalNbrZeroEnergyStates = EigenstateListFile.GetNbrLines();
	  for (int i = 0 ; i < TotalMaxLz; ++i)
	    LzDegeneracy[i] = 0;
	  LzDegeneracy[TotalMaxLz] = TotalNbrZeroEnergyStates;
	}
    }
  
  if (TotalNbrZeroEnergyStates == 0l)
    {
      cout << "Spectrum " << Manager.GetString("input-file") << " does not contain any zero energy state" << endl;      
      return -1;
    }

  char* OutputFileName = Manager.GetString("output");
  if (OutputFileName == 0)
    {
      if (Manager.GetString("input-file") != 0)
	OutputFileName = ReplaceExtensionToFileName (Manager.GetString("input-file"), "dat" , "root");
      else
	OutputFileName = ReplaceExtensionToFileName (Manager.GetString("eigenstate-list"), "dat" , "root");
    }
  ofstream File;
  File.open(OutputFileName, ios::out);
  cout << "total nbr of root partitions : " << TotalNbrZeroEnergyStates << endl;
  cout << "nbr of root partitions per Lz sector : " << endl;;
  File << "total nbr of root partitions : " << TotalNbrZeroEnergyStates << endl;
  File << "nbr of root partitions per Lz sector : " << endl;

  for (int i = MinLzValue; i <= MaxLzValue; i += 2)
    {
      cout << LzDegeneracy[(i - MinLzValue) >> 1] << " ";
      File << LzDegeneracy[(i - MinLzValue) >> 1] << " ";
    }
  File << endl << endl << "-----------------------------------------" << endl;
  cout << endl << endl << "-----------------------------------------" << endl;

  char* BaseFileName = 0;
  char* TmpFileName = 0;
  if (Manager.GetString("input-file") != 0)
    {
      BaseFileName = RemoveExtensionFromFileName(Manager.GetString("input-file"), "dat");
      TmpFileName = new char [strlen(BaseFileName) + 256]; 
      BaseFileName[strlen(BaseFileName) - 1] = '\0';
    }


  for (int TotalLz = MinLzValue; TotalLz <= MaxLzValue; TotalLz += 2)
    {
      if (FourDFlag == false)
	{
	  if ((TotalMaxLz & 1) != 0)
	    {
	      File << "Lz = " << TotalLz << "/2 : " << endl;
	    }
	  else
	    {
	      File << "Lz = " << (TotalLz  >> 1) << " : " << endl;
	    }
	}
      else
	{
	  int Shift = NbrParticles*NbrFluxQuanta*NbrParticles*NbrFluxQuanta ;
	  int TotalMaxJz = totalJz[TotalMaxLz + Shift];
	  int TotalMaxKz = totalKz[TotalMaxLz + Shift];
	  TotalJz = totalJz[TotalLz + Shift];
	  TotalKz = totalKz[TotalLz + Shift];
	  cout << TotalJz << "  " << TotalKz << endl;
	  if ((TotalMaxJz & 1) != 0)
	    {
	      File << "Jz = " << TotalJz << "/2 , ";
	      if ((TotalMaxKz & 1) != 0)
		{
		  File << "Kz = " << TotalKz << "/2 : " << endl;
		}
	      else 
		{
		  File << "Kz = " << (TotalKz  >> 1) << " : " << endl;
		}
	    }
	  else
	    {
	      File << "Jz = " << (TotalJz  >> 1) << " , ";
	      if ((TotalMaxKz & 1) != 0)
		{
		  File << "Kz = " << TotalKz << "/2 : " << endl;
		}
	      else 
		{
		  File << "Kz = " << (TotalKz  >> 1) << " : " << endl;
		}
	    }
	  
	}
      int TmpNbrStates = LzDegeneracy[(TotalLz  - MinLzValue) >> 1];
      RealVector* TmpVectors = 0;
      LongRationalVector* TmpRationalVectors = 0;
      if (Manager.GetBoolean("rational"))
	{
	  TmpRationalVectors = new LongRationalVector[TmpNbrStates];
	}
      else
	{
	  TmpVectors = new RealVector[TmpNbrStates];
	}
      if (Manager.GetString("input-file") != 0)
	{
	  for (int j = 0; j < TmpNbrStates; ++j)
	    {
	      sprintf (TmpFileName, "%s_%d.%d.vec", BaseFileName, TotalLz, j);
	      if (TmpRationalVectors == 0)
		{
		  if (TmpVectors[j].ReadVector(TmpFileName) == false)
		    {
		      cout << "error while reading " << TmpFileName << endl;
		      return -1;
		    }
		}
	      else
		{
		  if (TmpRationalVectors[j].ReadVector(TmpFileName) == false)
		    {
		      cout << "error while reading " << TmpFileName << endl;
		      return -1;
		    }
		}
	    }
	}
      else
	{
	  MultiColumnASCIIFile EigenstateListFile;
	  if (EigenstateListFile.Parse(Manager.GetString("eigenstate-list")) == false)
	    {
	      EigenstateListFile.DumpErrors(cout);
	      return -1;
	    }
	  for (int j = 0; j < TmpNbrStates; ++j)
	    {
	      if (TmpRationalVectors == 0)
		{
		  if (TmpVectors[j].ReadVector(EigenstateListFile(0, j)) == false)
		    {
		      cout << "error while reading " << EigenstateListFile(0, j) << endl;
		      return -1;
		    }
		}
	      else
		{
		  if (TmpRationalVectors[j].ReadVector(EigenstateListFile(0, j)) == false)
		    {
		      cout << "error while reading " << EigenstateListFile(0, j) << endl;
		      return -1;
		    }
		}
	    }
	}
      int* RootPositions = new int [TmpNbrStates];
      for (int j = 0; j < TmpNbrStates; ++j)
	RootPositions[j] = 0;

      int MinTmpPos = 0;
      if (TmpRationalVectors == 0)
	{
	  MinTmpPos = ReshuffleVectors(TmpVectors, TmpNbrStates, VectorError);
	}
      else
	{
	  MinTmpPos = ReshuffleVectors(TmpRationalVectors, TmpNbrStates);	  
	}
      RootPositions[0] = MinTmpPos;

      if (TmpRationalVectors == 0)
	{
	  for (int k = 1; k < TmpNbrStates; ++k)
	    {
	      for (int j = k; j < TmpNbrStates; ++j)
		TmpVectors[j].AddLinearCombination(-TmpVectors[j][MinTmpPos] / TmpVectors[k - 1][MinTmpPos], TmpVectors[k -1]);
	      MinTmpPos = ReshuffleVectors(TmpVectors + k, TmpNbrStates - k, VectorError);
	      RootPositions[k] = MinTmpPos;
 	      for (int j = 0; j < k; ++j)
 		{
 		  TmpVectors[j].AddLinearCombination(-TmpVectors[j][MinTmpPos] / TmpVectors[k][MinTmpPos], TmpVectors[k]);
 		}
	    }
	}
      else
	{
	  if (TmpRationalVectors[0].IsNullVector() == false)
	    {
	      LongRational Tmp = TmpRationalVectors[0][MinTmpPos];
	      TmpRationalVectors[0] /= Tmp;
	      for (int k = 1; k < TmpNbrStates; ++k)
		{
		  for (int j = k; j < TmpNbrStates; ++j)
		    {
		      TmpRationalVectors[j].AddLinearCombination(-TmpRationalVectors[j][MinTmpPos], TmpRationalVectors[k -1]);
		    }
		  MinTmpPos = ReshuffleVectors(TmpRationalVectors + k, TmpNbrStates - k);
		  RootPositions[k] = MinTmpPos;
		  if (TmpRationalVectors[k].IsNullVector() == false)
		    {
		      Tmp = TmpRationalVectors[k][MinTmpPos];
		      TmpRationalVectors[k] /= Tmp;
		      for (int j = 0; j < k; ++j)
			{
			  TmpRationalVectors[j].AddLinearCombination(-TmpRationalVectors[j][MinTmpPos], TmpRationalVectors[k]);
			}
		    }
		}
	    }
	}
      
      ParticleOnSphere* Space = 0;
      char* FilePrefix = new char [256];
      char* FileSuffix = new char [256];
      if (FermionFlag == false)
	{
	  if (SU2SpinFlag == false)
	    {
	      if (FourDFlag == false)
	      {
		Space = new BosonOnSphereShort(NbrParticles, TotalLz, NbrFluxQuanta);
		if (TmpRationalVectors == 0)
		  {
		    sprintf (FilePrefix, "bosons_");
		  }
		else
		  {
		    sprintf (FilePrefix, "bosons_rational_");
		  }
		sprintf (FileSuffix, "n_%d_2s_%d_lz_%d.", NbrParticles, NbrFluxQuanta, TotalLz);
	      }
	      else
	      {
		Space = new BosonOn4DSphere(NbrParticles, NbrFluxQuanta, TotalJz, TotalKz);
		sprintf (FilePrefix, "bosons_sphere4d_");
		sprintf (FileSuffix, "n_%d_2s_%d_jz_%d_kz_%d.", NbrParticles, NbrFluxQuanta, TotalJz, TotalKz);
	      }
	    }
	  else
	    {
	      Space = new BosonOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
	      sprintf (FilePrefix, "bosons_sphere_su2_");
	      sprintf (FileSuffix, "n_%d_2s_%d_sz_%d_lz_%d.", NbrParticles, NbrFluxQuanta, TotalSz, TotalLz);
	    }
	}
      else
	{
	  if ((SU2SpinFlag == false) && (SU3SpinFlag == false) && (SU4SpinFlag == false))
	    {
#ifdef __64_BITS__
	      if (NbrFluxQuanta <= 63)
		Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#else
	      if (NbrFluxQuanta <= 31)
		Space = new FermionOnSphere(NbrParticles, TotalLz, NbrFluxQuanta);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, NbrFluxQuanta);
#endif
	      if (TmpRationalVectors == 0)
		{
		  sprintf (FilePrefix, "fermions_");
		}
	      else
		{
		  sprintf (FilePrefix, "fermions_rational_");
		}
	      sprintf (FileSuffix, "n_%d_2s_%d_lz_%d.", NbrParticles, NbrFluxQuanta, TotalLz);
	    }
	  else
	    {
	      if (SU2SpinFlag == true)
		{
		  Space = new FermionOnSphereWithSpin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz);
		  if (TmpRationalVectors == 0)
		    {
		      sprintf (FilePrefix, "fermions_sphere_su2_");
		    }
		  else
		    {
		      sprintf (FilePrefix, "fermions_rational_sphere_su2_");
		    }
		  sprintf (FileSuffix, "n_%d_2s_%d_sz_%d_lz_%d.", NbrParticles, NbrFluxQuanta, TotalSz, TotalLz);
		}
	      else 
		{
		  if (SU3SpinFlag == true)
		    {
		      Space = new FermionOnSphereWithSU3Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalTz, TotalY);
		      if (TmpRationalVectors == 0)
			{
			  sprintf (FilePrefix, "fermions_sphere_su3_");
			}
		      else
			{
			  sprintf (FilePrefix, "fermions_rational_sphere_su3_");
			}
		      sprintf (FileSuffix, "n_%d_2s_%d_tz_%d_y_%d_lz_%d.", NbrParticles, NbrFluxQuanta, TotalTz, TotalY, TotalLz);
		    }
		  else
		    {
		      if (SU4SpinFlag == true)
			{
			  Space = new FermionOnSphereWithSU4Spin(NbrParticles, TotalLz, NbrFluxQuanta, TotalSz, TotalIz, TotalPz);	    
			  if (TmpRationalVectors == 0)
			    {
			      sprintf (FilePrefix, "fermions_sphere_su4_");
			    }
			  else
			    {
			      sprintf (FilePrefix, "fermions_rational_sphere_su4_");
			    }
			  sprintf (FileSuffix, "n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz_%d.", NbrParticles, NbrFluxQuanta, TotalSz, TotalIz, TotalPz, TotalLz);
			}
		    }
		}
	    }
	}
      
      int NbrSavedStates = 0;
      for (int j = 0; j < TmpNbrStates; ++j)
	{
	  if (TmpRationalVectors == 0)
	    {
	      Space->PrintState(File, RootPositions[j]) << endl;
	      if (SaveStates == true)
		{
		  char* FileName = new char[512];
		  sprintf (FileName, "%sroot_%s%d.vec", FilePrefix, FileSuffix, NbrSavedStates); 
		  TmpVectors[j].WriteVector(FileName);
		  delete[] FileName;
		}
	      ++NbrSavedStates;
	    }
	  else
	    {
	      if (TmpRationalVectors[j].IsNullVector() == false)
		{
		  Space->PrintState(File, RootPositions[j]) << endl;
		  if (SaveStates == true)
		    {
		      char* FileName = new char[512];
		      sprintf (FileName, "%sroot_%s%d.vec", FilePrefix, FileSuffix, NbrSavedStates); 
		      TmpRationalVectors[j].WriteVector(FileName);
		      delete[] FileName;
		    }
		  ++NbrSavedStates;
		}
	    }
	}
      File << "nbr states = " << NbrSavedStates << endl;
      File << "-----------------------------------------" << endl;

      delete[] FilePrefix;
      delete[] FileSuffix;
      delete Space;
      delete[] RootPositions;
      delete[] TmpVectors;
    }
  File.close();

  return 0;
}

// resuffle an array of vectors such that the first vector is the one with the first non-zero component having the smallest index
//
// vectors = array of vectors
// nbrVectors = numerb of vectors
// componentError = error threshold to assume a component is zero
// return value = index of the first non-zero component for the first vector 

int ReshuffleVectors (RealVector* vectors, int nbrVectors, double componentError)
{
  int MinTmpPos = vectors[0].GetVectorDimension();
  int TmpVectorPos = 0;
  for (int j = 0; j < nbrVectors; ++j)      
    {
      int TmpPos = 0;
      while ((TmpPos < vectors[j].GetVectorDimension()) && (fabs(vectors[j][TmpPos]) < componentError))
	{
	  TmpPos++;
	}
      if (TmpPos < MinTmpPos)
	{
	  TmpVectorPos = j;
	  MinTmpPos = TmpPos;
	}
    }
  if (TmpVectorPos != 0)
    {
      RealVector TmpVector = vectors[TmpVectorPos];
      vectors[TmpVectorPos] = vectors[0];
       vectors[0] = TmpVector;
    }
  TmpVectorPos = 0;
  double MaxComponent = fabs(vectors[0][MinTmpPos]);
  for (int j = 1; j < nbrVectors; ++j)      
    {
      double TmpComponent = fabs(vectors[j][MinTmpPos]);
      if (TmpComponent > MaxComponent)
	{
	  TmpVectorPos = j;
	  MaxComponent = TmpComponent;
	}
    }
  if (TmpVectorPos != 0)
    {
      RealVector TmpVector = vectors[TmpVectorPos];
      vectors[TmpVectorPos] = vectors[0];
       vectors[0] = TmpVector;
    }
  return MinTmpPos;
}

// resuffle an array of vectors such that the first vector is the one with the first non-zero component having the smallest index
//
// vectors = array of vectors
// nbrVectors = numerb of vectors
// return value = index of the first non-zero component for the first vector 

int ReshuffleVectors (LongRationalVector* vectors, int nbrVectors)
{
  int MinTmpPos = vectors[0].GetVectorDimension();
  int TmpVectorPos = 0;
  for (int j = 0; j < nbrVectors; ++j)      
    {
      int TmpPos = 0;
      while ((TmpPos < vectors[j].GetVectorDimension()) && (vectors[j][TmpPos].IsZero() == true))
	{
	  TmpPos++;
	}
      if (TmpPos < MinTmpPos)
	{
	  TmpVectorPos = j;
	  MinTmpPos = TmpPos;
	}
    }
  if (TmpVectorPos != 0)
    {
      LongRationalVector TmpVector = vectors[TmpVectorPos];
      vectors[TmpVectorPos] = vectors[0];
      vectors[0] = TmpVector;
    }
  return MinTmpPos;
}

// get the quantum numbers Jz, Kz of a N-particle 4D state 
//
//totalJz = array that gives the sector Jz for all values of the linearized index
//totalKz = array that gives the sector Kz for all values of the linearized index
void GetSectorsJzKzFromLinearizedIndex(int nbrParticles, int nbrFluxQuanta, int* totalJz, int* totalKz)
  {
    int Jzmax = nbrParticles*nbrFluxQuanta;
    int index = 0;
    int ShiftedIndex = 0;
    for (int Jz = -Jzmax; Jz <= Jzmax; ++Jz)
      {
	for (int Kz = -Jzmax + abs(Jz); Kz <= Jzmax - abs(Jz) ; ++Kz)
	  {
	      if (Jz >=0)
		  index = Kz + Jzmax + 2*Jzmax*Jz - Jz*(Jz - 1);
	      else
		  index = -(Kz + Jzmax + Jz + (2*Jzmax + 1 + Jz)*(-Jz - 1)) - 1;
		ShiftedIndex = index + Jzmax*Jzmax;
	      totalJz[ShiftedIndex] = Jz;
	      totalKz[ShiftedIndex] = Kz;
	      
// 	      cout << Jz << " " << Kz << " : " << index << " " << ShiftedIndex << endl;
	  }
	}
      
  }
