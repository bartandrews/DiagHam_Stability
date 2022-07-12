#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Options/Options.h"

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::cout;
using std::endl;
using std::hex;
using std::dec;
using std::ofstream;
using std::ios;


// build the binary mask for for a Hamiltonian Z term 
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 
unsigned long BuildHamiltonianZTermMask (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ);

// build the binary mask for for a Hamiltonian Z term  (unsigned long long version)
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 
ULONGLONG BuildHamiltonianZTermMaskLongLong (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ);

// build the linearized index from the spin coordinates
//
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// spinIndex = index of the spin we want to address (either 0 or 1)
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = linearized index
int GetHaahCodeLinearizedIndex (int x, int y, int z, int spinIndex, int nbrSiteX, int nbrSiteY, int nbrSiteZ);

// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
unsigned long GetSpin0Parity(unsigned long state);

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
unsigned long GetSpin1Parity(unsigned long state);

// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
ULONGLONG GetSpin0Parity(ULONGLONG state);

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)
ULONGLONG GetSpin1Parity(ULONGLONG state);


int main(int argc, char** argv)
{
  cout.precision(14); 

  // some running options and help
  OptionManager Manager ("HaahCodeEntropy" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('z', "nbr-sitez", "number of sites along the z direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbra-sitex", "number of sites along the x direction for the region A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbra-sitey", "number of sites along the y direction for the region A", 2);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbra-sitez", "number of sites along the z direction for the region A", 2);
  (*SystemGroup) += new SingleStringOption  ('\n', "kept-sites", "column-based text file that list sites/spins that have to be kept");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "gs-parity", "select the ground state parity sector", 0);
  (*SystemGroup) += new  BooleanOption ('\n', "low-memory", "use a slower but less memory demanding algorithm", 0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*OutputGroup) += new BooleanOption  ('\n', "export-entspectrum", "export the spectrum of reduced density matrix in an ASCII file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type HaahCodeEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrSitesX = Manager.GetInteger("nbr-sitex");
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSitesZ = Manager.GetInteger("nbr-sitez");
  int NbrSitesXA = Manager.GetInteger("nbra-sitex");
  int NbrSitesYA = Manager.GetInteger("nbra-sitey");
  int NbrSitesZA = Manager.GetInteger("nbra-sitez");

  int TotalNbrSites = NbrSitesX * NbrSitesY * NbrSitesZ;
  int TotalNbrSpins = 2 * TotalNbrSites;
  double* FullReducedDensityMatrixEigenvalues;
  long FullReducedDensityMatrixNbrEigenvalues;

  int MaximumNumberZTerms = NbrSitesX * NbrSitesY * NbrSitesZ;
  ULONGLONG ProductZTerm = 0x0ul;
  for (int x = 0; x < NbrSitesX; ++x)
    {
      for (int y = 0; y < NbrSitesY; ++y)
        {
          for (int z = 0; z < NbrSitesZ; ++z)
            {
              ProductZTerm ^= BuildHamiltonianZTermMaskLongLong(x, y, z, NbrSitesX, NbrSitesY, NbrSitesZ);
            }
        }
    }
  if (ProductZTerm == ((ULONGLONG) 0x0ul))
    {
      --MaximumNumberZTerms;
      cout << "product of all the z term is equal to the identity" << endl;
    }
  long GroundStateDimension = 1l << MaximumNumberZTerms;
  timeval StartingTime;
  timeval EndingTime;  


  if (Manager.GetBoolean("low-memory") == false)
    {
      if (TotalNbrSpins <= 64)
	{
	  gettimeofday (&(StartingTime), 0);
	  unsigned long* GroundState = new unsigned long[GroundStateDimension];
	  
	  GroundState[0l] = (unsigned long) Manager.GetInteger("gs-parity");
	  GroundStateDimension = 1l;
	  
	  for (int x = 0; (x < NbrSitesX) && (MaximumNumberZTerms > 0); ++x)
	    {
	      for (int y = 0; (y < NbrSitesY) && (MaximumNumberZTerms > 0); ++y)
		{
		  for (int z = 0; (z < NbrSitesZ) && (MaximumNumberZTerms > 0); ++z)
		    {
		      unsigned long TmpMask = BuildHamiltonianZTermMask(x, y, z, NbrSitesX, NbrSitesY, NbrSitesZ);
		      for (long i = 0l; i < GroundStateDimension; ++i)
			{
			  GroundState[i + GroundStateDimension] = GroundState[i] ^ TmpMask;
			}
		      GroundStateDimension *= 2l;
		      --MaximumNumberZTerms;
		    }	  
		}      
	    }
	  gettimeofday (&(EndingTime), 0);
	  double DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			       ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << GroundStateDimension << " intermediate components generated for the groundstate (done in " << DeltaTime << "s)" << endl;
	  gettimeofday (&(StartingTime), 0);
	  SortArrayUpOrderingAndRemoveDuplicates(GroundState, GroundStateDimension);
	  gettimeofday (&(EndingTime), 0);
	  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << GroundStateDimension << " distinct components generated for the groundstate (done in " << DeltaTime << "s)" << endl;      
	  
	  double GroundStateNormalizationCoefficient = 1.0 / sqrt((double) GroundStateDimension);
	  double GroundStateSqrNormalizationCoefficient = 1.0 / ((double) GroundStateDimension);
	  
	  gettimeofday (&(StartingTime), 0);
	  int RegionANbrSites = NbrSitesXA * NbrSitesYA * NbrSitesZA;
	  unsigned long RegionAMask = 0x0ul;  
	  
	  if (Manager.GetString("kept-sites") == 0)
	    {	  
	      for (int x = 0; x < NbrSitesXA; ++x)
		{
		  for (int y = 0; y < NbrSitesYA; ++y)
		    {
		      for (int z = 0; z < NbrSitesZA; ++z)
			{
			  RegionAMask |= 0x1ul << GetHaahCodeLinearizedIndex(x, y, z, 0, NbrSitesX, NbrSitesY, NbrSitesZ);
			  RegionAMask |= 0x1ul << GetHaahCodeLinearizedIndex(x, y, z, 1, NbrSitesX, NbrSitesY, NbrSitesZ);
			}
		    }
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile KeptOrbitalFile;
	      if (KeptOrbitalFile.Parse(Manager.GetString("kept-sites")) == false)
		{
		  KeptOrbitalFile.DumpErrors(cout);
		  return -1;
		}
	      if (KeptOrbitalFile.GetNbrColumns() < 3)
		{
		  cout << "error," << Manager.GetString("kept-sites") << " should contain at least 3 columns" << endl;
		  return 0;
		}
	      RegionANbrSites = KeptOrbitalFile.GetNbrLines();
	      int* TmpKeptSpinX = KeptOrbitalFile.GetAsIntegerArray(0);
	      int* TmpKeptSpinY = KeptOrbitalFile.GetAsIntegerArray(1);
	      int* TmpKeptSpinZ = KeptOrbitalFile.GetAsIntegerArray(2);
	      int* TmpKeptSpinIndex = KeptOrbitalFile.GetAsIntegerArray(3);
	      for (int i = 0; i < RegionANbrSites; ++i)
		{
		  RegionAMask |= 0x1ul << GetHaahCodeLinearizedIndex(TmpKeptSpinX[i], TmpKeptSpinY[i], TmpKeptSpinZ[i], TmpKeptSpinIndex[i], NbrSitesX, NbrSitesY, NbrSitesZ);	      
		}
	    }
	  
	  unsigned long RegionBMask = ((0x1ul << (2 * TotalNbrSites)) - 0x1ul) & (~RegionAMask);
	  
	  unsigned long* RegionAHilbertSpace = new unsigned long[GroundStateDimension];
	  long RegionAHilbertSpaceDimension = 0l;
	  for (long i = 0l; i < GroundStateDimension; ++i)
	    {
	      RegionAHilbertSpaceDimension += SearchInSortedArrayAndInsert(GroundState[i] & RegionAMask, RegionAHilbertSpace, RegionAHilbertSpaceDimension);
	    }
	  unsigned long* TmpArray = new unsigned long[RegionAHilbertSpaceDimension]; 
	  for (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
	    {
	      TmpArray[i] = RegionAHilbertSpace[i];
	    }
	  delete[] RegionAHilbertSpace;
	  RegionAHilbertSpace = TmpArray;
	  gettimeofday (&(EndingTime), 0);
	  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << "Hilbert space dimension for the A region = " << RegionAHilbertSpaceDimension << " (done in " << DeltaTime << "s)" << endl;
	  unsigned long* RegionBConfigurations = new unsigned long[GroundStateDimension];
	  int* RegionAIndices = new int[GroundStateDimension];
	  FullReducedDensityMatrixEigenvalues = new double[RegionAHilbertSpaceDimension];
	  FullReducedDensityMatrixNbrEigenvalues = 0;
	  unsigned long* RegionAHilbertSpaceFixedParities = new unsigned long[RegionAHilbertSpaceDimension];
	  
	  for (unsigned long TmpParity = 0x0ul; TmpParity <= 0x3ul; ++TmpParity)
	    {
	      gettimeofday (&(StartingTime), 0);
	      long RegionAHilbertSpaceDimensionFixedParities = 0l;
	      for  (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
		{
		  unsigned long Tmp = RegionAHilbertSpace[i];
		  if (((GetSpin1Parity(Tmp) << 1) | GetSpin0Parity(Tmp)) == TmpParity)
		    {
		      RegionAHilbertSpaceFixedParities[RegionAHilbertSpaceDimensionFixedParities] = Tmp;
		      ++RegionAHilbertSpaceDimensionFixedParities;
		    }
		}
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "Hilbert space dimension for the region A in  the parity sector = " << ((unsigned long) TmpParity) << " (done in " << DeltaTime << "s)"  << endl;
	      
	      if (RegionAHilbertSpaceDimensionFixedParities > 0)
		{
		  cout << "Building the entanglement matrix" << endl;
		  gettimeofday (&(StartingTime), 0);
		  long RegionBHilbertSpaceDimensionFixedParities = 0l;
		  for (long i = 0l; i < GroundStateDimension; ++i)
		    {
		      int Tmp = SearchInArray(GroundState[i] & RegionAMask, RegionAHilbertSpaceFixedParities, RegionAHilbertSpaceDimensionFixedParities);
		      if (Tmp >= 0)
			{
			  RegionBConfigurations[RegionBHilbertSpaceDimensionFixedParities] = GroundState[i] & RegionBMask;
			  RegionAIndices[RegionBHilbertSpaceDimensionFixedParities] = Tmp;
			  ++RegionBHilbertSpaceDimensionFixedParities;
			}
		    }
		  SortArrayDownOrdering(RegionBConfigurations, RegionAIndices, RegionBHilbertSpaceDimensionFixedParities);
		  
		  gettimeofday (&(EndingTime), 0);	      
		  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
				((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
		  cout << "done in " << DeltaTime << "s"  << endl;
		  gettimeofday (&(StartingTime), 0);
		  cout << "Building the reduced density matrix" << endl;
		  RealSymmetricMatrix ReducedDensityMatrix (RegionAHilbertSpaceDimensionFixedParities, true);
		  long TmpIndex = 0l;
		  while (TmpIndex < RegionBHilbertSpaceDimensionFixedParities)
		    {
		      long TmpIndex2 = TmpIndex + 1l;
		      while ((TmpIndex2 < RegionBHilbertSpaceDimensionFixedParities) && (RegionBConfigurations[TmpIndex] == RegionBConfigurations[TmpIndex2]))
			{
			  ++TmpIndex2;
			}
		      for (long i = TmpIndex; i < TmpIndex2; ++i)
			{
			  for (long j= TmpIndex; j < TmpIndex2; ++j)
			    {
			      if (RegionAIndices[i] <= RegionAIndices[j])
				{
				  ReducedDensityMatrix.AddToMatrixElement(RegionAIndices[i], RegionAIndices[j], GroundStateSqrNormalizationCoefficient);
				}
			    }
			}
		      TmpIndex = TmpIndex2;
		    }
		  gettimeofday (&(EndingTime), 0);	      
		  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
				((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
		  cout << "done in " << DeltaTime << "s"  << endl;
		  gettimeofday (&(StartingTime), 0);
		  cout << "Diagonalizing the reduced density matrix (" << RegionAHilbertSpaceDimensionFixedParities << "x" << RegionAHilbertSpaceDimensionFixedParities << ")" << endl;
		  RealDiagonalMatrix ReducedDensityMatrixEigenvalues(ReducedDensityMatrix.GetNbrRow(), true);
		  ReducedDensityMatrix.LapackDiagonalize(ReducedDensityMatrixEigenvalues);
		  gettimeofday (&(EndingTime), 0);	      
		  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
				((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
		  cout << "done in " << DeltaTime << "s"  << endl;
		  for (long i = 0l; i < ReducedDensityMatrixEigenvalues.GetNbrRow(); ++i)
		    {
		      FullReducedDensityMatrixEigenvalues[FullReducedDensityMatrixNbrEigenvalues] = ReducedDensityMatrixEigenvalues[i];
		      ++FullReducedDensityMatrixNbrEigenvalues;
		    }
		}
	    }
	  delete[] RegionBConfigurations;
	  delete[] RegionAIndices;
	  delete[] GroundState;
	}      
      else
	{
	  gettimeofday (&(StartingTime), 0);
	  ULONGLONG* GroundState = new ULONGLONG[GroundStateDimension];
	  
	  GroundState[0l] = ((ULONGLONG) Manager.GetInteger("gs-parity"));
	  GroundStateDimension = 1l;
	  
	  for (int x = 0; (x < NbrSitesX) && (MaximumNumberZTerms > 0); ++x)
	    {
	      for (int y = 0; (y < NbrSitesY) && (MaximumNumberZTerms > 0); ++y)
		{
		  for (int z = 0; (z < NbrSitesZ) && (MaximumNumberZTerms > 0); ++z)
		    {
		      ULONGLONG TmpMask = BuildHamiltonianZTermMaskLongLong(x, y, z, NbrSitesX, NbrSitesY, NbrSitesZ);
		      for (long i = 0l; i < GroundStateDimension; ++i)
			{
			  GroundState[i + GroundStateDimension] = GroundState[i] ^ TmpMask;
			}
		      GroundStateDimension *= 2l;
		      --MaximumNumberZTerms;
		    }	  
		}      
	    }
	  gettimeofday (&(EndingTime), 0);
	  double DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			       ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << GroundStateDimension << " intermediate components generated for the groundstate (done in " << DeltaTime << "s)" << endl;
	  gettimeofday (&(StartingTime), 0);
	  SortArrayUpOrderingAndRemoveDuplicates(GroundState, GroundStateDimension);
	  gettimeofday (&(EndingTime), 0);
	  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << GroundStateDimension << " distinct components generated for the groundstate (done in " << DeltaTime << "s)" << endl;
	  
	  double GroundStateNormalizationCoefficient = 1.0 / sqrt((double) GroundStateDimension);
	  double GroundStateSqrNormalizationCoefficient = 1.0 / ((double) GroundStateDimension);
	  
	  gettimeofday (&(StartingTime), 0);
	  int RegionANbrSites = NbrSitesXA * NbrSitesYA * NbrSitesZA;
	  ULONGLONG RegionAMask = ((ULONGLONG) 0x0ul);  
	  if (Manager.GetString("kept-sites") == 0)
	    {	  
	      for (int x = 0; x < NbrSitesXA; ++x)
		{
		  for (int y = 0; y < NbrSitesYA; ++y)
		    {
		      for (int z = 0; z < NbrSitesZA; ++z)
			{
			  RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x, y, z, 0, NbrSitesX, NbrSitesY, NbrSitesZ);
			  RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x, y, z, 1, NbrSitesX, NbrSitesY, NbrSitesZ);
			}
		    }
		}
	    }
	  else
	    {
	      MultiColumnASCIIFile KeptOrbitalFile;
	      if (KeptOrbitalFile.Parse(Manager.GetString("kept-sites")) == false)
		{
		  KeptOrbitalFile.DumpErrors(cout);
		  return -1;
		}
	      if (KeptOrbitalFile.GetNbrColumns() < 3)
		{
		  cout << "error," << Manager.GetString("kept-sites") << " should contain at least 3 columns" << endl;
		  return 0;
		}
	      RegionANbrSites = KeptOrbitalFile.GetNbrLines();
	      int* TmpKeptSpinX = KeptOrbitalFile.GetAsIntegerArray(0);
	      int* TmpKeptSpinY = KeptOrbitalFile.GetAsIntegerArray(1);
	      int* TmpKeptSpinZ = KeptOrbitalFile.GetAsIntegerArray(2);
	      int* TmpKeptSpinIndex = KeptOrbitalFile.GetAsIntegerArray(3);
	      for (int i = 0; i < RegionANbrSites; ++i)
		{
		  RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(TmpKeptSpinX[i], TmpKeptSpinY[i], TmpKeptSpinZ[i], TmpKeptSpinIndex[i], NbrSitesX, NbrSitesY, NbrSitesZ);	      
		}
	    }
	  
	  ULONGLONG RegionBMask = ((((ULONGLONG) 0x1ul) << (2 * TotalNbrSites)) - ((ULONGLONG) 0x1ul)) & (~RegionAMask);
	  
	  ULONGLONG* RegionAHilbertSpace = new ULONGLONG[GroundStateDimension];
	  long RegionAHilbertSpaceDimension = 0l;
	  for (long i = 0l; i < GroundStateDimension; ++i)
	    {
	      RegionAHilbertSpaceDimension += SearchInSortedArrayAndInsert(GroundState[i] & RegionAMask, RegionAHilbertSpace, RegionAHilbertSpaceDimension);
	    }
	  ULONGLONG* TmpArray = new ULONGLONG[RegionAHilbertSpaceDimension]; 
	  for (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
	    {
	      TmpArray[i] = RegionAHilbertSpace[i];
	    }
	  delete[] RegionAHilbertSpace;
	  RegionAHilbertSpace = TmpArray;
	  gettimeofday (&(EndingTime), 0);
	  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << "Hilbert space dimension for the A region = " << RegionAHilbertSpaceDimension << " (done in " << DeltaTime << "s)" << endl;
	  ULONGLONG* RegionBConfigurations = new ULONGLONG[GroundStateDimension];
	  int* RegionAIndices = new int[GroundStateDimension];
	  FullReducedDensityMatrixEigenvalues = new double[RegionAHilbertSpaceDimension];
	  FullReducedDensityMatrixNbrEigenvalues = 0;
	  ULONGLONG* RegionAHilbertSpaceFixedParities = new ULONGLONG[RegionAHilbertSpaceDimension];
	  
	  for (ULONGLONG TmpParity = ((ULONGLONG) 0x0ul); TmpParity <= ((ULONGLONG) 0x3ul); ++TmpParity)
	    {
	      gettimeofday (&(StartingTime), 0);
	      long RegionAHilbertSpaceDimensionFixedParities = 0l;
	      for  (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
		{
		  ULONGLONG Tmp = RegionAHilbertSpace[i];
		  if (((GetSpin1Parity(Tmp) << 1) | GetSpin0Parity(Tmp)) == TmpParity)
		    {
		      RegionAHilbertSpaceFixedParities[RegionAHilbertSpaceDimensionFixedParities] = Tmp;
		      ++RegionAHilbertSpaceDimensionFixedParities;
		    }
		}
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "Hilbert space dimension for the region A in  the parity sector = " << ((unsigned long) TmpParity) << " (done in " << DeltaTime << "s)"  << endl;
	      
	      if (RegionAHilbertSpaceDimensionFixedParities > 0)
		{
		  cout << "Building the entanglement matrix" << endl;
		  gettimeofday (&(StartingTime), 0);
		  long RegionBHilbertSpaceDimensionFixedParities = 0l;
		  for (long i = 0l; i < GroundStateDimension; ++i)
		    {
		      int Tmp = SearchInArray(GroundState[i] & RegionAMask, RegionAHilbertSpaceFixedParities, RegionAHilbertSpaceDimensionFixedParities);
		      if (Tmp >= 0)
			{
			  RegionBConfigurations[RegionBHilbertSpaceDimensionFixedParities] = GroundState[i] & RegionBMask;
			  RegionAIndices[RegionBHilbertSpaceDimensionFixedParities] = Tmp;
			  ++RegionBHilbertSpaceDimensionFixedParities;
			}
		    }
		  SortArrayDownOrdering(RegionBConfigurations, RegionAIndices, RegionBHilbertSpaceDimensionFixedParities);
		  gettimeofday (&(EndingTime), 0);	      
		  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
				((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
		  cout << "done in " << DeltaTime << "s"  << endl;
		  gettimeofday (&(StartingTime), 0);
		  cout << "Building the reduced density matrix" << endl;
		  
		  RealSymmetricMatrix ReducedDensityMatrix (RegionAHilbertSpaceDimensionFixedParities, true);
		  long TmpIndex = 0l;
		  while (TmpIndex < RegionBHilbertSpaceDimensionFixedParities)
		    {
		      long TmpIndex2 = TmpIndex + 1l;
		      while ((TmpIndex2 < RegionBHilbertSpaceDimensionFixedParities) && (RegionBConfigurations[TmpIndex] == RegionBConfigurations[TmpIndex2]))
			{
			  ++TmpIndex2;
			}
		      for (long i = TmpIndex; i < TmpIndex2; ++i)
			{
			  for (long j= TmpIndex; j < TmpIndex2; ++j)
			    {
			      if (RegionAIndices[i] <= RegionAIndices[j])
				{
				  ReducedDensityMatrix.AddToMatrixElement(RegionAIndices[i], RegionAIndices[j], GroundStateSqrNormalizationCoefficient);
				}
			    }
			}
		      TmpIndex = TmpIndex2;
		    }
		  gettimeofday (&(EndingTime), 0);	      
		  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
				((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
		  cout << "done in " << DeltaTime << "s"  << endl;
		  gettimeofday (&(StartingTime), 0);
		  cout << "Diagonalizing the reduced density matrix (" << RegionAHilbertSpaceDimensionFixedParities << "x" << RegionAHilbertSpaceDimensionFixedParities << ")" << endl;
		  RealDiagonalMatrix ReducedDensityMatrixEigenvalues(ReducedDensityMatrix.GetNbrRow(), true);
		  ReducedDensityMatrix.LapackDiagonalize(ReducedDensityMatrixEigenvalues);
		  gettimeofday (&(EndingTime), 0);	      
		  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
				((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
		  cout << "done in " << DeltaTime << "s"  << endl;
		  for (long i = 0l; i < ReducedDensityMatrixEigenvalues.GetNbrRow(); ++i)
		    {
		      FullReducedDensityMatrixEigenvalues[FullReducedDensityMatrixNbrEigenvalues] = ReducedDensityMatrixEigenvalues[i];
		      ++FullReducedDensityMatrixNbrEigenvalues;
		    }
		}
	    }
	  delete[] RegionBConfigurations;
	  delete[] RegionAIndices;
	  delete[] GroundState;
	}
    }
  else
    {
      gettimeofday (&(StartingTime), 0);
      ULONGLONG InitialGroundState = ((ULONGLONG) Manager.GetInteger("gs-parity"));
      double GroundStateNormalizationCoefficient = 1.0 / sqrt((double) GroundStateDimension);
      double GroundStateSqrNormalizationCoefficient = 1.0 / ((double) GroundStateDimension);
      
      ULONGLONG* HamiltonianZTermMasks = new ULONGLONG[MaximumNumberZTerms];
      int TmpIndex = 0;
      for (int x = 0; (x < NbrSitesX) && (TmpIndex < MaximumNumberZTerms); ++x)
	{
	  for (int y = 0; (y < NbrSitesY) && (TmpIndex < MaximumNumberZTerms); ++y)
	    {
	      for (int z = 0; (z < NbrSitesZ) && (TmpIndex < MaximumNumberZTerms); ++z)
		{
		  HamiltonianZTermMasks[TmpIndex] = BuildHamiltonianZTermMaskLongLong(x, y, z, NbrSitesX, NbrSitesY, NbrSitesZ);
		  ++TmpIndex;
		}	  
	    }      
	}
      
      int RegionANbrSites = NbrSitesXA * NbrSitesYA * NbrSitesZA;
      ULONGLONG RegionAMask = ((ULONGLONG) 0x0ul);  
      if (Manager.GetString("kept-sites") == 0)
	{	  
	  for (int x = 0; x < NbrSitesXA; ++x)
	    {
	      for (int y = 0; y < NbrSitesYA; ++y)
		{
		  for (int z = 0; z < NbrSitesZA; ++z)
		    {
		      RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x, y, z, 0, NbrSitesX, NbrSitesY, NbrSitesZ);
		      RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x, y, z, 1, NbrSitesX, NbrSitesY, NbrSitesZ);
		    }
		}
	    }
	}
      else
	{
	  MultiColumnASCIIFile KeptOrbitalFile;
	  if (KeptOrbitalFile.Parse(Manager.GetString("kept-sites")) == false)
	    {
	      KeptOrbitalFile.DumpErrors(cout);
	      return -1;
	    }
	  if (KeptOrbitalFile.GetNbrColumns() < 3)
	    {
	      cout << "error," << Manager.GetString("kept-sites") << " should contain at least 3 columns" << endl;
	      return 0;
	    }
	  RegionANbrSites = KeptOrbitalFile.GetNbrLines();
	  int* TmpKeptSpinX = KeptOrbitalFile.GetAsIntegerArray(0);
	  int* TmpKeptSpinY = KeptOrbitalFile.GetAsIntegerArray(1);
	  int* TmpKeptSpinZ = KeptOrbitalFile.GetAsIntegerArray(2);
	  int* TmpKeptSpinIndex = KeptOrbitalFile.GetAsIntegerArray(3);
	  for (int i = 0; i < RegionANbrSites; ++i)
	    {
	      RegionAMask |= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(TmpKeptSpinX[i], TmpKeptSpinY[i], TmpKeptSpinZ[i], TmpKeptSpinIndex[i], NbrSitesX, NbrSitesY, NbrSitesZ);	      
	    }
	}
      
      ULONGLONG RegionBMask = ((((ULONGLONG) 0x1ul) << (2 * TotalNbrSites)) - ((ULONGLONG) 0x1ul)) & (~RegionAMask);

      long RegionAHilbertSpaceDimension = 1l << (2 * RegionANbrSites);
      ULONGLONG* RegionAHilbertSpace = new ULONGLONG[RegionAHilbertSpaceDimension];
      RegionAHilbertSpaceDimension = 0l;

      ULONGLONG TmpConfiguration;
      long GrayCode = 0l;
      long ChangedBit;
      int ChangedBitIndex;      
      RegionAHilbertSpaceDimension += SearchInSortedArrayAndInsert(InitialGroundState & RegionAMask, RegionAHilbertSpace, 
								   RegionAHilbertSpaceDimension);
      TmpConfiguration = InitialGroundState;
      for (long i = 1l; i < GroundStateDimension; ++i)
	{
	  ChangedBit = (i ^ (i >> 1)) ^ GrayCode;
	  GrayCode = i ^ (i >> 1);
	  ChangedBitIndex  = ((ChangedBit & 0xffffffff00000000l) != 0x0ul) << 5;
	  ChangedBitIndex |= ((ChangedBit & 0xffff0000ffff0000l) != 0x0ul) << 4;
	  ChangedBitIndex |= ((ChangedBit & 0xff00ff00ff00ff00l) != 0x0ul) << 3;
	  ChangedBitIndex |= ((ChangedBit & 0xf0f0f0f0f0f0f0f0l) != 0x0ul) << 2;
	  ChangedBitIndex |= ((ChangedBit & 0xccccccccccccccccl) != 0x0ul) << 1;
	  ChangedBitIndex |= ((ChangedBit & 0xaaaaaaaaaaaaaaaal) != 0x0ul);
	  TmpConfiguration ^= HamiltonianZTermMasks[ChangedBitIndex];
	  RegionAHilbertSpaceDimension += SearchInSortedArrayAndInsert(TmpConfiguration & RegionAMask, RegionAHilbertSpace, RegionAHilbertSpaceDimension);
	}
      ULONGLONG* TmpArray = new ULONGLONG[RegionAHilbertSpaceDimension]; 
      for (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
	{
	  TmpArray[i] = RegionAHilbertSpace[i];
	}
      delete[] RegionAHilbertSpace;
      RegionAHilbertSpace = TmpArray;
      gettimeofday (&(EndingTime), 0);
      double DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			   ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
      cout << "Hilbert space dimension for the A region = " << RegionAHilbertSpaceDimension << " (done in " << DeltaTime << "s)" << endl;
      long MaxNbrRegionBConfigurations = GroundStateDimension / 2l;
      ULONGLONG* RegionBConfigurations = new ULONGLONG[MaxNbrRegionBConfigurations];
      int* RegionAIndices = new int[MaxNbrRegionBConfigurations];
      FullReducedDensityMatrixEigenvalues = new double[RegionAHilbertSpaceDimension];
      FullReducedDensityMatrixNbrEigenvalues = 0;
      ULONGLONG* RegionAHilbertSpaceFixedParities = new ULONGLONG[RegionAHilbertSpaceDimension];
	  
      for (ULONGLONG TmpParity = ((ULONGLONG) 0x0ul); TmpParity <= ((ULONGLONG) 0x3ul); ++TmpParity)
	{
	  gettimeofday (&(StartingTime), 0);
	  long RegionAHilbertSpaceDimensionFixedParities = 0l;
	  for  (long i = 0l; i < RegionAHilbertSpaceDimension; ++i)
	    {
	      ULONGLONG Tmp = RegionAHilbertSpace[i];
	      if (((GetSpin1Parity(Tmp) << 1) | GetSpin0Parity(Tmp)) == TmpParity)
		{
		  RegionAHilbertSpaceFixedParities[RegionAHilbertSpaceDimensionFixedParities] = Tmp;
		  ++RegionAHilbertSpaceDimensionFixedParities;
		}
	    }
	  gettimeofday (&(EndingTime), 0);	      
	  DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	  cout << "Hilbert space dimension for the region A in  the parity sector = " << ((unsigned long) TmpParity) << " (done in " << DeltaTime << "s)"  << endl;
	  
	  if (RegionAHilbertSpaceDimensionFixedParities > 0)
	    {
	      cout << "Building the entanglement matrix" << endl;
	      gettimeofday (&(StartingTime), 0);
	      long RegionBHilbertSpaceDimensionFixedParities = 0l;
	      int Tmp = SearchInArray(InitialGroundState & RegionAMask, RegionAHilbertSpaceFixedParities, RegionAHilbertSpaceDimensionFixedParities);
	      if (Tmp >= 0)
		{
		  RegionBConfigurations[RegionBHilbertSpaceDimensionFixedParities] = InitialGroundState & RegionBMask;
		  RegionAIndices[RegionBHilbertSpaceDimensionFixedParities] = Tmp;
		  ++RegionBHilbertSpaceDimensionFixedParities;
		}
	      GrayCode = 0l;
	      TmpConfiguration = InitialGroundState;
	      for (long i = 1l; i < GroundStateDimension; ++i)
		{
		  ChangedBit = (i ^ (i >> 1)) ^ GrayCode;
		  GrayCode = i ^ (i >> 1);
		  ChangedBitIndex  = ((ChangedBit & 0xffffffff00000000l) != 0x0ul) << 5;
		  ChangedBitIndex |= ((ChangedBit & 0xffff0000ffff0000l) != 0x0ul) << 4;
		  ChangedBitIndex |= ((ChangedBit & 0xff00ff00ff00ff00l) != 0x0ul) << 3;
		  ChangedBitIndex |= ((ChangedBit & 0xf0f0f0f0f0f0f0f0l) != 0x0ul) << 2;
		  ChangedBitIndex |= ((ChangedBit & 0xccccccccccccccccl) != 0x0ul) << 1;
		  ChangedBitIndex |= ((ChangedBit & 0xaaaaaaaaaaaaaaaal) != 0x0ul);
		  TmpConfiguration ^= HamiltonianZTermMasks[ChangedBitIndex];
		  int Tmp = SearchInArray(TmpConfiguration & RegionAMask, RegionAHilbertSpaceFixedParities, RegionAHilbertSpaceDimensionFixedParities);
		  if (Tmp >= 0)
		    {
		      RegionBConfigurations[RegionBHilbertSpaceDimensionFixedParities] = TmpConfiguration & RegionBMask;
		      RegionAIndices[RegionBHilbertSpaceDimensionFixedParities] = Tmp;
		      ++RegionBHilbertSpaceDimensionFixedParities;
		    }
		}
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "done in " << DeltaTime << "s"  << endl;
	      gettimeofday (&(StartingTime), 0);
	      cout << "sorting " << RegionBHilbertSpaceDimensionFixedParities << " elements" << endl;
	      SortArrayDownOrdering(RegionBConfigurations, RegionAIndices, RegionBHilbertSpaceDimensionFixedParities);
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "done in " << DeltaTime << "s"  << endl;
	      gettimeofday (&(StartingTime), 0);
	      cout << "Building the reduced density matrix" << endl;
	      
	      RealSymmetricMatrix ReducedDensityMatrix (RegionAHilbertSpaceDimensionFixedParities, true);
	      long TmpIndex = 0l;
	      while (TmpIndex < RegionBHilbertSpaceDimensionFixedParities)
		{
		  long TmpIndex2 = TmpIndex + 1l;
		  while ((TmpIndex2 < RegionBHilbertSpaceDimensionFixedParities) && (RegionBConfigurations[TmpIndex] == RegionBConfigurations[TmpIndex2]))
		    {
		      ++TmpIndex2;
		    }
		  for (long i = TmpIndex; i < TmpIndex2; ++i)
		    {
		      for (long j= TmpIndex; j < TmpIndex2; ++j)
			{
			  if (RegionAIndices[i] <= RegionAIndices[j])
			    {
			      ReducedDensityMatrix.AddToMatrixElement(RegionAIndices[i], RegionAIndices[j], GroundStateSqrNormalizationCoefficient);
			    }
			}
		    }
		  TmpIndex = TmpIndex2;
		}
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "done in " << DeltaTime << "s"  << endl;
	      gettimeofday (&(StartingTime), 0);
	      cout << "Diagonalizing the reduced density matrix (" << RegionAHilbertSpaceDimensionFixedParities << "x" << RegionAHilbertSpaceDimensionFixedParities << ")" << endl;
	      RealDiagonalMatrix ReducedDensityMatrixEigenvalues(ReducedDensityMatrix.GetNbrRow(), true);
	      ReducedDensityMatrix.LapackDiagonalize(ReducedDensityMatrixEigenvalues);
	      gettimeofday (&(EndingTime), 0);	      
	      DeltaTime  = ((double) (EndingTime.tv_sec - StartingTime.tv_sec) + 
			    ((EndingTime.tv_usec - StartingTime.tv_usec) / 1000000.0));
	      cout << "done in " << DeltaTime << "s"  << endl;
	      for (long i = 0l; i < ReducedDensityMatrixEigenvalues.GetNbrRow(); ++i)
		{
		  FullReducedDensityMatrixEigenvalues[FullReducedDensityMatrixNbrEigenvalues] = ReducedDensityMatrixEigenvalues[i];
		  ++FullReducedDensityMatrixNbrEigenvalues;
		}
	    }
	}
      delete[] RegionBConfigurations;
      delete[] RegionAIndices;
    }

  double ReducedDensityMatrixTrace = 0.0;
  double ReducedDensityMatrixEntanglementEntropy = 0.0;
  long ReducedDensityMatrixNbrNonZeroEigenvalues = 0l;
  for (long i = 0l; i < FullReducedDensityMatrixNbrEigenvalues; ++i)
    {
      ReducedDensityMatrixTrace += FullReducedDensityMatrixEigenvalues[i];
    }
  for (long i = 0l; i < FullReducedDensityMatrixNbrEigenvalues; ++i)
    {
      FullReducedDensityMatrixEigenvalues[i] /= ReducedDensityMatrixTrace;
    }
  for (long i = 0l; i < FullReducedDensityMatrixNbrEigenvalues; ++i)
    {
      if (FullReducedDensityMatrixEigenvalues[i] > 0.0)
	{
	  ReducedDensityMatrixEntanglementEntropy -= FullReducedDensityMatrixEigenvalues[i] * log(FullReducedDensityMatrixEigenvalues[i]);
	  ReducedDensityMatrixNbrNonZeroEigenvalues++;
	}
    }
  if (Manager.GetBoolean("export-entspectrum") == true)
    {
      char* ReducedDensityMatrixOutputFileName = 0;
      if (Manager.GetString("kept-sites") != 0)
	{
	  ReducedDensityMatrixOutputFileName = new char[256 + strlen(Manager.GetString("kept-sites"))];
	  sprintf(ReducedDensityMatrixOutputFileName, "haahcode_entspectrum_x_%d_y_%d_z_%d_%s.dat", NbrSitesX, NbrSitesY, NbrSitesZ, Manager.GetString("kept-sites"));
	}
      else
	{
	  ReducedDensityMatrixOutputFileName = new char[512];
	  sprintf(ReducedDensityMatrixOutputFileName, "haahcode_entspectrum_x_%d_y_%d_z_%d_xa_%d_ya_%d_za_%d.dat", NbrSitesX, NbrSitesY, NbrSitesZ, NbrSitesXA, NbrSitesYA, NbrSitesZA);
	}
      ofstream File;
      File.open(ReducedDensityMatrixOutputFileName, ios::binary | ios::out);
      File.precision(14);
      File << "# Trace of the reduced density matrix before normalization = " << ReducedDensityMatrixTrace << endl;
      File << "# Number of non zero eigenvalues for the reduced density matrix = " << ReducedDensityMatrixNbrNonZeroEigenvalues << endl;
      File << "# Entangement entropy = " << (ReducedDensityMatrixEntanglementEntropy / log(2.0)) << " * log 2" << endl;
      for (long i = 0l; i < FullReducedDensityMatrixNbrEigenvalues; ++i)
	{
	  if (FullReducedDensityMatrixEigenvalues[i] > 0.0)
	    {
	      File << FullReducedDensityMatrixEigenvalues[i] << endl;
	    }
	}
      File.close();
    }

  cout << "Trace of the reduced density matrix before normalization = " << ReducedDensityMatrixTrace << endl;
  cout << "Number of non zero eigenvalues for the reduced density matrix = " << ReducedDensityMatrixNbrNonZeroEigenvalues << endl;
  cout << "Entangement entropy = " << (ReducedDensityMatrixEntanglementEntropy / log(2.0)) << " * log 2" << endl;
   

  return 0;
}


// build the linearized index from the spin coordinates
//
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// spinIndex = index of the spin we want to address (either 0 or 1)
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = linearized index

inline int GetHaahCodeLinearizedIndex (int x, int y, int z, int spinIndex, int nbrSiteX, int nbrSiteY, int nbrSiteZ)
{
  return (spinIndex + 2 * ((z % nbrSiteZ) + nbrSiteZ * ((y % nbrSiteY) + nbrSiteY * (x % nbrSiteX))));  
}

// build the binary mask for for a Hamiltonian Z term 
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 

unsigned long BuildHamiltonianZTermMask (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ)
{
  unsigned long TmpMask = 0x0ul;
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z,     0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x,     y + 1, z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= 0x1ul << GetHaahCodeLinearizedIndex(x + 1, y + 1, z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  return TmpMask;
}

// build the binary mask for for a Hamiltonian Z term  (unsigned long long version)
// 
// x = x coordinate of the cube lower front leftmost corner
// y = y coordinate of the cube lower front leftmost corner
// z = z coordinate of the cube lower front leftmost corner
// nbrSiteX = number of sites along the x direction for the full system
// nbrSiteY = number of sites along the y direction for the full system
// nbrSiteZ = number of sites along the z direction for the full system
// return value = binary mask for for a Hamiltonian Z term 

ULONGLONG BuildHamiltonianZTermMaskLongLong (int x, int y, int z, int nbrSiteX, int nbrSiteY, int nbrSiteZ)
{
  ULONGLONG TmpMask = ((ULONGLONG) 0x0ul);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x,     y,     z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x + 1, y,     z,     0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x + 1, y + 1, z,     1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x,     y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x + 1, y,     z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x,     y + 1, z + 1, 1, nbrSiteX, nbrSiteY, nbrSiteZ);
  TmpMask ^= ((ULONGLONG) 0x1ul) << GetHaahCodeLinearizedIndex(x + 1, y + 1, z + 1, 0, nbrSiteX, nbrSiteY, nbrSiteZ);
  return TmpMask;
}

// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline unsigned long GetSpin0Parity(unsigned long state)
{
  state &= 0x5555555555555555ul;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state &0x1ul);
}

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline unsigned long GetSpin1Parity(unsigned long state)
{
  state &= 0xaaaaaaaaaaaaaaaaul;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state &0x1ul);
}


// get the parity of spin 0 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline ULONGLONG GetSpin0Parity(ULONGLONG state)
{
  state &= (((ULONGLONG) 0x5555555555555555ul) << 32) | ((ULONGLONG) 0x5555555555555555ul);
  state ^= state >> 64;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state & ((ULONGLONG) 0x1ul));
}

// get the parity of spin 1 for a given configuration
//
// state = binary encoded configuration
// return value = parity (either 0 or 1)

inline ULONGLONG GetSpin1Parity(ULONGLONG state)
{
  state &= (((ULONGLONG) 0xaaaaaaaaaaaaaaaaul) << 32) | ((ULONGLONG) 0xaaaaaaaaaaaaaaaaul);
  state ^= state >> 64;
  state ^= state >> 32;
  state ^= state >> 16;
  state ^= state >> 8;
  state ^= state >> 4;
  state ^= state >> 2;
  state ^= state >> 1;
  return (state & ((ULONGLONG) 0x1ul));
}

