#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/FermionOnDisk.h"
#include "HilbertSpace/FermionOnDiskUnlimited.h"
#include "HilbertSpace/FermionOnDiskHaldaneBasis.h"
#include "HilbertSpace/FermionOnDiskLong.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/FQHEOnDiskFileTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHEDiskFermionsRealSpaceEntanglementEntropy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  
  Manager += SystemGroup;
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\0', "ground-file", "name of the file corresponding to the ground state of the whole system");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "force-maxmomentum", "force the maximum single particle momentum to a particular value (negative from the number of particles and the state total angular momentum)", -1);
  (*SystemGroup) += new SingleIntegerOption ('\n',"min-na","minimum number of particlesin the region A whose entropy has to be evaluated", 0);
  (*SystemGroup) += new SingleIntegerOption ('\n',"max-na","maximum number of particlesin the region A whose entropy has to be evaluated", -1);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new SingleDoubleOption ('\n', "radius", "radius that defines the size of the real space parition", 1.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "shift-orbitals", "shift the angular momentum of all orbitals", 0);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with partent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");	
  
  if(Manager.ProceedOptions(argv, argc, cout) == false)
    {	
      cout << "see man page for option syntax or type FQHEDiskFermionsRealSpaceEntanglementEntropy -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereFermionEntanglementEntropyParticlePartition -h" << endl;
      return -1;
    }
  
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if ((Manager.GetString("degenerated-groundstate") != 0) && 
      (IsFile(Manager.GetString("degenerated-groundstate")) == false))
    {
      cout << "can't open file " << Manager.GetString("degenerated-groundstate") << endl;
      return -1;
    }
  
  
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int ForceMaxMomentum = Manager.GetInteger("force-maxmomentum");
  
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  unsigned long MemorySpace = Manager.GetInteger("fast-search") << 20;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  
  int* TotalLz = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphere** Spaces = 0;
  RealVector* GroundStates = 0;
  char** GroundStateFiles = 0;
  int* ReferenceState = 0;
  bool SVDFlag = Manager.GetBoolean("use-svd");
  double Radius = Manager.GetDouble ("radius");
  
  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("ground-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("ground-file"));      
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerated-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
      NbrSpaces = DegeneratedFile.GetNbrLines();
      GroundStateFiles = new char* [NbrSpaces];
      TotalLz = new int[NbrSpaces];
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	  strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	}
    }
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalLz[i] = 0;
      if (FQHEOnDiskFindSystemInfoFromFileName(GroundStateFiles[i], NbrParticles, ForceMaxMomentum, TotalLz[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == false)
	{
	  cout << GroundStateFiles[i] << " is not a fermionic state" << endl;
	  return -1;
	}
    }

  
  GroundStates = new RealVector [NbrSpaces];  
  for (int i = 0; i < NbrSpaces; ++i)
    if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[i] << endl;
	return -1;      
      }
  
  
  Spaces = new ParticleOnSphere* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
      int TmpMaxMomentum = (TotalLz[i] - (((NbrParticles - 1) * (NbrParticles - 2)) / 2));
      if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < TmpMaxMomentum))
	TmpMaxMomentum = ForceMaxMomentum;
      
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (TmpMaxMomentum <= 62)
#else
	    if (TmpMaxMomentum <= 30)
#endif
	      Spaces[i] = new FermionOnDisk(NbrParticles, TotalLz[i], TmpMaxMomentum, MemorySpace);
	    else
#ifdef __128_BIT_LONGLONG__
	      if (TmpMaxMomentum <= 126)
#else
		if (TmpMaxMomentum <= 62)
#endif
		  Spaces[i] = new FermionOnDiskLong(NbrParticles, TotalLz[i], TmpMaxMomentum, MemorySpace);
		else
		  Spaces[i] = new FermionOnDiskUnlimited(NbrParticles, TotalLz[i], TmpMaxMomentum, MemorySpace);
	}
      else
	{
#ifdef __64_BITS__
	  if (TmpMaxMomentum <= 62)
#else
	    if (TmpMaxMomentum <= 30)
#endif
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  Spaces[i] = new FermionOnDiskHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
		else
		  {
		    ConfigurationParser ReferenceStateDefinition;
		    if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
		      {	
			ReferenceStateDefinition.DumpErrors(cout) << endl;
			return -1;
		      }
		    if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
		      {
			cout << "NbrParticles is not defined or as a wrong value" << endl;
			return -1;
		      }
		    if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", ForceMaxMomentum) == false) || (ForceMaxMomentum <= 0))
		      {
			cout << "LzMax is not defined or as a wrong value" << endl;
			return -1;
		      }
		    int MaxNbrLz;
		    if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		      {
			cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
			return -1;     
		      }
		    if (MaxNbrLz != (ForceMaxMomentum + 1))
		      {
			cout << "wrong LzMax value in ReferenceState" << endl;
			return -1;     
		      }
		    TotalLz[i] = 0;
		    for (int k = 1; k <= ForceMaxMomentum; ++k)
		      TotalLz[i] += k * ReferenceState[k];
		    
		    Spaces[i] = new FermionOnDiskHaldaneBasis(NbrParticles, TotalLz[i],ForceMaxMomentum, ReferenceState, MemorySpace);
		  }
	      }
	}
      
      if (Spaces[i]->GetHilbertSpaceDimension() != GroundStates[i].GetVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    }
  
  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Lz    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }
  
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "partent");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFiles[0] << " file name" << endl;
	  return 0;
	}
      File.open(TmpFileName, ios::binary | ios::out);
      delete[] TmpFileName;
    }
  File.precision(14);
  cout.precision(14);
  
  
  int  MaxSubsystemNbrParticles = NbrParticles;
  if( Manager.GetInteger("max-na") != -1)
    {
      MaxSubsystemNbrParticles = Manager.GetInteger("max-na");
    }
  int SubsystemNbrParticles = Manager.GetInteger("min-na");
  
  int MaxMomentum = (TotalLz[0] - (((NbrParticles - 1) * (NbrParticles - 2)) / 2));
  
  if ((ForceMaxMomentum >= 0) && (ForceMaxMomentum < MaxMomentum))
    MaxMomentum = ForceMaxMomentum;
  
  double TotalTrace = 0.0;
  cout <<"MaxSubsystemNbrParticles = "<<MaxSubsystemNbrParticles<<endl;
  for (; SubsystemNbrParticles <= MaxSubsystemNbrParticles; ++SubsystemNbrParticles)
    {
      cout <<"SubsystemNbrParticles = " <<SubsystemNbrParticles<<endl;
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      
      int ComplementarySubsystemNbrParticles = NbrParticles - SubsystemNbrParticles;
      int SubsystemMaxTotalLz = SubsystemNbrParticles * MaxMomentum - ((SubsystemNbrParticles * (SubsystemNbrParticles - 1))>>1);
      int ComplementaryMaxTotalLz = ComplementarySubsystemNbrParticles * MaxMomentum - ((ComplementarySubsystemNbrParticles * (ComplementarySubsystemNbrParticles - 1))>>1);
      int ComplementaryMinTotalLz = ComplementarySubsystemNbrParticles * (ComplementarySubsystemNbrParticles - 1) >> 1;
      
      int SubsystemTotalLz = SubsystemNbrParticles * (SubsystemNbrParticles - 1)>>1; 
      
      while ( (TotalLz[0] - SubsystemMaxTotalLz ) < ComplementaryMinTotalLz )
	SubsystemMaxTotalLz--;
      
      while ((TotalLz[0]-SubsystemTotalLz)>ComplementaryMaxTotalLz)
	SubsystemTotalLz++;
      
      cout << "SubsystemMaxTotalLz = " << SubsystemMaxTotalLz << "    SubsystemTotalLz = " << SubsystemTotalLz << endl;
      for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz++)
	{
	  int ShiftedTotaLz = 2*SubsystemTotalLz - SubsystemNbrParticles * MaxMomentum;
	  cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
	  RealSymmetricMatrix PartialDensityMatrix;
	  RealMatrix PartialEntanglementMatrix;
	  
	  if (SVDFlag == false)
	    {
	      if(HaldaneBasisFlag == false )
		{	
		  PartialDensityMatrix = ((FermionOnDisk*) Spaces[0])->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, Radius, GroundStates[0],Manager.GetInteger("shift-orbitals"));
		}
	      else
		{
		  PartialDensityMatrix = ((FermionOnDiskHaldaneBasis*) Spaces[0])->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, Radius, GroundStates[0],Manager.GetInteger("shift-orbitals"));
		}
	    }
	  else
	    {
	      if(HaldaneBasisFlag == false)
		{
		  PartialEntanglementMatrix = ((FermionOnDisk*) Spaces[0])->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, ShiftedTotaLz, GroundStates[0] ,true);
		  if(PartialEntanglementMatrix.GetNbrRow() != 0)
		    {
		      ((FermionOnDisk*) Spaces[0])->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, Radius, Manager.GetInteger("shift-orbitals"), PartialEntanglementMatrix);
		    }
		}
	      else
		{
		  PartialEntanglementMatrix = ((FermionOnDiskHaldaneBasis*) Spaces[0])->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, ShiftedTotaLz, GroundStates[0] ,true);
		  
		  if(PartialEntanglementMatrix.GetNbrRow() != 0)
		    {
		      ((FermionOnDiskHaldaneBasis*) Spaces[0])->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, Radius, Manager.GetInteger("shift-orbitals"), PartialEntanglementMatrix);
		    }			
		}
	      
	    }
	  
	  for (int i = 1; i < NbrSpaces; ++i)
	    {
	      RealSymmetricMatrix TmpMatrix;
	      RealMatrix TmpEntanglementMatrix;
	      
	      if (SVDFlag == false)
		{	
		  if( HaldaneBasisFlag == false )
		    {	
		      TmpMatrix = ((FermionOnDisk*) Spaces[i])->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, Radius, GroundStates[i],Manager.GetInteger("shift-orbitals"));
		    }
		  else
		    {
		      TmpMatrix = ((FermionOnDiskHaldaneBasis*) Spaces[i])->EvaluatePartialDensityMatrixRealSpacePartition(SubsystemNbrParticles, SubsystemTotalLz, Radius, GroundStates[i],Manager.GetInteger("shift-orbitals"));
		    }
		  PartialDensityMatrix += TmpMatrix;
		}
	      else
		{
		  if(HaldaneBasisFlag == false)
		    {
		      TmpEntanglementMatrix = ((FermionOnDisk*) Spaces[0])->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, ShiftedTotaLz, GroundStates[i] ,true);
		      if(PartialEntanglementMatrix.GetNbrRow() != 0)
			{
			  ((FermionOnDisk*) Spaces[0])->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, Radius, Manager.GetInteger("shift-orbitals"), TmpEntanglementMatrix);
			}
		    }
		  else
		    {
		      TmpEntanglementMatrix = ((FermionOnDiskHaldaneBasis*) Spaces[0])->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, ShiftedTotaLz, GroundStates[i] ,true);
		      if(PartialEntanglementMatrix.GetNbrRow() != 0)
			{
			  ((FermionOnDiskHaldaneBasis*) Spaces[0])->EvaluateEntanglementMatrixRealSpacePartitionFromParticleEntanglementMatrix(SubsystemNbrParticles, SubsystemTotalLz, Radius, Manager.GetInteger("shift-orbitals"), TmpEntanglementMatrix);
			}			
		    }
		  PartialEntanglementMatrix += TmpEntanglementMatrix;
		}
	    }
	  
	  if (NbrSpaces > 1)
	    {
	      if (SVDFlag == false)
		PartialDensityMatrix /= ((double) NbrSpaces);
	      else
		PartialEntanglementMatrix /= sqrt((double) NbrSpaces);
	    }
	  if ((PartialDensityMatrix.GetNbrRow() > 1) || (PartialEntanglementMatrix.GetNbrRow() >= 1))
	    {
	      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
	      
	      if (SVDFlag == false)
		{
#ifdef __LAPACK__
		  if (LapackFlag == true)
		    PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		  else
		    PartialDensityMatrix.Diagonalize(TmpDiag);
#else
		  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
		}
	      else
		{
		  double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
		  int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
		  if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		    {
		      TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		    }
		  for (int i = 0; i < TmpDimension; ++i)
		    TmpValues[i] *= TmpValues[i];
		  TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
		}
	      TmpDiag.SortMatrixDownOrder();
	      if (DensityMatrixFileName != 0)
		{
		  ofstream DensityMatrixFile;
		  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		  DensityMatrixFile.precision(14);
		  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
		  DensityMatrixFile.close();
		}
	      
	      
	      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		{
		  if (TmpDiag[i] > 1e-14)
		    {
		      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		      DensitySum += TmpDiag[i];
		      TotalTrace += TmpDiag[i];
		    }
		}
	    }
	  else
	    if (PartialDensityMatrix.GetNbrRow() == 1)
	      {
		double TmpValue = PartialDensityMatrix(0,0);
		if (DensityMatrixFileName != 0)
		  {
		    ofstream DensityMatrixFile;
		    DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		    DensityMatrixFile.precision(14);
		    DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
		    DensityMatrixFile.close();
		  }		  
		if (TmpValue > 1e-14)
		  {
		    EntanglementEntropy += TmpValue * log(TmpValue);
		    DensitySum += TmpValue;
		    TotalTrace += TmpValue;
		  }
	      }
	}
      File << SubsystemNbrParticles << " " << (-EntanglementEntropy) << " " << DensitySum << " " << (1.0 - DensitySum) << " "<<TotalTrace<<endl;
      cout << "trace = " << DensitySum << endl;
    }
  File.close();
  cout <<"Total trace = " << TotalTrace<<endl;
  return 0;
}
