#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/BinomialCoefficients.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

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
  OptionManager Manager ("FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition" , "0.01");
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
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "na", "number of particles for the A part", 1);
  (*SystemGroup) += new SingleDoubleOption  ('c', "cut-off", "minus log of the cut-off to apply to the reduced density matrix eigenvalues", 32);

  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with trunc.vec extension)");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with full.parent extension) to store the reduced density matrix eigenvalues");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("ground-file") == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereFermionsTruncatedSchmidtDecompositionParticlePartition -h" << endl;
      return -1;
    }
  if ((Manager.GetString("ground-file") != 0) && 
      (IsFile(Manager.GetString("ground-file")) == false))
    {
      cout << "can't open file " << Manager.GetString("ground-file") << endl;
      return -1;
    }


  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
  unsigned long MemorySpace = Manager.GetInteger("fast-search") << 20;
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  int TotalLz = 0;
  bool Statistics = true;
  bool SVDFlag = Manager.GetBoolean("use-svd");
  char* DensityMatrixFileName = 0;
  double EigenvalueCutOff = exp(- Manager.GetDouble("cut-off"));
  double SVDEigenvalueCutOff = exp(- 0.5 * Manager.GetDouble("cut-off"));

  ParticleOnSphere* Space = 0;
  RealVector GroundState;
  char* GroundStateFile = 0;

  GroundStateFile = new char [strlen(Manager.GetString("ground-file")) + 1];
  strcpy (GroundStateFile, Manager.GetString("ground-file"));      

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFile,
						   NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << GroundStateFile << endl;
      return -1;
    }
  if (Statistics == false)
    {
      cout << GroundStateFile << " is not a fermionic state" << endl;
      return -1;
    }
  if (((NbrParticles * LzMax) & 1) != (TotalLz & 1))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << GroundStateFile << endl;
      return -1;
    }


  if (GroundState.ReadVector (GroundStateFile) == false)
    {
      cout << "can't open vector file " << GroundStateFile << endl;
      return -1;      
    }

  if (Manager.GetString("density-matrix") != 0)
    {
      DensityMatrixFileName = new char [strlen(Manager.GetString("density-matrix")) + 1];
      strcpy (DensityMatrixFileName, Manager.GetString("density-matrix"));      
    }
  else
    {
      DensityMatrixFileName = ReplaceExtensionToFileName(GroundStateFile, "vec", "full.parent");
      if (DensityMatrixFileName == 0)
	{
	  cout << "no vec extension was find in " << GroundStateFile << " file name" << endl;
	  return 0;
	}
    }


  if (Manager.GetBoolean("haldane") == false)
    {
#ifdef __64_BITS__
      if (LzMax <= 63)
#else
	if (LzMax <= 31)
#endif
	  {
	    Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax, MemorySpace);
	  }
	else
#ifdef __128_BIT_LONGLONG__
	  if (LzMax <= 126)
#else
	    if (LzMax <= 62)
#endif
	      {
		Space = new FermionOnSphereLong(NbrParticles, TotalLz, LzMax, MemorySpace);
	      }
	    else
	      Space = new FermionOnSphereUnlimited(NbrParticles, TotalLz, LzMax, MemorySpace);
	}
      else
	{
	  int* ReferenceState = 0;
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
	  if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax <= 0))
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
	  if (MaxNbrLz != (LzMax + 1))
	    {
	      cout << "wrong LzMax value in ReferenceState" << endl;
	      return -1;     
	    }
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	    if (LzMax <= 30)
#endif
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  Space = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
		else
		  Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
	      }
	    else
#ifdef __128_BIT_LONGLONG__
	      if (LzMax <= 126)
#else
		if (LzMax <= 62)
#endif
		  {
		    if (Manager.GetString("load-hilbert") != 0)
		      Space = new FermionOnSphereHaldaneBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
		    else
		      Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
		  } 
	}
  
  if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
    {
      cout << "dimension mismatch between Hilbert space and ground state" << endl;
      return 0;
    }

  RealVector SchmidtDecomposedState(GroundState.GetVectorDimension(), true);

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "#  N    Lz    lambda";
      DensityMatrixFile << endl;
      DensityMatrixFile.close();
    }

  cout.precision(14);
  
  int SubsystemNbrParticles = (NbrParticles >> 1) + (NbrParticles & 1);
  if (Manager.GetInteger("na") > 0)
    {
      SubsystemNbrParticles = Manager.GetInteger("na");
    }

  double EntanglementEntropy = 0.0;
  double DensitySum = 0.0;
  double TruncatedEntanglementEntropy = 0.0;
  double TruncatedDensitySum = 0.0;
  
  int ComplementarySubsystemNbrParticles = NbrParticles - SubsystemNbrParticles;
  int SubsystemMaxTotalLz = SubsystemNbrParticles * LzMax - (SubsystemNbrParticles * (SubsystemNbrParticles - 1));
  int ComplementaryMaxTotalLz = ComplementarySubsystemNbrParticles * LzMax - (ComplementarySubsystemNbrParticles * (ComplementarySubsystemNbrParticles - 1));
  cout << "SubsystemMaxTotalLz = " << SubsystemMaxTotalLz << "    ComplementaryMaxTotalLz = " << ComplementaryMaxTotalLz << endl;
  while (SubsystemMaxTotalLz > ComplementaryMaxTotalLz)
    SubsystemMaxTotalLz -= 2;
  int SubsystemTotalLz = -SubsystemMaxTotalLz;
  
  cout << "SubsystemMaxTotalLz = " << SubsystemMaxTotalLz << "    ComplementaryMaxTotalLz = " << ComplementaryMaxTotalLz << endl;
  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
    {
      cout << "processing subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
      RealSymmetricMatrix PartialDensityMatrix;
      RealMatrix PartialEntanglementMatrix;
      if (SVDFlag == false)
	{
	  PartialDensityMatrix = Space->EvaluatePartialDensityMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundState);
	}
      else
	{
	  PartialEntanglementMatrix = Space->EvaluatePartialEntanglementMatrixParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, GroundState,false);
	}      
      if ((PartialDensityMatrix.GetNbrRow() > 0) || (PartialEntanglementMatrix.GetNbrRow() > 0))
	{
	  RealDiagonalMatrix TmpDiag;
	  if (SVDFlag == false)
	    {
	      TmpDiag = RealDiagonalMatrix(PartialDensityMatrix.GetNbrRow());
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
	      RealMatrix AVectors;
	      RealMatrix BTVectors;	  
	      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition(AVectors, BTVectors);
	      int TmpDimension = PartialEntanglementMatrix.GetNbrColumn();
	      if (TmpDimension > PartialEntanglementMatrix.GetNbrRow())
		{
		  TmpDimension = PartialEntanglementMatrix.GetNbrRow();
		}
	      int NbrKeptEigenvalues = 0;
	      double* TmpTruncatedValues = new double[TmpDimension];
	      for (int i = 0; i < TmpDimension; ++i)
		{
		  if (TmpValues[i] > SVDEigenvalueCutOff)
		    {		      
		      ++NbrKeptEigenvalues;
		      TmpTruncatedValues[i] = TmpValues[i];
		    }
		  else
		    {
		      TmpTruncatedValues[i] = 0.0;
		    }
		}
	      if (NbrKeptEigenvalues > 0)
		{		  
		  Space->RebuildStateFromSchmidtDecompositionParticlePartition(SubsystemNbrParticles, SubsystemTotalLz, SchmidtDecomposedState, 
									       TmpDimension, TmpTruncatedValues, AVectors, BTVectors);
		}
	      cout << "number of kept eigenvalues = " << NbrKeptEigenvalues << endl;
	      for (int i = 0; i < TmpDimension; ++i)
		{
		  TmpValues[i] *= TmpValues[i];
		  if (TmpValues[i] > EigenvalueCutOff)
		    {
		      TruncatedEntanglementEntropy -= TmpValues[i] * log(TmpValues[i]);
		      TruncatedDensitySum += TmpValues[i];
		    }
		}
	      TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
	      delete[] TmpTruncatedValues;
	    }
	      
	  TmpDiag.SortMatrixDownOrder();
	  if (DensityMatrixFileName != 0)
	    {
	      ofstream DensityMatrixFile;
	      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
	      DensityMatrixFile.precision(14);
	      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
		{
		  DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
		  if (TmpDiag[i] > 0.0)
		    {
		      EntanglementEntropy -= TmpDiag[i] * log(TmpDiag[i]);
		      DensitySum += TmpDiag[i];
		    }
		}
	      DensityMatrixFile.close();
	    }
	}
//       else
// 	if (PartialDensityMatrix.GetNbrRow() == 1)
// 	  {
// 	    double TmpValue = PartialDensityMatrix(0,0);
// 	    if (DensityMatrixFileName != 0)
// 	      {
// 		ofstream DensityMatrixFile;
// 		DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
// 		DensityMatrixFile.precision(14);
// 		DensityMatrixFile << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
// 		DensityMatrixFile.close();
// 	      }		  
// 	    if (TmpValue > 0.0)
// 	      {
// 		EntanglementEntropy += TmpValue * log(TmpValue);
// 		DensitySum += TmpValue;
// 	      }
// 	    if (TmpValue > SVDEigenvalueCutOff)
// 	      {
// 		TruncatedEntanglementEntropy -= TmpValue * log(TmpValue);
// 		TruncatedDensitySum += TmpValue;
// 	      }
// 	  }
    }
  cout << "total trace = " << DensitySum << "   total entanglement entropy = " << EntanglementEntropy << endl;
  cout << "truncated trace = " << TruncatedDensitySum << "   truncated entanglement entropy = " << TruncatedEntanglementEntropy << endl;
  if (Manager.GetString("output-file") != 0)
    SchmidtDecomposedState.WriteVector(Manager.GetString("output-file"));
  else
    {
      char* TmpFileName;
      TmpFileName = ReplaceExtensionToFileName(Manager.GetString("ground-file"), "vec", "trunc.vec");
      if (TmpFileName == 0)
	{
	  cout << "no vec extension was find in " << Manager.GetString("ground-file") << " file name" << endl;
	  return 0;
	}
      SchmidtDecomposedState.WriteVector(TmpFileName);
      delete[] TmpFileName;
    }
  return 0;
}

