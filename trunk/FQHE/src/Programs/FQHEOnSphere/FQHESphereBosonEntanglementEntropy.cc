#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisLong.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"

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
  OptionManager Manager ("FQHESphereEntanglementEntropy" , "0.01");
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
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "shift-la", "index of the first orbital that is part of the subsystem whose entropy has to be evaluated (0 is the orbital at the north pole)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "stripe-subsystem", "use a stripe center around the equator as the subsystem");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*SystemGroup) += new BooleanOption  ('\n', "use-rational", "use rational numbers for the groundstates (only available in SVD mode)");
  (*ToolsGroup) += new BooleanOption  ('\n', "use-svd", "use singular value decomposition instead of diagonalization to compute the entropy");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices to a given file");
  (*OutputGroup) += new BooleanOption ('\n', "density-eigenstate", "compute the eigenstates of the reduced density matrix");
  (*OutputGroup) += new SingleIntegerOption  ('\n', "na-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed number of particles", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "lza-eigenstate", "compute the eigenstates of the reduced density matrix only for a subsystem with a fixed total Lz value", 0);
  (*OutputGroup) += new SingleIntegerOption  ('\n', "nbr-eigenstates", "number of reduced density matrix eigenstates to compute (0 if all)", 0);
  (*OutputGroup) += new SingleStringOption ('\n', "full-densitymatrix", "store full density matrices to a given file");
  (*OutputGroup) += new SingleStringOption ('\n', "rank-entanglementmatrix", "store the rank of each blook of the entanglement matrix (only available in SVD mode)");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "huge-memory", "maximum memory (in MBytes) that can allocated for precalculations when using huge mode", 100);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonEntanglementEntropy -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((Manager.GetString("ground-file") == 0) && (Manager.GetString("degenerated-groundstate") == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereEntanglementEntropy -h" << endl;
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


  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int LzMax = Manager.GetInteger("lzmax"); 
#ifdef __LAPACK__
  bool LapackFlag = Manager.GetBoolean("use-lapack");
#endif
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  char* RankEntanglementMatrixFileName = Manager.GetString("rank-entanglementmatrix");
  bool SVDFlag = Manager.GetBoolean("use-svd");
  if ((RankEntanglementMatrixFileName != 0) && (SVDFlag == false))
    {
      cout << "error, rank-entanglementmatrix is only available in SVD mode" << endl;
      return 0;
    }
  bool EigenstateFlag = Manager.GetBoolean("density-eigenstate");
  int FilterNa = Manager.GetInteger("na-eigenstate");
  int FilterLza = Manager.GetInteger("lza-eigenstate");
  int NbrEigenstates = Manager.GetInteger("nbr-eigenstates");
  int ShiftLa = Manager.GetInteger("shift-la");
  int* TotalLz = 0;
  int* NbrParticlesArray;
  int MaxNbrParticles;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphere** Spaces = 0;
  RealVector* GroundStates = 0;
  ComplexVector* ComplexGroundStates = 0;
  LongRationalVector* LongRationalGroundStates = 0;
  char** GroundStateFiles = 0;
  double* Weights =0;
  bool WeightFlag = false;
  bool ComplexFlag = false;
  bool RationalFlag = Manager.GetBoolean("use-rational");
  if ((RationalFlag == true)  && (SVDFlag == false))
    {
      cout << "error, rational mode requires SVD mode" << endl;
      return 0;
    }

  if (Manager.GetString("degenerated-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
      NbrParticlesArray = new int[1];
      Weights = new double[1];
      Weights[0] = 1.0;
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
      NbrParticlesArray = new int[NbrSpaces];
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	  strcpy (GroundStateFiles[i], DegeneratedFile(0, i));      	   
	}
      if (DegeneratedFile.GetNbrColumns() > 1)
	{
	  Weights = DegeneratedFile.GetAsDoubleArray(1);
	  WeightFlag = true;
	}
      else
	{
	  Weights = new double[NbrSpaces];
	  for (int i = 0; i < NbrSpaces; ++i)
	    Weights[i] = 1.0;
	}
    }
  
  
  for (int i = 0; i < NbrSpaces; ++i)
    {
      TotalLz[i] = 0;
      NbrParticlesArray[i] = NbrParticles;
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
						       NbrParticlesArray[i], LzMax, TotalLz[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == true)
	{
	  cout << GroundStateFiles[i] << " is not a bosonic state" << endl;
	  return -1;
	}
      if ((((NbrParticlesArray[i] * LzMax) & 1) != (TotalLz[i] & 1)) && (Manager.GetBoolean("haldane") == false))
	{
	  cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << GroundStateFiles[i] << endl;
	  return -1;
	}	
      if ( i == 0 ) 
	{
	  MaxNbrParticles = NbrParticlesArray[i];
	} 
      else if( NbrParticlesArray[i] > MaxNbrParticles )
	{
	  MaxNbrParticles = NbrParticlesArray[i];
	}
    }
  
  GroundStates = new RealVector [NbrSpaces];  
  ComplexGroundStates = new ComplexVector [NbrSpaces];  
  LongRationalGroundStates = new LongRationalVector [NbrSpaces];  
  if (RationalFlag == false)
    {
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (GroundStates[i].ReadVectorTest (GroundStateFiles[i]) == true)
	    {
	      if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
		{
		  cout << "can't open vector file " << GroundStateFiles[i] << endl;
		  return -1;      
		}
	    }
	  else
	    {
	      ComplexFlag = true;
	      if (ComplexGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
		{
		  cout << "can't open vector file " << GroundStateFiles[i] << endl;
		  return -1;
		}
	    }
	}
    }
  else
    {
      if (NbrSpaces > 1)
	{
	  cout << "error, rational mode cannot use more than one groundstate" << endl;
	  return 0;
	}
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  if (LongRationalGroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }
  
  
  Spaces = new ParticleOnSphere* [NbrSpaces];
  for (int i = 0; i < NbrSpaces; ++i)
    {
#ifdef  __64_BITS__
      if ((LzMax + NbrParticlesArray[i] - 1) < 63)
#else
	if ((LzMax + NbrParticlesArray[i] - 1) < 31)	
#endif
	  {
	    if (Manager.GetBoolean("huge-basis") == true)
	      {
		if (Manager.GetString("load-hilbert") == 0)
		  {
		    cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		    return -1;
		  }
		Spaces[i] = new  BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	      }
	    else
	      {
		if (Manager.GetBoolean("haldane") == false)
		  {
		    if ((SymmetrizedBasis == false) || (TotalLz != 0))
		      Spaces[i] = new BosonOnSphereShort (NbrParticlesArray[i], TotalLz[i], LzMax);
		    else
		      {
			Spaces[i] = new BosonOnSphereShort (NbrParticlesArray[i], TotalLz[i], LzMax);
			BosonOnSphereSymmetricBasisShort TmpSpace(NbrParticlesArray[i], LzMax);
			RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundStates[i], *((BosonOnSphereShort*) Spaces[i]));
			GroundStates[i] = OutputState;
		      }
		  }
		else
		  {
		    int* ReferenceState = 0;
		    if (Manager.GetString("reference-file") == 0)
		      {
			cout << "error, a reference file is needed" << endl;
			return 0;
		      }
		    ConfigurationParser ReferenceStateDefinition;
		    if (ReferenceStateDefinition.Parse(Manager.GetString("reference-file")) == false)
		      {
			ReferenceStateDefinition.DumpErrors(cout) << endl;
			return 0;
		      }
		    if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticlesArray[i]) == false) || (NbrParticlesArray[i] <= 0))
		      {
			cout << "NbrParticles is not defined or as a wrong value" << endl;
			return 0;
		      }
		    if ((ReferenceStateDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
		      {
			cout << "LzMax is not defined or as a wrong value" << endl;
			return 0;
		      }
		    int MaxNbrLz;
		    if (ReferenceStateDefinition.GetAsIntegerArray("ReferenceState", ' ', ReferenceState, MaxNbrLz) == false)
		      {
			cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
			return 0;     
		      }
		    if (MaxNbrLz != (LzMax + 1))
		      {
			cout << "wrong LzMax value in ReferenceState" << endl;
			return 0;     
		      }
		    if (Manager.GetString("load-hilbert") != 0)
		      Spaces[i] = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));
		    else
		      {
			Spaces[i] = new BosonOnSphereHaldaneBasisShort(NbrParticlesArray[i], TotalLz[i], LzMax, ReferenceState);	  
		      }
		  }
	      }
	  }
	else
	  {
	    if ((SymmetrizedBasis == false) || (TotalLz != 0))
	      Spaces[i] = new BosonOnSphereLong (NbrParticlesArray[i], TotalLz[i], LzMax);
	    else
	      {
		Spaces[i] = new BosonOnSphere (NbrParticlesArray[i], TotalLz[i], LzMax);
		BosonOnSphereSymmetricBasis TmpSpace(NbrParticlesArray[i], LzMax);
		RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundStates[i], *((BosonOnSphere*) Spaces[i]));
		GroundStates[i] = OutputState;
	      }
	  }
      
      if ((ComplexFlag == true) && (Spaces[i]->GetLargeHilbertSpaceDimension() != ComplexGroundStates[i].GetLargeVectorDimension()))
	{
	  cout << "dimension mismatch between Hilbert space (" << Spaces[i]->GetLargeHilbertSpaceDimension() 
	       << ") and ground state (" << ComplexGroundStates[i].GetLargeVectorDimension() << ")" << endl;
	  return 0;	  
	}
      if ((RationalFlag == true) && (Spaces[i]->GetLargeHilbertSpaceDimension() != LongRationalGroundStates[i].GetLargeVectorDimension()))
	{
	  cout << "dimension mismatch between Hilbert space (" << Spaces[i]->GetLargeHilbertSpaceDimension() 
	       << ") and ground state (" << LongRationalGroundStates[i].GetLargeVectorDimension() << ")" << endl;
	  return 0;	  
	}
      if ((RationalFlag == false) && (ComplexFlag == false) && 
	  (Spaces[i]->GetLargeHilbertSpaceDimension() != GroundStates[i].GetLargeVectorDimension()))
	{
	  cout << "dimension mismatch between Hilbert space (" << Spaces[i]->GetLargeHilbertSpaceDimension() 
	       << ") and ground state (" << GroundStates[i].GetLargeVectorDimension() << ")" << endl;
	  return 0;	  
	}
    }
  
  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# l_a    N    Lz    lambda" << endl;
      DensityMatrixFile.close();
    }
  if (RankEntanglementMatrixFileName != 0)
    {
      ofstream RankEntanglementMatrixFile;
      RankEntanglementMatrixFile.open(RankEntanglementMatrixFileName, ios::binary | ios::out); 
      RankEntanglementMatrixFile << "# l_a    N    Lz    rank    dimA    dimB" << endl;
      RankEntanglementMatrixFile.close();
    }
  ofstream File;
  if (Manager.GetString("output-file") != 0)
    File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
    {
      char* TmpFileName = ReplaceExtensionToFileName(GroundStateFiles[0], "vec", "ent");
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
  int MeanSubsystemSize = LzMax >> 1;
  if ((LzMax & 1) != 0)
    ++MeanSubsystemSize;
  if (Manager.GetInteger("max-la") > 0)
    {
      MeanSubsystemSize = Manager.GetInteger("max-la");
      if (MeanSubsystemSize > LzMax)
	MeanSubsystemSize = LzMax;
    }
  int SubsystemSize = Manager.GetInteger("min-la");
  if (SubsystemSize < 1)
    SubsystemSize = 1;
  BinomialCoefficients Coefs(MeanSubsystemSize + MaxNbrParticles - 1);
  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;      
      long MaximumSize = 0;
      for (int i = 0; i <= MaxNbrParticles; ++i)
	MaximumSize += Coefs(SubsystemSize + i - 1, i);
      double* TmpDensityMatrixEigenvalues = new double [MaximumSize];
      long TmpDensityMatrixEigenvaluePosition = 0;
      for (int SubsystemNbrParticles = 0; SubsystemNbrParticles <= MaxNbrParticles; ++SubsystemNbrParticles)
	{
	  int SubsystemTotalLz = 0;
	  int SubsystemLzMax = SubsystemSize - 1;
	  int SubsystemMaxTotalLz = SubsystemNbrParticles * SubsystemLzMax;
	  SubsystemTotalLz = -SubsystemMaxTotalLz;
	  RealSymmetricMatrix PartialDensityMatrix;
	  RealMatrix PartialEntanglementMatrix;
	  LongRationalMatrix PartialRationalEntanglementMatrix;
	  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
	    {
	      if (Manager.GetBoolean("stripe-subsystem") == true)
		{
		  ShiftLa = (LzMax - SubsystemSize + 1) >> 1;
		}
	      cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
	      
	      if ((ComplexFlag == false) && (RationalFlag == false))
		{
		  if (SVDFlag == false)
		    {
		      PartialDensityMatrix = Spaces[0]->EvaluateShiftedPartialDensityMatrix(SubsystemSize, ShiftLa, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0]);
		      if (WeightFlag == true)
			PartialDensityMatrix *= Weights[0];
		    }
		  else
		    {
		      PartialEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0]);
		      if (WeightFlag == true)
			PartialEntanglementMatrix *= sqrt(Weights[0]);
		    }
		  
		  for (int i = 1; i < NbrSpaces; ++i)
		    {
		      RealSymmetricMatrix TmpMatrix;	
		      RealMatrix TmpEntanglementMatrix;
		      
		      if (SVDFlag == false)
			{
			  TmpMatrix = Spaces[i]->EvaluateShiftedPartialDensityMatrix(SubsystemSize, ShiftLa, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i]);
			  if (WeightFlag == true)
			    TmpMatrix *= Weights[i];
			  if (PartialDensityMatrix.GetNbrRow() == 0)
			    PartialDensityMatrix = TmpMatrix;
			  else
			    PartialDensityMatrix += TmpMatrix;
			}
		      else
			{
			  TmpEntanglementMatrix = Spaces[i]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i]);
			  if (WeightFlag == true)
			    PartialEntanglementMatrix *= sqrt(Weights[i]);
			  PartialEntanglementMatrix += TmpEntanglementMatrix;
			}
		    }
		  
		  if (SVDFlag == false)
		    {
		      if ((NbrSpaces > 1) && (WeightFlag == false))
			PartialDensityMatrix /= ((double) NbrSpaces);
		    }
		  else
		    {
		      if ((NbrSpaces > 1) && (WeightFlag == false))
			PartialEntanglementMatrix /= sqrt((double) NbrSpaces);
		    }
		  
		  if ((Manager.GetString("full-densitymatrix") != 0) && (FilterNa == SubsystemNbrParticles) && (FilterLza = SubsystemTotalLz))
		    {
		      ofstream FullDensityMatrixFile;
		      FullDensityMatrixFile.open(Manager.GetString("full-densitymatrix"), ios::binary | ios::out | ios::app); 
		      FullDensityMatrixFile << PartialDensityMatrix;
		      FullDensityMatrixFile.close();
		    }
		  if ((PartialDensityMatrix.GetNbrRow() > 1) || ((PartialEntanglementMatrix.GetNbrRow() >= 1) && (PartialEntanglementMatrix.GetNbrColumn() >= 1)))
		    {
		      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
		      if (SVDFlag == false)
			{
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    {
			      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
				  && (FilterLza == SubsystemTotalLz ))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(), PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
						   MaxNbrParticles, LzMax, TotalLz[0], SubsystemSize,
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
				}
			    }
			  else
			    {
			      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
				  && (FilterLza == SubsystemTotalLz ))
				{
				  RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							    PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
						   MaxNbrParticles, LzMax, TotalLz[0], SubsystemSize,
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.Diagonalize(TmpDiag);
				}
			    }
#else
			  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
			}
		      else
			{  
			  if (RankEntanglementMatrixFileName != 0)
			    {
			      ofstream RankEntanglementMatrixFile;
			      RankEntanglementMatrixFile.open(RankEntanglementMatrixFileName, ios::binary | ios::out | ios::app); 
			      RankEntanglementMatrixFile.precision(14);
			      RankEntanglementMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << PartialEntanglementMatrix.Rank() << " "
							 <<  PartialEntanglementMatrix.GetNbrRow()<< " " << PartialEntanglementMatrix.GetNbrColumn() << endl;
			      RankEntanglementMatrixFile.close();
			      
			    }
			  if ((PartialEntanglementMatrix.GetNbrRow() > 1) && (PartialEntanglementMatrix.GetNbrColumn() > 1))
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
			  else
			    {
			      double TmpValue = 0.0;
			      if (PartialEntanglementMatrix.GetNbrRow() == 1)
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
				    TmpValue += PartialEntanglementMatrix[i][0] * PartialEntanglementMatrix[i][0];
				}
			      else
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
				    TmpValue += PartialEntanglementMatrix[0][i] * PartialEntanglementMatrix[0][i];				  
				}
			      TmpDiag = RealDiagonalMatrix(1, 1);
			      TmpDiag[0] = TmpValue;
			    }
			}
		      
		      TmpDiag.SortMatrixDownOrder();
		      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpDiag[i];
		      if (DensityMatrixFileName != 0)
			{
			  ofstream DensityMatrixFile;
			  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			  DensityMatrixFile.precision(14);
			  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			    DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			  DensityMatrixFile.close();
			}
		    }
		  else
		    {
		      if (PartialDensityMatrix.GetNbrRow() == 1)
			{
			  double TmpValue = PartialDensityMatrix(0,0);
			  TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpValue;
			  if (DensityMatrixFileName != 0)
			    {
			      ofstream DensityMatrixFile;
			      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			      DensityMatrixFile.precision(14);
			      DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
			      DensityMatrixFile.close();
			    }		  
			}
		    }
		}
	      else
		{
		  if (ComplexFlag == true)
		    {
		      HermitianMatrix PartialDensityMatrix = Spaces[0]->EvaluateShiftedPartialDensityMatrix(SubsystemSize, ShiftLa, SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[0]);
		      if (WeightFlag == true)
			PartialDensityMatrix *= Weights[0];
		      for (int i = 1; i < NbrSpaces; ++i)
			{
			  HermitianMatrix TmpMatrix = Spaces[i]->EvaluateShiftedPartialDensityMatrix(SubsystemSize, ShiftLa, SubsystemNbrParticles, SubsystemTotalLz, ComplexGroundStates[i]);
			  if ( TmpMatrix .GetNbrRow() > 0 )
			    {
			      if (WeightFlag == true)
				TmpMatrix *= Weights[i];		      
			      if (PartialDensityMatrix.GetNbrRow() == 0)
				PartialDensityMatrix = TmpMatrix;
			      else
				PartialDensityMatrix += TmpMatrix;
			    }
			}
		      if ((NbrSpaces > 1) && (WeightFlag == false))
			PartialDensityMatrix /= ((double) NbrSpaces);
		      if ((Manager.GetString("full-densitymatrix") != 0) && (FilterNa == SubsystemNbrParticles) && (FilterLza = SubsystemTotalLz))
			{
			  ofstream FullDensityMatrixFile;
			  FullDensityMatrixFile.open(Manager.GetString("full-densitymatrix"), ios::binary | ios::out | ios::app); 
			  FullDensityMatrixFile << PartialDensityMatrix;
			  FullDensityMatrixFile.close();
			}
		      if (PartialDensityMatrix.GetNbrRow() > 1)
			{
			  RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
			  if (LapackFlag == true)
			    {
			      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
				  && (FilterLza == SubsystemTotalLz ))
				{
				  ComplexMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							       PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
						   MaxNbrParticles, LzMax, TotalLz[0], SubsystemSize,
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.LapackDiagonalize(TmpDiag);
				}
			    }
			  else
			    {
			      if ((EigenstateFlag == true) && (FilterNa == SubsystemNbrParticles)
				  && (FilterLza == SubsystemTotalLz ))
				{
				  ComplexMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
							       PartialDensityMatrix.GetNbrRow(), true);
				  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				    TmpEigenstates[i][i] = 1.0;
				  PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates, Manager.GetDouble("diag-precision"));
				  TmpDiag.SortMatrixDownOrder(TmpEigenstates);
				  char* TmpEigenstateName = new char[512];
				  int MaxNbrEigenstates = NbrEigenstates;
				  if (NbrEigenstates == 0)
				    MaxNbrEigenstates = PartialDensityMatrix.GetNbrRow();
				  for (int i = 0; i < MaxNbrEigenstates; ++i)
				    {
				      if (TmpDiag[i] > 1e-14)
					{
					  sprintf (TmpEigenstateName,
						   "bosons_sphere_density_n_%d_2s_%d_lz_%d_la_%d_na_%d_lza_%d.%d.vec",
						   MaxNbrParticles, LzMax, TotalLz[0], SubsystemSize,
						   SubsystemNbrParticles, SubsystemTotalLz, i);
					  TmpEigenstates[i].WriteVector(TmpEigenstateName);
					}
				    }
				  delete[] TmpEigenstateName;
				}
			      else
				{
				  PartialDensityMatrix.Diagonalize(TmpDiag);
				}
			    }
#else
			  PartialDensityMatrix.Diagonalize(TmpDiag);
#endif		  
			  TmpDiag.SortMatrixDownOrder();
			  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
			    TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpDiag[i];
			  if (DensityMatrixFileName != 0)
			    {
			      ofstream DensityMatrixFile;
			      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			      DensityMatrixFile.precision(14);
			      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			      DensityMatrixFile.close();
			    }
			}
		      else
			if (PartialDensityMatrix.GetNbrRow() == 1)
			  {
			    double TmpValue = PartialDensityMatrix(0,0);
			    TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpValue;
			    if (DensityMatrixFileName != 0)
			      {
				ofstream DensityMatrixFile;
				DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
				DensityMatrixFile.precision(14);
				DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpValue << endl;
				DensityMatrixFile.close();
			      }		  
			  }
		    }
		  else
		    {
		      // rational case
		      PartialRationalEntanglementMatrix = Spaces[0]->EvaluatePartialEntanglementMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, LongRationalGroundStates[0]);
		      if ((PartialRationalEntanglementMatrix.GetNbrRow() >= 1) && (PartialRationalEntanglementMatrix.GetNbrColumn() >= 1))
			{
			  PartialEntanglementMatrix = RealMatrix(PartialRationalEntanglementMatrix);
			  if (RankEntanglementMatrixFileName != 0)
			    {
			      ofstream RankEntanglementMatrixFile;
			      RankEntanglementMatrixFile.open(RankEntanglementMatrixFileName, ios::binary | ios::out | ios::app); 
			      RankEntanglementMatrixFile.precision(14);
			      RankEntanglementMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << PartialRationalEntanglementMatrix.Rank() << " "
							 <<  PartialRationalEntanglementMatrix.GetNbrRow()<< " " << PartialRationalEntanglementMatrix.GetNbrColumn() << endl;
			      RankEntanglementMatrixFile.close();
			      
			    }
			  RealDiagonalMatrix TmpDiag;
			  if ((PartialRationalEntanglementMatrix.GetNbrRow() > 1) && (PartialRationalEntanglementMatrix.GetNbrColumn() > 1))
			    {
			      double* TmpValues = PartialEntanglementMatrix.SingularValueDecomposition();
			      int TmpDimension = PartialRationalEntanglementMatrix.GetNbrColumn();
			      if (TmpDimension > PartialRationalEntanglementMatrix.GetNbrRow())
				{
				  TmpDimension = PartialRationalEntanglementMatrix.GetNbrRow();
				}
			      for (int i = 0; i < TmpDimension; ++i)
				TmpValues[i] *= TmpValues[i];
			      TmpDiag = RealDiagonalMatrix(TmpValues, TmpDimension);
			    }
			  else
			    {
			      double TmpValue = 0.0;
			      if (PartialEntanglementMatrix.GetNbrRow() == 1)
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrColumn(); ++i)
				    TmpValue += PartialEntanglementMatrix[i][0] * PartialEntanglementMatrix[i][0];
				}
			      else
				{
				  for (int i = 0; i < PartialEntanglementMatrix.GetNbrRow(); ++i)
				    TmpValue += PartialEntanglementMatrix[0][i] * PartialEntanglementMatrix[0][i];				  
				}
			      TmpDiag = RealDiagonalMatrix(1, 1);
			      TmpDiag[0] = TmpValue;
			    }
			  TmpDiag.SortMatrixDownOrder();
			  for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
			    TmpDensityMatrixEigenvalues[TmpDensityMatrixEigenvaluePosition++] = TmpDiag[i];
			  if (DensityMatrixFileName != 0)
			    {
			      ofstream DensityMatrixFile;
			      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
			      DensityMatrixFile.precision(14);
			      for (int i = 0; i < TmpDiag.GetNbrRow(); ++i)
				DensityMatrixFile << SubsystemSize << " " << SubsystemNbrParticles << " " << SubsystemTotalLz << " " << TmpDiag[i] << endl;
			      DensityMatrixFile.close();
			    }
			}
		    }
		}
	    }
	}
      EntanglementEntropy = 0.0;
      DensitySum = 0.0;
      cout << "sorting density matrix eigenvalues and computing entanglement entropy" << endl;
      SortArrayDownOrdering(TmpDensityMatrixEigenvalues, TmpDensityMatrixEigenvaluePosition);
      unsigned TmpPos = 0;
      for (; (TmpPos < TmpDensityMatrixEigenvaluePosition) && (DensitySum < 1.0); ++TmpPos)
	{
	  if (TmpDensityMatrixEigenvalues[TmpPos] > 1e-14)
	    {
	      EntanglementEntropy += TmpDensityMatrixEigenvalues[TmpPos] * log(TmpDensityMatrixEigenvalues[TmpPos]);
	      DensitySum += TmpDensityMatrixEigenvalues[TmpPos];
	    }
	}
      double DensitySumError = 0.0;
      for (; TmpPos < TmpDensityMatrixEigenvaluePosition; ++TmpPos)
	DensitySumError += TmpDensityMatrixEigenvalues[TmpPos];
      delete[] TmpDensityMatrixEigenvalues;
      File << SubsystemSize << " " << (-EntanglementEntropy) << " " << DensitySum << " " << DensitySumError << endl;
    }
  File.close();
}

