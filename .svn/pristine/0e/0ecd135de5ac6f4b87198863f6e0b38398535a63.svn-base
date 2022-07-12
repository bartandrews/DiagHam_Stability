#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"

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
  (*SystemGroup) += new SingleIntegerOption  ('\n', "min-la", "minimum size of the subsystem whose entropy has to be evaluated", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "max-la", "maximum size of the subsystem whose entropy has to be evaluated (0 if equal to half the total system size)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerated-groundstate", "single column file describing a degenerated ground state");
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with ent extension");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the partial density matrices in the a given file");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonEntanglementEntropy -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if ((((SingleStringOption*) Manager["ground-file"])->GetString() == 0) && (((SingleStringOption*) Manager["degenerated-groundstate"])->GetString() == 0))
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereEntanglementEntropy -h" << endl;
      return -1;
    }
  if ((((SingleStringOption*) Manager["ground-file"])->GetString() != 0) && 
      (IsFile(((SingleStringOption*) Manager["ground-file"])->GetString()) == false))
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["ground-file"])->GetString() << endl;
      return -1;
    }
  if ((((SingleStringOption*) Manager["degenerated-groundstate"])->GetString() != 0) && 
      (IsFile(((SingleStringOption*) Manager["degenerated-groundstate"])->GetString()) == false))
    {
      cout << "can't open file " << ((SingleStringOption*) Manager["degenerated-groundstate"])->GetString() << endl;
      return -1;
    }


  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger(); 
#ifdef __LAPACK__
  bool LapackFlag = ((BooleanOption*) Manager["use-lapack"])->GetBoolean();
#endif
  char* DensityMatrixFileName = ((SingleStringOption*) Manager["density-matrix"])->GetString();
  int* TotalLz = 0;
  bool Statistics = true;
  int NbrSpaces = 1;
  ParticleOnSphere** Spaces = 0;
  RealVector* GroundStates = 0;
  char** GroundStateFiles = 0;

  if (((SingleStringOption*) Manager["degenerated-groundstate"])->GetString() == 0)
    {
      GroundStateFiles = new char* [1];
      TotalLz = new int[1];
      GroundStateFiles[0] = new char [strlen(((SingleStringOption*) Manager["ground-file"])->GetString()) + 1];
      strcpy (GroundStateFiles[0], ((SingleStringOption*) Manager["ground-file"])->GetString());      
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(((SingleStringOption*) Manager["degenerated-groundstate"])->GetString()) == false)
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
      if (FQHEOnSphereFindSystemInfoFromVectorFileName(GroundStateFiles[i],
						       NbrParticles, LzMax, TotalLz[i], Statistics) == false)
	{
	  cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	  return -1;
	}
      if (Statistics == true)
	{
	  cout << GroundStateFiles[i] << " is not a bosonic state" << endl;
	  return -1;
	}
      if (((NbrParticles * LzMax) & 1) != (TotalLz[i] & 1))
	{
	  cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << GroundStateFiles[i] << endl;
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
#ifdef  __64_BITS__
      if ((LzMax + NbrParticles - 1) < 63)
#else
	if ((LzMax + NbrParticles - 1) < 31)	
#endif
	  {
	    if ((SymmetrizedBasis == false) || (TotalLz != 0))
	      Spaces[i] = new BosonOnSphereShort (NbrParticles, TotalLz[i], LzMax);
	    else
	      {
		Spaces[i] = new BosonOnSphereShort (NbrParticles, TotalLz[i], LzMax);
		BosonOnSphereSymmetricBasisShort TmpSpace(NbrParticles, LzMax);
		RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundStates[i], *((BosonOnSphereShort*) Spaces[i]));
		GroundStates[i] = OutputState;
	      }
	    if (Manager.GetBoolean("haldane") == true)
	      {
		int* ReferenceState = 0;
		if (((SingleStringOption*) Manager["reference-file"])->GetString() == 0)
		  {
		    cout << "error, a reference file is needed" << endl;
		    return 0;
		  }
		ConfigurationParser ReferenceStateDefinition;
		if (ReferenceStateDefinition.Parse(((SingleStringOption*) Manager["reference-file"])->GetString()) == false)
		  {
		    ReferenceStateDefinition.DumpErrors(cout) << endl;
		    return 0;
		  }
		if ((ReferenceStateDefinition.GetAsSingleInteger("NbrParticles", NbrParticles) == false) || (NbrParticles <= 0))
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
		    cout << "error while parsing ReferenceState in " << ((SingleStringOption*) Manager["reference-file"])->GetString() << endl;
		    return 0;     
		  }
		if (MaxNbrLz != (LzMax + 1))
		  {
		    cout << "wrong LzMax value in ReferenceState" << endl;
		    return 0;     
		  }
		if (((SingleStringOption*) Manager["load-hilbert"])->GetString() != 0)
		  Spaces[i] = new BosonOnSphereHaldaneBasisShort(((SingleStringOption*) Manager["load-hilbert"])->GetString());
		else
		  {
		    Spaces[i] = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz[i], LzMax, ReferenceState);	  
		  }
	      }

	  }
	else
	  {
	    if ((SymmetrizedBasis == false) || (TotalLz != 0))
	      Spaces[i] = new BosonOnSphere (NbrParticles, TotalLz[i], LzMax);
	    else
	      {
		Spaces[i] = new BosonOnSphere (NbrParticles, TotalLz[i], LzMax);
		BosonOnSphereSymmetricBasis TmpSpace(NbrParticles, LzMax);
		RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundStates[i], *((BosonOnSphere*) Spaces[i]));
		GroundStates[i] = OutputState;
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
      DensityMatrixFile << "# l_a    N    Lz    lambda" << endl;
      DensityMatrixFile.close();
    }

  ofstream File;
  if (((SingleStringOption*) Manager["output-file"])->GetString() != 0)
    File.open(((SingleStringOption*) Manager["output-file"])->GetString(), ios::binary | ios::out);
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
  if (((SingleIntegerOption*) Manager["max-la"])->GetInteger() > 0)
    {
      MeanSubsystemSize = ((SingleIntegerOption*) Manager["max-la"])->GetInteger();
      if (MeanSubsystemSize > LzMax)
	MeanSubsystemSize = LzMax;
    }
  int SubsystemSize = ((SingleIntegerOption*) Manager["min-la"])->GetInteger();
  if (SubsystemSize < 1)
    SubsystemSize = 1;
  BinomialCoefficients Coefs(MeanSubsystemSize + NbrParticles - 1);
  for (; SubsystemSize <= MeanSubsystemSize; ++SubsystemSize)
    {
      double EntanglementEntropy = 0.0;
      double DensitySum = 0.0;
      int MaxSubsystemNbrParticles = NbrParticles;
      if (MaxSubsystemNbrParticles > SubsystemSize)
	MaxSubsystemNbrParticles = SubsystemSize;
      long MaximumSize = 0;
      for (int i = 0; i <= NbrParticles; ++i)
	MaximumSize += Coefs(SubsystemSize + i - 1, i);
      double* TmpDensityMatrixEigenvalues = new double [MaximumSize];
      long TmpDensityMatrixEigenvaluePosition = 0;
      for (int SubsystemNbrParticles = 0; SubsystemNbrParticles <= NbrParticles; ++SubsystemNbrParticles)
	{
	  int SubsystemTotalLz = 0;
	  int SubsystemLzMax = SubsystemSize - 1;
	  int SubsystemMaxTotalLz = SubsystemNbrParticles * SubsystemLzMax;
	  SubsystemTotalLz = -SubsystemMaxTotalLz; 
	  for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
// 	    if (((TotalLz[0] - (SubsystemTotalLz + ((LzMax - SubsystemSize + 1) * SubsystemNbrParticles))) <= (LzMax * (NbrParticles - SubsystemNbrParticles))) || (NbrSpaces > 1))
	      {
		cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
		RealSymmetricMatrix PartialDensityMatrix = Spaces[0]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[0]);
		for (int i = 1; i < NbrSpaces; ++i)
//		  if ((TotalLz[i] - (SubsystemTotalLz + ((LzMax - SubsystemSize + 1) * SubsystemNbrParticles))) <= (LzMax * (NbrParticles - SubsystemNbrParticles)))
		    {
		      RealSymmetricMatrix TmpMatrix = Spaces[i]->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundStates[i]);
		      PartialDensityMatrix += TmpMatrix;
		    }
		if (NbrSpaces > 1)
		  PartialDensityMatrix /= ((double) NbrSpaces);
		if (PartialDensityMatrix.GetNbrRow() > 1)
		  {
		    RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
#ifdef __LAPACK__
		    if (LapackFlag == true)
		      PartialDensityMatrix.LapackDiagonalize(TmpDiag);
		    else
		      PartialDensityMatrix.Diagonalize(TmpDiag);
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

