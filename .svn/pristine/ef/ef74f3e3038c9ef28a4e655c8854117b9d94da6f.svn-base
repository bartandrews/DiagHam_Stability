#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/RealMatrix.h"

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
#include "Options/SingleDoubleOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("FQHESphereBosonsTruncatedSchmidtDecomposition" , "0.01");
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
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "huge-basis", "use huge Hilbert space support");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('a', "subsystem-size", "subsystem size i.e. number of orbitals which defines the geometrical partition", 1);
  (*SystemGroup) += new SingleDoubleOption  ('c', "cut-off", "minus log of the cut-off to apply to the reduced density matrix eigenvalues", 32);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name instead of the one that can be deduced from the input file name (replacing the vec extension with trunc.vec extension");
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereBosonsTruncatedSchmidtDecomposition -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if (Manager.GetString("ground-file") == 0)
    {
      cout << "error, a ground state file should be provided. See man page for option syntax or type FQHESphereBosonsTruncatedSchmidtDecomposition -h" << endl;
      return -1;
    }

  bool HaldaneBasisFlag = ((BooleanOption*) Manager["haldane"])->GetBoolean();
  bool SymmetrizedBasis = ((BooleanOption*) Manager["symmetrized-basis"])->GetBoolean();
  int NbrParticles = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger(); 
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger(); 
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  double EigenvalueCutOff = exp(- Manager.GetDouble("cut-off"));
#ifdef __LAPACK__
  bool LapackFlag = ((BooleanOption*) Manager["use-lapack"])->GetBoolean();
#endif
  int TotalLz = 0;
  bool Statistics = true;
  ParticleOnSphere* Space = 0;

  if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("ground-file"),
						   NbrParticles, LzMax, TotalLz, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << Manager.GetString("ground-file") << endl;
      return -1;
    }
  if (Statistics == true)
    {
      cout <<Manager.GetString("ground-file")  << " is not a fermionic state" << endl;
      return -1;
    }
  if (((NbrParticles * LzMax) & 1) != (TotalLz & 1))
    {
      cout << "incompatible values for nbr-particles, nbr-flux and total-lz for ground state file " << Manager.GetString("ground-file") << endl;
      return -1;
    }


  RealVector GroundState;  
  if (GroundState.ReadVector (Manager.GetString("ground-file")) == false)
    {
      cout << "can't open vector file " << Manager.GetString("ground-file") << endl;
      return -1;      
    }
  RealVector SchmidtDecomposedState(GroundState.GetLargeVectorDimension(), true);

#ifdef  __64_BITS__
      if ((LzMax + NbrParticles - 1) < 63)
#else
	if ((LzMax + NbrParticles - 1) < 31)	
#endif
	  {
	    if (Manager.GetBoolean("huge-basis") == true)
	      {
		if (Manager.GetString("load-hilbert") == 0)
		  {
		    cout << "error : huge basis mode requires to save and load the Hilbert space" << endl;
		    return -1;
		  }
		Space = new  BosonOnSphereHaldaneHugeBasisShort (Manager.GetString("load-hilbert"), Manager.GetInteger("memory"));
	      }
	    else
	      {
		if (Manager.GetBoolean("haldane") == false)
		  {
		    if ((SymmetrizedBasis == false) || (TotalLz != 0))
		      Space = new BosonOnSphereShort (NbrParticles, TotalLz, LzMax);
		    else
		      {
			Space = new BosonOnSphereShort (NbrParticles, TotalLz, LzMax);
			BosonOnSphereSymmetricBasisShort TmpSpace(NbrParticles, LzMax);
			RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundState, *((BosonOnSphereShort*) Space));
			GroundState = OutputState;
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
			cout << "error while parsing ReferenceState in " << Manager.GetString("reference-file") << endl;
			return 0;     
		      }
		    if (MaxNbrLz != (LzMax + 1))
		      {
			cout << "wrong LzMax value in ReferenceState" << endl;
			return 0;     
		      }
		    if (Manager.GetString("load-hilbert") != 0)
		      Space = new BosonOnSphereHaldaneBasisShort(Manager.GetString("load-hilbert"));
		    else
		      {
			Space = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);	  
		      }
		  }
	      }
	  }
	else
	  {
	    if ((SymmetrizedBasis == false) || (TotalLz != 0))
	      Space = new BosonOnSphereLong (NbrParticles, TotalLz, LzMax);
	    else
	      {
		Space = new BosonOnSphere (NbrParticles, TotalLz, LzMax);
		BosonOnSphereSymmetricBasis TmpSpace(NbrParticles, LzMax);
		RealVector OutputState = TmpSpace.ConvertToNbodyBasis(GroundState, *((BosonOnSphere*) Space));
		GroundState = OutputState;
	      }
	  }

      if (Space->GetLargeHilbertSpaceDimension() != GroundState.GetLargeVectorDimension())
	{
	  cout << "dimension mismatch between Hilbert space and ground state" << endl;
	  return 0;
	}
    

  cout.precision(14);
  int SubsystemSize = Manager.GetInteger("subsystem-size");
  double EntanglementEntropy = 0.0;
  double DensitySum = 0.0;
  double TruncatedEntanglementEntropy = 0.0;
  double TruncatedDensitySum = 0.0;

  int MaxSubsystemNbrParticles = NbrParticles;
  if (MaxSubsystemNbrParticles > SubsystemSize)
    MaxSubsystemNbrParticles = SubsystemSize;

  for (int SubsystemNbrParticles = 0; SubsystemNbrParticles <= NbrParticles; ++SubsystemNbrParticles)
    {
      int SubsystemTotalLz = 0;
      int SubsystemLzMax = SubsystemSize - 1;
      int SubsystemMaxTotalLz = (SubsystemNbrParticles * SubsystemLzMax);
      SubsystemTotalLz = -SubsystemMaxTotalLz; 
      for (; SubsystemTotalLz <= SubsystemMaxTotalLz; SubsystemTotalLz += 2)
	{
            cout << "processing subsystem size=" << SubsystemSize << "  subsystem nbr of particles=" << SubsystemNbrParticles << " subsystem total Lz=" << SubsystemTotalLz << endl;
	    RealSymmetricMatrix PartialDensityMatrix = Space->EvaluatePartialDensityMatrix(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, GroundState);

//            if ((Manager.GetString("full-densitymatrix") != 0) && (FilterNa == SubsystemNbrParticles) && (FilterLza = SubsystemTotalLz))
//	      {
//		  ofstream FullDensityMatrixFile;
//		  FullDensityMatrixFile.open(Manager.GetString("full-densitymatrix"), ios::binary | ios::out | ios::app); 
//		  FullDensityMatrixFile << PartialDensityMatrix;
//		  FullDensityMatrixFile.close();
//	      }

	  if (PartialDensityMatrix.GetNbrRow() > 1)
	    {
	      RealDiagonalMatrix TmpDiag (PartialDensityMatrix.GetNbrRow());
	      RealMatrix TmpEigenstates(PartialDensityMatrix.GetNbrRow(),
					PartialDensityMatrix.GetNbrRow(), true);
	      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		TmpEigenstates[i][i] = 1.0;
#ifdef __LAPACK__
	      if (LapackFlag == true)
		PartialDensityMatrix.LapackDiagonalize(TmpDiag, TmpEigenstates);
	      else
#endif
		PartialDensityMatrix.Diagonalize(TmpDiag, TmpEigenstates);	      
	      TmpDiag.SortMatrixDownOrder(TmpEigenstates);
	      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		{
		  if (TmpDiag[i] > 1e-14)
		    {
		      EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		      DensitySum += TmpDiag[i];
		    }
		  if (TmpDiag[i] > EigenvalueCutOff)
		    {
		      TruncatedEntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		      TruncatedDensitySum += TmpDiag[i];			
		    }		    
		}
	      Space->EvaluatePartialSchmidtDecomposition(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, EigenvalueCutOff,
							 GroundState, SchmidtDecomposedState, TmpDiag, TmpEigenstates);
	    }
	  else
	    if (PartialDensityMatrix.GetNbrRow() == 1)
	      {
		double TmpValue = PartialDensityMatrix(0,0);
		if (TmpValue > 1e-14)
		  {
		    EntanglementEntropy += TmpValue * log(TmpValue);
		    DensitySum += TmpValue;
		  }
		if (TmpValue > EigenvalueCutOff)
		  {
		    TruncatedEntanglementEntropy += TmpValue * log(TmpValue);
		    TruncatedDensitySum += TmpValue;			
		  }		    
		RealDiagonalMatrix TmpDiag (1);
		RealMatrix TmpEigenstates(1, 1);
		TmpEigenstates(0, 0) = 1.0;
		TmpDiag[0] = TmpValue;
		Space->EvaluatePartialSchmidtDecomposition(SubsystemSize, SubsystemNbrParticles, SubsystemTotalLz, EigenvalueCutOff,
							   GroundState, SchmidtDecomposedState, TmpDiag, TmpEigenstates);
	      }
	}
    }

  cout << endl << "------------------------------------------------------" << endl << endl;
  cout << "reduced density matrix trace = " << DensitySum << " (truncated = " << TruncatedDensitySum << ")"<< endl;
  cout << "entanglement entropy = " << (-EntanglementEntropy) << " (truncated = " << (-TruncatedEntanglementEntropy) << ")"<< endl;
  double TmpNorm = SchmidtDecomposedState.Norm();
  cout << "truncated state normalization = " << TmpNorm << endl;
  SchmidtDecomposedState /= TmpNorm;

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

