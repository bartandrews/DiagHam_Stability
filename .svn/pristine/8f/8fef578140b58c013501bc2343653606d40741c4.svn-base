#include "config.h"

#include "Vector/RealVector.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Operator/ParticleOnSphereSquareTotalMomentumOperator.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzProjection.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include <cstring>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereLValue" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
 
  (*SystemGroup) += new SingleStringOption  ('\0', "state", "name of the file that contains the state whose average L value has to be evaluated");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total lz value of the system (0 if it has to be guessed from file name)", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "project-spin", "project onto single spin species of fermions with spin");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz-projection", "can be set to either +1 or -1", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "sz", "twice the total spin projection", 0);

  (*SystemGroup) += new SingleStringOption  ('s', "statistics", "particle statistics (bosons or fermions, try to guess it from file name if not defined)");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereLValue -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  if(Manager.GetString("state") == 0)
    {
      cout << "no input state" << endl << "see man page for option syntax or type FQHESphereLValue -h" << endl;
      return -1;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz = Manager.GetInteger("total-lz");
  bool FermionFlag = false;
  bool FixedLzFlag = true;
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  if (Manager.GetString("statistics") == 0)
    FermionFlag = true;
  if (NbrParticles==0)
    if (FQHEOnSphereFindSystemInfoFromVectorFileName(Manager.GetString("state"), NbrParticles, LzMax, TotalLz, FermionFlag) == false)
      {
	return -1;
      }
  if (Manager.GetString("statistics") != 0)
    {
      if ((strcmp ("fermions", Manager.GetString("statistics")) == 0))
	{
	  FermionFlag = true;
	}
      else
	{
	  if ((strcmp ("bosons", Manager.GetString("statistics")) == 0))
	    {
	      FermionFlag = false;
	    }
	  else
	    {
	      cout << Manager.GetString("statistics") << " is an undefined statistics" << endl;
	    }  
	}
    }
  int Parity = TotalLz & 1;
  if (Parity != ((NbrParticles * LzMax) & 1))
    {
      cout << "Lz and (NbrParticles * LzMax) must have the parity" << endl;
      return -1;           
    }

  char* StateFileName = Manager.GetString("state");
  if (IsFile(StateFileName) == false)
    {
      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
      return -1;           
    }
  RealVector State;
  if (State.ReadVector(StateFileName) == false)
    {
      cout << "error while reading " << StateFileName << endl;
      return -1;
    }


  long MemorySpace = 9l << 20;
  ParticleOnSphere* Space=0;
  if (FermionFlag == true)
    {
      if (Manager.GetBoolean("project-spin") )
	{
#ifdef __64_BITS__
	  if (LzMax <= 31)
#else
	  if (LzMax <= 15)
#endif
	    Space = new FermionOnSphereWithSpinSzProjection(NbrParticles, TotalLz, LzMax,
							    Manager.GetInteger("sz"),
							    Manager.GetInteger("sz-projection"), MemorySpace);
	  FixedLzFlag = false; // need to calculate value of Lz^2 also in this case!
	}
      else
	{
	  if (HaldaneBasisFlag == false)
	    {
#ifdef __64_BITS__
	      if (LzMax <= 62)
#else
		if (LzMax <= 30)
#endif
		  if ((SymmetrizedBasis == false) || (TotalLz != 0))
		    Space = new FermionOnSphere(NbrParticles, TotalLz, LzMax, MemorySpace);
		  else
		    Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
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
	      if (Manager.GetString("reference-file") == 0)
		{
		  ReferenceState = new int[LzMax + 1];
		  for (int i = 0; i <= LzMax; ++i)
		    ReferenceState[i] = 0;
		  if (strcasecmp(Manager.GetString("reference-state"), "laughlin") == 0)
		    for (int i = 0; i <= LzMax; i += 3)
		      ReferenceState[i] = 1;
		  else
		    if (strcasecmp(Manager.GetString("reference-state"), "pfaffian") == 0)
		      for (int i = 0; i <= LzMax; i += 4)
			{
			  ReferenceState[i] = 1;
			  ReferenceState[i + 1] = 1;
			}
		    else
		      if (strcasecmp(Manager.GetString("reference-state"), "readrezayi3") == 0)
			for (int i = 0; i <= LzMax; i += 5)
			  {
			    ReferenceState[i] = 1;
			    ReferenceState[i + 1] = 1;
			    ReferenceState[i + 2] = 1;
			  }
		      else
			{
			  cout << "unknown reference state " << Manager.GetString("reference-state") << endl;
			  return -1;
			}
		}
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
		  cout << "Reference =";
		  for (int i=0; i<=LzMax; ++i)
		    cout << " " << ReferenceState[i];
		  cout << endl;
		}
	      if (SymmetrizedBasis == false)
		{
		  if (Manager.GetString("load-hilbert") != 0)
		    Space = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
		  else
		    Space = new FermionOnSphereHaldaneBasis(NbrParticles, TotalLz, LzMax, ReferenceState, MemorySpace);
		  if (Manager.GetString("save-hilbert") != 0)
		    {
		      ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		      return 0;
		    }
		}
	      else
		{
		  if (Manager.GetString("load-hilbert") != 0)
		    Space = new FermionOnSphereHaldaneSymmetricBasis(Manager.GetString("load-hilbert"), MemorySpace);
		  else
		    Space = new FermionOnSphereHaldaneSymmetricBasis(NbrParticles, LzMax, ReferenceState, MemorySpace);
		  if (Manager.GetString("save-hilbert") != 0)
		    {
		      ((FermionOnSphereHaldaneSymmetricBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		      return 0;
		    }
		}
	    }
	}
    }
  else
    {
      if (HaldaneBasisFlag == false)
        {
#ifdef  __64_BITS__
	  if ((LzMax + NbrParticles - 1) < 63)
#else
	    if ((LzMax + NbrParticles - 1) < 31)	
#endif
	      {
		if ((SymmetrizedBasis == false) || (TotalLz != 0))
		  Space = new BosonOnSphereShort(NbrParticles, TotalLz, LzMax);
		else
		  Space = new BosonOnSphereSymmetricBasisShort(NbrParticles, LzMax);
	      }
	    else
	      {
		if ((SymmetrizedBasis == false) || (TotalLz != 0))
		  Space = new BosonOnSphere (NbrParticles, TotalLz, LzMax);
		else
		  Space = new BosonOnSphereSymmetricBasis(NbrParticles, LzMax);
	      }
	}
      else
	{
	  int* ReferenceState = 0;
          if (Manager.GetString("reference-file") == 0)
            {
              cout << "error, a reference file is needed for bosons in Haldane basis" << endl;
              return -1;
            }
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
#ifdef  __64_BITS__
          if ((LzMax + NbrParticles - 1) < 63)
#else
            if ((LzMax + NbrParticles - 1) < 31)
#endif
              Space = new BosonOnSphereHaldaneBasisShort(NbrParticles, TotalLz, LzMax, ReferenceState);
        }
    }
    
  
  if (Space->GetHilbertSpaceDimension() != State.GetVectorDimension())
    {
      cout << "dimension mismatch between the state (" << State.GetVectorDimension() << ") and the Hilbert space (" << Space->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }

  int TotalMaxLz = LzMax * NbrParticles;
  if (FermionFlag == true)
    {
      TotalMaxLz = (LzMax - NbrParticles + 1) * NbrParticles;
    }
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  ParticleOnSphereL2Hamiltonian Hamiltonian (Space, NbrParticles, LzMax, TotalLz, Architecture.GetArchitecture(), 1.0, 0, FixedLzFlag);
  RealVector TmpState(Space->GetHilbertSpaceDimension());
  VectorHamiltonianMultiplyOperation Operation (&Hamiltonian, &State, &TmpState);
  Operation.ApplyOperation(Architecture.GetArchitecture());
  double L2Value = TmpState * State;
  double RawTmpAngularMomentum = 0.5 * (sqrt ((4.0 * L2Value) + 1.0) - 1.0);
  cout << "<L^2> = " << L2Value << endl
       << "<L> = " << RawTmpAngularMomentum << endl;
  delete Space;
  return 0;
}

