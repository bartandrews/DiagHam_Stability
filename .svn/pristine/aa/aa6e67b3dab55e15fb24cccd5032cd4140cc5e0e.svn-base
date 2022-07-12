#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"

#include "Hamiltonian/ParticleOnSphereGenericThreeBodyHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsThreeBodyGeneric" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
  OptionGroup* LanczosGroup = Manager.GetOptionGroup("Lanczos options");

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 8);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new BooleanOption  ('\n', "haldane", "use Haldane basis instead of the usual n-body basis");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "reference state to start the Haldane algorithm from (can be laughlin, pfaffian or readrezayi3)", "laughlin");
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-file", "use a file as the definition of the reference state");
  (*SystemGroup) += new BooleanOption  ('\n', "symmetrized-basis", "use Lz <-> -Lz symmetrized version of the basis (only valid if total-lz=0)");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "l2-memory", "precalculation memory for L^2 operator",1000);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "l2-nbr-vectors", "number of states stored for L^2 projection",10);
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  (*LanczosGroup) += new  BooleanOption ('\n', "project-l2", "add a projector onto the L2 groundstate");
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-storage", "additional number of vectors in RAM when using projected Lanczos", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-iter-max", "maximum number of iterations for internal lanczos",100);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "projector-precision", "define Lanczos precision for projection (0 if automatically defined by the program)", 1e-14);
  (*LanczosGroup) += new  BooleanOption ('\n', "restart-projection", "allow lanczos projections to be restarted if full convergence not yet reached");


  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the Haldane basis)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the Haldane basis)",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsThreeBodyGeneric -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool GroundFlag = Manager.GetBoolean("ground");
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0l;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  bool SymmetrizedBasis = Manager.GetBoolean("symmetrized-basis");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");  
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  double* PseudoPotentials = 0;
  double* OneBodyPotentials = 0;
  double* ThreeBodyPotentials = 0;
  bool Normalize3Body = 0;
  int TmpNbrThreeBodyPseudoPotentials = 0;
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  return -1;
	}
      if (InteractionDefinition.GetAsDoubleArray("ThreebodyPseudopotentials", ' ', ThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials) == false)
	{
	  cout << "ThreebodyPseudopotentials are not defined or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	  return -1;
	}
      if (InteractionDefinition["NormalizeThreeBody"]!=NULL)
	InteractionDefinition.GetAsBoolean("NormalizeThreeBody",Normalize3Body);
      int TmpNbrPseudoPotentials;
      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', PseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials" << endl;
	      return -1;	  
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of onebody pseudo-potentials" << endl;
	      return -1;	  
	    }
	}
    }

  char* OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
  char* ExtraTerms = new char[50];
  ExtraTerms[0]='\0';
  if (Manager.GetBoolean("project-l2"))
    sprintf(ExtraTerms,"_Pl2");
  sprintf (OutputNameLz, "fermions_%s%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), ExtraTerms, NbrParticles, LzMax);

  int Max = ((LzMax - NbrParticles + 1) * NbrParticles);

  int  L = InitialLz;
  if (L < -Max)
    L = -Max;
  else
    if (L > Max)
      L = Max;
  if ((abs(Max) & 1) != (abs(InitialLz) & 1))
    L += 1;
  if (GroundFlag == true)
      Max = L;
  else
    {
      if (NbrLz > 0)
	{
	  Max = L + (2 * (NbrLz - 1));
	}
    }
  for (; L <= Max; L += 2)
    {
      ParticleOnSphere* Space = 0;
      if (HaldaneBasisFlag == false)
	{
#ifdef __64_BITS__
	  if (LzMax <= 62)
#else
	  if (LzMax <= 30)
#endif
	    if ((SymmetrizedBasis == false) || (L != 0))
	      Space = new FermionOnSphere(NbrParticles, L, LzMax, MemorySpace);
	    else
	      {
		if (Manager.GetString("load-hilbert") != 0)
		  Space = new FermionOnSphereSymmetricBasis(Manager.GetString("load-hilbert"), MemorySpace);
		else
		  Space = new FermionOnSphereSymmetricBasis(NbrParticles, LzMax, MemorySpace);
		if (Manager.GetString("save-hilbert") != 0)
		  {
		    ((FermionOnSphereSymmetricBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
		    return 0;
		  }
	      }
	  else
#ifdef __128_BIT_LONGLONG__
	    if (LzMax <= 126)
#else
	      if (LzMax <= 62)
#endif
		{
		  if ((SymmetrizedBasis == false) || (L != 0))
		    Space = new FermionOnSphereLong(NbrParticles, L, LzMax, MemorySpace);
		  else
		    {
		      if (Manager.GetString("load-hilbert") != 0)
			Space = new FermionOnSphereSymmetricBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
		      else
			Space = new FermionOnSphereSymmetricBasisLong(NbrParticles, LzMax, MemorySpace);
		      if (Manager.GetString("save-hilbert") != 0)
			{
			  ((FermionOnSphereSymmetricBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			  return 0;
			}
		    }
		}
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, L, LzMax, MemorySpace);
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
	    }
	  if (SymmetrizedBasis == false)
	     {
#ifdef __64_BITS__
	       if (LzMax <= 62)
#else
		 if (LzMax <= 30)
#endif
		   {
		     if (Manager.GetString("load-hilbert") != 0)
		       Space = new FermionOnSphereHaldaneBasis(Manager.GetString("load-hilbert"), MemorySpace);
		     else
		       Space = new FermionOnSphereHaldaneBasis(NbrParticles, L, LzMax, ReferenceState, MemorySpace);
		     if (Manager.GetString("save-hilbert") != 0)
		       {
			 ((FermionOnSphereHaldaneBasis*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			 return 0;
		       }
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
			 Space = new FermionOnSphereHaldaneBasisLong(NbrParticles, L, LzMax, ReferenceState, MemorySpace);
		       if (Manager.GetString("save-hilbert") != 0)
			 {
			   ((FermionOnSphereHaldaneBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			   return 0;
			 }
		     }	       
	     }
	  else
	    {
#ifdef __64_BITS__
	       if (LzMax <= 62)
#else
		 if (LzMax <= 30)
#endif
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
		 else
#ifdef __128_BIT_LONGLONG__
		   if (LzMax <= 126)
#else
		     if (LzMax <= 62)
#endif
		       {
			 if (Manager.GetString("load-hilbert") != 0)
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(Manager.GetString("load-hilbert"), MemorySpace);
			 else
			   Space = new FermionOnSphereHaldaneSymmetricBasisLong(NbrParticles, LzMax, ReferenceState, MemorySpace);
			 if (Manager.GetString("save-hilbert") != 0)
			   {
			     ((FermionOnSphereHaldaneSymmetricBasisLong*) Space)->WriteHilbertSpace(Manager.GetString("save-hilbert"));
			     return 0;
			   }
		       }
	    }
	}
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEOnSphereHamiltonian* Hamiltonian = 0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Hamiltonian = new ParticleOnSphereGenericThreeBodyHamiltonian(Space, NbrParticles, LzMax, ThreeBodyPotentials, TmpNbrThreeBodyPseudoPotentials - 1,
								    ((SingleDoubleOption*) Manager["l2-factor"])->GetDouble(), PseudoPotentials, OneBodyPotentials,
								      Architecture.GetArchitecture(), 
								    Memory, DiskCacheFlag,
								    LoadPrecalculationFileName, Normalize3Body);

      double Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
      if (Manager.GetString("energy-expectation") != 0 )
	{
	  char* StateFileName = Manager.GetString("energy-expectation");
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
	  if (State.GetVectorDimension()!=Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
	  RealVector TmpState(Space->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  double EnergyValue = State*TmpState;
	  cout << "< Energy > = "<<EnergyValue<<endl;
	  if (false)
	  {
	    AbstractQHEOnSphereHamiltonian* L2Operator = new ParticleOnSphereL2Hamiltonian(Space, NbrParticles, LzMax, L, 
											 Architecture.GetArchitecture(), 1.0,
											 0x250ul << 20,
											   false, false, NULL);
	    RealVector TmpState2(Space->GetHilbertSpaceDimension());
	    RealVector TmpStateA(Space->GetHilbertSpaceDimension());
	    RealVector TmpState2A(Space->GetHilbertSpaceDimension());
	    VectorHamiltonianMultiplyOperation Operation2 (L2Operator, &State, &TmpState2);
	    Operation2.ApplyOperation(Architecture.GetArchitecture());
	    VectorHamiltonianMultiplyOperation Operation3 (L2Operator, &TmpState, &TmpStateA);
	    Operation3.ApplyOperation(Architecture.GetArchitecture());
	    VectorHamiltonianMultiplyOperation Operation4 (Hamiltonian, &TmpState2, &TmpState2A);
	    Operation4.ApplyOperation(Architecture.GetArchitecture());
	    TmpStateA/=TmpStateA.Norm();
	    TmpState2A/=TmpState2A.Norm();
	    double Overlap = TmpState2A*TmpStateA;
	    cout << "< Overlap > = "<<Overlap<<endl;
	  }
	  cout << "< shifted energy > = "<<EnergyValue + Shift<<endl;
	  return 0;
	}

      Hamiltonian->ShiftHamiltonian(Shift);

            // add eventual projectors
      int NbrProjectors = 0;
      AbstractHamiltonian** Projectors = NULL;
      if (Manager.GetBoolean("project-l2")) ++NbrProjectors;
      Projectors = new AbstractHamiltonian*[NbrProjectors];
      NbrProjectors = 0;
      if (Manager.GetBoolean("project-l2"))
	{
	  AbstractHamiltonian* L2Projector =
	    new ParticleOnSphereL2Hamiltonian(Space, NbrParticles, LzMax, L, 
					      Architecture.GetArchitecture(), 1.0, ((long)Manager.GetInteger("l2-memory"))<<20);
	  L2Projector->ShiftHamiltonian(-0.25*(double)L*(L+2.0));
	  Projectors[NbrProjectors++]=L2Projector;
	}

      
      char* EigenvectorName = 0;
      if (Manager.GetBoolean("eigenstate") == true)	
	{
	  EigenvectorName = new char [64];
	  sprintf (EigenvectorName, "fermions_%s%s_n_%d_2s_%d_lz_%d", Manager.GetString("interaction-name"), ExtraTerms, NbrParticles, LzMax, L);
	}
      
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax, Projectors, NbrProjectors);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	}
      for (int p=0; p<NbrProjectors;++p)
	delete Projectors[p];
      delete [] Projectors;
      delete Hamiltonian;
      if (FirstRun == true)
	FirstRun = false;
    }
  delete [] OutputNameLz;
  delete [] ExtraTerms;
  return 0;
}
