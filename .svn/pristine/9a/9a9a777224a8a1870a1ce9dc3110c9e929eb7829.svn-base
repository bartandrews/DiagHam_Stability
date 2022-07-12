#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "Hamiltonian/ParticleOnSphereGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereL2Hamiltonian.h"
//#include "Hamiltonian/ProjectedQHEHamiltonian.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasis.h"
#include "HilbertSpace/FermionOnSphereLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasisLong.h"
#include "HilbertSpace/FermionOnSphereSymmetricBasisLong.h"
#include "HilbertSpace/FermionOnSphereHaldaneSymmetricBasisLong.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Tools/FQHEFiles/FQHESpherePseudopotentialTools.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cstring>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>

using std::ios;
using std::cout;
using std::endl;
using std::ofstream;



int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsTwoBodyGeneric" , "0.01");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);
  ParticleOnSphereManager ParticleManager(true, false, 1);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");  
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");

  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);  
  Manager += ToolsGroup;
  Manager += MiscGroup;
  OptionGroup* LanczosGroup = Manager.GetOptionGroup("Lanczos options");

  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new BooleanOption  ('\n', "general-onebody", "add a general one-body potential (of type c_m^+ c_n, for all combinations of m, n); should be used with option all-lz");
  (*SystemGroup) += new  SingleStringOption ('\n', "onebody-potential", "file specifying general one-body potential (each row gives m, n, U_{mn})", "unknown");
  (*SystemGroup) += new BooleanOption  ('g', "ground", "restrict to the largest subspace");
  (*SystemGroup) += new BooleanOption  ('\n', "all-lz", "use full Hilbert space (all Lz sectors)");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new BooleanOption ('\n', "l2-only", "compose Hamiltonian only of L2 terms");
  (*SystemGroup) += new BooleanOption  ('\n', "get-lvalue", "compute mean l value from <L^2> for each eigenvalue");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");
  
  (*LanczosGroup) += new  BooleanOption ('\n', "project-l2", "add a projector onto the L2 groundstate");
    (*LanczosGroup) += new SingleIntegerOption  ('\n', "l2-memory", "precalculation memory for L^2 operator",1000);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "l2-nbr-vectors", "number of states stored for L^2 projection",10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-storage", "additional number of vectors in RAM when using projected Lanczos", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-iter-max", "maximum number of iterations for internal lanczos",100);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "projector-precision", "define Lanczos precision for projection (0 if automatically defined by the program)", 1e-14);
  (*LanczosGroup) += new  BooleanOption ('\n', "restart-projection", "allow lanczos projections to be restarted if full convergence not yet reached");

  (*PrecalculationGroup) += new BooleanOption ('\n', "no-hermitian", "do not use hermitian symmetry of the hamiltonian");
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption('\n', "energy-variance", "in addition to energy expectation, also evaluate energy variance sqrt[<H^2>-<H>^2]");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHEFermionsTwoBodyGeneric -h" << endl;
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
  if (ULONG_MAX>>20 < (unsigned long)Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0l;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;
  double* PseudoPotentials = new double[LzMax + 1];
  double* OneBodyPotentials = 0;
  if (Manager.GetString("interaction-file") == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      if (FQHESphereGetPseudopotentials(Manager.GetString("interaction-file"), LzMax, PseudoPotentials, OneBodyPotentials) == false)
	return -1;
    }
  char* InteractionName = 0;

  if (Manager.GetBoolean("l2-only"))
    {      
      InteractionName=new char[3];
      sprintf (InteractionName, "l2");
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(Manager.GetString("interaction-file")) == false)
	{
	  InteractionDefinition.DumpErrors(cout) << endl;
	  exit(-1);
	}
      if (InteractionDefinition["Name"] != 0)
	{
	  InteractionName = new char[strlen(InteractionDefinition["Name"]) + 1];
	  strcpy(InteractionName, InteractionDefinition["Name"]);
	}
      else
	{
	  InteractionName = new char[32 + strlen(Manager.GetString("interaction-name"))];
	  strcpy(InteractionName, Manager.GetString("interaction-name"));
	  InteractionName=new char[15 + strlen(Manager.GetString("interaction-name"))];
	  if (Manager.GetDouble("l2-factor") != 0.0)
	    sprintf (InteractionName, "%s_l2_%g",Manager.GetString("interaction-name"), Manager.GetDouble("l2-factor"));
	  else
	    strcpy (InteractionName, Manager.GetString("interaction-name"));
	}
    }

  int* OneBodyMValues = 0;
  int* OneBodyNValues = 0;
  double* OneBodyPotentialValues = 0;
  int NbrGeneralOneBodyPotentials;
  bool HaveGeneralOneBodyPotentials = Manager.GetBoolean("general-onebody");

  if (HaveGeneralOneBodyPotentials)
   {
     cout << "Using general one-body pseudopotentials U_{mn} c_m^+ c_n " << endl;
     cout << "Warning: should be used with all-lz option because the potential does not conserve Lz!" << endl;
 
     std::ifstream infile(Manager.GetString("onebody-potential"));
     if (infile.is_open()) 
      {
        int m1, m2;
        double Um1m2Re, Um1m2Im;
        int NbrTerms = 0; 
        while (infile >> m1 >> m2 >> Um1m2Re >> Um1m2Im)
         {  
           if (fabs(Um1m2Re) > 1e-12)
             NbrTerms++; 
	   if (fabs(Um1m2Im) > 1e-12)
  	     {
               cout << "One body potential has imaginary part... exiting." << endl;
	       exit(1);
             }
         }
        //cout << "Number of one-body terms: " << NbrTerms << endl;

        OneBodyMValues = new int[NbrTerms];
        OneBodyNValues = new int[NbrTerms];
        OneBodyPotentialValues = new double[NbrTerms];
        NbrGeneralOneBodyPotentials = NbrTerms;

	NbrTerms = 0;
        infile.clear();
        infile.seekg(0, ios::beg);
        while (infile >> m1 >> m2 >> Um1m2Re >> Um1m2Im)
         {     
	   if (fabs(Um1m2Re) > 1e-12)
            {
              OneBodyMValues[NbrTerms] = m1;
              OneBodyNValues[NbrTerms] = m2;
              OneBodyPotentialValues[NbrTerms] = Um1m2Re;
              NbrTerms++; 
            }
	   if (fabs(Um1m2Im) > 1e-12)
  	     {
               cout << "One body potential has imaginary part... exiting." << endl;
	       exit(1);
             }
         }
   
        //for (int i = 0; i < NbrTerms; i++)
        //  cout << OneBodyMValues[i] << " " << OneBodyNValues[i] << " : " << OneBodyPotentialValues[i] << " ; ";
        //cout << endl; 
      }
    infile.close();
   }


  char* OutputNameLz = new char [256 + strlen(InteractionName)];
  char* ExtraTerms = new char[50];
  ExtraTerms[0]='\0';
  if (Manager.GetBoolean("project-l2"))
    sprintf(ExtraTerms,"_Pl2");
  sprintf (OutputNameLz, "fermions_%s%s_n_%d_2s_%d_lz.dat", InteractionName, ExtraTerms, NbrParticles, LzMax);

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

  if (Manager.GetBoolean("all-lz")) 
   {
     Max = L;
   }

  for (; L <= Max; L += 2)
    {
      cout << "Lz="<<L<<endl;
      ParticleOnSphere* Space;

      if (Manager.GetBoolean("all-lz")) 
          Space = new FermionOnSphereFull(NbrParticles, LzMax);
      else
          Space = (FermionOnSphere*) ParticleManager.GetHilbertSpace(L);

      if (Space==0) return 0; // happens if we wrote the Hilbert space to disk, for instance!
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      AbstractQHEHamiltonian* Hamiltonian = 0;
      double Shift=0.0;
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      if (Manager.GetBoolean("l2-only"))
	{
	  Hamiltonian = new ParticleOnSphereL2Hamiltonian(Space, NbrParticles, LzMax, L, 
								  Architecture.GetArchitecture(), 1.0,
							  ((unsigned long)Manager.GetInteger("l2-memory")) << 20,
							  false, false, NULL, !Manager.GetBoolean("no-hermitian"));
	}
      else
	{
	  if (HaveGeneralOneBodyPotentials)
              Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials, OneBodyPotentials,
                                                                 NbrGeneralOneBodyPotentials, OneBodyMValues, OneBodyNValues, OneBodyPotentialValues,
								 Manager.GetDouble("l2-factor"),
								 Architecture.GetArchitecture(), 
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName,
								 !Manager.GetBoolean("no-hermitian"));
          else
           { 
  	     if (OneBodyPotentials == 0)
	       Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials,
								 Manager.GetDouble("l2-factor"),
								 Architecture.GetArchitecture(), 
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName,
								 !Manager.GetBoolean("no-hermitian"));
	    else	
              Hamiltonian = new ParticleOnSphereGenericHamiltonian(Space, NbrParticles, LzMax, PseudoPotentials, OneBodyPotentials,
								 Manager.GetDouble("l2-factor"),
								 Architecture.GetArchitecture(), 
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName,
								 !Manager.GetBoolean("no-hermitian"));
          }	  
	  Shift = - 0.5 * ((double) (NbrParticles * NbrParticles)) / (0.5 * ((double) LzMax));
	}
      
    if ( (Manager.GetString("energy-expectation") != 0 ) || (Manager.GetBoolean("energy-variance") != 0 ) )
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
	  double EnergyValue = State * TmpState;
	  cout << "< Energy > = "<<EnergyValue<<endl;
	  cout << "< shifted energy > = "<<EnergyValue + Shift<<endl;

      if (Manager.GetBoolean("energy-variance") != 0 )
       {
   	     RealVector TmpState2(Space->GetHilbertSpaceDimension());
	     VectorHamiltonianMultiplyOperation Operation2 (Hamiltonian, &TmpState, &TmpState2);
	     Operation2.ApplyOperation(Architecture.GetArchitecture());
	     double varH = State * TmpState2 - EnergyValue * EnergyValue;
	     cout << "(varH)^2 = " << varH << endl;
       }   
	  return 0;
  	 }
      
      Hamiltonian->ShiftHamiltonian(Shift);
//       AbstractQHEHamiltonian *Projector=NULL;
//       if (Manager.GetBoolean("l2-project"))
// 	{
// 	  AbstractQHEHamiltonian *TmpHamiltonian = Hamiltonian;
// 	  Projector = new ParticleOnSphereL2Hamiltonian(Space, NbrParticles, LzMax, L,
// 							Architecture.GetArchitecture(), /* l2Factor */ 1.0, 
// 							((unsigned long)Manager.GetInteger("l2-memory")) << 20);
// 	  Hamiltonian = new ProjectedQHEHamiltonian (TmpHamiltonian, Projector, Architecture.GetArchitecture(),
// 						     Manager.GetInteger("l2-nbr-vectors"));//,  int maxIterProj);
// 	}

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
	  sprintf (EigenvectorName, "fermions_%s%s_n_%d_2s_%d_lz_%d", InteractionName, ExtraTerms, NbrParticles, LzMax, L);
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
      delete Space;
      if (FirstRun == true)
	FirstRun = false;
    }
  delete [] PseudoPotentials;
  if (HaveGeneralOneBodyPotentials)
   {
     delete[] OneBodyMValues;
     delete[] OneBodyNValues;
     delete[] OneBodyPotentialValues;
   }
  delete [] OutputNameLz;
  delete [] ExtraTerms;
  return 0;
}
