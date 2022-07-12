#include "HilbertSpace/BosonOnSphereWithSU4Spin.h"
#include "HilbertSpace/BosonOnSphereWithSU4SpinAllEntanglement.h"

#include "Hamiltonian/ParticleOnSphereWithSU4SpinS2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSU4SpinL2Hamiltonian.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include "LanczosAlgorithm/LanczosManager.h"
#include "HilbertSpace/ParticleOnSphereManager.h"

#include <iostream>
#include <stdlib.h>
#include <cstring>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


//#define __profiling__

#ifdef __profiling__
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#endif


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereWithSU4SpinL2Diagonalize" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(true, true, 4);
  ParticleManager.AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");  
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('z', "total-lz", "twice the total momentum projection for the system", 0);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "energy-shift", "if non zero, override energy shift using the indicated value ", 0.0);

  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag",
						"maximum Hilbert space dimension for which full diagonalization is applied",
						500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);  
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-vector", "maximum number of vector in RAM during Lanczos iteration", 10);
  (*LanczosGroup) += new BooleanOption  ('\n', "force-reorthogonalize", 
					 "force to use Lanczos algorithm with reorthogonalizion even if the number of eigenvalues to evaluate is 1", false);
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate", "evaluate eigenstates", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "eigenstate-convergence", "evaluate Lanczos convergence from eigenstate convergence", false);  
  (*LanczosGroup) += new BooleanOption  ('\n', "show-itertime", "show time spent for each Lanczos iteration", false); 
  (*LanczosGroup) += new SingleStringOption  ('\n', "initial-vector", "use file as the initial vector for the Lanczos algorithm" , 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "partial-lanczos", "only run a given number of Lanczos iterations" , false);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "lanczos-precision", "define Lanczos precision for eigenvalues (0 if automatically defined by the program)", 0);
  (*LanczosGroup) += new  BooleanOption ('\n', "fast-disk", "use disk storage to increase speed of ground state calculation and decrease memory footprint when using Lanczos algorithm");
  (*PrecalculationGroup) += new BooleanOption ('\n', "disk-cache", "use disk cache for fast multiplication", false);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-hilbert", "save Hilbert space description in the indicated file and exit (only available for the haldane or symmetrized bases)",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-hilbert", "load Hilbert space description from the indicated file (only available for the haldane or symmetrized bases)",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "show-basis", "show basis representation of the hilbert space");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);

#ifdef __profiling__
  // testing copying of arrays
  int MySize = 1<<24;
  int NbrCopies = 100;
  double *Data = new double[MySize];
  double *Copy = new double[MySize];
  NumRecRandomGenerator Generator;
  for (int i=0; i<MySize; ++i)
    Data[i] = Generator.GetRealRandomNumber();

  timeval TotalStartingTime2;
  timeval TotalEndingTime2;
  double Dt2;
  cout << "test on large arrays:"<<endl;
  cout << "------------------------------------------------------------------" << endl << endl;;
  cout << "start copy loops...";
  gettimeofday (&(TotalStartingTime2), 0);
  for (int k=0; k<NbrCopies; ++k)
    for (int i=0; i<MySize; ++i)
      Copy[i] = Data[i];
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "done" << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout << "------------------------------------------------------------------" << endl << endl;
  cout << "start memcpy...";
  gettimeofday (&(TotalStartingTime2), 0);
  for (int k=0; k<NbrCopies; ++k)
    memcpy (Copy, Data, MySize*sizeof(double));
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "done" << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout << "------------------------------------------------------------------" << endl << endl;
  cout << "start std::copy...";
  gettimeofday (&(TotalStartingTime2), 0);
  for (int k=0; k<NbrCopies; ++k)
    std::copy (&Data[0],&Data[MySize-1],Copy);
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "done" << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl<<endl;

    cout << "test on small arrays:"<<endl;
  cout << "------------------------------------------------------------------" << endl << endl;;
  cout << "start copy loops...";
  MySize = 10;
  NbrCopies = 2<<25;
  gettimeofday (&(TotalStartingTime2), 0);
  for (int k=0; k<NbrCopies; ++k)
    for (int i=0; i<MySize; ++i)
      Copy[i] = Data[i];
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "done" << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout << "------------------------------------------------------------------" << endl << endl;
  cout << "start memcpy...";
  gettimeofday (&(TotalStartingTime2), 0);
  for (int k=0; k<NbrCopies; ++k)
    memcpy (Copy, Data, MySize*sizeof(double));
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "done" << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;
  cout << "------------------------------------------------------------------" << endl << endl;
  cout << "start std::copy...";
  gettimeofday (&(TotalStartingTime2), 0);
  for (int k=0; k<NbrCopies; ++k)
    std::copy (&Data[0],&Data[MySize-1],Copy);
  gettimeofday (&(TotalEndingTime2), 0);
  cout << "done" << endl;
  Dt2 = (double) (TotalEndingTime2.tv_sec - TotalStartingTime2.tv_sec) + 
    ((TotalEndingTime2.tv_usec - TotalStartingTime2.tv_usec) / 1000000.0);
  cout << "time = " << Dt2 << endl;

  
  return 0;
#endif
  int NbrParticles = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int TotalLz  = Manager.GetInteger("total-lz");
  int TotalSz  = Manager.GetInteger("total-sz");
  int IsoSzTotal = Manager.GetInteger("total-isosz");
  int TotalEntanglement = Manager.GetInteger("total-entanglement");

//   bool LzSymmetrizedBasis = ((BooleanOption*) Manager["lzsymmetrized-basis"])->GetBoolean();
//   bool TzSymmetrizedBasis = ((BooleanOption*) Manager["tzsymmetrized-basis"])->GetBoolean();
//   bool Z3SymmetrizedBasis = ((BooleanOption*) Manager["z3symmetrized-basis"])->GetBoolean();
//   bool TzMinusParity = ((BooleanOption*) Manager["minus-tzparity"])->GetBoolean();
  unsigned long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  bool DiskCacheFlag = Manager.GetBoolean("disk-cache");
  bool FirstRun = true;  
  char* OutputNameLz, *InteractionName;


  if ((Manager.GetDouble("l2-factor")==0.0)&&(Manager.GetDouble("s2-factor")==0.0))
    {
      cout << "Either l2 or s2 need to be non-zero!"<<endl;
      return -1;
    }
  
  if (Manager.GetString("interaction-name")!=NULL)
    {
      OutputNameLz = new char [256 + strlen(Manager.GetString("interaction-name"))];
      InteractionName = Manager.GetString("interaction-name");
    }
  else
    {
      OutputNameLz = new char [369];
      InteractionName = new char[65];
      int lenFilePrefix=0;
      if (fabs(Manager.GetDouble("l2-factor"))>1e-12)
	{
	  if (fabs(Manager.GetDouble("l2-factor")-1.0)<1e-12)
	    lenFilePrefix += sprintf(InteractionName,"l2");
	  else
	    lenFilePrefix += sprintf(InteractionName,"l2_%g",Manager.GetDouble("l2-factor"));
	  if (fabs(Manager.GetDouble("s2-factor"))>1e-12)
	    lenFilePrefix += sprintf(lenFilePrefix+InteractionName,"_");
	}
      if (fabs(Manager.GetDouble("s2-factor"))>1e-12)
	{
	  if (fabs(Manager.GetDouble("s2-factor")-1.0)<1e-12)
	    sprintf(lenFilePrefix+InteractionName,"s2");
	  else
	    sprintf(lenFilePrefix+InteractionName,"s2_%g",Manager.GetDouble("s2-factor"));
	}
    }

  if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
    {  
      if (Manager.GetBoolean("use-entanglement"))
	sprintf (OutputNameLz, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz.dat", InteractionName, 
		 NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalEntanglement);
      else
	sprintf (OutputNameLz, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_lz.dat", InteractionName, 
		 NbrParticles, LzMax, TotalSz, IsoSzTotal);
    }
  else
    {
      if (Manager.GetBoolean("use-entanglement"))
	sprintf (OutputNameLz, "bosons_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz.dat", InteractionName, 
		 NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalEntanglement);
      else
	sprintf (OutputNameLz, "bosons_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_lz.dat", InteractionName, 
		 NbrParticles, LzMax, TotalSz, IsoSzTotal);
    }

  ParticleOnSphereWithSU4Spin* Space=(ParticleOnSphereWithSU4Spin*) ParticleManager.GetHilbertSpace(TotalLz);

  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    {
      cout << "Overriding memory by architecture"<<endl;
      Memory = Architecture.GetArchitecture()->GetLocalMemory();
    }
  AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian* Hamiltonian;
  if (Manager.GetDouble("l2-factor")!=0.0)
    {
      Hamiltonian = new ParticleOnSphereWithSU4SpinL2Hamiltonian(Space, NbrParticles, LzMax, TotalLz,
								 Architecture.GetArchitecture(), 
								 Manager.GetDouble("l2-factor"),
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName);
      
      if (Manager.GetDouble("s2-factor") != 0.0)
	Hamiltonian->AddS2(TotalLz, TotalSz, Manager.GetDouble("s2-factor")/Manager.GetDouble("l2-factor"), ((unsigned long)Manager.GetInteger("memory")) << 20);
      
    }
  else
    {
      Hamiltonian = new ParticleOnSphereWithSU4SpinS2Hamiltonian(Space, NbrParticles, LzMax, TotalLz, TotalSz,
								 Architecture.GetArchitecture(), 
								 Manager.GetDouble("s2-factor"),
								 Memory, DiskCacheFlag,
								 LoadPrecalculationFileName);
    }
       


  double Shift = Manager.GetDouble("energy-shift");
  Hamiltonian->ShiftHamiltonian(Shift);
  char* EigenvectorName = 0;


  if (Manager.GetBoolean("eigenstate") == true)	
    {
      /*EigenvectorName = new char[strlen(OutputNameLz)+10];
	char *TmpC=RemoveExtensionFromFileName(OutputNameLz,".dat")
	sprintf(EigenvectorName,"%s_%d",TmpC,TotalLz);
	delete [] TmpC;
       */
      
      EigenvectorName =  new char [128 + strlen(InteractionName)];
      if (strcmp ("fermions", Manager.GetString("statistics")) == 0)
	{  
	  if (Manager.GetBoolean("use-entanglement"))
	    sprintf (EigenvectorName, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz_%d", InteractionName, 
		     NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalEntanglement, TotalLz);
	  else
	    sprintf (EigenvectorName, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_lz_%d", InteractionName, 
		     NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalLz);
	}
      else
	{
	  if (Manager.GetBoolean("use-entanglement"))
	    sprintf (EigenvectorName, "bosons_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz_%d", InteractionName, 
		     NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalEntanglement, TotalLz);
	  else
	    sprintf (EigenvectorName, "bosons_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_lz_%d", InteractionName, 
		     NbrParticles, LzMax, TotalSz, IsoSzTotal, TotalLz);
	}
    }

  QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, TotalLz, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  if (EigenvectorName != 0)
    delete[] EigenvectorName;
  delete Hamiltonian;
  delete Space;
  if (FirstRun == true)
    FirstRun = false;

  return 0;
}
