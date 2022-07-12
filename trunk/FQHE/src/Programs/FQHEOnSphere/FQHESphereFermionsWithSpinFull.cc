#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/FermionOnSphereWithSpinAllSz.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSqueezedBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"

#include "Hamiltonian/ParticleOnSphereWithSpinFullGenericHamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinL2Hamiltonian.h"
#include "Hamiltonian/ParticleOnSphereWithSpinS2Hamiltonian.h"
#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "MainTask/QHEOnSphereMainTask.h"

#include "Options/Options.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstdlib>
#include <climits>
#include <cmath>
#include <cstring>
#include <sys/time.h>
#include <stdio.h>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::ios;
using std::cout;
using std::endl;

using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  const int NbrPseudopotentialsCoeffs = 10;

  OptionManager Manager ("FQHESphereFermionsWithSpinFull" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(true, false, -2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleDoubleOption ('\n', "subband", "manual entry for the subband gap (overrides OneBodyPseudopotentials)",-1.0);
  
  (*LanczosGroup) += new SingleIntegerOption  ('n', "nbr-eigen", "number of eigenvalues", 30);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "full-diag", 
                                                "maximum Hilbert space dimension for which full diagonalization is applied", 
                                                500, true, 100);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "iter-max", "maximum number of lanczos iteration", 3000);
  (*LanczosGroup)  += new BooleanOption  ('\n', "block-lanczos", "use block Lanczos algorithm", false);
  (*LanczosGroup)  += new SingleIntegerOption  ('\n', "block-size", "size of the block used in the block Lanczos algorithm", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "limit-time", "use limit in time instead of a number of lanczos iteration (0 if none, time in seconds)", 0);
  (*LanczosGroup)  += new BooleanOption  ('d', "disk", "enable disk resume capabilities", false);
  (*LanczosGroup) += new BooleanOption  ('r', "resume", "resume from disk datas", false);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "nbr-iter", "number of lanczos iteration (for the current run)", 10);
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
  (*LanczosGroup) += new  BooleanOption ('\n', "resume-fastdisk", "resume the fast-disk mode Lanczos algorithm from a stopped one (for example due to computer crash)");
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "s2-memory", "amount of memory that can be allocated for fast multiplication of s2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "l2-memory", "amount of memory that can be allocated for fast multiplication of l2 term (in Mbytes)", 500);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  
  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinFull -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  double SubbandGap = Manager.GetDouble("subband");
  bool ManualSubbandGap = (SubbandGap != -1.0 );
  if (LONG_MAX>>20 < Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  long Memory = (Manager.GetInteger("memory")) << 20;  
  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0l;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();
  bool onDiskCacheFlag = ((BooleanOption*) Manager["allow-disk-storage"])->GetBoolean();
  bool FirstRun = true;
  double** PseudoPotentials  = new double*[NbrPseudopotentialsCoeffs];
  for (int i = 0; i < NbrPseudopotentialsCoeffs; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
      for (int j = 0; j <= LzMax; ++j)
        PseudoPotentials[i][j] = 0.0;
    };
  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;
  double* OneBodyPotentialUpDown = 0;
  double* OneBodyPotentialDownUp = 0;

  if (((SingleStringOption*) Manager["interaction-file"])->GetString() == 0)
    {
      cout << "an interaction file has to be provided" << endl;
      return -1;
    }
  else
    {
      ConfigurationParser InteractionDefinition;
      if (InteractionDefinition.Parse(((SingleStringOption*) Manager["interaction-file"])->GetString()) == false)
        {
          InteractionDefinition.DumpErrors(cout) << endl;
          return -1;
        }
      int TmpNbrPseudoPotentials;
      double* TmpPseudoPotentials;
      bool Flag = false;

      if (InteractionDefinition.GetAsDoubleArray("Pseudopotentials", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in Pseudopotentials" << endl;
              return -1;          
            }
          for (int i = 0; i < NbrPseudopotentialsCoeffs; ++i)
            for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
              PseudoPotentials[i][j] = TmpPseudoPotentials[j];
          delete []TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["Pseudopotentials"] != 0)
          {
            cout << "Pseudopotentials has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }


// UpUpUpUp Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpUpUpUp", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpUpUpUp" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[0][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsUpUpUpUp"] != 0)
          {
            cout << "PseudopotentialsUpUpUpUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// UpDownUpUp Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpDownUpUp", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpDownUpUp" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[1][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsUpDownUpUp"] != 0)
          {
            cout << "PseudopotentialsUpDownUpUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// DownDownUpUp Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownDownUpUp", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsDownDownUpUp" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[2][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsDownDownUpUp"] != 0)
          {
            cout << "PseudopotentialsDownDownUpUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// UpUpUpDown Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpUpUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpUpUpDown" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[3][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsUpUpUpDown"] != 0)
          {
            cout << "PseudopotentialsUpUpUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// UpDownUpDown Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpDownUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpDownUpDown" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[4][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsUpDownUpDown"] != 0)
          {
            cout << "PseudopotentialsUpDownUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// DownUpUpDown Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownUpUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsDownUpUpDown" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[5][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsDownUpUpDown"] != 0)
          {
            cout << "PseudopotentialsDownUpUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// DownDownUpDown Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownDownUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsDownDownUpDown" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[6][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsDownDownUpDown"] != 0)
          {
            cout << "PseudopotentialsDownDownUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// UpUpDownDown Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpUpDownDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpUpDownDown" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[7][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsUpUpDownDown"] != 0)
          {
            cout << "PseudopotentialsUpUpDownDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// UpDownDownDown Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpDownDownDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpDownDownDown" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[8][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsUpDownDownDown"] != 0)
          {
            cout << "PseudopotentialsUpDownDownDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
// DownDownDownDown Sector
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownDownDownDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsDownDownDownDown" << endl;
              return -1;          
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[9][j] = TmpPseudoPotentials[j];
          delete [] TmpPseudoPotentials;
        }
      else
        if (InteractionDefinition["PseudopotentialsDownDownDownDown"] != 0)
          {
            cout << "PseudopotentialsDownDownDownDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }

// One body potentials
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', OneBodyPotentialUpUp, TmpNbrPseudoPotentials) == true)
        {
          if (TmpNbrPseudoPotentials != (LzMax + 1))
            {
              cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
              return -1;
            }
        }
      else OneBodyPotentialUpUp=NULL;
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', OneBodyPotentialDownDown, TmpNbrPseudoPotentials) == true)
        {
          if (TmpNbrPseudoPotentials != (LzMax + 1))
            {
              cout << "OneBodyPotentialDownDown has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
              return -1;
            }
        }
      else OneBodyPotentialDownDown=NULL;
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpDown", ' ', OneBodyPotentialUpDown, TmpNbrPseudoPotentials) == true)
        {
          if (TmpNbrPseudoPotentials != (LzMax + 1))
            {
              cout << "OneBodyPotentialUpDown has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
              return -1;
            }
        }
      else OneBodyPotentialUpDown=NULL;
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownUp", ' ', OneBodyPotentialDownUp, TmpNbrPseudoPotentials) == true)
        {
          if (TmpNbrPseudoPotentials != (LzMax + 1))
            {
              cout << "OneBodyPotentialDownUp has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
              return -1;
            }
        }
      else OneBodyPotentialDownUp=NULL;

      double *OneBodyPotentials;
      if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
        {
          if (TmpNbrPseudoPotentials != (LzMax + 1))
            {
              cout << "Onebodypotentials has a wrong number of components or has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
              return -1;
            }
          if (OneBodyPotentialUpUp==NULL)
            {
              OneBodyPotentialUpUp = new double [LzMax+1];
              for (int i=0; i<=LzMax; ++i)
                OneBodyPotentialUpUp[i]=OneBodyPotentials[i];
            }
          if (OneBodyPotentialDownDown==NULL)
            OneBodyPotentialDownDown = OneBodyPotentials;
          else delete [] OneBodyPotentials;
        }
    }
//===========================================================================================
  char* OutputNameLz = new char [512 + strlen(Manager.GetString("interaction-name"))];
  char* ExtraTerms = new char[50];
  ExtraTerms[0]='\0';

    if (ManualSubbandGap)
      {
	OneBodyPotentialUpUp = new double [LzMax+1];
	for (int i=0; i<=LzMax; ++i)
	  OneBodyPotentialUpUp[i]=0;
	OneBodyPotentialDownDown = new double [LzMax+1];
	for (int i=0; i<=LzMax; ++i)
	  OneBodyPotentialDownDown[i]=SubbandGap;
	OneBodyPotentialUpDown = new double [LzMax+1];
	for (int i=0; i<=LzMax; ++i)
	  OneBodyPotentialUpDown[i]=0;
	OneBodyPotentialDownUp = new double [LzMax+1];
	for (int i=0; i<=LzMax; ++i)
	  OneBodyPotentialDownUp[i]=0;
      sprintf(ExtraTerms,"_subband_%g",SubbandGap);
      }
  
  // Maximum value of the total angular momentum depends on the parity of NbrFermions
  int Max;
  if( NbrFermions % 2 )
    {
      Max=(LzMax-(NbrFermions-1)/2+1)*(NbrFermions-1)/2+(LzMax-(NbrFermions+1)/2+1)*(NbrFermions+1)/2; // NbrFermions odd
    }
  else 
    {
      Max=((LzMax-NbrFermions/2+1)*NbrFermions); // NbrFermions even
    }

  cout << "maximum Lz value = " << Max << endl;

  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz >= 0)
    {
      L = InitialLz;
      if ((abs(Max) & 1) != 0)
        L |= 1;
      else
        L &= ~0x1;
    }
  if (NbrLz > 0)
    {
      if (L + (2 * (NbrLz - 1)) < Max)
        Max = L + (2 * (NbrLz - 1));
    }

  if (NbrLz==1) 
    sprintf (OutputNameLz, "fermions_sphere_su2_%s%s_n_%d_2s_%d_lz_%d.dat", Manager.GetString("interaction-name"), ExtraTerms, NbrFermions, LzMax, L);
  else
    sprintf (OutputNameLz, "fermions_sphere_su2_%s%s_n_%d_2s_%d_lz.dat", Manager.GetString("interaction-name"), ExtraTerms, NbrFermions, LzMax);
 
  for (; L <= Max; L += 2)
    {
      double Shift = -10.0;

      ParticleOnSphereWithSpin* Space = (ParticleOnSphereWithSpin*)ParticleManager.GetHilbertSpace(L);
      
      cout << "l=" <<  L << endl;
      
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEOnSphereWithSpinFullHamiltonian* Hamiltonian;      
      Hamiltonian = new ParticleOnSphereWithSpinFullGenericHamiltonian(Space, NbrFermions, LzMax, PseudoPotentials, OneBodyPotentialUpUp, OneBodyPotentialDownDown, OneBodyPotentialUpDown, OneBodyPotentialDownUp , Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);

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
          cout << "< shifted energy > = "<<EnergyValue + Shift<<endl;
          return 0;
        }
      // add eventual projectors
      int NbrProjectors = 0;
      AbstractHamiltonian** Projectors = NULL;
      
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
        {
          Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
        }
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)       
        {
          EigenvectorName = new char [512];
          sprintf (EigenvectorName, "fermions_sphere_su2_%s%s_n_%d_2s_%d_lz_%d",
                   ((SingleStringOption*) Manager["interaction-name"])->GetString(), ExtraTerms,
                   NbrFermions, LzMax, L);
        }
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax, Projectors, NbrProjectors);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      for (int p=0; p<NbrProjectors;++p)
        delete Projectors[p];
      delete [] Projectors;
      delete Hamiltonian;
      delete Space;      
      if (EigenvectorName != 0)
        {
          delete[] EigenvectorName;
          EigenvectorName = 0;
        }
      if (FirstRun == true)
        FirstRun = false;
      //if (HaldaneBasisFlag) return 0; // only one subspace defined...
    }
  delete[] OutputNameLz;
  delete[] ExtraTerms;
  for (int i = 0; i < NbrPseudopotentialsCoeffs; ++i)
    delete [] PseudoPotentials[i];
  delete [] PseudoPotentials;
  return 0;
}
