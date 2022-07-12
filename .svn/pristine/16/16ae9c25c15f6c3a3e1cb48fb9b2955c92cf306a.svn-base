#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

#include "Hamiltonian/ParticleOnSphereWithSU4SpinGenericHamiltonian.h"
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

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#ifdef __MPI__
#include <mpi.h>
#endif


using std::ios;
using std::cout;
using std::endl;

using std::ofstream;

// find populations knowing total spin, isospin and entanglement projection
//
// nbrParticles = number of particles
// szTotal = spin total projection
// isoSzTotal = isospin total projection
// pzTotal = entanglement total projection
// nUpPlus = reference on the number of particles with spin up and isospin plus
// nUpMinus = reference on the number of particles with spin up and isospin minus
// nDownPlus = reference on the number of particles with spin down and isospin plus
// nDownPMinus = reference on the number of particles with spin down and isospin minus
// return value = true if values of spin, isospin and entanglement lead to valid populations
bool GetPopulations(int nbrParticles, int szTotal, int isoSzTotal, int pzTotal, int& nUpPlus, int& nUpMinus, int& nDownPlus, int& nDownMinus);


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHESphereFermionsGraphene" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += PrecalculationGroup;
  Manager += MiscGroup;
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

  (*LanczosGroup) += new  BooleanOption ('\n', "project-s2", "add a projector onto the spin S2 groundstate");
  (*LanczosGroup) += new  BooleanOption ('\n', "project-p2", "add a projector onto the isospin P2 groundstate");
  (*LanczosGroup) += new  BooleanOption ('\n', "project-l2", "add a projector onto the angular momentum L2 groundstate");
  (*LanczosGroup) += new  BooleanOption ('\n', "project-l2-s2", "add a projector onto the common groundstate of L2+S2");
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-storage", "additional number of vectors in RAM when using projected Lanczos", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-iter-max", "maximum number of iterations for internal lanczos",100);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "projector-precision", "define Lanczos precision for projection (0 if automatically defined by the program)", 1e-14);
  (*LanczosGroup) += new  BooleanOption ('\n', "restart-projection", "allow lanczos projections to be restarted if full convergence not yet reached");

  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "lzmax", "twice the maximum momentum for a single particle", 15);
  (*SystemGroup) += new SingleIntegerOption  ('s', "total-sz", "twice the z component of the total spin of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('i', "total-isosz", "twice the z component of the total isospin (i.e valley SU(2) degeneracy) of the system", 0);
  (*SystemGroup) += new BooleanOption  ('\n', "use-entanglement", "use a define value for the spin-isopsin entanglement of the system");
  (*SystemGroup) += new SingleIntegerOption  ('e', "total-entanglement", "twice the projection of the total spin-isopsin entanglement of the system", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-lz", "number of lz value to evaluate", -1);
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-file", "file describing the 2-body interaction in terms of the pseudo-potential");
  (*SystemGroup) += new  SingleStringOption ('\n', "interaction-name", "interaction name (as it should appear in output files)", "unknown");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "zeeman-field", "value for the Zeeman field that couples to the Sz quantum number", 0.0);
  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");

  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 0);
  (*PrecalculationGroup) += new BooleanOption  ('\n', "allow-disk-storage", "expand memory for fast multiplication using disk storage",false);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "load-precalculation", "load precalculation from a file",0);
  (*PrecalculationGroup) += new SingleStringOption  ('\n', "save-precalculation", "save precalculation in a file",0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEFermionsSphereWithSpin -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int LzMax = ((SingleIntegerOption*) Manager["lzmax"])->GetInteger();
  int SzTotal = ((SingleIntegerOption*) Manager["total-sz"])->GetInteger();
  int IsoSzTotal = ((SingleIntegerOption*) Manager["total-isosz"])->GetInteger();
  int TotalEntanglement = ((SingleIntegerOption*) Manager["total-entanglement"])->GetInteger();
  double Zeeman = ((SingleDoubleOption*) Manager["zeeman-field"])->GetDouble();

  long Memory = ((unsigned long) ((SingleIntegerOption*) Manager["memory"])->GetInteger()) << 20;
  unsigned long MemorySpace = ((unsigned long) ((SingleIntegerOption*) Manager["fast-search"])->GetInteger()) << 20;
  int InitialLz = ((SingleIntegerOption*) Manager["initial-lz"])->GetInteger();
  int NbrLz = ((SingleIntegerOption*) Manager["nbr-lz"])->GetInteger();
  char* LoadPrecalculationFileName = ((SingleStringOption*) Manager["load-precalculation"])->GetString();
  char* SavePrecalculationFileName = ((SingleStringOption*) Manager["save-precalculation"])->GetString();
  bool onDiskCacheFlag = ((BooleanOption*) Manager["allow-disk-storage"])->GetBoolean();
  bool FirstRun = true;
  double** PseudoPotentials  = new double*[10];
  for (int i = 0; i < 10; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
      for (int j = 0; j <= LzMax; ++j)
	PseudoPotentials[i][j] = 0.0;
    };

  int NbrUp = (NbrFermions + SzTotal) >> 1;
  int NbrDown = (NbrFermions - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  int NbrPlus = (NbrFermions + IsoSzTotal) >> 1;
  int NbrMinus = (NbrFermions - IsoSzTotal) >> 1;
  if ((NbrPlus < 0) || (NbrMinus < 0))
    {
      cout << "This value of the isospin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

  int NbrEntanglementPlus = (NbrFermions + TotalEntanglement) >> 1;
  int NbrEntanglementMinus = (NbrFermions - TotalEntanglement) >> 1;
  if ((((BooleanOption*) Manager["use-entanglement"])->GetBoolean()) && 
      ((NbrEntanglementPlus < 0) || (NbrEntanglementMinus < 0)))
    {
      cout << "This value of the entanglement projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
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
	  for (int i = 0; i < 10; ++i)
	    for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	      PseudoPotentials[i][j] = TmpPseudoPotentials[j];
	}
      else
	if (InteractionDefinition["Pseudopotentials"] != 0)
	  {
	    cout << "Pseudopotentials has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpUp", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsUpUp" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    {
	      PseudoPotentials[0][j] = TmpPseudoPotentials[j];
	      PseudoPotentials[1][j] = TmpPseudoPotentials[j];
	      PseudoPotentials[4][j] = TmpPseudoPotentials[j];
	    }
	}
      else
	if (InteractionDefinition["PseudopotentialsUpUp"] != 0)
	  {
	    cout << "PseudopotentialsUpUp has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsDownDown" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    {
	      PseudoPotentials[7][j] = TmpPseudoPotentials[j];
	      PseudoPotentials[8][j] = TmpPseudoPotentials[j];
	      PseudoPotentials[9][j] = TmpPseudoPotentials[j];
	    }
	}
      else
	if (InteractionDefinition["PseudopotentialsDownDown"] != 0)
	  {
	    cout << "PseudopotentialsDownDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpDown", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
	{
	  Flag = true;
	  if (TmpNbrPseudoPotentials != (LzMax +1))
	    {
	      cout << "Invalid number of pseudo-potentials in PseudopotentialsUpDown" << endl;
	      return -1;	  
	    }
	  for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	    {
	      PseudoPotentials[2][j] = TmpPseudoPotentials[j];
	      PseudoPotentials[3][j] = TmpPseudoPotentials[j];
	      PseudoPotentials[5][j] = TmpPseudoPotentials[j];
	      PseudoPotentials[6][j] = TmpPseudoPotentials[j];
	    }
	}
      else
	if (InteractionDefinition["PseudopotentialsUpDown"] != 0)
	  {
	    cout << "PseudopotentialsUpDown has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusUpPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpPlusUpPlus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[0][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsUpPlusUpPlus"] != 0)
          {
            cout << "PseudopotentialsUpPlusUpPlus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusUpMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpPlusUpMinus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[1][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsUpPlusUpMinus"] != 0)
          {
            cout << "PseudopotentialsUpPlusUpMinus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpPlusDownPlus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[2][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsUpPlusDownPlus"] != 0)
          {
            cout << "PseudopotentialsUpPlusDownPlus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpPlusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpPlusDownMinus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[3][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsUpPlusDownMinus"] != 0)
          {
            cout << "PseudopotentialsUpPlusDownMinus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusUpMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpMinusUpMinus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[4][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsUpMinusUpMinus"] != 0)
          {
            cout << "PseudopotentialsUpMinusUpMinus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpMinusDownPlus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[5][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsUpMinusDownPlus"] != 0)
          {
            cout << "PseudopotentialsUpMinusDownPlus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsUpMinusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsUpMinusDownMinus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[6][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsUpMinusDownMinus"] != 0)
          {
            cout << "PseudopotentialsUpMinusDownMinus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownPlusDownPlus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsDownPlusDownPlus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[7][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsDownPlusDownPlus"] != 0)
          {
            cout << "PseudopotentialsDownPlusDownPlus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownPlusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsDownPlusDownMinus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[8][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsDownPlusDownMinus"] != 0)
          {
            cout << "PseudopotentialsDownPlusDownMinus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
      if (InteractionDefinition.GetAsDoubleArray("PseudopotentialsDownMinusDownMinus", ' ', TmpPseudoPotentials, TmpNbrPseudoPotentials) == true)
        {
          Flag = true;
          if (TmpNbrPseudoPotentials != (LzMax +1))
            {
              cout << "Invalid number of pseudo-potentials in PseudopotentialsDownMinusDownMinus" << endl;
              return -1;
            }
          for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
            PseudoPotentials[9][j] = TmpPseudoPotentials[j];
        }
      else
        if (InteractionDefinition["PseudopotentialsDownMinusDownMinus"] != 0)
          {
            cout << "PseudopotentialsDownMinusDownMinus has a wrong value in " << ((SingleStringOption*) Manager["interaction-file"])->GetString() << endl;
            return -1;
          }
     delete[] TmpPseudoPotentials;
    }

  char* OutputNameLz = new char [512 + strlen(Manager.GetString("interaction-name"))];
  char* ExtraTerms = new char[50];
  ExtraTerms[0]='\0';
  if (Manager.GetBoolean("project-s2"))
    {
      sprintf(ExtraTerms,"%s_Ps2", ExtraTerms);
    }
  if (Manager.GetBoolean("project-p2"))
    {
      sprintf(ExtraTerms,"%s_Pp2", ExtraTerms);
    }
  if (Manager.GetBoolean("project-l2"))
    {
      sprintf(ExtraTerms,"%s_Pl2", ExtraTerms);
    }


  if (((BooleanOption*) Manager["use-entanglement"])->GetBoolean())
    sprintf (OutputNameLz, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
	     NbrFermions, LzMax, SzTotal, IsoSzTotal, TotalEntanglement);
  else
    sprintf (OutputNameLz, "fermions_sphere_su2su2_%s_n_%d_2s_%d_sz_%d_iz_%d_lz.dat", ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
	     NbrFermions, LzMax, SzTotal, IsoSzTotal);
  int Max = 0;

  int NbrUpPlus;
  int NbrDownPlus;
  int NbrUpMinus;
  int NbrDownMinus;
  if (((BooleanOption*) Manager["use-entanglement"])->GetBoolean())
    {
      if (GetPopulations(NbrFermions, SzTotal, IsoSzTotal, TotalEntanglement, NbrUpPlus, NbrUpMinus, NbrDownPlus, NbrDownMinus) == false)
	{
	  cout << "incompatible values of spin, isospin and entanglement" << endl;
	  return -1;
	}    
      Max  = (((LzMax - NbrUpPlus + 1) * NbrUpPlus) + ((LzMax - NbrUpMinus + 1) * NbrUpMinus) + 
	      ((LzMax - NbrDownPlus+ 1) * NbrDownPlus) + ((LzMax - NbrDownMinus + 1) * NbrDownMinus));
      cout << "populations : NUpPlus = " << NbrUpPlus << " NUpMinus  = " << NbrUpMinus << "  NDownPlus = " << NbrDownPlus << "  NDownMinus = " << NbrDownMinus << endl;
    }
  else
    {
      for (int i = -NbrFermions; i <= NbrFermions; i += 2)
	{
	  if (GetPopulations(NbrFermions, SzTotal, IsoSzTotal, i, NbrUpPlus, NbrUpMinus, NbrDownPlus, NbrDownMinus) == true)
	    {
	      int TmpMax = (((LzMax - NbrUpPlus + 1) * NbrUpPlus) + ((LzMax - NbrUpMinus + 1) * NbrUpMinus) + 
			    ((LzMax - NbrDownPlus+ 1) * NbrDownPlus) + ((LzMax - NbrDownMinus + 1) * NbrDownMinus));
	      if (TmpMax > Max)
		Max = TmpMax;
	    }
	}
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
  for (; L <= Max; L += 2)
    {
      double Shift = -10.0;
      ParticleOnSphereWithSU4Spin* Space;
#ifdef __64_BITS__
      if (LzMax <= 15)
#else
      if (LzMax <= 7)
#endif
        {
	  if (((BooleanOption*) Manager["use-entanglement"])->GetBoolean())
	    Space = new FermionOnSphereWithSU4Spin(NbrFermions, L, LzMax, SzTotal, IsoSzTotal, TotalEntanglement, MemorySpace);
	  else
	    Space = new FermionOnSphereWithSU4Spin(NbrFermions, L, LzMax, SzTotal, IsoSzTotal, MemorySpace);
        }
      else
	{
	  cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	  return -1;
	}	
      if ((((BooleanOption*) Manager["use-entanglement"])->GetBoolean()) && (Space->GetHilbertSpaceDimension() == 0))
	{
	  cout << "zero dimension Hilbert space" << endl;
	  return -1;	  
	}
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEHamiltonian* Hamiltonian;      
      Hamiltonian = new ParticleOnSphereWithSU4SpinGenericHamiltonian(Space, NbrFermions, LzMax, PseudoPotentials, Zeeman,
								       Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);      
      Hamiltonian->ShiftHamiltonian(Shift);
      if (SavePrecalculationFileName != 0)
	{
	  Hamiltonian->SavePrecalculation(SavePrecalculationFileName);
	}
      // add eventual projectors
      int NbrProjectors = 0;
      AbstractHamiltonian** Projectors = NULL;
      if (Manager.GetBoolean("project-s2")) ++NbrProjectors;
      if (Manager.GetBoolean("project-p2")) ++NbrProjectors;
      if (Manager.GetBoolean("project-l2")) ++NbrProjectors;
      Projectors = new AbstractHamiltonian*[NbrProjectors];
      NbrProjectors = 0;
      if (Manager.GetBoolean("project-s2"))
	{
	  AbstractHamiltonian* S2Projector =
	    new ParticleOnSphereWithSU4SpinS2Hamiltonian(Space, NbrFermions, LzMax, L, SzTotal,
						      Architecture.GetArchitecture(), 1.0,
						      ((unsigned long)Manager.GetInteger("s2-memory")) << 20,
						      onDiskCacheFlag);
	  S2Projector->ShiftHamiltonian(-0.25*(double)SzTotal*(SzTotal+2.0));
	  Projectors[NbrProjectors++]=S2Projector;
	}
      if (Manager.GetBoolean("project-l2"))
	{
	  AbstractHamiltonian* L2Projector =
	    new ParticleOnSphereWithSU4SpinL2Hamiltonian(Space, NbrFermions, LzMax, L,
							 Architecture.GetArchitecture(), 1.0,
							 ((unsigned long)Manager.GetInteger("l2-memory")) << 20,
							 onDiskCacheFlag);
	  L2Projector->ShiftHamiltonian(-0.25*(double)L*(L+2.0));
	  Projectors[NbrProjectors++]=L2Projector;
	}      
      if (Manager.GetBoolean("project-l2-s2"))
	{
	  AbstractQHEOnSphereWithSU4SpinCasimirHamiltonian* L2S2Projector =
	    new ParticleOnSphereWithSU4SpinL2Hamiltonian(Space, NbrFermions, LzMax, L,
							 Architecture.GetArchitecture(), 1.0,
							 ((unsigned long)Manager.GetInteger("l2-memory")) << 20,
							 onDiskCacheFlag);
	  if (Manager.GetDouble("s2-factor") != 0.0)
	    L2S2Projector->AddS2(L, SzTotal, Manager.GetDouble("s2-factor")/Manager.GetDouble("l2-factor"), ((unsigned long)Manager.GetInteger("l2-memory")) << 20);

	  L2S2Projector->ShiftHamiltonian(-0.25*(double)L*(L+2.0)-0.25*(double)SzTotal*(SzTotal+2.0));
	  Projectors[NbrProjectors++]=L2S2Projector;
	}
      char* EigenvectorName = 0;
      if (((BooleanOption*) Manager["eigenstate"])->GetBoolean() == true)	
	{
	  EigenvectorName = new char [120];
	  if (((BooleanOption*) Manager["use-entanglement"])->GetBoolean())
	    sprintf (EigenvectorName, "fermions_sphere_su4_%s_n_%d_2s_%d_sz_%d_iz_%d_pz_%d_lz_%d",
		     ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
		     NbrFermions, LzMax, SzTotal, IsoSzTotal, TotalEntanglement, L);
	  else
	    sprintf (EigenvectorName, "fermions_sphere_su2su2_%s_n_%d_2s_%d_sz_%d_iz_%d_lz_%d",
		     ((SingleStringOption*) Manager["interaction-name"])->GetString(), 
		     NbrFermions, LzMax, SzTotal, IsoSzTotal, L);
	}
      QHEOnSphereMainTask Task (&Manager, Space, Hamiltonian, L, Shift, OutputNameLz, FirstRun, EigenvectorName, LzMax);
      MainTaskOperation TaskOperation (&Task);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete Space;
      if (EigenvectorName != 0)
	{
	  delete[] EigenvectorName;
	  EigenvectorName = 0;
	}
      if (FirstRun == true)
	FirstRun = false; 
    }
  delete[] OutputNameLz;
  return 0;
}


// find populations knowing total spin, isospin and entanglement projection
//
// nbrParticles = number of particles
// szTotal = spin total projection
// isoSzTotal = isospin total projection
// pzTotal = entanglement total projection
// nUpPlus = reference on the number of particles with spin up and isospin plus
// nUpMinus = reference on the number of particles with spin up and isospin minus
// nDownPlus = reference on the number of particles with spin down and isospin plus
// nDownPMinus = reference on the number of particles with spin down and isospin minus
// return value = true if values of spin, isospin and entanglement lead to valid populations

bool GetPopulations(int nbrParticles, int szTotal, int isoSzTotal, int pzTotal, int& nUpPlus, int& nUpMinus, int& nDownPlus, int& nDownMinus)
{
  nUpPlus = nbrParticles + szTotal + isoSzTotal + pzTotal;
  nUpMinus = nbrParticles + szTotal - isoSzTotal - pzTotal;
  nDownPlus = nbrParticles - szTotal + isoSzTotal - pzTotal;
  nDownMinus = nbrParticles - szTotal - isoSzTotal + pzTotal;
  if ((nUpPlus < 0) || ((nUpPlus & 0x3) != 0) || (nUpMinus < 0) || ((nUpMinus & 0x3) != 0) ||
      (nDownPlus < 0) || ((nDownPlus & 0x3) != 0) || (nDownMinus < 0) || ((nDownMinus & 0x3) != 0))
    return false;
  nUpPlus >>= 2;
  nUpMinus >>= 2;
  nDownPlus >>= 2;
  nDownMinus >>= 2;
  return true;
}
