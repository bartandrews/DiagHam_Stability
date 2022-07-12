#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSqueezedBasis.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinSzSymmetryLong.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetry.h"
#include "HilbertSpace/FermionOnSphereWithSpinLzSymmetryLong.h"

#include "Hamiltonian/ParticleOnSphereWithSpinGenericHamiltonian.h"
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


  OptionManager Manager ("FQHESphereFermionsWithSpin" , "0.01");
  OptionGroup* LanczosGroup  = new OptionGroup ("Lanczos options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  ParticleOnSphereManager ParticleManager(true, false, 2);
  ParticleManager.AddOptionGroup(&Manager);
  OptionGroup* SystemGroup = Manager.GetOptionGroup("system options");
  OptionGroup* PrecalculationGroup = Manager.GetOptionGroup("precalculation options");
  
  Architecture.AddOptionGroup(&Manager);
  Manager += LanczosGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbr-chebyshev", "number of times H acts on a vector)", 0);
  (*SystemGroup) += new SingleStringOption  ('\n', "input-file", "file describing the input vector");
  (*SystemGroup) += new SingleDoubleOption  ('\n', "energy-min", "lowest energy in the spectrum", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "energy-max", "highest energy in the spectrum", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "chemical-potential", "chemical potential shift, H -> H + mu", 0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "chebyshev-epsilon", "epsilon parameter for Chebyshev", 0);

  (*SystemGroup) += new SingleDoubleOption ('\n', "l2-factor", "multiplicative factor in front of an optional L^2 operator than can be added to the Hamiltonian", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "s2-factor", "multiplicative factor in front of an optional S^2 operator than can be added to the Hamiltonian", 0.0);
  
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

  (*LanczosGroup) += new  BooleanOption ('\n', "project-l2", "add a projector onto the L2 groundstate");
  (*LanczosGroup) += new  BooleanOption ('\n', "project-s2", "add a projector onto the S2 groundstate");
  (*LanczosGroup) += new  BooleanOption ('\n', "project-l2-s2", "add a projector onto the common groundstate of L2+S2");
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-storage", "additional number of vectors in RAM when using projected Lanczos", 2);
  (*LanczosGroup) += new SingleIntegerOption  ('\n', "projector-iter-max", "maximum number of iterations for internal lanczos",100);
  (*LanczosGroup) += new SingleDoubleOption ('\n', "projector-precision", "define Lanczos precision for projection (0 if automatically defined by the program)", 1e-14);
  (*LanczosGroup) += new  BooleanOption ('\n', "restart-projection", "allow lanczos projections to be restarted if full convergence not yet reached");
  
  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "initial-lz", "twice the inital momentum projection for the system", 0);
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
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpin -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int NbrFermions = Manager.GetInteger("nbr-particles");
  int LzMax = Manager.GetInteger("lzmax");
  int SzTotal = Manager.GetInteger("total-sz");
  bool HaldaneBasisFlag = Manager.GetBoolean("haldane");
  if (LONG_MAX>>20 < Manager.GetInteger("memory"))
    cout << "Warning: integer overflow in memory request - you might want to use 64 bit code."<<endl;
  long Memory = (Manager.GetInteger("memory")) << 20;  
  if (Manager.GetString("energy-expectation") != 0 ) Memory = 0x0l;
  int InitialLz = Manager.GetInteger("initial-lz");
  int NbrLz = Manager.GetInteger("nbr-lz");
  char* LoadPrecalculationFileName = Manager.GetString("load-precalculation");
  char* SavePrecalculationFileName = Manager.GetString("save-precalculation");
  bool onDiskCacheFlag = Manager.GetBoolean("allow-disk-storage");
  bool FirstRun = true;
  double** PseudoPotentials  = new double*[10];
  for (int i = 0; i < 3; ++i)
    {
      PseudoPotentials[i] = new double[LzMax + 1];
      for (int j = 0; j <= LzMax; ++j)
	PseudoPotentials[i][j] = 0.0;
    };
  double* OneBodyPotentialUpUp = 0;
  double* OneBodyPotentialDownDown = 0;

  int NbrUp = (NbrFermions + SzTotal) >> 1;
  int NbrDown = (NbrFermions - SzTotal) >> 1;
  if ((NbrUp < 0 ) || (NbrDown < 0 ))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

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
	  for (int i = 0; i < 3; ++i)
	    for (int j = 0; j < TmpNbrPseudoPotentials; ++j)
	      PseudoPotentials[i][j] = TmpPseudoPotentials[j];
	  delete []TmpPseudoPotentials;
	}
      else
	if (InteractionDefinition["Pseudopotentials"] != 0)
	  {
	    cout << "Pseudopotentials has a wrong value in " << Manager.GetString("interaction-file") << endl;
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
	    PseudoPotentials[0][j] = TmpPseudoPotentials[j];
	  delete [] TmpPseudoPotentials;
	}
      else
	if (InteractionDefinition["PseudopotentialsUpUp"] != 0)
	  {
	    cout << "PseudopotentialsUpUp has a wrong value in " << Manager.GetString("interaction-file") << endl;
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
	    PseudoPotentials[1][j] = TmpPseudoPotentials[j];
	  delete [] TmpPseudoPotentials;
	}
      else
	if (InteractionDefinition["PseudopotentialsDownDown"] != 0)
	  {
	    cout << "PseudopotentialsDownDown has a wrong value in " << Manager.GetString("interaction-file") << endl;
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
	    PseudoPotentials[2][j] = TmpPseudoPotentials[j];
	  delete [] TmpPseudoPotentials;
	}
      else
	if (InteractionDefinition["PseudopotentialsUpDown"] != 0)
	  {
	    cout << "PseudopotentialsUpDown has a wrong value in " << Manager.GetString("interaction-file") << endl;
	    return -1;
	  }
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialUpUp", ' ', OneBodyPotentialUpUp, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialUpUp has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
      if (InteractionDefinition.GetAsDoubleArray("OneBodyPotentialDownDown", ' ', OneBodyPotentialDownDown, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "OneBodyPotentialDownDown has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
	      return -1;
	    }
	}
      double *OneBodyPotentials;
      if (InteractionDefinition.GetAsDoubleArray("Onebodypotentials", ' ', OneBodyPotentials, TmpNbrPseudoPotentials) == true)
	{
	  if (TmpNbrPseudoPotentials != (LzMax + 1))
	    {
	      cout << "Onebodypotentials has a wrong number of components or has a wrong value in " << Manager.GetString("interaction-file") << endl;
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

  char* OutputNameLz = new char [512 + strlen(Manager.GetString("interaction-name"))];
  char* ExtraTerms = new char[50];
  ExtraTerms[0]='\0';
  if (Manager.GetBoolean("project-l2-s2"))
    {
      sprintf(ExtraTerms,"_Pl2-s2_%g",Manager.GetDouble("s2-factor")/Manager.GetDouble("l2-factor"));
    }
  else
    {      
      if (Manager.GetDouble("l2-factor") != 0.0)
	{
	  if (Manager.GetDouble("s2-factor") != 0.0)
	    sprintf(ExtraTerms,"_l2_%g_s2_%g",Manager.GetDouble("l2-factor"), Manager.GetDouble("s2-factor"));
	  else
	    sprintf(ExtraTerms,"_l2_%g",Manager.GetDouble("l2-factor"));
	}
      else
	if (Manager.GetDouble("s2-factor") != 0.0)
	  sprintf(ExtraTerms,"_s2_%g",Manager.GetDouble("s2-factor"));
      if (Manager.GetBoolean("project-l2"))
	{
	  sprintf(ExtraTerms,"%s_Pl2", ExtraTerms);
	}
      if (Manager.GetBoolean("project-s2"))
	{
	  sprintf(ExtraTerms,"%s_Ps2", ExtraTerms);
	}
    }
  int Max = (((LzMax - NbrUp + 1) * NbrUp) + ((LzMax - NbrDown + 1) * NbrDown));
  cout << "maximum Lz value = " << Max << endl;

  int  L = 0;
  if ((abs(Max) & 1) != 0)
     L = 1;
  if (InitialLz != 0)
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
    sprintf (OutputNameLz, "fermions_sphere_su2_%s%s_n_%d_2s_%d_sz_%d_lz_%d.dat", Manager.GetString("interaction-name"), ExtraTerms, NbrFermions, LzMax, SzTotal, L);
  else
    sprintf (OutputNameLz, "fermions_sphere_su2_%s%s_n_%d_2s_%d_sz_%d_lz.dat", Manager.GetString("interaction-name"), ExtraTerms, NbrFermions, LzMax, SzTotal);
  
  for (; L <= Max; L += 2)
    {
      double Shift = -10.0;

      ParticleOnSphereWithSpin* Space = (ParticleOnSphereWithSpin*)ParticleManager.GetHilbertSpace(L);
      
      cout << "l=" <<  L << endl;
      
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
        Memory = Architecture.GetArchitecture()->GetLocalMemory();

      AbstractQHEOnSphereWithSpinHamiltonian* Hamiltonian;      
      Hamiltonian = new ParticleOnSphereWithSpinGenericHamiltonian(Space, NbrFermions, LzMax, PseudoPotentials, OneBodyPotentialUpUp, OneBodyPotentialDownDown, NULL, 
								   Architecture.GetArchitecture(), Memory, onDiskCacheFlag, LoadPrecalculationFileName);

      if (Manager.GetDouble("s2-factor") != 0.0)
	Hamiltonian->AddS2(L, SzTotal, Manager.GetDouble("s2-factor"), ((unsigned long)Manager.GetInteger("s2-memory")) << 20);
      if (Manager.GetDouble("l2-factor") != 0.0)
	Hamiltonian->AddL2(L, SzTotal, Manager.GetDouble("l2-factor"), ((unsigned long)Manager.GetInteger("l2-memory")) << 20);



  RealVector InitialVector; 
  if (InitialVector.ReadVector(Manager.GetString("input-file")) == false)
    {
      cout << "error while reading " << Manager.GetString("input-file") << endl;
      return -1;
    }
	
   if (InitialVector.GetVectorDimension()!=Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
    cout << "Read in vector " << InitialVector.GetVectorDimension() << endl;

   //Computes <v|T_n(H')|v> for n <  2 N
    //Uses rescaled H' = (H-b)/a if provided
    //Returns [<v|T_0(H)|v>, <v|T_1(H)|v>, <v|T_2(H)|v>, ...]

   int NbrChebyshev = Manager.GetInteger("nbr-chebyshev");
   double Emin = Manager.GetDouble("energy-min");
   double Emax = Manager.GetDouble("energy-max");
   double epsilon = Manager.GetDouble("chebyshev-epsilon");
   if (epsilon == 0)
     epsilon = 16.0/NbrChebyshev;
   double ChemicalPotential = Manager.GetDouble("chemical-potential");
  
    if (ChemicalPotential != 0)
      Hamiltonian->ShiftHamiltonian(ChemicalPotential);


    double a = (Emax-Emin)/(2.0 - epsilon);
    double b = (Emax+Emin)/2.0;

    cout << "Rescaled a= " << a << " b= " << b << endl;

    RealVector alpha = InitialVector; 
    RealVector beta(Space->GetHilbertSpaceDimension(), true); 
    VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &alpha, &beta);
    Operation.ApplyOperation(Architecture.GetArchitecture());

    if (b != 0.)
        beta.AddLinearCombination(-b, alpha);
    if ((a != 1.) && (a > 1e-14))
        beta *= (1.0/a);

    double* mu = new double[NbrChebyshev];

    mu[0] = InitialVector * alpha;
    mu[1] = InitialVector * beta;


    for (int n = 2; n < NbrChebyshev; n += 2)
    	{
    		RealVector y(Space->GetHilbertSpaceDimension(), true);

    		VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &beta, &y);
		Operation.ApplyOperation(Architecture.GetArchitecture());
    		if (b !=0.)
        		y.AddLinearCombination(-b, beta);
    		if (a !=1.)
        		y *= (1.0/a);
        	y *= 2.0;
        	y.AddLinearCombination(-1.0, alpha);

    	    	alpha = beta;
        	beta = y;
        	mu[n] = 2.0 * (alpha * alpha) - mu[0]; 
        	mu[n+1] = 2.0 * (beta * alpha) - mu[1];
	    }

      for (int i = 0; i < NbrChebyshev; i++)
    	cout << "mu " << i << " : " << mu[i] << endl;
	  
      delete[] mu;  

      return 0;


    //Computes <v|T_n(H')|v> for n <  2 N
    //Uses rescaled H' = (H-b)/a if provided
    //Returns [<v|T_0(H)|v>, 2 <u|T_1(H)|v>, 2 <u|T_2(H)|v>, ...]

  //int NbrChebyshev = Manager.GetInteger("nbr-chebyshev");
  //double Emin = Manager.GetDouble("energy-min");
  //double Emax = Manager.GetDouble("energy-max");
  //double epsilon = 4.0/NbrChebyshev;
  //double ChemicalPotential = Manager.GetDouble("chemical-potential");
  
  //if (ChemicalPotential != 0)
  //  Hamiltonian->ShiftHamiltonian(ChemicalPotential);


  //  double a = (Emax-Emin)/(2.0 - 4.0 * epsilon);
  //  double b = (Emax+Emin)/2.0;

  //  cout << "Rescaled a= " << a << " b= " << b << endl;

  //  RealVector TestVector(Space->GetHilbertSpaceDimension(), true); 
  //  VectorHamiltonianMultiplyOperation TestOperation (Hamiltonian, &InitialVector, &TestVector);
  //  TestOperation.ApplyOperation(Architecture.GetArchitecture());
  //  cout << "<psi_0|H|psi_0>= " << InitialVector * TestVector << endl;	

  //  RealVector alpha = InitialVector; 
  //  RealVector beta(Space->GetHilbertSpaceDimension(), true); 
  //  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &alpha, &beta);
  //  Operation.ApplyOperation(Architecture.GetArchitecture());

  //  if (b != 0.0)
  //      beta.AddLinearCombination(-b, alpha);

  //  if (a > 1e-14)
  //      beta *= (1.0/a);
  //  else
  //     {
//	 cout << "a= " << a << " Warning: division with zero " << endl;
  //       exit(2);
    //   }   

 //   cout << "Norm of alpha= " << alpha*alpha << " Norm of beta= " << beta*beta << endl;	
 //   double* mu = new double[NbrChebyshev];

//    mu[0] = InitialVector * alpha;
//    mu[1] = InitialVector * beta;
 
  // cout << mu[0] << " " << mu[1] << endl;

    //for (int n = 2; n < NbrChebyshev; n += 2)
    	//{
    	//	RealVector y(Space->GetHilbertSpaceDimension(), true);

    	//	VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &beta, &y);
	//		Operation.ApplyOperation(Architecture.GetArchitecture());
    	//	if (b !=0.)
        //		y.AddLinearCombination(-b, beta);
    	//	if (a !=1.)
        //		y *= (1.0/a);
        //	y *= 2.0;
        //	y.AddLinearCombination(-1.0, alpha);

    	  //  alpha = beta;
        	//beta = y;
        //	mu[n] = 2.0 * (alpha * alpha) - mu[0]; 
        //	mu[n+1] = 2.0 * (beta * alpha) - mu[1];
	 //  }
//
 //   for (int i = 1; i < NbrChebyshev; i++)	    
 //	   mu[i] *= 2.0;

//    for (int i = 0; i < NbrChebyshev; i++)
  //  	cout << "mu " << i << " : " << mu[i] << endl;
	//  
//    delete[] mu;  

//    return 0;


      delete Hamiltonian;
      delete Space;
      if (FirstRun == true)
	FirstRun = false;

      if (HaldaneBasisFlag) return 0; // only one subspace defined...
    }
  delete[] OutputNameLz;
  delete[] ExtraTerms;
  for (int i = 0; i < 3; ++i)
    delete [] PseudoPotentials[i];
  delete [] PseudoPotentials;
  return 0;
}


