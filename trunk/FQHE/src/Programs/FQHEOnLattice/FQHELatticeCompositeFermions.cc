#include "HilbertSpace/BosonOnLattice.h"
#include "HilbertSpace/BosonOnLatticeLong.h"
#include "HilbertSpace/FermionOnLattice.h"
#include "HilbertSpace/HardCoreBosonOnLattice.h"
#include "HilbertSpace/HardCoreBosonLong.h"
#include "Hamiltonian/ParticleOnLatticeDeltaHamiltonian.h"
#include "Operator/ParticleOnLatticeTranslationOperator.h"

#include "Matrix/ComplexMatrix.h"
#include "Matrix/ComplexLapackDeterminant.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "MainTask/QHEOnLatticeMainTask.h"

#include "MathTools/IntegerAlgebraTools.h"
#include "MathTools/JacobiThetaFunction.h"

#include "GeneralTools/FilenameTools.h"

#include "Options/Options.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <cmath>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


#ifdef HAVE_GSL
#include "gsl/gsl_multimin.h"
#endif


// some global variables to avoid too many parameters in function calls
int NbrBosons;
int Lx;
int Ly;
ArchitectureManager Architecture;

// Manager = OptionManager indicating all contents
// Space = total space to work in
// solenoidCF_X = solenoid flux for CF in x-direction
// solenoidCF_y = solenoid flux for CF in y-direction
// writeState = flag indicating whether final state shall be written to predefined output name
ComplexVector* GetTrialState(OptionManager &Manager, ParticleOnLattice* Space, double solenoidCF_x, double solenoidCF_y, bool verbose=false, bool writeState=false);

double GetTrialStateSqrOverlap(ComplexVector &reference, OptionManager &Manager, ParticleOnLattice* Space, double solenoidCF_x, double solenoidCF_y, bool verbose=false, bool writeState=false)
{
  ComplexVector *TmpV = GetTrialState(Manager, Space, solenoidCF_x, solenoidCF_y, verbose, writeState);
  double Overlap = SqrNorm(reference* (*TmpV));
  delete TmpV;
  return Overlap;
}

struct OptParams
{
  ComplexVector *ReferenceVector;
  OptionManager *Manager;
  ParticleOnLattice* Space;
  double Offset;
  int EvaluationsF;
  int EvaluationsdF;
};

#ifdef HAVE_GSL

double TargetFunction_F (const gsl_vector * X, void * PARAMS)
{
  OptParams *p=(OptParams*)PARAMS;
  double rst = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0), gsl_vector_get(X,1));
  ++(p->EvaluationsF);
  return rst;
}

void TargetFunction_dF (const gsl_vector * X, void * PARAMS, gsl_vector * G)
{
  OptParams *p=(OptParams*)PARAMS;
  double h = p->Offset;
  double Oplusx = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0)+h, gsl_vector_get(X,1));
  double Ominusx = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0)-h, gsl_vector_get(X,1));
  double Oplusy =  -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0), gsl_vector_get(X,1)+h);
  double Ominusy = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0), gsl_vector_get(X,1)-h);
  gsl_vector_set (G, 0, (Oplusx-Ominusx)/(2.0*h));
  gsl_vector_set (G, 1, (Oplusy-Ominusy)/(2.0*h));
  ++(p->EvaluationsdF);
}

  
void TargetFunction_FdF (const gsl_vector * X, void * PARAMS, double * f, gsl_vector * G)
{
  OptParams *p=(OptParams*)PARAMS;
  double h = p->Offset;
  *f = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0), gsl_vector_get(X,1));
  double Oplusx = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0)+h, gsl_vector_get(X,1));
  double Ominusx = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0)-h, gsl_vector_get(X,1));
  double Oplusy =  -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0), gsl_vector_get(X,1)+h);
  double Ominusy = -GetTrialStateSqrOverlap(*(p->ReferenceVector),*(p->Manager), p->Space, gsl_vector_get(X,0), gsl_vector_get(X,1)-h);
  ++(p->EvaluationsF);
  ++(p->EvaluationsdF);
  gsl_vector_set (G, 0, (Oplusx-Ominusx)/(2.0*h));
  gsl_vector_set (G, 1, (Oplusy-Ominusy)/(2.0*h));
}

#endif // HAVE_GSL

// values = eigenvalues
// vectors = complex matrix with many vectors to be analyzed and overwritten with basis diagonal in k
// NbrFlux = number of flux quanta to be considered
// start, end = indices of range of eigenvalues to be considered
ComplexMatrix& DiagonalizeMomentaInSubspace(RealDiagonalMatrix &values, ComplexMatrix &vectors, ParticleOnLatticeTranslationOperator *translationOperator, int NbrFlux, int start, int end);

void GetTranslationMatrix(ParticleOnLatticeTranslationOperator *Operator, int NbrVectors,
			  ComplexVector *Vectors, ComplexMatrix &MatrixRepresentation,
			  ComplexVector &TmpState, ArchitectureManager &Architecture);


int main(int argc, char** argv)
{
  cout.precision(14);

  OptionManager Manager ("FQHELatticeCompositeFermions" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OptimizationGroup = new OptionGroup ("optimization options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");  

  Manager += SystemGroup;
  Manager += OptimizationGroup;
  Architecture.AddOptionGroup(&Manager);  
  Manager += PrecalculationGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 3);
  (*SystemGroup) += new SingleIntegerOption  ('x', "lx", "length in x-direction of given lattice", 6);
  (*SystemGroup) += new SingleIntegerOption  ('y', "ly", "length in y-direction of given lattice", 9);
  (*SystemGroup) += new SingleIntegerOption  ('q', "flux", "number of flux quanta piercing the lattice ", 6);
  (*SystemGroup) += new SingleIntegerOption  ('f', "flux-per-CF", "number of flux attached to each boson (allowed values: +/-1)", 1);
  (*SystemGroup) += new BooleanOption('c',"hard-core","Use Hilbert-space of hard-core bosons");

  (*SystemGroup) += new MultipleDoubleOption  ('\n', "solenoid-flux", "twist in periodic boundary conditions for total wavefunction phi_x[,phi_y])",',');
  (*SystemGroup) += new MultipleDoubleOption  ('s', "solenoid-CF", "twist in periodic boundary conditions for CF part phi_x[,phi_y])",',');
  (*SystemGroup) += new SingleStringInternalOption('\n',"CF","externally supply single particle states for CF basis (base)",NULL, true);
  (*SystemGroup) += new SingleStringInternalOption('\n',"all-CF","externally supply single particle states (list of full file-names)",NULL, true);
  (*SystemGroup) +=
    new SingleStringInternalOption('\n',"J","externally supply single particle states for Jastrow basis (base)",NULL, true);

  (*OptimizationGroup) += new SingleStringOption  ('\n', "optimize", "vector file to optimize CF solenoid flux against",NULL);

#ifdef HAVE_GSL
  (*OptimizationGroup) += new SingleDoubleOption  ('\n', "offset", "offset h to be used for numerical derivatives",0.005);
  (*OptimizationGroup) += new BooleanOption  ('\n', "opt-gradient", "use gradient-method instead of simplex algorithm");
  
#endif
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*MiscGroup) += new BooleanOption  ('d', "omit-diag", "omit diagonalizing in momentum basis");
  (*MiscGroup) += new BooleanOption  ('a', "analytic", "also generate the analytic wavefunctions");
  (*MiscGroup) += new BooleanOption  ('\n', "write-basis", "write the single particle basis states that were used");
  (*MiscGroup) += new BooleanOption  ('\n', "write-product", "write the product states of pairs of basis states");
  (*MiscGroup) += new BooleanOption  ('\n', "write-slater", "write the slater determinant part of the wavefunction");
  (*MiscGroup) += new BooleanOption  ('\n', "write-jastrow", "write the jastrow factor part of the wavefunction");
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  NbrBosons = Manager.GetInteger("nbr-particles");
  Lx = Manager.GetInteger("lx");
  Ly = Manager.GetInteger("ly");
  int NbrFluxQuanta = Manager.GetInteger("flux");  
  bool HardCore = Manager.GetBoolean("hard-core");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;


  double SolenoidCF_X=0.0, SolenoidCF_Y=0.0;
  {
    int tmpI;
    double *Fluxes=Manager.GetDoubles("solenoid-CF", tmpI);
    if (tmpI>0) SolenoidCF_X=Fluxes[0];
    if (tmpI>1) SolenoidCF_Y=Fluxes[1];
    if (tmpI>0) delete [] Fluxes;	
  }


  cout << "* Full Hilbert-space: N="<<NbrBosons<<" bosons in "<<Lx<<" x "<<Ly<<" cells at N_phi="<<NbrFluxQuanta<<endl;
  ParticleOnLattice* Space;
  if (HardCore)
    {
      if (Lx*Ly + NbrBosons - 1 < 63)
	{
	  Space =new HardCoreBosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
	}	
      else		
	{			
	  cout <<"Not implemented Hilberspace"<<endl;
	  return -1;
	}
    }	
  else 
    {
      if (Lx*Ly + NbrBosons - 1 < 63)
	Space = new BosonOnLattice(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
      else	
	Space = new BosonOnLatticeLong(NbrBosons, Lx, Ly, NbrFluxQuanta, MemorySpace);
    }	
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());


  if (Manager.GetString("optimize")!=NULL)
    {
      if ((Manager.GetString("CF")!=NULL)||(Manager.GetString("all-CF")!=NULL)||(Manager.GetString("J")!=NULL))
	{
	  cout << "Cannot optimize with externally supplied basis, ignoring the given basis."<<endl;
	  ((SingleStringInternalOption*)Manager["CF"])->SetString(NULL);
	  ((SingleStringInternalOption*)Manager["all-CF"])->SetString(NULL);
	  ((SingleStringInternalOption*)Manager["J"])->SetString(NULL);
	}
      ComplexVector ReferenceState;

      if (ReferenceState.ReadVector(Manager.GetString("optimize"))==false)
	{
	  cout << "error reading reference state "<<Manager.GetString("optimize")<<endl;
	  exit(-1);
	}

      if (ReferenceState.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	{
	  cout << "Reference state "<<Manager.GetString("optimize")<< " does not match the dimension of the given Hilbert-space"<<endl;
	  exit(-1);
	}

      double Onod = GetTrialStateSqrOverlap(ReferenceState, Manager, Space, SolenoidCF_X, SolenoidCF_Y);
      cout << "Checking initial O="<<Onod<<endl;

#ifdef HAVE_GSL
      // common assignments
      OptParams MyParameters;
      MyParameters.ReferenceVector=&ReferenceState;
      MyParameters.Manager = &Manager;
      MyParameters.Space = Space;
      MyParameters.Offset = Manager.GetDouble("offset");
      MyParameters.EvaluationsF=0;
      MyParameters.EvaluationsdF=0;

      if (Manager.GetBoolean("opt-gradient"))
	{
	  size_t iter = 0;
	  int status;

	  const gsl_multimin_fdfminimizer_type *T;
	  gsl_multimin_fdfminimizer *s;
      
	  gsl_vector *x;
	  gsl_multimin_function_fdf my_func;

	  my_func.n = 2;
	  my_func.f = &TargetFunction_F;
	  my_func.df = &TargetFunction_dF;
	  my_func.fdf = &TargetFunction_FdF;
	  my_func.params = &MyParameters;

	  /* Starting point, x = (5,7) */
	  x = gsl_vector_alloc (2);
	  gsl_vector_set (x, 0, SolenoidCF_X);
	  gsl_vector_set (x, 1, SolenoidCF_Y);
      
	  T = gsl_multimin_fdfminimizer_conjugate_fr;
	  s = gsl_multimin_fdfminimizer_alloc (T, 2);

	  gsl_multimin_fdfminimizer_set (s, &my_func, x, 0.01, 1e-4);

	  do
	    {
	      iter++;
	      status = gsl_multimin_fdfminimizer_iterate (s);

	      if (status)
		break;

	      status = gsl_multimin_test_gradient (s->gradient, 1e-3);

	      if (status == GSL_SUCCESS)
		printf ("Minimum found at:\n");

	      printf ("%5d %.5f %.5f %10.5f\n", (int)iter,
		      gsl_vector_get (s->x, 0),
		      gsl_vector_get (s->x, 1),
		      s->f);
	    }
	  while (status == GSL_CONTINUE && iter < 100);

	  SolenoidCF_X = gsl_vector_get(s->x, 0);
	  SolenoidCF_Y = gsl_vector_get (s->x, 1);
	  cout << "Final overlap at theta=( "<<SolenoidCF_X<<","<<SolenoidCF_Y<<" ), O="<<-s->f<<endl;
	  cout << "Derivatives =( "<<SolenoidCF_X<<","<<SolenoidCF_Y<<" ), O="<<-s->f<<endl;
	  cout << "Total number of function evaluations:        "<<MyParameters.EvaluationsF<<endl;
	  cout << "Total number of evaluations for derivative : "<<MyParameters.EvaluationsdF<<endl;
	  cout << "Total points where overlap was calculated :  "<<MyParameters.EvaluationsF+4*MyParameters.EvaluationsdF<<endl;
	  gsl_multimin_fdfminimizer_free (s);
	  gsl_vector_free (x);
	}
      else // use simplex minimizer
	{
	  const gsl_multimin_fminimizer_type *T =
	    gsl_multimin_fminimizer_nmsimplex;
	  gsl_multimin_fminimizer *s = NULL;
	  gsl_vector *ss, *x;
	  gsl_multimin_function minex_func;

	  size_t iter = 0;
	  int status;
	  double size;

	  /* Starting point */
	  x = gsl_vector_alloc (2);
	  gsl_vector_set (x, 0, SolenoidCF_X);
	  gsl_vector_set (x, 1, SolenoidCF_Y);

	  /* Set initial step sizes to 1 */
	  ss = gsl_vector_alloc (2);
	  gsl_vector_set_all (ss, 1.0);

	  /* Initialize method and iterate */
	  minex_func.n = 2;
	  minex_func.f = &TargetFunction_F;
	  minex_func.params = &MyParameters;
	  
	  s = gsl_multimin_fminimizer_alloc (T, 2);
	  gsl_multimin_fminimizer_set (s, &minex_func, x, ss);

	  do
	    {
	      iter++;
	      status = gsl_multimin_fminimizer_iterate(s);

	      if (status)
		break;

	      size = gsl_multimin_fminimizer_size (s);
	      status = gsl_multimin_test_size (size, 5e-3);

	      if (status == GSL_SUCCESS)
		{
		  printf ("converged to minimum at\n");
		}

	      printf ("%5d %10.3e %10.3ef f() = %7.3f size = %.3f\n",
		      (int)iter,
		      gsl_vector_get (s->x, 0),
		      gsl_vector_get (s->x, 1),
		      s->fval, size);
	    }
	  while (status == GSL_CONTINUE && iter < 100);
	  
	  SolenoidCF_X = gsl_vector_get(s->x, 0);
	  SolenoidCF_Y = gsl_vector_get (s->x, 1);
	  cout << "Final overlap at theta=( "<<SolenoidCF_X<<","<<SolenoidCF_Y<<" ), O="<<-s->fval<<endl;
	  cout << "Total number of function evaluations:        "<<MyParameters.EvaluationsF<<endl;
	  
	  gsl_vector_free(x);
	  gsl_vector_free(ss);
	  gsl_multimin_fminimizer_free (s);

	}
      
#else
      

      cout << "Attention, this minimizer still need debugging!"<<endl;
      double h = Manager.GetDouble("offset");
      bool Convergence=false;
      int iteration = 0;
      while (Convergence==false)
	{
	  cout << "Calculating for theta=("<<SolenoidCF_X<<", "<<SolenoidCF_Y<<")"<<endl;
	  double Onod = GetTrialStateSqrOverlap(ReferenceState, Manager, Space, SolenoidCF_X, SolenoidCF_Y);
	  double Oplusx = GetTrialStateSqrOverlap(ReferenceState, Manager, Space, SolenoidCF_X+h, SolenoidCF_Y);
	  double Ominusx = GetTrialStateSqrOverlap(ReferenceState, Manager, Space, SolenoidCF_X-h, SolenoidCF_Y);
	  double Oplusy =  GetTrialStateSqrOverlap(ReferenceState, Manager, Space, SolenoidCF_X, SolenoidCF_Y+h);
	  double Ominusy = GetTrialStateSqrOverlap(ReferenceState, Manager, Space, SolenoidCF_X, SolenoidCF_Y-h);

	  cout << "Iteration "<<iteration++<<": theta=("<<SolenoidCF_X<<", "<<SolenoidCF_Y<<"), O="<<Onod<<endl;
	  
	  double dOx = (Oplusx-Ominusx)/(2.0*h);
	  double dOy = (Oplusy-Ominusy)/(2.0*h);
	  double d2Ox = (Oplusx+Ominusx-2.0*Onod)/(2.0*h*h);
	  double d2Oy = (Oplusy+Ominusy-2.0*Onod)/(2.0*h*h);
	  Complex Direction = Complex(dOx,dOy);
	  double Angle = Arg(Direction);
	  double alpha = 0.05; // damping parameter 0<=alpha<1
	  double Radius=0.0;
	  double cs = cos(Angle), sn = sin(Angle);	  

	  cout << "O="<<Onod<<" O+x="<<Oplusx<<" O-x="<<Ominusx<<" O+y="<<Oplusy<<" O-y="<<Ominusy<<endl;
	  cout << "dOx="<<dOx<<" dOy="<<dOy<<" d2Ox="<<d2Ox<<" d2Oy="<<d2Oy<<" Angle="<<Angle<<endl;
	  Convergence = ((sqrt(dOx*dOx+dOy*dOy)<1e-6)&&(d2Ox<0.0)&&(d2Oy<0.0));
	  
	  if (Convergence == false)
	    {
	      if ((sqrt(dOx*dOx+dOy*dOy)>1e-6))
		{
		  double a = 0.5*(cs*cs*d2Ox+sn*sn*d2Oy);
		  double b = cs*dOx + sn *dOy;
		  double c = (1.0-alpha)*(Onod-1.0);
		  cout << "a="<<a<<" b="<<b<<" c="<<c<<endl;
		  if (b<0) cout << " wrong assumptions..."<<endl;
		  if ((fabs(a)>1e-6) && (b*b-4.0*a*c>0)) // can use quadratic relation?
		    {
		      if (a>0)
			Radius = (-b+sqrt(b*b-4.0*a*c))/(2.0*a);
		      else
			Radius = (-b-sqrt(b*b-4.0*a*c))/(2.0*a);
		    }
		  else // use linear interpolation
		    { 
		      Radius = c/b;
		    }
		  cout << "Angle=" <<Angle<<" Radius = "<<Radius<<endl;
		}
	      else
		{
		  if ((d2Ox>0.0)&&(d2Oy>0.0))
		    {
		      Direction = Complex(d2Ox,d2Oy);
		      Angle = Arg(Direction);
		      cs = cos(Angle);
		      sn = sin(Angle);
		      double a = 0.5*(cs*cs*d2Ox+sn*sn*d2Oy);
		      double c = (1.0-alpha)*(Onod-1.0);		      
		      Radius = sqrt(-c/a);
		    }
		  else if (d2Ox<0.0)
		    {
		      cs=0.0; sn=1.0;
		      Radius = sqrt(-(1.0-alpha)*(Onod-1.0)/d2Oy);
		    }
		  else if (d2Oy<0.0)
		    {
		      cs=1.0; sn=0.0;
		      Radius = sqrt(-(1.0-alpha)*(Onod-1.0)/d2Oy);
		    }
		  else
		    {
		      cout << "This case should not happen!"<<endl;
		    }
		}
	    }
	  
	  SolenoidCF_X += Radius * cs;
	  SolenoidCF_Y += Radius * sn;

	}

#endif
      cout << "Writing final state to disk for : theta = ( "<<SolenoidCF_X<<","<<SolenoidCF_Y<<" )"<<endl;
      GetTrialState(Manager, Space, SolenoidCF_X, SolenoidCF_Y, /*verbose*/ true, /*writeState*/ true);
      
    }
  else GetTrialState(Manager, Space, SolenoidCF_X, SolenoidCF_Y, /*verbose*/ true, /*writeState*/ true);

  delete Space;
  
  return 0;
}



ComplexVector* GetTrialState(OptionManager &Manager, ParticleOnLattice* Space, double solenoidCF_X, double solenoidCF_Y, bool verbose, bool writeState)
{
  int NbrFluxQuanta = Manager.GetInteger("flux");
  int CFFlux = Manager.GetInteger("flux-per-CF");
  bool HardCore = Manager.GetBoolean("hard-core");
  bool NoMomentumDiagonalize = Manager.GetBoolean("omit-diag");
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;

  char boundaryCdStr[30]="";
  double SolenoidX=0.0, SolenoidY=0.0;
  {
    int tmpI;
    double *Fluxes=Manager.GetDoubles("solenoid-flux", tmpI);
    if (tmpI>0) SolenoidX=Fluxes[0];
    if (tmpI>1) SolenoidY=Fluxes[1];
    
    if (tmpI>0)
      {
	delete [] Fluxes;
	sprintf(boundaryCdStr,"_s_%g_%g",SolenoidX,SolenoidY);
      }
  }
  
  double SolenoidCF_X=solenoidCF_X, SolenoidCF_Y=solenoidCF_Y;
  if ((SolenoidCF_X!=0.0)||(SolenoidCF_Y!=0.0))
    {
      sprintf(boundaryCdStr,"%s_s_%g_%g",boundaryCdStr,SolenoidCF_X,SolenoidCF_Y);
    }


  char* OutputName;
  char interactionStr[20]="";
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [300];
      if (HardCore)
	sprintf(interactionStr,"_hardcore");      
      sprintf (OutputName, "bosons_lattice_CF_n_%d_x_%d_y_%d%s_q_%d_p_%d%s", NbrBosons, Lx, Ly, interactionStr, NbrFluxQuanta, CFFlux,boundaryCdStr);
    }
  char *TmpC = new char[strlen(OutputName)+20];

  
  int AttachedFlux = CFFlux * NbrBosons;  
  ParticleOnLatticeTranslationOperator *TranslationOperator;
  
  // constructing 1P states:
  
  if (verbose) cout << "* CF states contribute "<<NbrFluxQuanta-AttachedFlux<<" flux"<<endl;
  // space in which CF's live (statistics doesn't matter as we consider single particle physics!)
  //BosonOnLattice *CFSpace = new BosonOnLattice(/*NbrParticles*/ 1, Lx, Ly, NbrFluxQuanta-AttachedFlux, MemorySpace);
  FermionOnLattice *CFSpace = new FermionOnLattice(/*NbrParticles*/ 1, Lx, Ly, NbrFluxQuanta-AttachedFlux, MemorySpace, SolenoidCF_X, SolenoidCF_Y, /*silent*/ verbose);
  TranslationOperator = new ParticleOnLatticeTranslationOperator(CFSpace);
  
  // corresponding Hamiltonians
  AbstractQHEOnLatticeHamiltonian* CFHamiltonian = new ParticleOnLatticeDeltaHamiltonian(CFSpace, /*NbrParticles*/ 1, Lx, Ly, NbrFluxQuanta-AttachedFlux, /* U */ 0.0 , /*ReverseHopping*/ false, /* Delta */ 0.0, /* Random */ 0.0, Architecture.GetArchitecture(), 0, 0x0ul);    
  
  HermitianMatrix HCF(CFHamiltonian->GetHilbertSpaceDimension(), true);
  ComplexMatrix CFEigenVecs(CFHamiltonian->GetHilbertSpaceDimension(), CFHamiltonian->GetHilbertSpaceDimension());
  CFHamiltonian->GetHamiltonian(HCF);
  RealDiagonalMatrix CFEigenVals(CFHamiltonian->GetHilbertSpaceDimension());
  HCF.Diagonalize(CFEigenVals, CFEigenVecs, /* error */ 1e-10 , /* maxIter */ 250);
  ParticleOnLatticeTranslationOperator *CFTranslationOperator= new ParticleOnLatticeTranslationOperator(CFSpace);
  if (!NoMomentumDiagonalize)
    DiagonalizeMomentaInSubspace(CFEigenVals, CFEigenVecs,  CFTranslationOperator, NbrFluxQuanta-AttachedFlux, 0, NbrBosons);
  delete CFTranslationOperator;

  if (Manager.GetString("CF")!=NULL)
    {
      char *InputBase=Manager.GetString("CF");
      for (int i=0; i<NbrBosons; ++i)
	{
	  sprintf(TmpC,"%s.%d.vec",InputBase,i);
	  CFEigenVecs[i].ReadVector(TmpC);
	}
      cout << "Read CF basis from vectors "<<InputBase<<".?.vec"<<endl;
    }
  else
    if (Manager.GetString("all-CF")!=NULL)
    {
      int NbrF;
      char **InputFiles=Manager.GetStrings("all-CF", NbrF);
      if (NbrF!=NbrBosons)
	{
	  cout << "Wrong number of states for CF basis!"<<endl;
	  exit(1);
	}
      for (int i=0; i<NbrBosons; ++i)
	CFEigenVecs[i].ReadVector(InputFiles[i]);
    }

  for (int i=0; i<NbrBosons; ++i)
    if (verbose) cout << "E_CF["<<i<<"]="<<CFEigenVals[i]<<" norm of EVec: "<<CFEigenVecs[i].Norm()<<endl;
  if (verbose) cout << "E_other["<<NbrBosons<<"]="<<CFEigenVals[NbrBosons]<<endl;

  if (Manager.GetBoolean("write-basis"))
    {
      for (int i=0; i<NbrBosons; ++i)
	{
	  sprintf(TmpC,"%s.CF.%d.vec",OutputName,i);
	  CFEigenVecs[i].WriteVector(TmpC);
	}
    }


  if (verbose) cout << "* LLL states for Jastrow-factor contribute "<<AttachedFlux<<" flux"<<endl;  
  
  
  // calculate states required to build Jastrow factor:
  // BosonOnLattice *JastrowSpace = new BosonOnLattice(/*NbrParticles*/ 1, Lx, Ly, AttachedFlux, MemorySpace);
  FermionOnLattice *JastrowSpace = new FermionOnLattice(/*NbrParticles*/ 1, Lx, Ly, AttachedFlux, MemorySpace, SolenoidX-SolenoidCF_X, SolenoidY-SolenoidCF_Y, /*silent*/ verbose);

  AbstractQHEOnLatticeHamiltonian* JastrowHamiltonian = new ParticleOnLatticeDeltaHamiltonian(JastrowSpace, /*NbrParticles*/ 1, Lx, Ly, AttachedFlux, /* U */ 0.0 , /*ReverseHopping*/ false, /* Delta */ 0.0, /* Random */ 0.0, Architecture.GetArchitecture(), 0, 0x0ul);
  delete TranslationOperator;
  
  HermitianMatrix HJastrow(JastrowHamiltonian->GetHilbertSpaceDimension(), true);
  ComplexMatrix JastrowEigenVecs(JastrowHamiltonian->GetHilbertSpaceDimension(),
				 JastrowHamiltonian->GetHilbertSpaceDimension());
  JastrowHamiltonian->GetHamiltonian(HJastrow);
  RealDiagonalMatrix JastrowEigenVals(JastrowHamiltonian->GetHilbertSpaceDimension());
  HJastrow.Diagonalize(JastrowEigenVals, JastrowEigenVecs, /* error */ 1e-10 , /* maxIter */ 250);
  ParticleOnLatticeTranslationOperator *JastrowTranslationOperator= new ParticleOnLatticeTranslationOperator(JastrowSpace);
  if (!NoMomentumDiagonalize)
    DiagonalizeMomentaInSubspace(JastrowEigenVals, JastrowEigenVecs,  JastrowTranslationOperator, AttachedFlux, 0, NbrBosons);
  delete JastrowTranslationOperator;

  for (int i=0; i<NbrBosons; ++i)
    if (verbose) cout << "E_Jastrow["<<i<<"]="<<JastrowEigenVals[i]<<" norm of EVec: "<<JastrowEigenVecs[i].Norm()<<endl;
  if (verbose) cout << "E_other["<<NbrBosons<<"]="<<JastrowEigenVals[NbrBosons]<<endl;

  if (Manager.GetString("J")!=NULL)
    {
      char *InputBase=Manager.GetString("J");
      for (int i=0; i<NbrBosons; ++i)
	{
	  sprintf(TmpC,"%s.%d.vec",InputBase,i);
	  JastrowEigenVecs[i].ReadVector(TmpC);
	}
    }
  
  if (Manager.GetBoolean("write-basis"))
    {
      for (int i=0; i<NbrBosons; ++i)
	{
	  sprintf(TmpC,"%s.Jastrow.%d.vec",OutputName,i);
	  JastrowEigenVecs[i].WriteVector(TmpC);
	}
    }


  if (Manager.GetBoolean("write-product")) 
    {
      // build some product vectors, mainly for testing
      ComplexVector TmpVector(JastrowHamiltonian->GetHilbertSpaceDimension());
      for (int i=0; i<NbrBosons; ++i)
	for (int j=0; j<NbrBosons; ++j)
	  {
	    for (int k=0; k<JastrowHamiltonian->GetHilbertSpaceDimension();++k)
	      {
		TmpVector[k]=JastrowEigenVecs[i][k] * CFEigenVecs[j][k];	    
	      }
	    TmpVector/=TmpVector.Norm();
	    sprintf(TmpC,"%s.prod_J%d_C%d.vec",OutputName,i,j);
	    TmpVector.WriteVector(TmpC);
	  }
    }
    
  
  // cycle through all configurations of the Hilbert-space and calculate the corresponding Slater determinants

  int *QuantumNumbers = new int[NbrBosons];
  double Multiplicity;

#ifdef __LAPACK__
  ComplexLapackDeterminant SlaterCF(NbrBosons);
  ComplexLapackDeterminant SlaterJastrow(NbrBosons);
#else
  ComplexMatrix SlaterCF(NbrBosons, NbrBosons);
  ComplexMatrix SlaterJastrow(NbrBosons, NbrBosons);
#endif
  Complex Value;
  Complex Tmp;

  ComplexVector CFState(Space->GetHilbertSpaceDimension(), true);
  ComplexVector JastrowState(Space->GetHilbertSpaceDimension(), true);

  
  ComplexVector AnalyticJastrowX(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticRelativeX(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticCM1X(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticCM2X(Space->GetHilbertSpaceDimension(), true);
  ComplexVector Analytic1X(Space->GetHilbertSpaceDimension(), true);
  ComplexVector Analytic2X(Space->GetHilbertSpaceDimension(), true);

  ComplexVector AnalyticJastrowY(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticRelativeY(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticCM1Y(Space->GetHilbertSpaceDimension(), true);
  ComplexVector AnalyticCM2Y(Space->GetHilbertSpaceDimension(), true);
  ComplexVector Analytic1Y(Space->GetHilbertSpaceDimension(), true);
  ComplexVector Analytic2Y(Space->GetHilbertSpaceDimension(), true);
  

  // for calculation of analytic Laughlin state:
  // with Landau-gauge along x-axis:
  JacobiThetaFunction ThetaRelX(0.5,0.5,Complex(0.0,((double)Ly)/Lx));
  JacobiThetaFunction ThetaCM1X(1.0/2.0+(NbrFluxQuanta-2.0)/2.0,(NbrFluxQuanta-2.0)/2.0,Complex(0.0,(2.0*(double)Ly)/Lx));
  JacobiThetaFunction ThetaCM2X((NbrFluxQuanta-2.0)/2.0,-(2.0-NbrFluxQuanta)/2.0,Complex(0.0,(2.0*(double)Ly)/Lx));

  // with Landau-gauge along y-axis:
  JacobiThetaFunction ThetaRelY(0.5,0.5,Complex(0.0,((double)Lx)/Ly));
  JacobiThetaFunction ThetaCM1Y(1.0/2.0+(NbrFluxQuanta-2.0)/2.0,(NbrFluxQuanta-2.0)/2.0,Complex(0.0,2.0*((double) Lx)/((double)Ly)));
  JacobiThetaFunction ThetaCM2Y((NbrFluxQuanta-2.0)/2.0,-(2.0-NbrFluxQuanta)/2.0,Complex(0.0,2.0*((double)Lx)/((double) Ly)));


  int *PosX = new int[NbrBosons];
  int *PosY = new int[NbrBosons];
  int Subl;
  double SumX, SumY, SumSqX, SumSqY;
  Complex FRelX, FRelY;  
//   double XSpacing = sqrt(((double)2.0*M_PI*NbrFluxQuanta)/((double)Lx*Ly));
//   double YSpacing = sqrt(((double)2.0*M_PI*NbrFluxQuanta)/((double)Ly*Lx));
  double AreaScaling = 2.0*M_PI*NbrFluxQuanta/((double)(Lx*Ly));

  bool EvaluateAnalytic=Manager.GetBoolean("analytic");


  
  for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
    {
      Space->ListQuantumNumbers(i,QuantumNumbers,Multiplicity);
      //cout << "Q=["<<QuantumNumbers[0];
      //for (int k=1; k<NbrBosons; ++k) cout<<", "<<QuantumNumbers[k];
      //cout << "]"<<endl;
      if (Multiplicity<1.01) // ignore state with multiple occupation
	{	  
	  for (int p = 0; p < NbrBosons; ++p)
	    {
// 	      ComplexVector TmpEigenFctCF = CFEigenVecs[p];
// 	      ComplexVector TmpEigenFctJastrow = JastrowEigenVecs[p];
	      for (int q = 0; q < NbrBosons; ++q)
		{
		  // need to consider proper ordering of matrix elements
		  // in Hilbert-space, largest quantum number q corresponds to position 0!
		  SlaterCF.SetMatrixElement(p,q,CFEigenVecs[p][CFHamiltonian->GetHilbertSpaceDimension()-1-QuantumNumbers[q]]);
		  SlaterJastrow.SetMatrixElement(p,q,JastrowEigenVecs[p][JastrowHamiltonian->GetHilbertSpaceDimension()-1-QuantumNumbers[q]]);
		}	      
	    }
	  CFState[i] = SlaterCF.Determinant();
	  JastrowState[i] = SlaterJastrow.Determinant();
	  
	  // analytic states:
	  if (EvaluateAnalytic)
	    {

	      SumX=0;
	      SumY=0;
	      SumSqX=0;
	      SumSqY=0;
	      for (int q = 0; q < NbrBosons; ++q)
		{
		  Space->DecodeQuantumNumber(QuantumNumbers[q],PosX[q],PosY[q],Subl);
		  //cout << "q="<<QuantumNumbers[q]<<" => [x="<<PosX[q]<<", y="<<PosY[q]<<"]"<<endl;
		  SumX+=PosX[q];
		  SumY+=PosY[q];
		  SumSqX+=PosX[q]*PosX[q]*AreaScaling; // XSpacing*XSpacing;
		  SumSqY+=PosY[q]*PosY[q]*AreaScaling; // YSpacing*YSpacing;
		}
	  
	      FRelX=1.0;
	      for (int bi = 1; bi < NbrBosons; ++bi)
		for (int bj = 0; bj < bi; ++bj)
		  FRelX*=ThetaRelX.GetValue(Complex( ((double)(PosX[bi]-PosX[bj]))/Lx,-((double)(PosY[bi]-PosY[bj]))/Lx));
	  
	      AnalyticJastrowX[i] = FRelX * exp(-0.5*SumSqY);
	      AnalyticRelativeX[i] = FRelX*FRelX * exp(-0.5*SumSqY);
	      AnalyticCM1X[i] = ThetaCM1X.GetValue(2.0*Complex(((double)SumX)/(Lx),-((double)SumY)/(Lx)));
	      AnalyticCM2X[i] = ThetaCM2X.GetValue(2.0*Complex(((double)SumX)/(Lx),-((double)SumY)/(Lx)));
	      //  	  ThetaCM1X.PrintValue(cout,Complex((2.0*(double)SumX)/Lx,(2.0*(double)SumY)/Lx))<<endl;
	      // 	  ThetaCM2X.PrintValue(cout,Complex((2.0*(double)SumX)/Lx,(2.0*(double)SumY)/Lx))<<endl;
	      Analytic1X[i] = AnalyticRelativeX[i] * AnalyticCM1X[i];
	      Analytic2X[i] = AnalyticRelativeX[i] * AnalyticCM2X[i];

	      FRelY=1.0;
	      for (int bi = 1; bi < NbrBosons; ++bi)
		for (int bj = 0; bj < bi; ++bj)
		  FRelY*=ThetaRelY.GetValue(Complex( ((double)(PosY[bi]-PosY[bj]))/Ly,((double)(PosX[bi]-PosX[bj]))/Ly));
	  
	      AnalyticJastrowY[i] = FRelY * exp(-0.5*SumSqX);
	      AnalyticRelativeY[i] = FRelY*FRelY * exp(-0.5*SumSqX);
	      AnalyticCM1Y[i] = ThetaCM1Y.GetValue(2.0*Complex(((double)SumY)/Ly,((double)SumX)/Ly));
	      AnalyticCM2Y[i] = ThetaCM2Y.GetValue(2.0*Complex(((double)SumY)/Ly,((double)SumX)/Ly));
	      //  	  ThetaCM1Y.PrintValue(cout,Complex(((double)SumY)/Ly,((double)SumX)/Ly))<<endl;
	      // 	  ThetaCM2Y.PrintValue(cout,Complex((2.0*(double)SumY)/Ly,(2.0*(double)SumX)/Ly))<<endl;
	      Analytic1Y[i] = AnalyticRelativeY[i] * AnalyticCM1Y[i];
	      Analytic2Y[i] = AnalyticRelativeY[i] * AnalyticCM2Y[i];
	    }
	}
    }
  
  if (Manager.GetBoolean("write-slater"))
    {
      CFState /= CFState.Norm();  
      sprintf(TmpC,"%s.CF.vec",OutputName);
      CFState.WriteVector(TmpC);
    }
  
  if (Manager.GetBoolean("write-jastrow"))
    {
      JastrowState /= JastrowState.Norm();
      sprintf(TmpC,"%s.J.vec",OutputName);
      JastrowState.WriteVector(TmpC);
    }
  
  sprintf(TmpC,"%s.vec",OutputName);
  for (int i=0; i<Space->GetHilbertSpaceDimension(); ++i)
    CFState[i] *= JastrowState[i];
  CFState /= CFState.Norm();
  if (writeState)
    CFState.WriteVector(TmpC);


  
  // writing Laughlin state in Landau gauge along x-direction
  if (Manager.GetBoolean("analytic"))
    {
      AnalyticJastrowX /= AnalyticJastrowX.Norm();
      sprintf(TmpC,"%s.AJX.vec",OutputName);
      AnalyticJastrowX.WriteVector(TmpC);

      AnalyticRelativeX /= AnalyticRelativeX.Norm();
      sprintf(TmpC,"%s.ACFX.vec",OutputName);
      AnalyticRelativeX.WriteVector(TmpC);

      AnalyticCM1X /= AnalyticCM1X.Norm();
      sprintf(TmpC,"%s.CM1X.vec",OutputName);
      AnalyticCM1X.WriteVector(TmpC);

      AnalyticCM2X /= AnalyticCM2X.Norm();
      sprintf(TmpC,"%s.CM2X.vec",OutputName);
      AnalyticCM2X.WriteVector(TmpC);

      Analytic1X /= Analytic1X.Norm();
      sprintf(TmpC,"%s.A1X.vec",OutputName);
      Analytic1X.WriteVector(TmpC);

      Analytic2X /= Analytic2X.Norm();
      sprintf(TmpC,"%s.A2X.vec",OutputName);
      Analytic2X.WriteVector(TmpC);


      // writing Laughlin state in Landau gauge along y-direction
      AnalyticJastrowY /= AnalyticJastrowY.Norm();
      sprintf(TmpC,"%s.AJY.vec",OutputName);
      AnalyticJastrowY.WriteVector(TmpC);

      AnalyticRelativeY /= AnalyticRelativeY.Norm();
      sprintf(TmpC,"%s.ACFY.vec",OutputName);
      AnalyticRelativeY.WriteVector(TmpC);

      AnalyticCM1Y /= AnalyticCM1Y.Norm();
      sprintf(TmpC,"%s.CM1Y.vec",OutputName);
      AnalyticCM1Y.WriteVector(TmpC);

      AnalyticCM2Y /= AnalyticCM2Y.Norm();
      sprintf(TmpC,"%s.CM2Y.vec",OutputName);
      AnalyticCM2Y.WriteVector(TmpC);

      Analytic1Y /= Analytic1Y.Norm();
      sprintf(TmpC,"%s.A1Y.vec",OutputName);
      Analytic1Y.WriteVector(TmpC);

      Analytic2Y /= Analytic2Y.Norm();
      sprintf(TmpC,"%s.A2Y.vec",OutputName);
      Analytic2Y.WriteVector(TmpC);
    }
  
  
  delete [] TmpC;
  delete [] QuantumNumbers;

  delete [] PosX;
  delete [] PosY;

  delete [] OutputName;

  delete CFHamiltonian;
  delete CFSpace;  

  delete JastrowHamiltonian;
  delete JastrowSpace;

  return new ComplexVector(CFState);
}


// values = eigenvalues
// vectors = complex matrix with many vectors to be analyzed and overwritten with basis diagonal in k
// NbrFlux = number of flux quanta to be considered
// start, end = indices of range of eigenvalues of the same eigenvalue
// based on code from FQHELatticeDensityMatrix
ComplexMatrix& DiagonalizeMomentaAtEnergy(ComplexMatrix &vectors, ParticleOnLatticeTranslationOperator *translationOperator, int nbrFlux, int start, int end)
{
  int NbrVectors = end - start;
  ComplexMatrix XTranslationMatrix(NbrVectors, NbrVectors);
  ComplexMatrix YTranslationMatrix(NbrVectors, NbrVectors);
  ComplexVector TmpState(vectors.GetNbrRow());

  ComplexVector TmpState3(vectors.GetNbrRow());

  ComplexVector *Vectors = new ComplexVector[NbrVectors];
  for (int i=start, j=0; i<end; ++i,++j)
    Vectors[j]=ComplexVector(vectors[i], true);

  int Degeneracy=1;
  int n1=1, n2=1;
  int FluxModulo = FindGCD(nbrFlux, Lx*Ly);
  int r=nbrFlux/FluxModulo;
  int t=Lx*Ly/FluxModulo;

  while ((((Ly*n1)%t)!=0) && (n1<Lx)) ++n1;
  while ((((Lx*n2)%t)!=0) && (n2<Ly)) ++n2;

  while ((r*NbrBosons*n1*n2*Degeneracy)%t != 0) ++Degeneracy;
  
//   cout << "N_phi = "<<r<<"/"<<t<<endl;
//   cout << "n1="<<n1<<", n2="<<n2<<", global degeneracy: "<<Degeneracy<<endl;

  ComplexMatrix EVecX(NbrVectors, NbrVectors);
  ComplexMatrix EVecY(NbrVectors, NbrVectors);  
  ComplexDiagonalMatrix EValX(NbrVectors, NbrVectors);
  ComplexDiagonalMatrix EValY(NbrVectors, NbrVectors);
  
  translationOperator->SetTranslationComponents(n1,0);
  GetTranslationMatrix(translationOperator, NbrVectors, Vectors, XTranslationMatrix, TmpState, Architecture);

  translationOperator->SetTranslationComponents(0,n2);
  GetTranslationMatrix(translationOperator, NbrVectors, Vectors, YTranslationMatrix, TmpState, Architecture);

  ComplexMatrix EVecXY(NbrVectors, NbrVectors);


  // form linear superposition of Tx and Ty to diagonalize:
  ComplexMatrix Z((Matrix&)XTranslationMatrix);
  Z*=log(91.0); // scale with some random number > 1
  Z+=YTranslationMatrix;
  Z.Diagonalize(EValX,EVecXY);
  ComplexMatrix QH=EVecXY.GetAdjoint();

  for (int i=0;i<NbrVectors;++i)
    {
      vectors[start+i].ClearVector();
      for (int j=0; j<NbrVectors;++j)
	vectors[start+i].AddLinearCombination(Conj(EVecXY[i][j]),Vectors[j]);      
    }
    
  delete [] Vectors;
  return vectors;
}


// values = eigenvalues
// vectors = complex matrix with many vectors to be analyzed and overwritten with basis diagonal in k
// nbrFlux = number of flux quanta to be considered
// start, end = indices of range of eigenvalues to be considered
ComplexMatrix& DiagonalizeMomentaInSubspace(RealDiagonalMatrix &values, ComplexMatrix &vectors, ParticleOnLatticeTranslationOperator *translationOperator, int nbrFlux, int start, int end)

{
  // make sure we get entire multiplets:
  while ((start > 0)&&(fabs(values[start]-values[start-1])<1e-12)) --start;
  while ((end < values.GetNbrRow()-1)&&(fabs(values[end]-values[end-1])<1e-12))
    {
      ++end;
    }
  
  int startSeg=start, endSeg=start+1;

  while (endSeg <= end)
    {
      while ((endSeg<end)&&(fabs(values[endSeg]-values[endSeg-1])<1e-12))
	++endSeg;      
      if (endSeg-startSeg > 1)
	{
	  // cout << "Multiplet ["<<startSeg<<", "<<endSeg-1<<"]"<<endl;
	  DiagonalizeMomentaAtEnergy(vectors, translationOperator, nbrFlux, startSeg, endSeg);
	}
      startSeg=endSeg;
      endSeg=startSeg+1;
    }
  return vectors;
}


void GetTranslationMatrix(ParticleOnLatticeTranslationOperator *Operator, int NbrVectors,
			  ComplexVector *Vectors, ComplexMatrix &MatrixRepresentation,
			  ComplexVector &TmpState, ArchitectureManager &Architecture)
{
  Complex Tmp;
  for (int i=0; i<NbrVectors; ++i)
    {
      VectorOperatorMultiplyOperation Operation (Operator, &(Vectors[i]), &TmpState);      
      Operation.ApplyOperation(Architecture.GetArchitecture());           
      for (int j=0; j<NbrVectors; ++j)
	{
	  Tmp = Vectors[j] * TmpState;
	  MatrixRepresentation.SetMatrixElement(i,j,Tmp);
	}
    }
}
