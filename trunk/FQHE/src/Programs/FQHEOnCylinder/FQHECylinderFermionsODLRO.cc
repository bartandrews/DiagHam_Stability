#include "Vector/RealVector.h"
#include "Vector/ComplexVector.h"

#include "HilbertSpace/FermionOnTorus.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereFull.h"
#include "HilbertSpace/FermionOnSphereUnlimited.h"
#include "HilbertSpace/FermionOnSphereLong.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "GeneralTools/ConfigurationParser.h"

#include "Hamiltonian/ParticleOnCylinderOrbitalProjection.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
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
  OptionManager Manager ("FQHECylinderFermionsODLRO" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('e', "eigenstate", "name of the file containing the eigenstate");
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "ky-max", "twice the maximum momentum for a single particle (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('y', "total-y", "twice the total momentum projection for the system (override autodetection from input file name if greater or equal to zero)", 0);
  (*SystemGroup) += new SingleDoubleOption  ('r', "ratio", "aspect ratio of the cylinder", 1.0);
  (*SystemGroup) += new SingleDoubleOption ('a', "anisotropy", "shape (anisotropy) parameter", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "x0-start", "starting x-coordinate of the center of the projected orbital", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "y0-start", "starting y-coordinate of the center of the projected orbital", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "x0-end", "ending x-coordinate of the center of the projected orbital", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "y0-end", "ending y-coordinate of the center of the projected orbital", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrpoints-x", "calculate ODLRO at nbrpoints between x0start and x0end", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "nbrpoints-y", "calculate ODLRO at nbrpoints between y0start and y0end", 1);
  (*SystemGroup) += new BooleanOption  ('\n', "project-only", "project state with operator [1-n_2][1-n_1]n_0 and exit", false);
  (*OutputGroup) += new SingleStringOption ('o', "output-file", "use this file name for output");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsCorrelation -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = Manager.GetInteger("nbr-particles");
  int KyMax = Manager.GetInteger("ky-max");
  int TotalKy = Manager.GetInteger("total-y");
  double XRatio = Manager.GetDouble("ratio");
  double X0Start = Manager.GetDouble("x0-start");
  double Y0Start = Manager.GetDouble("y0-start");
  double X0End = Manager.GetDouble("x0-end");
  double Y0End = Manager.GetDouble("y0-end");
  double Anisotropy = ((SingleDoubleOption*) Manager["anisotropy"])->GetDouble();
  cout<<"Orbital shape: anisotropy = "<<Anisotropy<<endl;


  int NbrPointsX = Manager.GetInteger("nbrpoints-x");
  int NbrPointsY = Manager.GetInteger("nbrpoints-y");

  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  if (Manager.GetString("eigenstate") == 0)
    {
      cout << "FQHECylinderFermionsCorrelation requires a state" << endl;
      return -1;
    }

  if (IsFile(Manager.GetString("eigenstate")) == false)
    {
      cout << "can't find vector file " << Manager.GetString("eigenstate") << endl;
      return -1;      
    }

  ParticleOnSphere* Space = 0;
#ifdef __64_BITS__
	  if (TotalKy <= 62)
#else
	  if (TotalKy <= 30)
#endif
  	    Space = new FermionOnSphere(NbrParticles, TotalKy, KyMax);

	  else
#ifdef __128_BIT_LONGLONG__
	    if (TotalKy <= 126)
#else
	      if (TotalKy <= 62)
#endif
 	        Space = new FermionOnSphereLong(NbrParticles, TotalKy, KyMax);
	      else
		Space = new FermionOnSphereUnlimited(NbrParticles, TotalKy, KyMax);


  cout << " Hilbert space dimension = " << Space->GetHilbertSpaceDimension() << endl;

  ParticleOnSphere* FullSpace;
  FullSpace = new FermionOnSphereFull(NbrParticles, KyMax);

  cout << " full Hilbert space dimension = " << FullSpace->GetHilbertSpaceDimension() << endl;


  ofstream File;
  File.precision(14);
  if (Manager.GetString("output-file") != 0)
     File.open(Manager.GetString("output-file"), ios::binary | ios::out);
  else
   {
    cout << "Enter output file! " << endl;
    exit(1);
   }
  
  ComplexVector State;

  if (State.ReadVector (Manager.GetString("eigenstate")) == false)
   {
     cout << "can't open vector file " << Manager.GetString("eigenstate") << endl;
     return -1;      
   }

  char* OutputNameLz = new char [512];
  sprintf (OutputNameLz, "fermions_allky_n_%d_2s_%d.0.vec", NbrParticles, KyMax);

  ComplexVector FullState(FullSpace->GetHilbertSpaceDimension(), true);

  FermionOnSphereFull* FermionFullSpace = (FermionOnSphereFull*) FullSpace;

  FermionFullSpace->ConvertToAllLz(State, Space, FullState);

  //FullState.WriteVector(OutputNameLz);

  double Length = sqrt(2.0 * M_PI * XRatio * (KyMax + 1));
  double Height = sqrt(2.0 * M_PI * (KyMax + 1) / XRatio);
  cout<<"Length = "<<Length<<" Circumference= "<<Height<<endl;
  if ((fabs(X0End - X0Start) >= Length) || (fabs(Y0End - Y0Start) >= Height))
   {
     cout<<"Point is at the boundary or beyond! " <<endl;
     exit(-1);
   }
  double StepX = (X0End-X0Start)/NbrPointsX;
  double StepY = (Y0End-Y0Start)/NbrPointsY;

  //*************************************************************************
  //********* Act with n_0(0)(1-n_1(0))(1-n_2(0)) ***************************
  //*************************************************************************

  ComplexVector StateAt0(FullSpace->GetHilbertSpaceDimension(), true);
  ComplexVector StateAt01(FullSpace->GetHilbertSpaceDimension(), true);
  ComplexVector StateAt012(FullSpace->GetHilbertSpaceDimension(), true);

  AbstractHamiltonian* Hamiltonian0 = new ParticleOnCylinderOrbitalProjection (FullSpace, NbrParticles, KyMax, XRatio, 0, Anisotropy, X0Start, Y0Start, Architecture.GetArchitecture(), Memory);
  VectorHamiltonianMultiplyOperation Operation0 (Hamiltonian0, &FullState, &StateAt0);
  Operation0.ApplyOperation(Architecture.GetArchitecture());
  cout<<"Completed projecting 0 orbital at (0,0); norm= "<<StateAt0.Norm()<<endl;
  delete Hamiltonian0;

  AbstractHamiltonian* Hamiltonian01 = new ParticleOnCylinderOrbitalProjection (FullSpace, NbrParticles, KyMax, XRatio, 1, Anisotropy, X0Start, Y0Start, Architecture.GetArchitecture(), Memory);
  VectorHamiltonianMultiplyOperation Operation01 (Hamiltonian01, &StateAt0, &StateAt01);
  Operation01.ApplyOperation(Architecture.GetArchitecture());
  cout<<"Completed projecting 1 orbital at (0,0); norm= "<< StateAt01.Norm() << endl;
  delete Hamiltonian01;

  StateAt01 = StateAt0 - StateAt01;

  AbstractHamiltonian* Hamiltonian012 = new ParticleOnCylinderOrbitalProjection (FullSpace, NbrParticles, KyMax, XRatio, 2, Anisotropy, X0Start, Y0Start, Architecture.GetArchitecture(), Memory);
  VectorHamiltonianMultiplyOperation Operation012 (Hamiltonian012, &StateAt01, &StateAt012);
  Operation012.ApplyOperation(Architecture.GetArchitecture());
  cout<<"Completed projecting 2 orbital at (0,0); norm= "<< StateAt012.Norm() << endl;
  delete Hamiltonian012;

  StateAt012 = StateAt01 - StateAt012;
  StateAt012 /= StateAt012.Norm();
  
  //*************************************************************************
  if (((BooleanOption*) Manager["project-only"])->GetBoolean() == true)
    {
      StateAt012.WriteVector(OutputNameLz); 
      cout<<"Overlap with initial state: "<<(FullState * StateAt012) << endl;
      StateAt0 /= StateAt0.Norm();
      cout<<"Overlap between n0 Psi and (1-n2)(1-n1)n0 Psi: "<<(StateAt0 * StateAt012) << endl;
      return 0;
    }
  //*************************************************************************

  int counter=0;
  for (int i = 0; i < NbrPointsX; i++)
    for (int j = 0; j < NbrPointsY; j++)
    {
     double X = X0Start + i * StepX;
     double Y = Y0Start + j * StepY;
     cout<<"---------------Step "<<counter<<" out of "<<(NbrPointsX * NbrPointsY)<<" X= " << X<<" Y= "<<Y<<"---------"<<endl;    
  
    //*************************************************************************
    //********* Act with n_0(x0,y0)(1-n_1(x0,y0))(1-n_2(x0,y0)) ***************************
    //*************************************************************************

    ComplexVector RStateAt0(FullSpace->GetHilbertSpaceDimension(), true);
    ComplexVector RStateAt01(FullSpace->GetHilbertSpaceDimension(), true);
    ComplexVector RStateAt012(FullSpace->GetHilbertSpaceDimension(), true);

    AbstractHamiltonian* Hamiltonian0R = new ParticleOnCylinderOrbitalProjection (FullSpace, NbrParticles, KyMax, XRatio, 0, Anisotropy, X, Y, Architecture.GetArchitecture(), Memory);
    VectorHamiltonianMultiplyOperation Operation0R (Hamiltonian0R, &FullState, &RStateAt0);
    Operation0R.ApplyOperation(Architecture.GetArchitecture());
    cout<<"Completed projecting 0 orbital at (x0,y0); norm= "<< RStateAt0.Norm() <<endl;
    delete Hamiltonian0R;

    AbstractHamiltonian* Hamiltonian01R = new ParticleOnCylinderOrbitalProjection (FullSpace, NbrParticles, KyMax, XRatio, 1, Anisotropy, X, Y, Architecture.GetArchitecture(), Memory);
    VectorHamiltonianMultiplyOperation Operation01R (Hamiltonian01R, &RStateAt0, &RStateAt01);
    Operation01R.ApplyOperation(Architecture.GetArchitecture());
    cout<<"Completed projecting 1 orbital at (x0,y0); norm= "<< RStateAt01.Norm() <<endl;
    delete Hamiltonian01R;

    RStateAt01 = RStateAt0 - RStateAt01;

    AbstractHamiltonian* Hamiltonian012R = new ParticleOnCylinderOrbitalProjection (FullSpace, NbrParticles, KyMax, XRatio, 2, Anisotropy, X, Y, Architecture.GetArchitecture(), Memory);
    VectorHamiltonianMultiplyOperation Operation012R (Hamiltonian012R, &RStateAt01, &RStateAt012);
    Operation012R.ApplyOperation(Architecture.GetArchitecture());
    cout<<"Completed projecting 2 orbital at (x0,y0); norm= "<< RStateAt012.Norm() <<endl;
    delete Hamiltonian012R;

    RStateAt012 = RStateAt01 - RStateAt012;
    RStateAt012 /= RStateAt012.Norm();

    //*************************************************************************

    Complex Overlap = StateAt012 * RStateAt012;
    cout<<"<Psi_0|Psi_R>="<<endl;
    cout<<Overlap.Re<<" "<<Overlap.Im<<endl;
    File << X << " " << Y << " " << Overlap.Re << " "<<Overlap.Im << endl;  
    counter++;
  }


  File.close();
 
  return 0;
}
