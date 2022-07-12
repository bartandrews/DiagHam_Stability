#include "Options/Options.h"

#include "HilbertSpace/FermionOnLattice.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeLong.h"

#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"

#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"

#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
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

void  FindMagneticCell(const int nbrFluxQuanta,const int lx, const int ly, int & nxZero,int & nyZero);

int main(int argc, char** argv)
{
  cout.precision(14);
  
  OptionManager Manager ("FCIHofstadterModelCompositeFermions" , "0.01");  
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OptimizationGroup = new OptionGroup ("optimization options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");  
  ArchitectureManager Architecture;
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
  (*SystemGroup) += new BooleanOption('\n',"no-translation","Use Hilbert-space of hard-core bosons");
  
  (*SystemGroup) += new MultipleDoubleOption  ('\n', "solenoid-flux", "twist in periodic boundary conditions for total wavefunction phi_x[,phi_y])",',');
  (*SystemGroup) += new MultipleDoubleOption  ('s', "solenoid-CF", "twist in periodic boundary conditions for CF part phi_x[,phi_y])",',');
  
  
  (*PrecalculationGroup) += new SingleIntegerOption  ('\n', "fast-search", "amount of memory that can be allocated for fast state search (in Mbytes)", 9);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
  
  (*MiscGroup) += new BooleanOption  ('a', "analytic", "also generate the analytic wavefunctions");
  (*MiscGroup) += new BooleanOption  ('\n', "write-basis", "write the single particle basis states that were used");
  (*MiscGroup) += new BooleanOption  ('\n', "write-product", "write the product states of pairs of basis states");
  (*MiscGroup) += new BooleanOption  ('\n', "write-slater", "write the slater determinant part of the wavefunction");
  (*MiscGroup) += new BooleanOption  ('\n', "write-jastrow", "write the jastrow factor part of the wavefunction");
  (*MiscGroup) += new SingleStringOption  ('o', "output-file", "redirect output to this file",NULL);
  
  (*MiscGroup) += new BooleanOption  ('\n', "debug", "display debug information");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  Manager.StandardProceedings(argv, argc, cout);
  
  int NbrBosons = Manager.GetInteger("nbr-particles");
  int Lx = Manager.GetInteger("lx");
  int Ly = Manager.GetInteger("ly");
  int NbrSites = Lx*Ly;
  int NbrFluxQuanta = Manager.GetInteger("flux");  
  bool HardCore = Manager.GetBoolean("hard-core");
  int CFFlux = Manager.GetInteger("flux-per-CF");
  char Axis ='y';
  bool EmbeddingFlag = true;
  bool NoTranslationFlag = Manager.GetBoolean("no-translation");
  
  unsigned long MemorySpace = ((unsigned long) Manager.GetInteger("fast-search")) << 20;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  
  double SolenoidCF_X=0.0, SolenoidCF_Y=0.0;
  {
    int tmpI;
    double *Fluxes=Manager.GetDoubles("solenoid-CF", tmpI);
    if (tmpI>0) SolenoidCF_X=Fluxes[0];
    if (tmpI>1) SolenoidCF_Y=Fluxes[1];
    if (tmpI>0) delete [] Fluxes;	
  }
  
  int NxZero,NyZero;
  FindMagneticCell(NbrFluxQuanta,Lx,Ly, NxZero,NyZero);
  int MaxMomentumX =  Lx/ NxZero;
  int MaxMomentumY =  Ly/ NyZero;

  cout <<  "NxZero = "<<NxZero<<" NyZero = " <<NyZero <<" Flux Per unit cell " << (((double) NbrFluxQuanta )* ((double) NxZero*(double) NyZero)/((double )NbrSites )) <<endl;
  
  int XMomentum=0;
  int YMomentum=0;
  
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
  
  char* OutputName;
  if ( (OutputName = Manager.GetString("output-file")) == NULL)
    {
      OutputName = new char [300];
      
      sprintf (OutputName, "bosons_lattice_CF_n_%d_X_%d_Y_%d_x_%d_y_%d_q_%d_p_%d_kx_%d_ky_%d%s.vec", NbrBosons, NxZero, NyZero, MaxMomentumX ,MaxMomentumY , NbrFluxQuanta, CFFlux, XMomentum , YMomentum ,boundaryCdStr);
    }
  
  int AttachedFlux = CFFlux * NbrBosons;
  
  int TotalFluxCF =  (NbrFluxQuanta-AttachedFlux);
  int NxZeroCF,NyZeroCF;
  FindMagneticCell(TotalFluxCF,Lx,Ly, NxZeroCF,NyZeroCF);
  int MaxMomentumXCF =  Lx/NxZeroCF;
  int MaxMomentumYCF =  Ly/NyZeroCF;
  
  int TotalFluxJastrow =  AttachedFlux;
  int NxZeroJastrow,NyZeroJastrow;
  FindMagneticCell(TotalFluxJastrow,Lx,Ly, NxZeroJastrow,NyZeroJastrow);
  int MaxMomentumXJastrow =  Lx/NxZeroJastrow;
  int MaxMomentumYJastrow =  Ly/NyZeroJastrow;
  
  int FluxPerCellCF =  TotalFluxCF/( MaxMomentumXCF* MaxMomentumYCF);
  int FluxPerCellJastrow =  TotalFluxJastrow/( MaxMomentumXJastrow* MaxMomentumYJastrow);  
  
  // constructing 1P states:
  bool verbose = Manager.GetBoolean("debug");
  if (verbose) cout << "* CF states contribute "<<NbrFluxQuanta-AttachedFlux<<" flux"<<endl;
  
  // corresponding Hamiltonians


  TightBindingModelHofstadterSquare TightBindingModelCF (MaxMomentumXCF, MaxMomentumYCF, NxZeroCF, NyZeroCF, FluxPerCellCF, Axis, SolenoidCF_X, SolenoidCF_Y, Architecture.GetArchitecture());
  
  char* OutputNameCFEnergy;
  OutputNameCFEnergy = new char [300];
  sprintf (OutputNameCFEnergy, "HofstadterModel_bandstructure_X_%d_Y_%d_x_%d_y_%d_q_%d.dat", NxZeroCF,NyZeroCF, MaxMomentumXCF,MaxMomentumYCF, FluxPerCellCF);
  TightBindingModelCF.WriteAsciiSpectrum(OutputNameCFEnergy);
  
  if(NoTranslationFlag == true)
    {
      NxZeroCF = NxZeroCF* MaxMomentumXCF;
      NyZeroCF = NyZeroCF* MaxMomentumYCF;
      FluxPerCellCF = TotalFluxCF;
      MaxMomentumXCF = 1;
      MaxMomentumYCF = 1;
      EmbeddingFlag = false;  
    }
  
  TightBindingModelHofstadterSquare TightBindingModelCF2 (MaxMomentumXCF, MaxMomentumYCF, NxZeroCF, NyZeroCF, FluxPerCellCF, Axis, SolenoidCF_X, SolenoidCF_Y, Architecture.GetArchitecture(),true, EmbeddingFlag);
  
  cout <<"Building CF factor with MaxMomentumXCF =" << MaxMomentumXCF << "  MaxMomentumYCF =" << MaxMomentumYCF << " NxZeroCF  = "<<  NxZeroCF<< " NyZeroCF = " <<NyZeroCF <<" FluxPerCellCF = "<< FluxPerCellCF <<endl;
  
  sprintf (OutputNameCFEnergy, "HofstadterModel_bandstructureWithEmbedding_X_%d_Y_%d_x_%d_y_%d_q_%d.dat", NxZeroCF, NyZeroCF, MaxMomentumXCF,MaxMomentumYCF, FluxPerCellCF);
  TightBindingModelCF2.WriteAsciiSpectrum(OutputNameCFEnergy);
  
 

  ComplexMatrix CFEigenVecs =  TightBindingModelCF2.GetRealSpaceTightBindingEigenstates();

  /*
  HermitianMatrix TmpHamCF =   TightBindingModelCF2.GetRealSpaceTightBindingHamiltonian();

  ComplexMatrix CFEigenVecs(Lx* Ly, Lx* Ly,true);
  CFEigenVecs.SetToIdentity();
  RealDiagonalMatrix TmpDiagCF;
  TmpHamCF.LapackDiagonalize(TmpDiagCF, CFEigenVecs);*/
  
  if (verbose) cout << "* LLL states for Jastrow-factor contribute "<<AttachedFlux<<" flux"<<endl;  
  
  
if(NoTranslationFlag == true)
    {
      NxZeroJastrow *= MaxMomentumXJastrow;
      NyZeroJastrow *= MaxMomentumYJastrow;
      FluxPerCellJastrow = TotalFluxJastrow;
      MaxMomentumXJastrow= 1;
      MaxMomentumYJastrow= 1;
      EmbeddingFlag =  false;
    }
 
 
 TightBindingModelHofstadterSquare  JastrowTightBindingModel (MaxMomentumXJastrow, MaxMomentumYJastrow, NxZeroJastrow, NyZeroJastrow, FluxPerCellJastrow, Axis, SolenoidCF_X, SolenoidCF_Y, Architecture.GetArchitecture(),true,EmbeddingFlag);
 
 ComplexMatrix JastrowEigenVecs =  JastrowTightBindingModel.GetRealSpaceTightBindingEigenstates();
 
 cout <<"Building Jastrow factor with MaxMomentumXJastrow =" << MaxMomentumXJastrow << "  MaxMomentumYJastrow =" << MaxMomentumYJastrow << " NxZeroJastrow  = "<<  NxZeroJastrow<< " NyZeroJastrow = " <<NyZeroJastrow <<" FluxPerCellJastrow = "<< FluxPerCellJastrow <<endl;
  
 /*HermitianMatrix TmpHamJastrow = JastrowTightBindingModel.GetRealSpaceTightBindingHamiltonian();
 ComplexMatrix JastrowEigenVecs( Lx* Ly, Lx* Ly,true);
 JastrowEigenVecs.SetToIdentity();
 RealDiagonalMatrix TmpDiagJastrow;
 TmpHamJastrow.LapackDiagonalize(TmpDiagJastrow, JastrowEigenVecs);*/
 
  
  
  if(NoTranslationFlag == false)
    {
      BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation Space (NbrBosons, Lx, Ly , XMomentum, MaxMomentumX, YMomentum, MaxMomentumY);
      cout << "Using Hilbert space: BosonOnLatticeGutzwillerProjectionRealSpaceOneOrbitalPerSiteAnd2DTranslation"<<endl;
      char* OutputName;
      if ( (OutputName = Manager.GetString("output-file")) == NULL)
	{
	  OutputName = new char [300];
	  sprintf (OutputName, "bosons_lattice_CF_n_%d_X_%d_Y_%d_x_%d_y_%d_q_%d_p_%d_kx_%d_ky_%d%s.vec", NbrBosons, NxZero, NyZero, MaxMomentumX ,MaxMomentumY , NbrFluxQuanta, CFFlux, XMomentum , YMomentum ,boundaryCdStr);
	}
      
      ComplexVector TrialState(Space.GetHilbertSpaceDimension(),true);
      double PhaseTranslationX = -2.0*M_PI*NbrFluxQuanta/(Lx*Ly) * NxZero;
      HermitianMatrix TmpHam = JastrowTightBindingModel.GetRealSpaceTightBindingHamiltonian();
      ComplexMatrix TmpMatrix( Lx* Ly, Lx* Ly,true);
      TmpMatrix.SetToIdentity();
      RealDiagonalMatrix TmpDiag;
      TmpHam.LapackDiagonalize(TmpDiag, TmpMatrix);
      
      Space.GetCompositeFermionWavefunction(TrialState, JastrowEigenVecs, CFEigenVecs,PhaseTranslationX );
      cout <<"State Norm " << TrialState.Norm()<<endl;
      TrialState/= TrialState.Norm();
      TrialState.WriteVector(OutputName);
    }
  else
    {
      BosonOnLatticeGutzwillerProjectionRealSpace Space1 (NbrBosons,Lx*Ly);
      cout << "Using Hilbert space: BosonOnLatticeGutzwillerProjectionRealSpace"<<endl;
      ComplexVector TrialState1(Space1.GetHilbertSpaceDimension(),true);
      Space1.GetCompositeFermionWavefunction(TrialState1, JastrowEigenVecs, CFEigenVecs);
      
      char* OutputName1;
      if ( (OutputName1 = Manager.GetString("output-file")) == NULL)
	{
	  OutputName1 = new char [300];
	  sprintf (OutputName1, "bosons_lattice_CF_n_%d_x_%d_y_%d_q_%d_p_%d%s.vec", NbrBosons, Lx , Ly, NbrFluxQuanta, CFFlux, boundaryCdStr);
	}
      
      cout <<"State Norm " << TrialState1.Norm()<<endl;
      TrialState1/= TrialState1.Norm();
      TrialState1.WriteVector(OutputName1);
    }
  return 0;
}




void  FindMagneticCell(const int nbrFluxQuanta,const int lx, const int ly, int & nxZero,int & nyZero)
{
  int FluxModulo = FindGCD(nbrFluxQuanta,  lx*ly);
  int p = nbrFluxQuanta /  FluxModulo;
  int q =  lx*ly /  FluxModulo;
  nxZero = lx;
  nyZero = ly;
  double OptimizeAspectRatio=0;
  cout <<"p= "<< p << " q = " <<q<<endl;

  for(int i =1 ; i <=sqrt(q) ; i++)
  {
 	if (q%i == 0 )
	{
           int TmpNxZero;
	   int TmpNyZero;
           bool ChangeFlag = false;
           if(( lx%i==0 ) && ( ly% (q/i) == 0))
	   {	
	     TmpNxZero = i;	
	     TmpNyZero = q/i;
	     ChangeFlag=true;
            }
	    else
	   {
            	 if(( ly%i == 0 ) && ( lx% (q/i) == 0))
		 {
	            TmpNxZero = q/i;	
         	    TmpNyZero = i;
		    ChangeFlag=true;
		 }
	   }
	   if (ChangeFlag)
	     {
	       double TmpAspectRatio =  ((double) TmpNyZero)/ ((double) TmpNxZero);
	       if(  TmpAspectRatio > 1  ) 
		 {
		   TmpAspectRatio = 1.0/TmpAspectRatio;
		 }
	       if ( TmpAspectRatio > OptimizeAspectRatio)
		 {
		   OptimizeAspectRatio = TmpAspectRatio;
		   nxZero = TmpNxZero;
		   nyZero = TmpNyZero;
		 }
	     }
	}
  }
}

