#include "Options/Options.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd1DTranslation.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"


#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DMagneticTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"


#include "Tools/FTITightBinding/TightBindingModelHofstadterSquare.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterTriangularQuarter.h"
#include "Tools/FTITightBinding/TightBindingModelHofstadterFiniteCylinder.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"
#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"

#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Operator/ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>

using std::cout;
using std::cin;
using std::endl;
using std::ios;
using std::ofstream;

void  ComputeDensity(int NbrFermions, ComplexMatrix & eigenVectors, int X1, int Y1,  int X2, int Y2, TightBindingModelHofstadterFiniteCylinder * tightBindingModel);
void  ComputeCurrent(FermionOnLatticeRealSpaceAnd1DTranslation * space, double * tunnelElementX, double * tunnelElementY, double fluxDensity, ComplexVector & groundState);


int main(int argc, char** argv)
{
  OptionManager Manager ("FCIHofstadterModelFiniteCylinder" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(true);  
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Lanczos.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 5);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbrsitex", "number of unit cells along the x direction", 5);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbrsitey", "number of sites in the y direction", 7);
  (*SystemGroup) += new SingleIntegerOption  ('q', "total-flux", "number of flux quanta per unit cell", 1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "compute a single sector of momentum", -1);
  
  (*SystemGroup) += new SingleDoubleOption  ('\n', "flux-inserted", "flux insert in the x direction", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive onsite(boson) or NN (fermion) potential strength", 1.0);

  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-x", "boundary condition twisting angle along x (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "gamma-y", "boundary condition twisting angle along y (in 2 Pi unit)", 0.0);
  (*SystemGroup) += new BooleanOption  ('\n', "torus", "compute mean value of the Hamiltonian against each eigenstate",false);

  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-spectrum", "only compute the one body spectrum");
  (*SystemGroup) += new BooleanOption  ('\n', "singleparticle-chernnumber", "compute the chern number (only in singleparticle-spectrum mode)");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebody", "export the one-body information (band structure and eigenstates) in a binary file");
  (*SystemGroup) += new BooleanOption  ('\n', "export-onebodytext", "export the one-body information (band structure and eigenstates) in an ASCII text file");
  (*SystemGroup) += new SingleStringOption  ('\n', "export-onebodyname", "optional file name for the one-body information output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenvalue-file", "filename for eigenvalues output");
  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "filename for eigenstates output; to be appended by _kx_#_ky_#.#.vec");
  (*SystemGroup) += new BooleanOption  ('\n', "get-hvalue", "compute mean value of the Hamiltonian against each eigenstate");

  (*SystemGroup) += new BooleanOption  ('\n', "no-translation", "use the real space representation when considering the system with all bandswithout the translations");
  (*SystemGroup) += new BooleanOption  ('\n', "synthetic-dimension", "use synthetic dimension coupling");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-density", "");

  (*SystemGroup) += new  SingleStringOption ('\n', "use-hilbert", "name of the file that contains the vector files used to describe the reduced Hilbert space (replace the n-body basis)");
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);

  (*SystemGroup) += new BooleanOption  ('\n', "atomic-limit", "use atomic limit tight-binding model to test interaction terms");

  (*SystemGroup) += new BooleanOption  ('\n', "compute-current", " ",false);
  (*SystemGroup) += new SingleStringOption  ('\n', "groundstate-file", "filename for the groundstate whose current will be computed");

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");
  (*ToolsGroup) += new SingleDoubleOption  ('\n',"testhermitian-error", "precision of the hermeticity test",0);

  (*MiscGroup) += new SingleStringOption('\n', "energy-expectation", "name of the file containing the state vector, whose energy expectation value shall be calculated");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);

  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSiteX = Manager.GetInteger("nbrsitex"); 
  int NbrSiteY = Manager.GetInteger("nbrsitey");
  int Flux = Manager.GetInteger("total-flux");
  char Axis ='y';

  double * TunnelElementX = new double[NbrSiteY];
  double * TunnelElementY = new double[NbrSiteY];
  
  cout <<"Give Tunneling element"<<endl;
/*  for (int i = 0; i <NbrSiteY; i++)
    cin >> TunnelElementX[i];*/

  TunnelElementX[0] = 1.0374; 
  TunnelElementX[1] = 1.0212;
  TunnelElementX[2] = 1.0094;
  TunnelElementX[3] = 1.0022;
  TunnelElementX[4] = 1.0;
  TunnelElementX[5] = 1.0022;
  TunnelElementX[6] = 1.0094;
  TunnelElementX[7] = 1.0212;
  TunnelElementX[8] = 1.0374;


/*  for (int i = 0; i <NbrSiteY; i++)
    cin >> TunnelElementY[i];*/

  for (int i = 0; i <NbrSiteY; i++)
    TunnelElementY[i] = 1.0;

  if( Manager.GetBoolean("compute-current") == true)
{
   TightBindingModelHofstadterFiniteCylinder  TightBindingModel (NbrSiteX, NbrSiteY, Flux, TunnelElementX, TunnelElementY, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), Manager.GetDouble("flux-inserted"),true,Manager.GetBoolean("torus"));

   FermionOnLatticeRealSpaceAnd1DTranslation  Space (NbrParticles,NbrSiteX*NbrSiteY, Manager.GetInteger("only-kx"),NbrSiteX);
   char* StateFileName = Manager.GetString("groundstate-file");
   if (IsFile(StateFileName) == false)
    {
       cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
       return -1;           
    }
  ComplexVector State;
  if (State.ReadVector(StateFileName) == false)
   {
      cout << "error while reading " << StateFileName << endl;
      return -1;
   }
  if (State.GetVectorDimension() != Space.GetHilbertSpaceDimension())
    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
   ComputeCurrent(&Space, TunnelElementX, TunnelElementY,TightBindingModel.GetFluxDensity() ,State);
   return 0;

}



  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;

  char* StatisticPrefix = new char [16];
  sprintf (StatisticPrefix, "fermions");
  
/*  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }
  */
  
  char* FilePrefix = new char [512];
  int lenFilePrefix=0;


if (Manager.GetBoolean("synthetic-dimension") == false)
      lenFilePrefix += sprintf (FilePrefix, "%s_realspace_hofstadter_q_%d", StatisticPrefix,Flux);
else
      lenFilePrefix += sprintf (FilePrefix, "%s_realspace_synth_hofstadter_q_%d", StatisticPrefix, Flux);

  // common naming options:
  lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_n_%d_x_%d_y_%d", NbrParticles, NbrSiteX, NbrSiteY);

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx");
  char* EigenvalueOutputFile = new char [512];
  if (Manager.GetString("eigenvalue-file")!=0)
    {
      strcpy(EigenvalueOutputFile, Manager.GetString("eigenvalue-file"));
      delete [] FilePrefix;
       FilePrefix = RemoveExtensionFromFileName(EigenvalueOutputFile,".dat");
      if (FilePrefix == 0)
	strcpy(FilePrefix, EigenvalueOutputFile);
    }
  else
    {
      lenFilePrefix += sprintf (FilePrefix+lenFilePrefix, "_u_%g_gx_%g",Manager.GetDouble("u-potential"),Manager.GetDouble("flux-inserted"));
      sprintf (EigenvalueOutputFile,"%s.dat",FilePrefix);
    }
  


  if (Manager.GetBoolean("singleparticle-spectrum") == true)
    {
      bool ExportOneBody = false;
      if ((Manager.GetBoolean("export-onebody") == true) || (Manager.GetBoolean("export-onebodytext") == true) || (Manager.GetBoolean("singleparticle-chernnumber") == true))
	ExportOneBody = true;

      Abstract2DTightBindingModel * TightBindingModel;
      TightBindingModel = new TightBindingModelHofstadterFiniteCylinder(NbrSiteX, NbrSiteY, Flux, TunnelElementX, TunnelElementY, Axis,Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(), Manager.GetDouble("flux-inserted"), ExportOneBody,Manager.GetBoolean("torus"));
      
      
    TightBindingModel->WriteAsciiSpectrum(EigenvalueOutputFile);

    HermitianMatrix TmpHam (TightBindingModel->GetRealSpaceTightBindingHamiltonian());
    RealDiagonalMatrix TmpHam2(TmpHam.GetNbrRow());
    ComplexMatrix Q(NbrSiteX, NbrSiteY, true); 
#ifdef __LAPACK__     
    TmpHam.LapackDiagonalize(TmpHam2, Q);
#else
    TmpHam.Diagonalize(TmpHam2, Q);
#endif
    if (Manager.GetBoolean("compute-density") == true)
      {
	for(int X1 = 0; X1 <NbrSiteX ; X1++)
	  {
	    for(int Y1 = 0; Y1 <NbrSiteY ; Y1++)
	      {
		for(int X2 = 0; X2 <NbrSiteX ; X2++)
		  {
		    for(int Y2 = 0; Y2 <NbrSiteY ; Y2++)
		      {
			ComputeDensity(NbrParticles, Q,X1,Y1, X2,Y2,(TightBindingModelHofstadterFiniteCylinder *) TightBindingModel);
		      }
		  }
	      }
	  }
      }
    
    delete TightBindingModel;
    return 0;
    }
  

  int MinKx = 0;
  int MaxKx = NbrSiteX - 1;

  if(Manager.GetBoolean("no-translation") == true)
    {  
      MaxKx = 0;
    }

   if(Manager.GetInteger("only-kx") != -1)
   {
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = Manager.GetInteger("only-kx");
   }

  TightBindingModel2DAtomicLimitLattice * TightBindingModel1 = 0;
  double * ChemicalPotential= new double[NbrSiteX* NbrSiteY];
  for(int i = 0 ; i <NbrSiteX* NbrSiteY ; i++)
   ChemicalPotential[i] =0.0 ;
  if(Manager.GetBoolean("atomic-limit") == true)
	{
  double * ChemicalPotential= new double[NbrSiteX* NbrSiteY];
  for(int i = 0 ; i <NbrSiteX* NbrSiteY ; i++)
   ChemicalPotential[i] =0.0 ;
 TightBindingModel1 = new  TightBindingModel2DAtomicLimitLattice(NbrSiteX, 1,NbrSiteY, ChemicalPotential,0,0, Architecture.GetArchitecture());
}
 
  TightBindingModelHofstadterFiniteCylinder  * TightBindingModel  = new TightBindingModelHofstadterFiniteCylinder(NbrSiteX, NbrSiteY, Flux, TunnelElementX, TunnelElementY, Axis, Manager.GetDouble("gamma-x"), Manager.GetDouble("gamma-y"), Architecture.GetArchitecture(),Manager.GetDouble("flux-inserted"), true,Manager.GetBoolean("torus"));

  HermitianMatrix TightBindingMatrix = 0;

  if(Manager.GetBoolean("atomic-limit") == true)
     TightBindingMatrix = TightBindingModel1->GetRealSpaceTightBindingHamiltonian();  
  else
     TightBindingMatrix = TightBindingModel->GetRealSpaceTightBindingHamiltonian();  

  bool FirstRunFlag = true;
 RealDiagonalMatrix TmpHam2(TightBindingMatrix.GetNbrRow());
#ifdef __LAPACK__     
  TightBindingMatrix.LapackDiagonalize(TmpHam2);
#else
  TightBindingMatrix.Diagonalize(TmpHam2);
#endif   
/*
  for (int i = 0; i < TmpHam2.GetNbrRow(); ++i)
  {
    cout << i << " : " << TmpHam2[i] << endl;
  }
*/

  RealSymmetricMatrix DensityDensityInteraction(TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), true);
  double UPotential = Manager.GetDouble("u-potential");

  if ( UPotential!= 0.0)
	{
	  if(Manager.GetBoolean("synthetic-dimension") ==false)
	      {
  			for (int x = 0; x <  NbrSiteX; ++x)
			{
			  int y = 0;
			  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x+1, y), UPotential);
        	 	  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y+1), UPotential);
			  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-1, y), UPotential);
			  for (y=1; y <  NbrSiteY-1; y++)
			    {
				  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x+1, y), UPotential);
				  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y+1), UPotential);
				  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-1, y), UPotential);
				  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y-1), UPotential);
		  	    }
                	  y = NbrSiteY - 1;
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x+1, y), UPotential);
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x-1, y), UPotential);
		  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y-1), UPotential);
		}
	}
	else
	{
		for (int x = 0; x <  NbrSiteX; ++x)
		{
		  for (int y = 0; y <  NbrSiteY; ++y)
		    {
                       for (int yint = 0; yint <  NbrSiteY; ++yint)
		    {
	        	  DensityDensityInteraction.SetMatrixElement(TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x, y),  TightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(x,  yint), UPotential);
	  	    }
}
		}
	
	}
//cout <<"interaction Matrix " <<endl<<DensityDensityInteraction <<endl;
}


  for (int i = MinKx; i <= MaxKx; ++i)
    {
	  cout << "(kx=" << i << ") : " << endl;
	  ParticleOnSphere* Space = 0;
	  AbstractQHEHamiltonian* Hamiltonian = 0;

         if(Manager.GetBoolean("no-translation") == true)
   {
       Space = new FermionOnLatticeRealSpace(NbrParticles,TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand());


 if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
     Memory = Architecture.GetArchitecture()->GetLocalMemory();

 Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());


       Hamiltonian = new ParticleOnLatticeRealSpaceHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(), TightBindingMatrix, DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
}
   else
{
 Space = new FermionOnLatticeRealSpaceAnd1DTranslation(NbrParticles,NbrSiteX*NbrSiteY,i,NbrSiteX);

 if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
     Memory = Architecture.GetArchitecture()->GetLocalMemory();

 Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());

 Hamiltonian = new ParticleOnLatticeRealSpaceAnd1DTranslationHamiltonian (Space, NbrParticles, TightBindingModel->GetNbrBands() * TightBindingModel->GetNbrStatePerBand(),  i,  NbrSiteX, TightBindingMatrix, DensityDensityInteraction, Architecture.GetArchitecture(), Memory);
}

  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
  double Shift = 0.0;
  Hamiltonian->ShiftHamiltonian(Shift);


 if (Manager.GetString("energy-expectation") != 0 )
	{
	  char* StateFileName = Manager.GetString("energy-expectation");
	  if (IsFile(StateFileName) == false)
	    {
	      cout << "state " << StateFileName << " does not exist or can't be opened" << endl;
	      return -1;           
	    }
	  ComplexVector State;
	  if (State.ReadVector(StateFileName) == false)
	    {
	      cout << "error while reading " << StateFileName << endl;
	      return -1;
	    }
	  if (State.GetVectorDimension() != Space->GetHilbertSpaceDimension())
	    {
	      cout << "error: vector and Hilbert-space have unequal dimensions"<<endl;
	      return -1;
	    }
	  ComplexVector TmpState(Space->GetHilbertSpaceDimension());
	  VectorHamiltonianMultiplyOperation Operation (Hamiltonian, &State, &TmpState);
	  Operation.ApplyOperation(Architecture.GetArchitecture());
	  Complex EnergyValue = State * TmpState;
	  cout << "< Energy > = " << (EnergyValue.Re - Shift) << " " << EnergyValue.Im << endl;
	  return 0; 
	}

	  
	  char* ContentPrefix = new char[256];
	  sprintf (ContentPrefix, "%d", i);
	  char* EigenstateOutputFile = new char [512];

	  sprintf (EigenstateOutputFile,"%s_kx_%d",FilePrefix, i);
	    
	  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  cout << "------------------------------------" << endl;
	  delete Hamiltonian;
	  delete Space;
	  delete[] EigenstateOutputFile;
	  delete[] ContentPrefix;
    }
  delete TightBindingModel;
  return 0;
}


void  ComputeDensity(int NbrFermions, ComplexMatrix & eigenVectors, int X1, int Y1,  int X2, int Y2, TightBindingModelHofstadterFiniteCylinder * tightBindingModel)
{
  int Index1 = tightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(X1, Y1);
  int Index2 = tightBindingModel->GetRealSpaceTightBindingLinearizedIndexSafe(X2, Y2);
  Complex Result = 0.0;
  for(int i =0; i < NbrFermions;i++)
  {
     for(int j =0; j < NbrFermions;j++)
     {
        Result+= eigenVectors.GetMatrixElement(i,Index1)*Conj(eigenVectors.GetMatrixElement(i,Index1)) * eigenVectors.GetMatrixElement(j,Index2)*Conj(eigenVectors.GetMatrixElement(j,Index2));
     }
  }
 cout <<Index1 << " " << Index2<< " " << Result<<endl;
}


void  ComputeCurrent(FermionOnLatticeRealSpaceAnd1DTranslation * space, double * tunnelElementX, double * tunnelElementY, double fluxDensity, ComplexVector & groundState)
{
 Complex Test = 0.0;
 int Ly = space->GetNbrSites() / space->GetMaxXMomentum();
cout <<"Flux density : "<<fluxDensity<<endl;
 for(int y =0; y < Ly;y++)
  {
 for(int x =0; x < space->GetMaxXMomentum();x++)
  { 
   ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator Operator(space, space->GetLinearizedIndexSafe(x+1, y), space->GetLinearizedIndexSafe(x, y));

   Complex Result = Operator.PartialMatrixElement(groundState,groundState,0,space->GetHilbertSpaceDimension());

   Result *= tunnelElementX[y]*Phase(2.0*M_PI*fluxDensity*((double) y));
   Test += 2.0*Result.Re;

   cout <<x<<" " <<y <<" " <<Result<<endl;
   }
 }

 for(int x =0; x < space->GetMaxXMomentum();x++)
  { 
 for(int y =0; y < Ly - 1;y++)
  {
   ParticleOnLatticeRealSpaceAnd1DTranslationOneBodyOperator Operator(space, space->GetLinearizedIndexSafe(x, y+1), space->GetLinearizedIndexSafe(x, y));
   Complex Result = tunnelElementY[y]*Operator.PartialMatrixElement(groundState,groundState,0,space->GetHilbertSpaceDimension());
   cout <<x<<" " <<y <<" " <<Result<<endl;
   Test += 2.0*Result.Re;
   }
 }

cout <<"Estimated Energy = " << Test<<endl;

}
