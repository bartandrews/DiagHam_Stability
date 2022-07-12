#include "Options/Options.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"

#include "Architecture/MonoProcessorArchitecture.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeRealSpaceHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionSingleBandHamiltonian.h"
#include "Hamiltonian/ParticleOnLatticeGenericDensityDensityInteractionTwoBandHamiltonian.h"


#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceAnd2DTranslation.h"


#include "Tools/FTITightBinding/TightBindingModelCheckerboardLattice.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"
#include "Tools/FTITightBinding/TightBindingModel2DAtomicLimitLattice.h"


#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"


#include "Vector/ComplexVector.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/StringTools.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>



// compute the description of the density-density interaction for the unit cell at the origin
//
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// bosonFlag = true if we are dealing with bosons
// uPotential = nearest neighbor (for fermions) or on-site (for bosons) interaction amplitude
// vPotential = next nearest neighbor (for fermions) or nearest neighbor (for bosons) interaction amplitude
void FCICheckerboardLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double uPotential, double vPotential);


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;

int main(int argc, char** argv)
{
  OptionManager Manager ("FCIOFLNOrbitalTriangularLatticeModel" , "0.01");
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
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 4);
  (*SystemGroup) += new SingleIntegerOption  ('x', "nbr-sitex", "number of sites along the x direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('y', "nbr-sitey", "number of sites along the y direction", 3);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "kx", "total x momentum", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "ky", "total y momentum", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "band-index", 0);


  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "nearest neighbor hoping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "next nearest neighbor hoping amplitude", 1.0 - 0.5 * M_SQRT2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "tpp", "second next nearest neighbor hoping amplitude", 0.5 * (M_SQRT2 - 1.0));
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);


  (*SystemGroup) += new SingleIntegerOption  ('\n', "number-point-x", "number of points in the discretization of theta_x", 10);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "number-point-y", "number of points in the discretization of theta_y", 10);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");

  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive next nearest neighbor potential strength", 0.0);
  (*PrecalculationGroup) += new SingleIntegerOption  ('m', "memory", "amount of memory that can be allocated for fast multiplication (in Mbytes)", 500);
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
#ifdef __SCALAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-scalapack", "use SCALAPACK libraries instead of DiagHam or LAPACK libraries");
#endif
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIOFLNOrbitalTriangularLatticeModel -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
    
  int NbrParticles = Manager.GetInteger("nbr-particles"); 
  int NbrSitesX = Manager.GetInteger("nbr-sitex"); 
  int NbrSitesY = Manager.GetInteger("nbr-sitey");
  int NbrSites = 2*NbrSitesX * NbrSitesY;
  int TotalKx = Manager.GetInteger("kx"); 
  int TotalKy = Manager.GetInteger("ky");
  int NbrPointX = Manager.GetInteger("number-point-x");
  int NbrPointY = Manager.GetInteger("number-point-y");
  int IncNbrPointX =  NbrPointX + 3;
  int IncNbrPointY =  NbrPointY + 3;
  long Memory = ((unsigned long) Manager.GetInteger("memory")) << 20;
  int BandIndex =  Manager.GetInteger("band-index");

  bool ExportOneBody = true; 
  bool FirstRunFlag = true;

  char* StatisticPrefix = new char [16];
  if (Manager.GetBoolean("boson") == false)
    {
      sprintf (StatisticPrefix, "fermions");
    }
  else
    {
      sprintf (StatisticPrefix, "bosons");
    }

  char* FilePrefix = new char [256];    
  sprintf (FilePrefix, "%s_realspace_checkerboardlattice_n_%d_ns_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSites, NbrSitesX, NbrSitesY);

  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t1_%f_t2_%f_tpp_%f", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"));
  
  char* ContentPrefix = new char[256];
  sprintf (ContentPrefix, "%d %d",  TotalKx,  TotalKy);
  
  double TrueGammaX =   0.0;
  double TrueGammaY =   0.0;
  
  TightBindingModelCheckerboardLattice  TightBindingModel (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"),  TrueGammaX,  TrueGammaY, Architecture.GetArchitecture(), true, true);
  
  
  int* NbrInteractingOrbitals;
  int** InteractingOrbitalsOrbitalIndices;
  int** InteractingOrbitalsSpatialIndices;
  double** InteractingOrbitalsPotentials;
  FCICheckerboardLatticeModelComputeInteractingOrbitals(NbrInteractingOrbitals, InteractingOrbitalsOrbitalIndices, 
							InteractingOrbitalsSpatialIndices, InteractingOrbitalsPotentials,
							Manager.GetBoolean("boson"), Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"));
  

  RealSymmetricMatrix DensityDensityInteraction(TightBindingModel.GetNbrBands() * TightBindingModel.GetNbrStatePerBand(), true);
  for (int x = 0; x < NbrSitesX; ++x)
    {
      for (int y = 0; y < NbrSitesY; ++y)
	{
	  for (int OrbitalIndex = 0; OrbitalIndex < TightBindingModel.GetNbrBands(); ++OrbitalIndex)
	    {
	      for (int k = 0; k < NbrInteractingOrbitals[OrbitalIndex]; ++k)
		{
		  DensityDensityInteraction.AddToMatrixElement(TightBindingModel.GetRealSpaceTightBindingLinearizedIndexSafe(x, y, OrbitalIndex),  TightBindingModel.GetRealSpaceTightBindingLinearizedIndexSafe(x + InteractingOrbitalsSpatialIndices[OrbitalIndex][2 * k], y + InteractingOrbitalsSpatialIndices[OrbitalIndex][(2 * k) + 1], InteractingOrbitalsOrbitalIndices[OrbitalIndex][k]),  InteractingOrbitalsPotentials[OrbitalIndex][k]);				  
		}
	    }
	}
    }
  
  
  ComplexVector * ManyBodyState = new ComplexVector [IncNbrPointX * IncNbrPointY];
  AbstractHamiltonian *  Hamiltonian  = 0;
  
  
  FermionOnLatticeRealSpaceAnd2DTranslation *  Space = new FermionOnLatticeRealSpaceAnd2DTranslation(NbrParticles, TightBindingModel.GetNbrBands() * TightBindingModel.GetNbrStatePerBand(), TotalKx, NbrSitesX, TotalKy, NbrSitesY);
  cout << "dim = " << Space->GetHilbertSpaceDimension()  << endl;
  
  
  MonoProcessorArchitecture TmpArchitecture;
  
  for (int GammaX = 0; GammaX < IncNbrPointX; GammaX++)
    {
      for (int GammaY = 0; GammaY <  IncNbrPointY ; GammaY++)
	{
	  cout <<GammaX<< " " <<GammaY<<endl;
	  TrueGammaX =   ( ( (double) GammaX - 1) / ( (double) NbrPointX));
	  TrueGammaY =   ( ( (double) GammaY - 1) / ( (double) NbrPointY));
	  
	  TightBindingModelCheckerboardLattice  TightBindingModel2 (NbrSitesX, NbrSitesY, Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("tpp"), Manager.GetDouble("mu-s"),  TrueGammaX,  TrueGammaY, &TmpArchitecture, true, true);
	  
	  
	  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	    Memory = Architecture.GetArchitecture()->GetLocalMemory();
	  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	  
	  
	  HermitianMatrix TightBindingMatrix = TightBindingModel2.GetRealSpaceTightBindingHamiltonian();
	  Hamiltonian = new ParticleOnLatticeRealSpaceAnd2DTranslationHamiltonian (Space, NbrParticles, TightBindingModel2.GetNbrBands() * TightBindingModel2.GetNbrStatePerBand(),   TotalKx, NbrSitesX, TotalKy, NbrSitesY,								       TightBindingMatrix, DensityDensityInteraction,  Architecture.GetArchitecture(), Memory);
	  char *  EigenstateOutputFile = new char [512];
	  char * EigenvalueOutputFile = new char [512];
	  
	  if (Manager.GetDouble("mu-s") == 0.0)
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f_kx_%d_ky_%d.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString,   TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
	      sprintf (EigenstateOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f_kx_%d_ky_%d",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString,   TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
	    }
	  else
	    {
	      sprintf (EigenvalueOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d.dat",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString,   TrueGammaX,  TrueGammaY, Manager.GetDouble("mu-s"), TotalKx, TotalKy);
	      sprintf (EigenstateOutputFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString,   TrueGammaX,  TrueGammaY, Manager.GetDouble("mu-s"), TotalKx, TotalKy);
	    }
	  
       
	  FirstRunFlag = true;
	  GenericComplexMainTask Task1(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateOutputFile);
	  FirstRunFlag = false;
	  MainTaskOperation TaskOperation (&Task1);
	  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
	  
	  delete Hamiltonian;
	  delete[] EigenstateOutputFile;
	  delete[] EigenvalueOutputFile;
	}
    }
  
  cout <<"reading states"<<endl;
  char* EigenstateFile = new char [512];
  
  for (int GammaX = 0; GammaX <   IncNbrPointX; GammaX++)
    {
      for (int GammaY = 0; GammaY <   IncNbrPointY ; GammaY++)
	{
	  
	  double TrueGammaX =   ( ( (double) GammaX - 1) / ( (double) NbrPointX));
	  double TrueGammaY =   ( ( (double) GammaY - 1) / ( (double) NbrPointY));
	  
	  
	  if (Manager.GetDouble("mu-s") == 0.0)
	    {
	      sprintf (EigenstateFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f_kx_%d_ky_%d.0.vec",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString,   TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
	    }
	  else
	    {
	      sprintf (EigenstateFile, "%s_u_%f_v_%f_%s_gx_%f_gy_%f_mus_%f_kx_%d_ky_%d.0.vec",FilePrefix, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), FileParameterString,   TrueGammaX,  TrueGammaY, Manager.GetDouble("mu-s"), TotalKx, TotalKy);
	    }
	  
	  
	  if ( ( ManyBodyState[GammaX*IncNbrPointY+GammaY]).ReadVector(EigenstateFile) == false)
	    {
	      cout << "error while reading " << EigenstateFile << endl;
	      return -1;
	    }
	}
    }
  
  for (int GammaX = 0; GammaX <  IncNbrPointX-1; GammaX++)
    {
      Complex Tmp =  ManyBodyState[GammaX*IncNbrPointY] * ManyBodyState[(GammaX+1)*IncNbrPointY];
      ManyBodyState[(GammaX+1)*IncNbrPointY] *= (Conj(Tmp)/Norm(Tmp));
    }
  
  for (int GammaY = 0; GammaY <  IncNbrPointY-1; GammaY++)
    {
      Complex Tmp =  ManyBodyState[GammaY] * ManyBodyState[GammaY+1];
      ManyBodyState[GammaY+1] *= (Conj(Tmp)/Norm(Tmp));
    }
  for (int GammaY = 1; GammaY <  IncNbrPointY ; GammaY++)
    {
      for (int GammaX = 0; GammaX < IncNbrPointX-1; GammaX++)
	{
	  Complex Tmp =  ManyBodyState[GammaX*IncNbrPointY+GammaY] * ManyBodyState[(GammaX+1)*IncNbrPointY+GammaY];
	  ManyBodyState[(GammaX+1)*IncNbrPointY+GammaY] *= (Conj(Tmp)/Norm(Tmp));
	}
    }
  
  Complex ManyBodyChernNumber (0.0);
  
  for (int ShiftedGammaX = 1; ShiftedGammaX <  NbrPointX +1; ShiftedGammaX++)
    {
      for (int ShiftedGammaY = 1; ShiftedGammaY <  NbrPointY + 1; ShiftedGammaY++)
	{
	  int GammaX = ShiftedGammaX - 1; 
	  int GammaY = ShiftedGammaY - 1;
	  
	  
	  int GammaPlusXModulo = (GammaX + 1);
	  int ShiftedGammaPlusXModulo = GammaPlusXModulo+1;
	  
	  int GammaPlusYModulo = (GammaY + 1);
	  int ShiftedGammaPlusYModulo = GammaPlusYModulo+1;
	  
	  int GammaMinusXModulo = (GammaX - 1);
	  int ShiftedGammaMinusXModulo = GammaMinusXModulo+1;
	  
	  int GammaMinusYModulo = (GammaY - 1);
	  int ShiftedGammaMinusYModulo = GammaMinusYModulo+1;
	  
	  
	  Complex Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusYModulo] * ManyBodyState[ShiftedGammaPlusXModulo*IncNbrPointY + ShiftedGammaY];
	  ManyBodyChernNumber-= Tmp;
	  
	  cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaMinusYModulo << " ) -> ( " << ShiftedGammaPlusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;
	  
	  Tmp =  ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusYModulo] *  ManyBodyState[ShiftedGammaMinusXModulo*IncNbrPointY + ShiftedGammaY];
	  
	  cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaPlusYModulo << " ) -> ( " << ShiftedGammaMinusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;
	  
	  ManyBodyChernNumber-= Tmp;   
	  
	  Tmp =  ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusYModulo] *  ManyBodyState[ShiftedGammaPlusXModulo*IncNbrPointY + ShiftedGammaY];
	  ManyBodyChernNumber+= Tmp;
	  cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaPlusYModulo << " ) -> ( " << ShiftedGammaPlusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;
	  
	  Tmp =  ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusYModulo] *  ManyBodyState[ShiftedGammaMinusXModulo*IncNbrPointY + ShiftedGammaY];
	  ManyBodyChernNumber+= Tmp;
	  
	  cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaMinusYModulo << " ) -> ( " << ShiftedGammaMinusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;
	  
	}
    }
  
  cout.precision(8);
  cout << "ManyBodyChernNumber = "<<ManyBodyChernNumber/(4.0 * M_PI)<<endl; 
  
  
  
  Complex ManyBodyChernNumberPathIntegral (0.0);
  
  for (int ShiftedGammaX = 1; ShiftedGammaX <=  NbrPointX + 1; ShiftedGammaX++)
    {
      int ShiftedGammaY = 1;
      
      int GammaX = ShiftedGammaX - 1; 
      int GammaY = ShiftedGammaY - 1;
      
      
      int GammaPlusXModulo = (GammaX + 1);
      int ShiftedGammaPlusX = GammaPlusXModulo+1;
      
      int GammaMinusXModulo = (GammaX - 1);
      int ShiftedGammaMinusX = GammaMinusXModulo+1;
      
      
      Complex Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaPlusX*IncNbrPointY + ShiftedGammaY];
      ManyBodyChernNumberPathIntegral+= Tmp;
      
      
      Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaMinusX*IncNbrPointY + ShiftedGammaY];
      
      ManyBodyChernNumberPathIntegral-= Tmp;   
    }
  
  
  
  
  
  for (int ShiftedGammaY = 1; ShiftedGammaY <=  NbrPointY + 1; ShiftedGammaY++)
    {
      int ShiftedGammaX =  NbrPointX + 1;
      
      int GammaX = ShiftedGammaX - 1; 
      int GammaY = ShiftedGammaY - 1;
      
      
      int GammaPlusYModulo = (GammaY + 1);
      int ShiftedGammaPlusY = GammaPlusYModulo+1;
      
      int GammaMinusYModulo = (GammaY - 1);
      int ShiftedGammaMinusY = GammaMinusYModulo+1;
      
      Complex Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusY];
      ManyBodyChernNumberPathIntegral+= Tmp;
      
      Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusY];
      
      ManyBodyChernNumberPathIntegral-= Tmp;   
    }

 
  
  for (int ShiftedGammaX = 1; ShiftedGammaX <=  NbrPointX + 1; ShiftedGammaX++)
    {
      int ShiftedGammaY =  NbrPointY + 1;
      
      int GammaX = ShiftedGammaX - 1; 
      int GammaY = ShiftedGammaY - 1;
      
      
      int GammaPlusXModulo = (GammaX + 1);
      int ShiftedGammaPlusX = GammaPlusXModulo+1;
      
      int GammaMinusXModulo = (GammaX - 1);
      int ShiftedGammaMinusX = GammaMinusXModulo+1;
      
      Complex Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaPlusX*IncNbrPointY + ShiftedGammaY];
      ManyBodyChernNumberPathIntegral-= Tmp;
      
      
      Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaMinusX*IncNbrPointY + ShiftedGammaY];
      
      ManyBodyChernNumberPathIntegral+= Tmp;   
    }
  
  
  
  
  for (int ShiftedGammaY = 1; ShiftedGammaY <=  NbrPointY + 1; ShiftedGammaY++)
    {
      int ShiftedGammaX =  1;
      
      int GammaX = ShiftedGammaX - 1; 
      int GammaY = ShiftedGammaY - 1;
      
      
      int GammaPlusYModulo = (GammaY + 1);
      int ShiftedGammaPlusY = GammaPlusYModulo+1;
      
      int GammaMinusYModulo = (GammaY - 1);
      int ShiftedGammaMinusY = GammaMinusYModulo+1;
      
      Complex Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusY];
      ManyBodyChernNumberPathIntegral-= Tmp;

      
      Tmp = ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY] * ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusY];
      
      ManyBodyChernNumberPathIntegral+= Tmp;   
      
    }
  cout << "ManyBodyChernNumberPathIntegral = "<<ManyBodyChernNumberPathIntegral/(4.0 * M_PI)<<endl; 
}



// compute the description of the density-density interaction for the unit cell at the origin
//
// nbrInteractingOrbitals = number of orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsOrbitalIndices = orbital indices of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsSpatialIndices = spatial indices (sorted as 2 consecutive integers) of the orbitals interacting with each orbital within the unit cell at the origin through a density-density term
// interactingOrbitalsPotentials = intensity of each density-density term 
// bosonFlag = true if we are dealing with bosons
// uPotential = nearest neighbor (for fermions) or on-site (for bosons) interaction amplitude
// vPotential = next nearest neighbor (for fermions) or nearest neighbor (for bosons) interaction amplitude

void FCICheckerboardLatticeModelComputeInteractingOrbitals(int*& nbrInteractingOrbitals, int**& interactingOrbitalsOrbitalIndices,
							   int**& interactingOrbitalsSpatialIndices, double**& interactingOrbitalsPotentials,
							   bool bosonFlag, double uPotential, double vPotential)
{
  nbrInteractingOrbitals = new int[2];
  interactingOrbitalsOrbitalIndices = new int*[2];
  interactingOrbitalsSpatialIndices = new int*[2];
  interactingOrbitalsPotentials = new double*[2];
  if (bosonFlag == false)
    {
      nbrInteractingOrbitals[0] = 1;
      nbrInteractingOrbitals[1] = 3;
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 2;
	  nbrInteractingOrbitals[1] += 2;
	}
      interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 2];
      interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 2];
      interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];

      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 1;
      interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[0][Index] = uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[0][Index] = 0;
	  interactingOrbitalsSpatialIndices[0][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[0][Index] = vPotential;
	  ++Index;	  
	  interactingOrbitalsOrbitalIndices[0][Index] = 0;
	  interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[0][Index] = vPotential;
	  ++Index;	  
	}
      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      interactingOrbitalsOrbitalIndices[1][Index] = 0;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
      interactingOrbitalsPotentials[1][Index] = uPotential;		  
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[1][Index] = 1;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[1][Index] = vPotential;
	  ++Index;	  
	  interactingOrbitalsOrbitalIndices[1][Index] = 1;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[1][Index] = vPotential;
	  ++Index;	  
	}
    }
  else
    {
      nbrInteractingOrbitals[0] = 1;
      nbrInteractingOrbitals[1] = 1;
      if (vPotential != 0.0)
	{
	  nbrInteractingOrbitals[0] += 1;
	  nbrInteractingOrbitals[1] += 3;
	}
      interactingOrbitalsOrbitalIndices[0] = new int[nbrInteractingOrbitals[0]];
      interactingOrbitalsSpatialIndices[0] = new int[nbrInteractingOrbitals[0] * 2];
      interactingOrbitalsPotentials[0] = new double[nbrInteractingOrbitals[0]];
      interactingOrbitalsOrbitalIndices[1] = new int[nbrInteractingOrbitals[1]];
      interactingOrbitalsSpatialIndices[1] = new int[nbrInteractingOrbitals[1] * 2];
      interactingOrbitalsPotentials[1] = new double[nbrInteractingOrbitals[1]];
  
      int Index = 0;
      interactingOrbitalsOrbitalIndices[0][Index] = 0;
      interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[0][Index] = 0.5 * uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[0][Index] = 1;
	  interactingOrbitalsSpatialIndices[0][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[0][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[0][Index] = vPotential;	  
	}
      Index = 0;
      interactingOrbitalsOrbitalIndices[1][Index] = 1;
      interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
      interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
      interactingOrbitalsPotentials[1][Index] = 0.5 * uPotential;
      ++Index;
      if (vPotential != 0.0)
	{
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 0;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 0;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	  interactingOrbitalsOrbitalIndices[1][Index] = 0;
	  interactingOrbitalsSpatialIndices[1][2 * Index] = 1;
	  interactingOrbitalsSpatialIndices[1][(2 * Index) + 1] = 1;
	  interactingOrbitalsPotentials[1][Index] = vPotential;		  
	  ++Index;
	}
    }
}


