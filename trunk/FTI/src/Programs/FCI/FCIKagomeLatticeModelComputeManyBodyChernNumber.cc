#include "Options/Options.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU2SpinMomentumSpace.h"



#include "Hamiltonian/ParticleOnLatticeKagomeLatticeSingleBandHamiltonian.h"
#include "Tools/FTITightBinding/TightBindingModelKagomeLattice.h"
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
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t1", "real part of the nearest neighbor hopping amplitude", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "t2", "real part of the next nearest neighbor hopping amplitude", -0.3);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l1", "imaginary part of the nearest neighbor hopping amplitude", 0.28);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "l2", "imaginary part of the next nearest neighbor hopping amplitude", 0.2);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "mu-s", "sublattice staggered chemical potential", 0.0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "number-point-x", "number of points in the discretization of theta_x", 10);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "number-point-y", "number of points in the discretization of theta_y", 10);
  (*SystemGroup) += new BooleanOption  ('\n', "boson", "use bosonic statistics instead of fermionic statistics");

  (*SystemGroup) += new SingleDoubleOption  ('\n', "u-potential", "repulsive nearest neighbor potential strength (or on-site density-density potential for bosons)", 1.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "v-potential", "repulsive next nearest neighbor potential strength (or nearest neighbor density-density potential for bosons)", 0.0);
  (*SystemGroup) += new SingleDoubleOption  ('\n', "w-potential", "repulsive next next nearest neighbor potential strength (or next nearest neighbor density-density potential for bosons)", 0.0);
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
  
  int TotalKx = Manager.GetInteger("kx"); 
  int TotalKy = Manager.GetInteger("ky");
  int NbrPointX = Manager.GetInteger("number-point-x");
  int NbrPointY = Manager.GetInteger("number-point-y");
  
  BosonOnSquareLatticeMomentumSpace * Space = new BosonOnSquareLatticeMomentumSpace (NbrParticles, NbrSitesX, NbrSitesY,  TotalKx, TotalKy);
  
  int IncNbrPointX =  NbrPointX+3;
  int IncNbrPointY =  NbrPointY+3;
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
    
  
  char* StateFilePrefix = new char [256];
  sprintf (StateFilePrefix, "%s_singleband_kagomelattice_n_%d_x_%d_y_%d", StatisticPrefix, NbrParticles, NbrSitesX, NbrSitesY);
  
  char* FileParameterString = new char [256];
  sprintf (FileParameterString, "t1_%g_t2_%g_l1_%g_l2_%g", Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"));
  
  char* InteractionName = new char [256];
  if ((Manager.GetDouble("v-potential") == 0.0) && (Manager.GetDouble("w-potential") == 0.0))
    {
      sprintf(InteractionName, "u_%g", Manager.GetDouble("u-potential"));
    }
  else
    {
      if (Manager.GetDouble("w-potential") == 0.0)
	{
	  sprintf(InteractionName, "u_%g_v_%g", Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"));
	}
      else
	{
	  sprintf(InteractionName, "u_%g_v_%g_w_%g", Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"));
	}
   }
 
  char* EigenstateFile = new char [512 + strlen(StateFilePrefix) + strlen(FileParameterString) + strlen(InteractionName)];
  char* EigenvalueOutputFile = new char [512 + strlen(StateFilePrefix) + strlen(FileParameterString) + strlen(InteractionName)];
  char* CommentLine = new char [256];
  sprintf (CommentLine, "eigenvalues\n# kx ky ");
  
  char* FilePrefix = new char [256];
  char* ContentPrefix = new char[256];
  sprintf (ContentPrefix, "%d %d",  TotalKx,  TotalKy);
  
  
  ComplexVector ** BandEigenvectors = new ComplexVector * [NbrSitesX*NbrSitesY];
  for(int i = 0; i <NbrSitesX*NbrSitesY; i++) 
    BandEigenvectors[i] = new ComplexVector[IncNbrPointX * IncNbrPointY];
  
  
  ComplexVector * ManyBodyState = new ComplexVector [IncNbrPointX * IncNbrPointY];
  
  
  AbstractHamiltonian *  Hamiltonian  = 0;
  
  double TrueGammaX =   ( ( - 1.0) / ( (double) NbrPointX));
  double TrueGammaY =   ( ( - 1.0) / ( (double) NbrPointY));
  
  TightBindingModelKagomeLattice TightBindingModel1 (NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
						     TrueGammaX,  TrueGammaY , Architecture.GetArchitecture(), ExportOneBody);
  
  for (int Kx = 0; Kx < NbrSitesX; Kx++)
    for (int Ky = 0; Ky < NbrSitesY; Ky++)
      {
	int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	BandEigenvectors[MomentumIndex][0] =  TightBindingModel1.GetOneBodyMatrix(MomentumIndex)[BandIndex];
      }
  
  if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
    Memory = Architecture.GetArchitecture()->GetLocalMemory();
  Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
  Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"),
									&TightBindingModel1, true, Architecture.GetArchitecture(), Memory);

  sprintf (EigenvalueOutputFile, "%s_%s_%s_gx_%g_gy_%g_kx_%d_ky_%d.dat", StateFilePrefix, FileParameterString, InteractionName,  TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
  sprintf (EigenstateFile, "%s_%s_%s_gx_%g_gy_%g_kx_%d_ky_%d", StateFilePrefix, FileParameterString, InteractionName, TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
  

  FirstRunFlag = true;
  GenericComplexMainTask Task(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateFile);
  FirstRunFlag = false;
  MainTaskOperation TaskOperation (&Task);
  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
  delete Hamiltonian;
  delete[] EigenstateFile;
  delete[] EigenvalueOutputFile;
  
  for (int GammaX = 1; GammaX < IncNbrPointX; GammaX++)
    {
      TrueGammaX =   ( ( (double) GammaX - 1.0) / ( (double) NbrPointX));
      TrueGammaY =   ( (  - 1.0) / ( (double) NbrPointY));
      
      TightBindingModelKagomeLattice TightBindingModel2 (NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
							 TrueGammaX,  TrueGammaY , Architecture.GetArchitecture(), ExportOneBody);
      
      for (int Kx = 0; Kx < NbrSitesX; Kx++)
	for (int Ky = 0; Ky < NbrSitesY; Ky++)
	  {
	    int MomentumIndex = TightBindingModel2.GetLinearizedMomentumIndexSafe(Kx,Ky);
	    BandEigenvectors[MomentumIndex][GammaX*IncNbrPointY] =  TightBindingModel2.GetOneBodyMatrix(MomentumIndex)[BandIndex];
	  }
      
      for (int LinearizedMomentum = 0; LinearizedMomentum < NbrSitesX * NbrSitesY ; LinearizedMomentum++)
	{
	  Complex Tmp =  BandEigenvectors[LinearizedMomentum][(GammaX-1)*IncNbrPointY] * BandEigenvectors[LinearizedMomentum][GammaX *IncNbrPointY];
	  BandEigenvectors[LinearizedMomentum][GammaX*IncNbrPointY] *= (Conj(Tmp)/Norm(Tmp));
	}
      
      if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	Memory = Architecture.GetArchitecture()->GetLocalMemory();
      Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
      Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
									    Manager.GetDouble("w-potential"), &TightBindingModel2, true, Architecture.GetArchitecture(), Memory);
      EigenstateFile = new char [512];
      EigenvalueOutputFile = new char [512];
      
      
      sprintf (EigenvalueOutputFile, "%s_%s_%s_gx_%g_gy_%g_kx_%d_ky_%d.dat",StateFilePrefix, FileParameterString, InteractionName,  TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
      sprintf (EigenstateFile, "%s_%s_%s_gx_%g_gy_%g_kx_%d_ky_%d", StateFilePrefix, FileParameterString, InteractionName, TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
      
      FirstRunFlag = true;
      GenericComplexMainTask Task1(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateFile);
      FirstRunFlag = false;
      MainTaskOperation TaskOperation (&Task1);
      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
      delete Hamiltonian;
      delete[] EigenstateFile;
      delete[] EigenvalueOutputFile;   
    }


    for (int GammaY = 1; GammaY <  IncNbrPointY; GammaY++)
     {
       TrueGammaX =   ( ( -1.0) / ( (double) NbrPointX));
       TrueGammaY =   ( ( (double)  GammaY - 1.0) / ( (double) NbrPointY));

       TightBindingModelKagomeLattice TightBindingModel2 (NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
							  TrueGammaX,  TrueGammaY , Architecture.GetArchitecture(), ExportOneBody);
       
       for (int Kx = 0; Kx < NbrSitesX; Kx++)
	 for (int Ky = 0; Ky < NbrSitesY; Ky++)
	   {
	     int MomentumIndex = TightBindingModel2.GetLinearizedMomentumIndexSafe(Kx,Ky);
	     BandEigenvectors[MomentumIndex][GammaY] =  TightBindingModel2.GetOneBodyMatrix(MomentumIndex)[BandIndex];
	   }
       
       for (int LinearizedMomentum = 0; LinearizedMomentum < NbrSitesX * NbrSitesY ; LinearizedMomentum++)
	 {
	   Complex Tmp =  BandEigenvectors[LinearizedMomentum][GammaY-1]  * BandEigenvectors[LinearizedMomentum][GammaY];
	   BandEigenvectors[LinearizedMomentum][GammaY] *= (Conj(Tmp)/Norm(Tmp));
	 }
       
       if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	 Memory = Architecture.GetArchitecture()->GetLocalMemory();
       Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
       Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), 
									     Manager.GetDouble("w-potential"), &TightBindingModel2, true, Architecture.GetArchitecture(), Memory);
       EigenstateFile = new char [512];
       EigenvalueOutputFile = new char [512];
       
       sprintf (EigenvalueOutputFile, "%s_%s_%s_gx_%g_gy_%g_kx_%d_ky_%d.dat", StateFilePrefix, FileParameterString, InteractionName, TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
       sprintf (EigenstateFile, "%s_%s_%s_gx_%g_gy_%g_kx_%d_ky_%d",StateFilePrefix, FileParameterString, InteractionName, TrueGammaX, TrueGammaY, TotalKx, TotalKy);
       
       FirstRunFlag = true;
       GenericComplexMainTask Task1(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateFile);
       FirstRunFlag = false;
       MainTaskOperation TaskOperation1 (&Task1);
       TaskOperation1.ApplyOperation(Architecture.GetArchitecture());
       delete Hamiltonian;
       delete[] EigenstateFile;
       delete[] EigenvalueOutputFile;   
     }
    
    
    for (int GammaY = 1; GammaY <  IncNbrPointY ; GammaY++)
      {
	for (int  GammaX = 1; GammaX < IncNbrPointX; GammaX++)
	  {
	    
	    cout << GammaY << " " << GammaX<<endl;
	    TrueGammaX =   ( ( (double) GammaX - 1) / ( (double) NbrPointX));
	    TrueGammaY =   ( ( (double) GammaY - 1) / ( (double) NbrPointY));
	    
	    TightBindingModelKagomeLattice TightBindingModel2 (NbrSitesX, NbrSitesY,  Manager.GetDouble("t1"), Manager.GetDouble("t2"), Manager.GetDouble("l1"), Manager.GetDouble("l2"), Manager.GetDouble("mu-s"), 
							       TrueGammaX,  TrueGammaY , Architecture.GetArchitecture(), ExportOneBody);
	    
	    for (int Kx = 0; Kx < NbrSitesX; Kx++)
	      for (int Ky = 0; Ky < NbrSitesY; Ky++)
		{
		  int MomentumIndex = TightBindingModel2.GetLinearizedMomentumIndexSafe(Kx,Ky);
		  BandEigenvectors[MomentumIndex][GammaX*IncNbrPointY+ GammaY] =  TightBindingModel2.GetOneBodyMatrix(MomentumIndex)[BandIndex];
		}
	    
	    for (int LinearizedMomentum = 0; LinearizedMomentum < NbrSitesX * NbrSitesY ; LinearizedMomentum++)
	      {
		Complex Tmp =  BandEigenvectors[LinearizedMomentum][(GammaX-1)*IncNbrPointY + GammaY]  * BandEigenvectors[LinearizedMomentum][GammaX*IncNbrPointY+GammaY];
		BandEigenvectors[LinearizedMomentum][GammaX*IncNbrPointY + GammaY] *= (Conj(Tmp)/Norm(Tmp));
	      }
	    
	    if (Architecture.GetArchitecture()->GetLocalMemory() > 0)
	      Memory = Architecture.GetArchitecture()->GetLocalMemory();
	    Architecture.GetArchitecture()->SetDimension(Space->GetHilbertSpaceDimension());	
	    Hamiltonian = new ParticleOnLatticeKagomeLatticeSingleBandHamiltonian(Space, NbrParticles, NbrSitesX, NbrSitesY,Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"), Manager.GetDouble("w-potential"), 
										  &TightBindingModel2, true, Architecture.GetArchitecture(), Memory);
	    EigenstateFile = new char [512];
	    EigenvalueOutputFile = new char [512];
	    
	    sprintf (EigenvalueOutputFile, "%s_%s_%sg_gx_%g_gy_%g_kx_%d_ky_%d.dat", StateFilePrefix, FileParameterString, InteractionName,  TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
	    sprintf (EigenstateFile, "%s_%s_%s_gx_%g_gy_%g_kx_%d_ky_%d", StateFilePrefix, FileParameterString, InteractionName, TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
	    
	    FirstRunFlag = true;
	    GenericComplexMainTask Task1(&Manager, Hamiltonian->GetHilbertSpace(), &Lanczos, Hamiltonian, ContentPrefix, CommentLine, 0.0,  EigenvalueOutputFile, FirstRunFlag, EigenstateFile);
	    FirstRunFlag = false;
	    MainTaskOperation TaskOperation1 (&Task1);
	    TaskOperation1.ApplyOperation(Architecture.GetArchitecture());
	    delete Hamiltonian;
	    delete[] EigenstateFile;
	    delete[] EigenvalueOutputFile;   
	  }
      }
    
    
    EigenstateFile = new char [512];

    for (int GammaX = 0; GammaX <   IncNbrPointX; GammaX++)
      {
	for (int GammaY = 0; GammaY <   IncNbrPointY ; GammaY++)
	  {
	    
	    double TrueGammaX =   ( ( (double) GammaX - 1) / ( (double) NbrPointX));
	    double TrueGammaY =   ( ( (double) GammaY - 1) / ( (double) NbrPointY));
	    
	    if (Manager.GetDouble("v-potential") == 0.0)
	      {
        	sprintf (EigenstateFile, "%s_%s_u_%g_gx_%g_gy_%g_kx_%d_ky_%d.0.vec",StateFilePrefix, FileParameterString, Manager.GetDouble("u-potential"), TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
		
	      }
	    else
	      {
        	sprintf (EigenstateFile, "%s_%s_u_%g_v_%g_gx_%g_gy_%g_kx_%d_ky_%d.0.vec",StateFilePrefix, FileParameterString, Manager.GetDouble("u-potential"), Manager.GetDouble("v-potential"),  TrueGammaX,  TrueGammaY, TotalKx, TotalKy);
	      }
	    
	    if ( ( ManyBodyState[GammaX*IncNbrPointY+GammaY]).ReadVector(EigenstateFile) == false)
	      {
		cout << "error while reading " << EigenstateFile << endl;
		return -1;
	      }
	  }
      }
    
    Complex * OverlapMatrix = 0;
    
    for (int GammaX = 0; GammaX <  IncNbrPointX-1; GammaX++)
      {
	OverlapMatrix = new Complex [NbrSitesX*NbrSitesY];
	for (int Kx = 0; Kx < NbrSitesX; Kx++)
	  for (int Ky = 0; Ky < NbrSitesY; Ky++)
	    {
	      int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	      OverlapMatrix[MomentumIndex] = BandEigenvectors[MomentumIndex][GammaX*IncNbrPointY] * BandEigenvectors[MomentumIndex][(GammaX+1)*IncNbrPointY];
	    }
	Complex Tmp = ManyBodyState[GammaX*IncNbrPointY] *  ManyBodyState[(GammaX+1)*IncNbrPointY];
	ManyBodyState[(GammaX+1)*IncNbrPointY] *= (Conj(Tmp)/Norm(Tmp));
	delete [] OverlapMatrix ;
      }
    for (int GammaY = 0; GammaY <  IncNbrPointY-1; GammaY++)
      {
	OverlapMatrix = new Complex [NbrSitesX*NbrSitesY];
	for (int Kx = 0; Kx < NbrSitesX; Kx++)
	  for (int Ky = 0; Ky < NbrSitesY; Ky++)
	    {
	      int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	      OverlapMatrix[MomentumIndex] = BandEigenvectors[MomentumIndex][GammaY] * BandEigenvectors[MomentumIndex][GammaY+1];
	    }
	Complex Tmp = ManyBodyState[GammaY] * ManyBodyState[GammaY+1];
	ManyBodyState[GammaY+1] *= (Conj(Tmp)/Norm(Tmp));
	delete [] OverlapMatrix ;
     }

    for (int GammaY = 1; GammaY <  IncNbrPointY ; GammaY++)
     {
       for (int GammaX = 0; GammaX < IncNbrPointX-1; GammaX++)
       {
       OverlapMatrix = new Complex [NbrSitesX*NbrSitesY];
       for (int Kx = 0; Kx < NbrSitesX; Kx++)
        for (int Ky = 0; Ky < NbrSitesY; Ky++)
	   {
	     int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	     OverlapMatrix[MomentumIndex] = BandEigenvectors[MomentumIndex][GammaX*IncNbrPointY+GammaY] * BandEigenvectors[MomentumIndex][(GammaX+1)*IncNbrPointY+GammaY];
	   }
       Complex Tmp = ManyBodyState[GammaX*IncNbrPointY+GammaY] * ManyBodyState[(GammaX+1)*IncNbrPointY+GammaY];
       ManyBodyState[(GammaX+1)*IncNbrPointY+GammaY] *= (Conj(Tmp)/Norm(Tmp));
       delete [] OverlapMatrix ;
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
       
       
       Complex * OverlapMatrixPlusPlus = new Complex [NbrSitesX*NbrSitesY];
       Complex * OverlapMatrixMinusMinus = new Complex [NbrSitesX*NbrSitesY];
       Complex * OverlapMatrixPlusMinus = new Complex [NbrSitesX*NbrSitesY];
       Complex * OverlapMatrixMinusPlus = new Complex [NbrSitesX*NbrSitesY];
       
       for (int Kx = 0; Kx < NbrSitesX; Kx++)
	 for (int Ky = 0; Ky < NbrSitesY; Ky++)
	   {
	     int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	     OverlapMatrixPlusPlus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaPlusXModulo*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusYModulo];
	     OverlapMatrixMinusMinus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaMinusXModulo*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusYModulo];
	     OverlapMatrixPlusMinus[MomentumIndex] =  BandEigenvectors[MomentumIndex][ShiftedGammaPlusXModulo*IncNbrPointY + ShiftedGammaY] *  BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusYModulo];
	     OverlapMatrixMinusPlus[MomentumIndex] =  BandEigenvectors[MomentumIndex][ShiftedGammaMinusXModulo*IncNbrPointY + ShiftedGammaY] *  BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusYModulo];
	   }
       
       
       Complex Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusYModulo], ManyBodyState[ShiftedGammaPlusXModulo*IncNbrPointY + ShiftedGammaY],OverlapMatrixPlusMinus);
       ManyBodyChernNumber-= Tmp;
       cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaMinusYModulo << " ) -> ( " << ShiftedGammaPlusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;

       Tmp =  Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusYModulo],  ManyBodyState[ShiftedGammaMinusXModulo*IncNbrPointY + ShiftedGammaY], OverlapMatrixMinusPlus);
       cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaPlusYModulo << " ) -> ( " << ShiftedGammaMinusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;
       ManyBodyChernNumber-= Tmp;   
       
       Tmp =  Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusYModulo],  ManyBodyState[ShiftedGammaPlusXModulo*IncNbrPointY + ShiftedGammaY],OverlapMatrixPlusPlus);
       ManyBodyChernNumber+= Tmp;
       cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaPlusYModulo << " ) -> ( " << ShiftedGammaPlusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;
       
       Tmp =   Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusYModulo],  ManyBodyState[ShiftedGammaMinusXModulo*IncNbrPointY + ShiftedGammaY],OverlapMatrixMinusMinus);
       ManyBodyChernNumber+= Tmp;
       cout << " ( " << ShiftedGammaX<< " , " << ShiftedGammaMinusYModulo << " ) -> ( " << ShiftedGammaMinusXModulo << " , " <<ShiftedGammaY<< " ) "<< Tmp<<endl;
       
       
   delete [] OverlapMatrixPlusPlus;
   delete [] OverlapMatrixMinusMinus;
   delete [] OverlapMatrixMinusPlus;
   delete [] OverlapMatrixPlusMinus;
   
    }
  }
  
  
  
  
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
      
      
      Complex * OverlapMatrixPlus = new Complex [NbrSitesX*NbrSitesY];
      Complex * OverlapMatrixMinus = new Complex [NbrSitesX*NbrSitesY];
      
      for (int Kx = 0; Kx < NbrSitesX; Kx++)
	for (int Ky = 0; Ky < NbrSitesY; Ky++)
	  {
	    int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	    OverlapMatrixPlus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaPlusX*IncNbrPointY + ShiftedGammaY];
	    OverlapMatrixMinus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaMinusX*IncNbrPointY + ShiftedGammaY];
	  }
      
      Complex Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaPlusX*IncNbrPointY + ShiftedGammaY],OverlapMatrixPlus);
      ManyBodyChernNumberPathIntegral+= Tmp;
      
      
      Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaMinusX*IncNbrPointY + ShiftedGammaY],OverlapMatrixMinus);
      
      ManyBodyChernNumberPathIntegral-= Tmp;   
      
      
      delete [] OverlapMatrixPlus;
      delete [] OverlapMatrixMinus;
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
     
     
     Complex * OverlapMatrixPlus = new Complex [NbrSitesX*NbrSitesY];
     Complex * OverlapMatrixMinus = new Complex [NbrSitesX*NbrSitesY];
     
     for (int Kx = 0; Kx < NbrSitesX; Kx++)
       for (int Ky = 0; Ky < NbrSitesY; Ky++)
	 {
	   int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	   OverlapMatrixPlus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusY];
	   OverlapMatrixMinus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusY];
	 }
     
     Complex Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusY],OverlapMatrixPlus);
     ManyBodyChernNumberPathIntegral+= Tmp;
     
     Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusY],OverlapMatrixMinus);
     
     ManyBodyChernNumberPathIntegral-= Tmp;   
     
     
     delete [] OverlapMatrixPlus;
     delete [] OverlapMatrixMinus;
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
      
      
      Complex * OverlapMatrixPlus = new Complex [NbrSitesX*NbrSitesY];
      Complex * OverlapMatrixMinus = new Complex [NbrSitesX*NbrSitesY];
      
      for (int Kx = 0; Kx < NbrSitesX; Kx++)
	for (int Ky = 0; Ky < NbrSitesY; Ky++)
	  {
	    int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	    OverlapMatrixPlus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaPlusX*IncNbrPointY + ShiftedGammaY];
	    OverlapMatrixMinus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaMinusX*IncNbrPointY + ShiftedGammaY];
	  }
      
      Complex Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaPlusX*IncNbrPointY + ShiftedGammaY],OverlapMatrixPlus);
      ManyBodyChernNumberPathIntegral-= Tmp;
      
      
      Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaMinusX*IncNbrPointY + ShiftedGammaY],OverlapMatrixMinus);
      
      ManyBodyChernNumberPathIntegral+= Tmp;   
      
      
      delete [] OverlapMatrixPlus;
      delete [] OverlapMatrixMinus;
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
      
      
      Complex * OverlapMatrixPlus = new Complex [NbrSitesX*NbrSitesY];
      Complex * OverlapMatrixMinus = new Complex [NbrSitesX*NbrSitesY];
      
      for (int Kx = 0; Kx < NbrSitesX; Kx++)
	for (int Ky = 0; Ky < NbrSitesY; Ky++)
	  {
	    int MomentumIndex = TightBindingModel1.GetLinearizedMomentumIndexSafe(Kx,Ky);
	    OverlapMatrixPlus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusY];
	    OverlapMatrixMinus[MomentumIndex] = BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaY] * BandEigenvectors[MomentumIndex][ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusY];
	  }
      
      Complex Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaPlusY],OverlapMatrixPlus);
      ManyBodyChernNumberPathIntegral-= Tmp;
      
      
      Tmp = Space->ComputeOverlapWaveFunctionsWithDifferentGamma ( ManyBodyState[ShiftedGammaX*IncNbrPointY +  ShiftedGammaY], ManyBodyState[ShiftedGammaX*IncNbrPointY + ShiftedGammaMinusY],OverlapMatrixMinus);
      
      ManyBodyChernNumberPathIntegral+= Tmp;   
      
      
      delete [] OverlapMatrixPlus;
      delete [] OverlapMatrixMinus;
    }
  

  
  cout.precision(8);
  ManyBodyChernNumber /=(4.0 * M_PI);
  ManyBodyChernNumberPathIntegral/=(4.0 * M_PI);
  cout << "ManyBodyChernNumber = "<< ManyBodyChernNumber <<" " << Norm(ManyBodyChernNumber)<<" "<< Arg(ManyBodyChernNumber)/ M_PI <<endl; 
  
  cout << "ManyBodyChernNumberPathIntegral = "<<ManyBodyChernNumberPathIntegral<<" " << Norm(ManyBodyChernNumberPathIntegral)<<" "<< Arg(ManyBodyChernNumberPathIntegral)/ M_PI <<endl;  
}

