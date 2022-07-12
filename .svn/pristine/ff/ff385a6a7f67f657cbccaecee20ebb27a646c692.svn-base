#include "Options/Options.h"

#include "HilbertSpace/FermionOnSquareLatticeWithSU3SpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/FermionOnSquareLatticeWithSpinMomentumSpaceLong.h"
#include "HilbertSpace/FermionOnSquareLatticeMomentumSpaceLong.h"
#include "HilbertSpace/BosonOnSquareLatticeMomentumSpace.h"
#include "HilbertSpace/BosonOnSquareLatticeWithSU3SpinMomentumSpace.h"

#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "Operator/AbstractOperator.h"
#include "Operator/ParticleOnSphereDensityOperator.h"


#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "Architecture/ArchitectureOperation/VectorOperatorMultiplyOperation.h"

#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "MainTask/GenericComplexMainTask.h"
#include "GeneralTools/FilenameTools.h"

#include "Tools/FQHEFiles/FQHEOnSquareLatticeFileTools.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

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
  OptionManager Manager ("FCIGenerateSMA" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* PrecalculationGroup = new OptionGroup ("precalculation options");

  ArchitectureManager Architecture;
  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += PrecalculationGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('\n', "eigenstate-file", "name of the vector file to which the SMA should be applied");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-groundstate", "name of the file that gives the vector files to which the SMA should be applied (in all-bilinear or sma mode)");	
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-kx", "only evalute a given x momentum sector (negative if all kx sectors have to be computed)", -1);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-ky", "only evalute a given y momentum sector (negative if all ky sectors have to be computed)", -1);  
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-gx", "only evalute a given x momentum sector (negative if only the first Brillouin zone has to be computed)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "only-gy", "only evalute a given y momentum sector (negative if only the first Brillouin zone has to be computed)", 0);  
  (*SystemGroup) += new BooleanOption  ('\n', "all-bilinear", "apply all bilinear operators to the ground state, without summing on the sector");
  (*SystemGroup) += new BooleanOption  ('\n', "compute-allsma", "compute all SMA");
  (*SystemGroup) += new BooleanOption  ('\n', "quantum-distance", "include the quantum distance when computing the SMA");
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a binary file");
  (*SystemGroup) += new SingleStringOption  ('\n', "embedding-file", "import embedding information from a text file (override information from onebody file)");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if ((Manager.GetString("eigenstate-file") == 0) && (Manager.GetString("degenerate-groundstate") == 0))
    {
      cout << "error, an eigenstate state file should be provided. See man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
  if (Manager.GetString("import-onebody") == 0)
    {
      cout << "error, a file giving the one-body information should be provided. See man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
 
  
  bool BilinearFlag = Manager.GetBoolean("all-bilinear");
  int NbrParticles = 0;
  int NbrSiteX = 0;
  int NbrSiteY = 0;
  int* TotalKx = 0;
  int* TotalKy = 0;
  int NbrSpaces = 1;
  char** GroundStateFiles = 0;
  bool Statistics = true;
  ComplexVector* GroundStates;
  
  if (Manager.GetString("degenerate-groundstate") == 0)
    {
      GroundStateFiles = new char* [1];
      TotalKx = new int[1];
      TotalKy = new int[1];
      GroundStates = new ComplexVector[1];
      GroundStateFiles[0] = new char [strlen(Manager.GetString("eigenstate-file")) + 1];
      strcpy (GroundStateFiles[0], Manager.GetString("eigenstate-file"));
      if (GroundStates[0].ReadVector (GroundStateFiles[0]) == false)
      {
	cout << "can't open vector file " << GroundStateFiles[0] << endl;
	return -1;      
      }	
      if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[0], NbrParticles, NbrSiteX, NbrSiteY, TotalKx[0], TotalKy[0], Statistics) == false)
      {
	cout << "error while retrieving system parameters from file name " << GroundStateFiles[0] << endl;
	return -1;
      }
    }
  else
    {
      MultiColumnASCIIFile DegeneratedFile;
      if (DegeneratedFile.Parse(Manager.GetString("degenerate-groundstate")) == false)
	{
	  DegeneratedFile.DumpErrors(cout);
	  return -1;
	}
      NbrSpaces = DegeneratedFile.GetNbrLines();
      GroundStateFiles = new char* [NbrSpaces];
      TotalKx = new int[NbrSpaces];
      TotalKy = new int[NbrSpaces];
      GroundStates = new ComplexVector[NbrSpaces];
      for (int i = 0; i < NbrSpaces; ++i)
	{
	  GroundStateFiles[i] = new char [strlen(DegeneratedFile(0, i)) + 1];
	  strcpy (GroundStateFiles[i], DegeneratedFile(0, i));	
	  if (GroundStates[i].ReadVector (GroundStateFiles[i]) == false)
	    {
	      cout << "can't open vector file " << GroundStateFiles[i] << endl;
	      return -1;      
	    }	
	  if (FQHEOnSquareLatticeFindSystemInfoFromVectorFileName(GroundStateFiles[i], NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i], Statistics) == false)
	    {
	      cout << "error while retrieving system parameters from file name " << GroundStateFiles[i] << endl;
	      return -1;
	    }
	}
    }
 
  
  int MinKx = 0;
  int MaxKx = NbrSiteX - 1;
  if (Manager.GetInteger("only-kx") >= 0)
    {						
      MinKx = Manager.GetInteger("only-kx");
      MaxKx = MinKx;
    }
  int MinKy = 0;
  int MaxKy = NbrSiteY - 1;
  if (Manager.GetInteger("only-ky") >= 0)
    {						
      MinKy = Manager.GetInteger("only-ky");
      MaxKy = MinKy;
    }
    
  
  ParticleOnSphere* SpaceSource = 0;
  Generic2DTightBindingModel TightBindingModel(Manager.GetString("import-onebody"));
  if (Manager.GetString("embedding-file") != 0 && TightBindingModel.SetEmbeddingFromAsciiFile(Manager.GetString("embedding-file")) == false )
    return 0;
	
  ParticleOnSphere* SpaceDestination = 0;

  if (Manager.GetBoolean("compute-allsma") == false)
    {
      for (int Qx0 = MinKx; Qx0 <= MaxKx; ++Qx0)
	    {
	      for (int Qy0 = MinKy; Qy0 <= MaxKy; ++Qy0)
		{
		  if (Statistics == true)
		    SpaceDestination = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, Qx0, Qy0);
		  else
		    SpaceDestination = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, Qx0, Qy0);
		  for (int i = 0; i < NbrSpaces; ++i)
		    {
		      cout << "N = " << NbrParticles << " Nx = " << NbrSiteX << " Ny = " << NbrSiteY << " Kx = "<< TotalKx[i] << " Ky = "<< TotalKy[i] <<  endl;
		      if (Statistics == true)
			SpaceSource = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		      else
			SpaceSource = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[i], TotalKy[i]);
		      
		      ComplexVector EigenstateOutput(SpaceDestination->GetHilbertSpaceDimension(), true);
		      SpaceSource->SetTargetSpace(SpaceDestination);
		      char* EigenstateOutputFile;
		      char* TmpExtention = new char [512];
		      for (int kx = 0; kx < NbrSiteX; ++kx)
			{
			  for (int ky = 0; ky < NbrSiteY; ++ky)
			    {
			      int KxQx = (kx + Qx0 - TotalKx[i]) % NbrSiteX;
			      if (KxQx < 0)
				KxQx += NbrSiteX;
			      int KyQy = (ky + Qy0 - TotalKy[i]) % NbrSiteY;
			      if (KyQy < 0)
				KyQy += NbrSiteY;
			      int indexDagger = TightBindingModel.GetLinearizedMomentumIndexSafe(KxQx, KyQy);
			      int index = TightBindingModel.GetLinearizedMomentumIndexSafe(kx, ky);
			      cout << "computing c^+_{"<< KxQx << "," << KyQy << "} c_{"<< kx << "," << ky << "} |Psi_" << i << ">" << endl;
			      ParticleOnSphereDensityOperator Projector(SpaceSource, indexDagger, index);
			      VectorOperatorMultiplyOperation Operation(&Projector, &(GroundStates[i]), &EigenstateOutput);
			      Operation.ApplyOperation(Architecture.GetArchitecture());
			      sprintf (TmpExtention, "_bilinear_kx_%d_ky_%d_qx0_%d_qy0_%d.vec", kx, ky, Qx0, Qy0);
			      EigenstateOutputFile = ReplaceExtensionToFileName(GroundStateFiles[i], ".vec", TmpExtention);
			      // 		EigenstateOutput.Normalize();
			      EigenstateOutput.WriteVector(EigenstateOutputFile);
			    }
			}
		      delete[] EigenstateOutputFile;
		    }   
		  delete SpaceDestination;
		  delete SpaceSource;
		}
	    }
	}
  else
    {
      int Qx = MinKx;
      int Qy = MinKy;
      int MinGx = 0;
      int MaxGx = 0;
      if (Manager.GetInteger("only-gx") != 0)
	{						
	  MinGx = Manager.GetInteger("only-gx");
	  MaxGx = MinGx;
	}
      int MinGy = 0;
      int MaxGy = 0;
      if (Manager.GetInteger("only-gy") != 0)
	{						
	  MinGy = Manager.GetInteger("only-gy");
	  MaxGy = MinGy;
	}
      if (Statistics == true)
	SpaceDestination = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, Qx, Qy);
      else
	SpaceDestination = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, Qx, Qy);
      if (Statistics == true)
	SpaceSource = new FermionOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[0], TotalKy[0]);
      else
	SpaceSource = new BosonOnSquareLatticeMomentumSpace(NbrParticles, NbrSiteX, NbrSiteY, TotalKx[0], TotalKy[0]);
      SpaceSource->SetTargetSpace(SpaceDestination);
      ComplexVector EigenstateOutput(SpaceDestination->GetHilbertSpaceDimension());
      ComplexVector TmpEigenstateOutput(SpaceDestination->GetHilbertSpaceDimension(), true);
      for (int Gx = MinGx; Gx <= MaxGx; ++Gx)
	{
	  for (int Gy = MinGy; Gy <= MaxGy; ++Gy)
	    {
	      EigenstateOutput.ClearVector();
	      char* EigenstateOutputFile;
	      char* TmpExtention = new char [512];
	      cout << "computing rho_{"<< Qx << "," << Qy << "}^{"<< Gx << "," << Gy << "} |Psi_0>" << endl;
	      for (int kx = 0; kx < NbrSiteX; ++kx)
		{
		  for (int ky = 0; ky < NbrSiteY; ++ky)
		    {
		      int KxQx = (kx + Qx - TotalKx[0]) % NbrSiteX;
		      if (KxQx < 0)
			KxQx += NbrSiteX;
		      int KyQy = (ky + Qy - TotalKy[0]) % NbrSiteY;
		      if (KyQy < 0)
			KyQy += NbrSiteY;
		      int indexDagger = TightBindingModel.GetLinearizedMomentumIndexSafe(KxQx, KyQy);
		      int index = TightBindingModel.GetLinearizedMomentumIndexSafe(kx, ky);
		      ParticleOnSphereDensityOperator Projector(SpaceSource, indexDagger, index);
		      VectorOperatorMultiplyOperation Operation(&Projector, &(GroundStates[0]), &TmpEigenstateOutput);
		      Operation.ApplyOperation(Architecture.GetArchitecture());
		      Complex TmpFactor = 0.0;
// 		      if (Manager.GetBoolean("quantum-distance"))
// 			TmpFactor = Conj(TightBindingModel.GetAbelianConnectionQuantumDistance(kx, ky, Qx + Gx * NbrSiteX, Qy + Gy * NbrSiteY, 0));
// 		      else
// 			TmpFactor = Conj(TightBindingModel.GetAbelianConnection(kx, ky, Qx + Gx * NbrSiteX, Qy + Gy * NbrSiteY, 0));
		      
		      if (Manager.GetBoolean("quantum-distance"))
			TmpFactor = Conj(TightBindingModel.GetAbelianConnectionQuantumDistance(kx, ky, Qx + Gx * NbrSiteX - TotalKx[0], Qy + Gy * NbrSiteY - TotalKy[0], 0));
		      else
			TmpFactor = Conj(TightBindingModel.GetAbelianConnection(kx, ky, Qx + Gx * NbrSiteX- TotalKx[0], Qy + Gy * NbrSiteY- TotalKy[0], 0));
			

		      EigenstateOutput. AddLinearCombination(TmpFactor, TmpEigenstateOutput);
		    }
		}
	      double TmpNorm = EigenstateOutput.Norm();
	      cout << "norm = " << TmpNorm << endl;
	      EigenstateOutput /= TmpNorm;
	      sprintf (TmpExtention, "_sma_qx_%d_qy_%d_gx_%d_gy_%d.vec", Qx, Qy, Gx, Gy);
	      EigenstateOutputFile = ReplaceExtensionToFileName(GroundStateFiles[0], ".vec", TmpExtention);
	      EigenstateOutput.WriteVector(EigenstateOutputFile);
	    }
	}
      delete SpaceDestination;
      delete SpaceSource;
    }  
  return 0;
}
