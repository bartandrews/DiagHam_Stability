#include "Vector/ComplexVector.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/OperatorMatrixElementOperation.h"

#include "Tools/FTIFiles/FTIHubbardModelFileTools.h"
#include "Tools/FTITightBinding/Abstract2DTightBindingModel.h"
#include "Tools/FTITightBinding/Generic2DTightBindingModel.h"

#include "HilbertSpace/FermionOnLatticeRealSpace.h"
#include "HilbertSpace/FermionOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeRealSpace.h"
#include "HilbertSpace/BosonOnLatticeRealSpaceAnd2DTranslation.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpace.h"
#include "HilbertSpace/BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation.h"

#include "Operator/ParticleOnSphereDensityOperator.h"
#include "Operator/ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <sys/time.h>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);
  OptionManager Manager ("FCIRealSpaceDensity" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleStringOption  ('i', "input-state", "name of the file ");
  (*SystemGroup) += new SingleStringOption  ('\n', "degenerate-states", "single column file describing a set of degenerate states");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "band-index", "index of the band whose occupation has to be computed, -1 if all have to be computed", -1);
  (*SystemGroup) += new SingleStringOption('\n', "import-onebody", "import information on the tight binding model from a binary file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FCIRealSpaceDensity -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrParticles = 0;
  int NbrSites = 0;
  bool Statistics = true;
  int NbrStates = 0;
  int MomentumFlag = false;
  int KxMomentum = 0;
  int XPeriodicity = 0;
  int KyMomentum = 0;
  int YPeriodicity = 0;
  bool GutzwillerFlag = false;

  if ((Manager.GetString("input-state") == 0) && (Manager.GetString("degenerate-states") == 0))
    {
      cout << "error, a state file should be provided. See man page for option syntax or type FTIRealSpaceComputeS2 -h" << endl;
      return -1;
    }
    
  if (Manager.GetString("import-onebody") == 0)
    {
      cout << "error, a file giving the one-body information should be provided. See man page for option syntax or type FCIGenerateSMA -h" << endl;
      return -1;
    }
    
    
  if (Manager.GetString("input-state") != 0)
    {
      NbrStates = 1;
      if (IsFile(Manager.GetString("input-state")) == false)
	{
	  cout << "can't open file " << Manager.GetString("input-state") << endl;
	  return -1;
	}
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites,
									   KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									   Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelFindSystemInfoFromVectorFileName(Manager.GetString("input-state"), NbrParticles, NbrSites, 
									   Statistics, GutzwillerFlag) == false)
	    {
	      
		cout << "error while retrieving system parameters from file name " << Manager.GetString("input-state") << endl;
		return -1;
	    }
	}
	else
	{
	  MomentumFlag = true;
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      NbrStates = DegenerateFile.GetNbrLines();
      if (FTIHubbardModelWith2DTranslationFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites,
									   KxMomentum, KyMomentum, XPeriodicity, YPeriodicity, 
									   Statistics, GutzwillerFlag) == false)
	{
	  if (FTIHubbardModelFindSystemInfoFromVectorFileName(DegenerateFile(0, 0), NbrParticles, NbrSites,
									       Statistics, GutzwillerFlag) == false)
	    {
		cout << "error while retrieving system parameters from file name " << DegenerateFile(0, 0) << endl;
		return -1;
	    }
	}
      else
	{
	  MomentumFlag = true;
	}
    }

      
  ComplexVector* InputStates = new ComplexVector[NbrStates];
  char** InputStateNames = new char*[NbrStates];
  if (Manager.GetString("input-state") != 0)
    {
      InputStateNames[0] = new char[strlen(Manager.GetString("input-state")) + 1];
      strcpy (InputStateNames[0], Manager.GetString("input-state"));
      if (InputStates[0].ReadVector (Manager.GetString("input-state")) == false)
	{
	  cout << "can't open vector file " << Manager.GetString("input-state") << endl;
	  return -1;      
	}
    }
  else
    {
      MultiColumnASCIIFile DegenerateFile;
      if (DegenerateFile.Parse(Manager.GetString("degenerate-states")) == false)
	{
	  DegenerateFile.DumpErrors(cout);
	  return -1;
	}
      if (InputStates[0].ReadVector (DegenerateFile(0, 0)) == false)
	{
	  cout << "can't open vector file " << DegenerateFile(0, 0) << endl;
	  return -1;      
	}	  
      InputStateNames[0] = new char[strlen(DegenerateFile(0, 0)) + 1];
      strcpy (InputStateNames[0], DegenerateFile(0, 0));
      for (int i = 1; i < NbrStates; ++i)
	{
	  InputStateNames[i] = new char[strlen(DegenerateFile(0, i)) + 1];
	  strcpy (InputStateNames[i], DegenerateFile(0, i));
	  if (InputStates[i].ReadVector (DegenerateFile(0, i)) == false)
	    {
	      cout << "can't open vector file " << DegenerateFile(0, i) << endl;
	      return -1;      
	    }	  
	  if (InputStates[0].GetVectorDimension() != InputStates[i].GetVectorDimension())
	    {
	      cout << "error, " << DegenerateFile(0, 0) << " and " <<  DegenerateFile(0, i) << "don't have the same  dimension (" << InputStates[0].GetVectorDimension() << " and " << InputStates[i].GetVectorDimension()<< ")" << endl;
	      return -1;
	    }
	}
    }
  ParticleOnSphere* InputSpace = 0;
  if (MomentumFlag == false)
    {
      if (Statistics == true)
	  InputSpace = new FermionOnLatticeRealSpace (NbrParticles, NbrSites);
      else
	{
	  if (GutzwillerFlag == false)
	    InputSpace = new BosonOnLatticeRealSpace (NbrParticles, NbrSites);
	  else
	    InputSpace = new BosonOnLatticeGutzwillerProjectionRealSpace (NbrParticles, NbrSites);
	}
    }
  else
    {
      if (Statistics == true)
	  InputSpace = new FermionOnLatticeRealSpaceAnd2DTranslation (NbrParticles, NbrSites, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
      else
      {
	if (GutzwillerFlag == false)
	  InputSpace = new BosonOnLatticeRealSpaceAnd2DTranslation (NbrParticles, NbrSites, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
	else
	  InputSpace = new BosonOnLatticeGutzwillerProjectionRealSpaceAnd2DTranslation (NbrParticles, NbrSites, KxMomentum, XPeriodicity, KyMomentum, YPeriodicity);
      }
    }


  if (InputSpace->GetHilbertSpaceDimension() != InputStates[0].GetVectorDimension())
    {
      cout << "error, " << Manager.GetString("input-state")  << " has a wrong dimension (" <<InputStates[0].GetVectorDimension() << ", should be " << InputSpace->GetHilbertSpaceDimension() << ")" << endl;
      return -1;
    }
  

  ComplexMatrix* RealSpaceDensityMatrixElements = new ComplexMatrix[NbrStates];
  for (int i = 0; i < NbrStates; ++i)
    RealSpaceDensityMatrixElements[i] = ComplexMatrix (NbrSites, NbrSites, true);
  Generic2DTightBindingModel TightBindingModel(Manager.GetString("import-onebody"));
  int TmpIndex1;
  int TmpIndex2;
  
  ParticleOnSphereDensityOperator* Operator;
  	      
  for (int pos1 = 0; pos1 < NbrSites; ++pos1)
  {
    for (int pos2 = 0; pos2 < NbrSites; ++pos2)
    {
      if (MomentumFlag == false)
	Operator = new ParticleOnSphereDensityOperator (InputSpace, pos1, pos2);
      else
	Operator = new ParticleOnLatticeRealSpaceAnd2DTranslationOneBodyOperator (InputSpace, pos1, pos2);
      
      for (int i = 0; i < NbrStates; ++i)
      {
	OperatorMatrixElementOperation Operation(Operator, InputStates[i], InputStates[i], InputStates[i].GetLargeVectorDimension());		Operation.ApplyOperation(Architecture.GetArchitecture());
	RealSpaceDensityMatrixElements[i].SetMatrixElement(pos1, pos2, Operation.GetScalar());
	if (pos1 == pos2)
	  RealSpaceDensityMatrixElements[i].SetMatrixElement(pos1, pos2, Operation.GetScalar());
	else
	  RealSpaceDensityMatrixElements[i].SetMatrixElement(pos1, pos2, (-Operation.GetScalar()));
      }
    }
  }
   
//   cout << RealSpaceDensityMatrixElements[0] << endl;
  
  int MinBandIndex = 0;
  int MaxBandIndex = TightBindingModel.GetNbrBands();
  if (Manager.GetInteger("band-index") > 0)
  {
    MinBandIndex = Manager.GetInteger("band-index");
    MaxBandIndex = MinBandIndex;
  }
  Complex* NbrParticlesPerBand = new Complex[NbrStates];
  Complex TmpDensity;
  int TmpMomentum;
  int pos1; 
  int posX1; 
  int posY1;
  int pos2;
  int posX2;
  int posY2;
  int index1;
  int index2;
  
  ComplexVector bra;
  ComplexVector ket;
  Complex inner;
  
  int NbrSiteX = TightBindingModel.GetNbrSiteX();
  int NbrSiteY = TightBindingModel.GetNbrSiteY();
  double Normalization = 1.0 / (TightBindingModel.GetNbrStatePerBand());
  
  Complex TotalNumber = 0.0;
  
  char* TmpOutputName = ReplaceExtensionToFileName(InputStateNames[0], "vec", "band_occupation.dat");
  ofstream File;
  File.precision(14);
  File.open(TmpOutputName, ios::binary | ios::out);
  File << "# bandIndex kx ky rho" << endl;
  
  for (int bandIndex = MinBandIndex; bandIndex < MaxBandIndex; ++bandIndex)
  {
    for (int i = 0; i < NbrStates; ++i)
      NbrParticlesPerBand[i] = 0.0;
    for (int momentumKx = 0; momentumKx < NbrSiteX; ++momentumKx)
    {
      for (int momentumKy = 0; momentumKy < NbrSiteY; ++momentumKy)
      {
  	    File << bandIndex << " " << momentumKx << " " << momentumKy << " " ;
	    TmpMomentum = TightBindingModel.GetLinearizedMomentumIndex(momentumKx, momentumKy);
	    
	    bra = TightBindingModel.GetOneBodyMatrix(TmpMomentum)[bandIndex];
	    ket = TightBindingModel.GetOneBodyMatrix(TmpMomentum)[bandIndex];
	    
	    for (int i = 0; i < NbrStates; ++i)
	      {
		TmpDensity = 0.0;
		
		for (int pos1 = 0; pos1 < NbrSites; ++pos1)
		{
		  for (int pos2 = 0; pos2 < NbrSites; ++pos2)
		  {		    
		    TightBindingModel.GetRealSpaceTightBindingLinearizedIndex(pos1, posX1, posY1, index1);
		    TightBindingModel.GetRealSpaceTightBindingLinearizedIndex(pos2, posX2, posY2, index2);
		    
		    	    
			    
// 		    TmpDensity += ((Phase(-2.0 * M_PI * ((double) ((posX1 - posX2) * momentumKx)) / ((double) NbrSiteX) 
// 				     - 2.0 * M_PI * ((double) ((posY1 - posY2) * momentumKy)) / ((double) NbrSiteY)) )
// 				    * (Conj(TightBindingModel.GetOneBodyMatrix(TmpMomentum)[bandIndex][index2]) * (TightBindingModel.GetOneBodyMatrix(TmpMomentum)[bandIndex][index1]))
// 			       	       * Normalization * RealSpaceDensityMatrixElements[i][pos2][pos1]);
		    
		    cout << posX1 << " " << posY1 << " " << index1 << " " << posX2 << " " << posY2 << " " << index2 << " " << ((Phase(-2.0 * M_PI * ((double) ((posX1 - posX2) * momentumKx)) / ((double) NbrSiteX)				     - 2.0 * M_PI * ((double) ((posY1 - posY2) * momentumKy)) / ((double) NbrSiteY))) * Normalization * Conj(ket[index1]) * (bra[index2]) * RealSpaceDensityMatrixElements[i][pos1][pos2]) << endl;
		    
		    TmpDensity += ((Phase(-2.0 * M_PI * ((double) ((posX1 - posX2) * momentumKx)) / ((double) NbrSiteX) 
				     - 2.0 * M_PI * ((double) ((posY1 - posY2) * momentumKy)) / ((double) NbrSiteY))) * Normalization * Conj(ket[index1]) * (bra[index2]) * RealSpaceDensityMatrixElements[i][pos1][pos2]);
		  }
		}
		File << TmpDensity.Re << " " ;
		NbrParticlesPerBand[i] += TmpDensity;
		TotalNumber += TmpDensity;
	      }
	      File << endl;
	      
	      }      
    }
    cout << "rho_" << bandIndex  << " = ";
    for (int i = 0; i < NbrStates; ++i)
      cout << NbrParticlesPerBand[i] << " " ;
    cout << endl;
  }
  cout << TotalNumber << endl;
  File.close();      
  
  delete[] NbrParticlesPerBand;
  delete[] RealSpaceDensityMatrixElements;
  delete InputSpace;
    
  return 0;
}

