#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"


#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers.h"


#include "MPSObjects/AbstractTransfertMatrixPBC.h"
#include "MPSObjects/TransfertMatrixPBCWithTranslationsFromFile.h"
#include "MPSObjects/AbstractPEPSTransfertMatrixPBC.h"
#include "MPSObjects/ComplexPEPSTransfertMatrixPBC.h"
#include "MPSObjects/ComplexPEPSTransfertMatrixPBCWithTranslations.h"
#include "MPSObjects/ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations.h"

#include "MainTask/GenericNonSymmetricMainTask.h"
#include "MainTask/GenericNonHermitianMainTask.h"
#include "MainTask/GenericRealMainTask.h"
#include "MainTask/GenericComplexMainTask.h"

#include "Matrix/RealDiagonalMatrix.h"

#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "LanczosAlgorithm/LanczosManager.h"

#include "LanczosAlgorithm/BasicArnoldiAlgorithm.h"
#include "LanczosAlgorithm/BasicArnoldiAlgorithmWithDiskStorage.h"


#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>

using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);
  
  // some running options and help
  OptionManager Manager ("EDTransfertMatrixSU2InvariantWithSublatticeQuantumNumbersPEPS" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  OptionGroup* ArnoldiGroup  = new OptionGroup ("Arnoldi options");
  LanczosManager Lanczos(false);
  
  Manager += ArnoldiGroup;
  ArchitectureManager Architecture;
  Lanczos.AddOptionGroup(&Manager);
  Architecture.AddOptionGroup(&Manager);
  Manager += SystemGroup;
  Manager += MiscGroup;
  Manager += ToolsGroup;

  (*SystemGroup) += new  SingleIntegerOption ('L', "length", "length of the spin chain", 4);
  (*SystemGroup) += new SingleStringOption  ('\n', "tensor-file", "name of the file containing the eigenstate to be displayed");
  (*SystemGroup) += new SingleStringOption  ('\n', "second-tensor", "name of the file containing the second PEPS (use the sublattice mode)");
  (*SystemGroup) += new SingleStringOption  ('\n', "peps-name", "name of the peps used to form the output file name");
  (*SystemGroup) += new BooleanOption ('\n', "translation", "use translation symmetry");
  (*SystemGroup) += new BooleanOption ('\n', "left", "compute left eigenvalue");
  (*SystemGroup) += new BooleanOption ('\n', "vison", "add a vison line in both layers");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sz", "consider a specific value of sz", -1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "k", "consider a specific value of k", -1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sublatticezerobra","consider a specific value of DeltaNo Bra" , -1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sublatticezeroket","consider a specific value of DeltaNo Ket" , -1);

  (*SystemGroup) += new  SingleIntegerOption ('\n', "translation-step", "", 1);
  
#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif
  (*ToolsGroup) += new BooleanOption  ('\n', "show-hamiltonian", "show matrix representation of the hamiltonian");
  (*ToolsGroup) += new BooleanOption  ('\n', "test-hermitian", "test if the hamiltonian is hermitian");

  (*MiscGroup) += new  BooleanOption ('\n', "print-tensor", "print the tensor elements", false);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type EDTransfertMatrixSU2InvariantWithSublatticeQuantumNumbersPEPS -h" << endl;
      return -1;
    }
  
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  bool TranslationFlag = Manager.GetBoolean("translation");
  bool LeftFlag = Manager.GetBoolean("left");
  MultiColumnASCIIFile TensorsElementsDefinition;
  MultiColumnASCIIFile TensorsElementsDefinitionB;
  
  if (TensorsElementsDefinition.Parse(Manager.GetString("tensor-file")) == false)
    {
      TensorsElementsDefinition.DumpErrors(cout) << endl;
      return -1;
    } 
  
  if (TensorsElementsDefinition.GetNbrColumns() != 6)
    {
      cout <<" The tensor file should have 6 columnns"<<endl;
    }
  

  bool SublatticeMode = false;
  if (Manager.GetString("second-tensor")!= 0 ) 
    {
      SublatticeMode = true;
      if (TensorsElementsDefinitionB.Parse(Manager.GetString("second-tensor")) == false)
	{
	  TensorsElementsDefinitionB.DumpErrors(cout) << endl;
	  return -1;
	} 
    }
    
  bool FirstRunFlag = true;   
  AbstractHilbertSpace *  Space = 0;  
  //  AbstractTransfertMatrixPBC * TransferMatrix = 0;
  int NbrSites = Manager.GetInteger("length");
  int  TranslationStep  = Manager.GetInteger("translation-step");  

  RealDiagonalMatrix BoundaryConditions(9,true);
  BoundaryConditions.SetToIdentity();
  if(Manager.GetBoolean("vison") == true)
    {
      for(int i =0;i <9; i++)
	{
	  double  Tmp[3];
	  Tmp[0] = -1.0;Tmp[1] = -1.0;Tmp[2] = 1.0;
	  BoundaryConditions.SetMatrixElement(i,i,Tmp[i%3] *Tmp[i/3] ); 
      	}
    }

  AbstractTransfertMatrixPBC *  TransferMatrix =0;
  
  if (SublatticeMode == true)
    {
      TransferMatrix = new  ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations(TensorsElementsDefinition,TensorsElementsDefinitionB,Architecture.GetArchitecture());
    }
  else
    {  
      if(TranslationFlag == true)
	{
	  TransferMatrix = new  ComplexPEPSTransfertMatrixPBCWithTranslations(TensorsElementsDefinition,Architecture.GetArchitecture());
	}
      else
	{
	  TransferMatrix = new ComplexPEPSTransfertMatrixPBC(TensorsElementsDefinition,&BoundaryConditions,Architecture.GetArchitecture()); 
	}
    }
	  
  
  if(Manager.GetBoolean("print-tensor") == true)
    {
      TransferMatrix->PrintTensorElements();
      return 0;
    }            
  int MinKx = 0;
  int MaxKx = NbrSites / TranslationStep  - 1;
  
  if (Manager.GetInteger("k") != -1 )
    {
      MinKx = Manager.GetInteger("k");
      MaxKx = MinKx;
    }
    
  int SzMin = -2*NbrSites;
  int SzMax = 2*NbrSites;

  char * SubspaceLegend = new char [200];
  char * TmpSzString = new char [200];
  if (Manager.GetInteger("sz") != -1 )
    {
      SzMin = Manager.GetInteger("sz");
      SzMax = SzMin;
    }

  
  if (TranslationFlag == false) 
    { 
      MaxKx=0;
    }
  
  if (TranslationFlag)
    {
      sprintf(SubspaceLegend,"Sz Kx SubLatticeZeroBra  SubLatticeZeroKet  SubLatticeZeroProduct"); 
    }
  else
    {
      sprintf(SubspaceLegend,"Sz SubLatticeZeroBra  SubLatticeZeroKet  SubLatticeZeroProduct"); 
    }
  
  char * FullOutputFileName = new char [200];
  sprintf(FullOutputFileName,"TransfertMatrix_%s_l_%d.dat",Manager.GetString("peps-name"),NbrSites); 

  int SubLatticeZeroBraMax=NbrSites/2;
  int SubLatticeZeroKetMax=NbrSites/2;  
  int SubLatticeZeroBraMin = 0;
  int SubLatticeZeroKetMin = 0;
  int SubLatticeZeroProductIncrement = 2;
  if (Manager.GetInteger("sublatticezerobra") != -1 )
    {
      SubLatticeZeroBraMax =Manager.GetInteger("sublatticezerobra");
      SubLatticeZeroBraMin =Manager.GetInteger("sublatticezerobra");
    }
  
  if (Manager.GetInteger("sublatticezeroket") != -1 )
    {
      SubLatticeZeroKetMax =Manager.GetInteger("sublatticezeroket");
      SubLatticeZeroKetMin =Manager.GetInteger("sublatticezeroket");
    }

  if (SublatticeMode == true)
    {
      SubLatticeZeroBraMin=-SubLatticeZeroKetMax;
      SubLatticeZeroKetMin=-SubLatticeZeroKetMax;
      SubLatticeZeroProductIncrement=3;
    }
  

  for(int  SubLatticeZeroBra =  SubLatticeZeroBraMin; SubLatticeZeroBra<= SubLatticeZeroBraMax; SubLatticeZeroBra+=1)
    {
      for(int SubLatticeZeroKet = SubLatticeZeroKetMin ; SubLatticeZeroKet <= SubLatticeZeroKetMax; SubLatticeZeroKet+=1)
	{
	  SzMax = 2*NbrSites - abs(SubLatticeZeroBra)-abs(SubLatticeZeroKet);
	  SzMin = - SzMax;
	  if (Manager.GetInteger("sz") != -1 )
	    {
	      if( ( Manager.GetInteger("sz") > SzMax) || ( Manager.GetInteger("sz") <  SzMin ))
		{
		  SzMin = 2;
		  SzMax = 0;
		}
	      else
		{
		  SzMin =  Manager.GetInteger("sz");
		  SzMax = SzMin;
		}
	    }
	  for(int Sz = SzMin; Sz<= SzMax ;Sz+=2)
	    {
	      cout <<"Sz = "<<Sz<<endl;
	      for (int i = MinKx; i <= MaxKx; ++i)
		{
		  cout <<" K = "<<i<<endl;
		  
		  
		  for(int SubLatticeZeroProduct = -1; SubLatticeZeroProduct<=1; SubLatticeZeroProduct+=SubLatticeZeroProductIncrement)
		    {
		      if  ( SubLatticeZeroBra  * SubLatticeZeroKet == 0 )
			SubLatticeZeroProduct = 0;
		      if (SublatticeMode == true)
			{
			  Space = new  DoubledSpin0_1_2_ChainWithSquareTranslationsAndSublatticeQuantumNumbers(NbrSites,i, Sz,SubLatticeZeroKet, SubLatticeZeroBra, 100000,100000);
			  cout << SubLatticeZeroBra << " " <<SubLatticeZeroKet<<endl;
			}
		      else
			{
			  if (TranslationFlag) 
			    {
			      Space = new DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers (NbrSites,i,TranslationStep ,Sz,SubLatticeZeroKet*SubLatticeZeroKet,SubLatticeZeroBra*SubLatticeZeroBra, SubLatticeZeroProduct*SubLatticeZeroKet*SubLatticeZeroBra, 100000,100000);
			    }
			  else
			    Space = new  DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetryAndSublatticeQuantumNumbers (NbrSites,Sz,SubLatticeZeroKet*SubLatticeZeroKet,SubLatticeZeroBra*SubLatticeZeroBra, SubLatticeZeroProduct*SubLatticeZeroKet*SubLatticeZeroBra, 100000,100000); 
			  cout << SubLatticeZeroBra << " " <<SubLatticeZeroKet<<" "<<SubLatticeZeroProduct<<endl;
			}
		      
		      TransferMatrix->SetHilbertSpace(Space);	  
		      
		      char * TmpEigenstateString = new char [200] ;
		      
		      if (Space->GetHilbertSpaceDimension() > 0 ) 
			{
			  cout <<"Hilbert Space dimension = "<<Space->GetHilbertSpaceDimension()<<endl;
		      if (SublatticeMode == true)
			{
			  sprintf(TmpSzString,"%d %d %d %d",Sz,i,SubLatticeZeroBra, SubLatticeZeroKet);
			  sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_sz_%d_k_%d_sbra_%d_sket_%d",Manager.GetString("peps-name"),NbrSites,Sz,i,SubLatticeZeroBra, SubLatticeZeroKet);
			}
		      else
			{			  
			  if (TranslationFlag)
			    {
			      sprintf(TmpSzString,"%d %d %d %d %d",Sz,i,SubLatticeZeroBra, SubLatticeZeroKet,SubLatticeZeroProduct);
			      sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_sz_%d_k_%d_sbra_%d_sket_%d_sprod_%d",Manager.GetString("peps-name"),NbrSites,Sz,i,SubLatticeZeroBra, SubLatticeZeroKet,SubLatticeZeroProduct);
			    }
			  else
			    {
			      sprintf(TmpSzString,"%d %d %d %d",Sz, SubLatticeZeroBra, SubLatticeZeroKet,SubLatticeZeroProduct);
			      sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_sz_%d_sbra_%d_sket_%d_sprod_%d",Manager.GetString("peps-name"),NbrSites,Sz, SubLatticeZeroBra, SubLatticeZeroKet,SubLatticeZeroProduct);
			    }
			}
			  Lanczos.SetComplexAlgorithms();
			  
			  int NbrEigenvalues = Manager.GetInteger("nbr-eigen" );
			  GenericNonHermitianMainTask Task (&Manager,  TransferMatrix, NbrEigenvalues, Manager.GetBoolean("eigenstate"), LeftFlag, 1e-12, TmpSzString, SubspaceLegend,0.0,FirstRunFlag, FullOutputFileName,TmpEigenstateString);
			  FirstRunFlag = false;
			  MainTaskOperation TaskOperation (&Task);
			  TaskOperation.ApplyOperation(Architecture.GetArchitecture());
			  
			  cout << "------------------------------------" << endl;
			}
		      delete Space;
		      delete [] TmpEigenstateString;
		    }
		}
	    }
	}
    }
  
  delete TransferMatrix;
  delete []  TmpSzString;
  delete [] SubspaceLegend;
  delete [] FullOutputFileName;
  return 0;
}
