#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"

#include "HilbertSpace/UndescribedHilbertSpace.h"

#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslations.h"
#include "HilbertSpace/VirtualSpaceTransferMatrixWithTranslations.h"
#include "HilbertSpace/Spin1_2ChainNew.h"
#include "HilbertSpace/Spin1_2ChainWithTranslations.h"
#include "HilbertSpace/DoubledSpin1_2_Chain.h"
#include "HilbertSpace/DoubledSpin1_2_ChainWithTranslations.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsStaggered.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry.h"
#include "HilbertSpace/DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry.h"

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
  OptionManager Manager ("EDTransfertMatrixMPOGivenByInputFile" , "0.01");
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
  (*SystemGroup) += new SingleStringOption  ('\n', "tensor-file", "name of the file containing the PEPS");
  (*SystemGroup) += new SingleStringOption  ('\n', "second-tensor", "name of the file containing the second PEPS (use the sublattice mode)");
  (*SystemGroup) += new SingleStringOption  ('\n', "peps-name", "name of the peps used to form the output file name");
  (*SystemGroup) += new BooleanOption ('\n', "no-spin", "use undescribed hilbert space ");

  (*SystemGroup) += new BooleanOption ('\n', "spin-half", "use spin-half virtual space");
  (*SystemGroup) += new BooleanOption ('c', "complex", "use complex version of the code");
  (*SystemGroup) += new BooleanOption ('\n', "doubled", "use double version of the code");
  (*SystemGroup) += new BooleanOption ('\n', "translation", "use translation symmetry");
  (*SystemGroup) += new BooleanOption ('\n', "zz-symmetry", "use the ZZ symmetry");
  (*SystemGroup) += new BooleanOption ('\n', "staggered", "use the stagerred Hilbert Space");
  (*SystemGroup) += new BooleanOption ('\n', "left", "compute left eigenvalue");
  (*SystemGroup) += new BooleanOption ('\n', "vison", "add a vison line in both layers");
  (*SystemGroup) += new  SingleIntegerOption ('\n', "sz", "consider a specific value of sz", -1);
  (*SystemGroup) += new  SingleIntegerOption ('d', "bond-dimension", "enter the bond dimension in mode no spin", -1);
  (*SystemGroup) += new  SingleIntegerOption ('\n', "k", "consider a specific value of k", -1);
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
      cout << "see man page for option syntax or type GenericOverlap -h" << endl;
      return -1;
    }

  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  int  TranslationStep  = Manager.GetInteger("translation-step");
  bool StaggeredFlag = Manager.GetBoolean("staggered");
  bool ComplexFlag = Manager.GetBoolean("complex");
  bool DoubledFlag = Manager.GetBoolean("doubled");
  bool TranslationFlag = Manager.GetBoolean("translation");
  bool SymmetryFlag = Manager.GetBoolean("zz-symmetry");
  bool UndescribedHilbertSpaceFlag =  Manager.GetBoolean("no-spin");
  bool LeftFlag = Manager.GetBoolean("left");
  MultiColumnASCIIFile TensorsElementsDefinition;
  MultiColumnASCIIFile TensorsElementsDefinitionB;
  int BondDimension = 3;
  if (TensorsElementsDefinition.Parse(Manager.GetString("tensor-file")) == false)
    {
      TensorsElementsDefinition.DumpErrors(cout) << endl;
      return -1;
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

  if (DoubledFlag)
    {
      if (TensorsElementsDefinition.GetNbrColumns() != 6)
	{
	  cout <<" The tensor file should have 6 columnns"<<endl;
	}
    }
  else
    {
      if (TensorsElementsDefinition.GetNbrColumns() != 5)
	{
	  cout <<" The tensor file should have 5 columnns"<<endl;
	}
    }

  bool FirstRunFlag = true;   
  AbstractHilbertSpace *  Space = 0;  
  int NbrSites = Manager.GetInteger("length");
  
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
  if (SublatticeMode) 
    {
      if(TranslationFlag == true)
	{
	  if(ComplexFlag == true)
	    {
	      if (DoubledFlag)
		{
		  TransferMatrix = new  ComplexPEPSABSublatticeTransfertMatrixPBCWithTranslations(TensorsElementsDefinition,TensorsElementsDefinition,Architecture.GetArchitecture());
		}
	      else
		{
		  cout <<"Work only for transfer matrix"<<endl;
		}
	    }
	  else
	    {
	        cout <<"Work only for complex PEPS"<<endl;
	    }
	}
      else
	{
	  cout <<"Work only with translation"<<endl;
	}
    }
  else
    {
      if(TranslationFlag == true)
	{
	  if(ComplexFlag == true)
	    {
	      if (DoubledFlag)
		{
		  TransferMatrix = new  ComplexPEPSTransfertMatrixPBCWithTranslations(TensorsElementsDefinition,Architecture.GetArchitecture());
		}
	      else
		{
		  cout <<" Complex transfert Matrix with translations not yet implemented"<<endl;
		  exit(-1);
		}
	    }
	  else
	    TransferMatrix = new TransfertMatrixPBCWithTranslationsFromFile(TensorsElementsDefinition,Architecture.GetArchitecture());
	}
      else
	{
	  if (ComplexFlag)
	    {
	      if (DoubledFlag)
		{
		  TransferMatrix = new ComplexPEPSTransfertMatrixPBC(TensorsElementsDefinition,&BoundaryConditions,Architecture.GetArchitecture()); 
		}
	      else
		{
		  cout <<" Complex transfert Matrix not yet implemented"<<endl;
		  exit(-1);
		}
	    }
	  else
	    {
	      if (DoubledFlag)
		{ 
		  TransferMatrix = new AbstractPEPSTransfertMatrixPBC(TensorsElementsDefinition,Architecture.GetArchitecture()); 
		}
	      else
		TransferMatrix = new AbstractTransfertMatrixPBC(TensorsElementsDefinition,Architecture.GetArchitecture());
	    }
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
  int SzMin = -NbrSites;
  int SzMax = NbrSites;
  if (DoubledFlag)
    {
      SzMax*=2;
      SzMin*=2;
    }
  char * SubspaceLegend = new char [50];
  char * TmpSzString = new char [30];
  if (Manager.GetInteger("sz") != -1 )
    {
      SzMin = Manager.GetInteger("sz");
      SzMax = SzMin;
    }
  int ZvalueMax = 0;
  if(SymmetryFlag == true)
    ZvalueMax = 1;
  
  if (TranslationFlag == false) 
    { 
      MaxKx=0;
    }
  
  if (UndescribedHilbertSpaceFlag ==true)
    {
      SzMin = 0;
      SzMax = 0;
      ZvalueMax = 0;
    }
  
  if (TranslationFlag)
    {
      if(SymmetryFlag)
	sprintf(SubspaceLegend,"Sz Kx ZBra ZKet"); 
      else
	sprintf(SubspaceLegend,"Sz Kx"); 
    }
  else
    {
      if(SymmetryFlag)
	sprintf(SubspaceLegend,"Sz ZBra ZKet"); 
      else
	sprintf(SubspaceLegend,"Sz");
    }
  
  char * FullOutputFileName = new char [200];
  sprintf(FullOutputFileName,"TransfertMatrix_%s_l_%d.dat",Manager.GetString("peps-name"),NbrSites); 
  int ZvalueKet = 0;
  for(int Sz = SzMin; Sz<= SzMax ;Sz+=1)
    {
      cout <<"Sz = "<<Sz<<endl;
      for (int i = MinKx; i <= MaxKx; ++i)
	{
	  cout <<" K = "<<i<<endl;
	  for(int ZvalueBra = 0 ; ZvalueBra <= ZvalueMax;ZvalueBra++)
	    {
	      
	      if (Sz %2 == 0)
		ZvalueKet =  ZvalueBra;
	      else
		{
		  ZvalueKet = 1 - ZvalueBra;
		}
	      
	      if (UndescribedHilbertSpaceFlag)
		{
int BondDimension = Manager.GetInteger("bond-dimension");
		  if (TranslationFlag) 
		    {
		      Space = new VirtualSpaceTransferMatrixWithTranslations(NbrSites, BondDimension,i,TranslationStep,100000,100000);
		    }
		  else
		    {
		      Space = new VirtualSpaceTransferMatrixWithTranslations(NbrSites, BondDimension,100000,100000);
		    }
		}
	      else	      
		{
		  if(Manager.GetBoolean("spin-half") == true)
		    {
		      if (DoubledFlag)
			{
			  if (TranslationFlag) 
			    {
			      Space = new DoubledSpin1_2_ChainWithTranslations (NbrSites,i,Sz,100000,100000);
			    }
			  else
			    Space =  new DoubledSpin1_2_Chain(NbrSites,Sz,100000,100000);
			}
		      else
			{
			  if (TranslationFlag) 
			    {
			      Space = new  Spin1_2ChainWithTranslations (NbrSites,i,1,Sz,1000000,1000000);
			    }
			  else
			    Space = new  Spin1_2ChainNew (NbrSites,Sz,100000);
			}
		    }
		  else
		    {
		      if (TranslationFlag) 
			{
			  if (SymmetryFlag)
			    {
			      if(StaggeredFlag)
				Space = new DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry (NbrSites,i,Sz, ZvalueBra, ZvalueKet,100000,100000);
			      else
				Space = new DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry (NbrSites,i,TranslationStep, Sz, ZvalueBra, ZvalueKet,100000,100000);
			    }
			  else
			    {
			      if(StaggeredFlag)
				Space = new DoubledSpin0_1_2_ChainWithTranslationsStaggered (NbrSites,i,Sz,100000,100000);
			      else
				Space = new DoubledSpin0_1_2_ChainWithTranslations (NbrSites,i,TranslationStep,Sz,100000,100000);
			    }
			}
		      else
			{
			  if (SymmetryFlag)
			    {
			      if(StaggeredFlag)
				Space = new DoubledSpin0_1_2_ChainWithTranslationsStaggeredAndZZSymmetry (NbrSites,Sz, ZvalueBra, ZvalueKet,100000,100000);
			      else
				Space = new DoubledSpin0_1_2_ChainWithTranslationsAndZZSymmetry (NbrSites,Sz, ZvalueBra, ZvalueKet,100000,100000);
			    }
			  else
			    {
			      if(StaggeredFlag)
				Space = new DoubledSpin0_1_2_ChainWithTranslationsStaggered (NbrSites,Sz, 100000,100000);
			      else
				{
				  cout <<" Space = new DoubledSpin0_1_2_ChainWithTranslations (NbrSites,Sz, 100000,100000);" <<endl;
				  Space = new DoubledSpin0_1_2_ChainWithTranslations (NbrSites,Sz, 100000,100000);
				}
			    }
			}
		    }
		}

	      if (SymmetryFlag)
		{
		  cout <<"Symmetry sector = "<< ZvalueBra<<" "<< ZvalueKet<<endl;
		}
	      TransferMatrix->SetHilbertSpace(Space);	  
	      
	      char * TmpEigenstateString = new char [200] ;
	      if ( (TranslationFlag) || (ComplexFlag) )
		{ 
		  
		  if (Space->GetHilbertSpaceDimension() > 0 ) 
		    {
		      cout <<"Hilbert Space dimension = "<<Space->GetHilbertSpaceDimension()<<endl;
		      if (UndescribedHilbertSpaceFlag)
			{
			  if (TranslationFlag)
			    {
			      sprintf(TmpSzString,"%d",i);			  
			      sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_d_%d_k_%d",Manager.GetString("peps-name"),NbrSites,BondDimension,i);
			    }
			  else
			    {
			      sprintf(TmpSzString,"");			  
			      sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_d_%d",Manager.GetString("peps-name"),BondDimension,NbrSites);
			    }
			}
		      else
			{
			  if (TranslationFlag)
			    {
			      if(SymmetryFlag)
				{
				  sprintf(TmpSzString,"%d %d %d %d",Sz,i, ZvalueBra, ZvalueKet);
				  sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_sz_%d_k_%d_zbra_%d_zket_%d",Manager.GetString("peps-name"),NbrSites,Sz,i, ZvalueBra, ZvalueKet);
				}
			      else
				{
 				  sprintf(TmpSzString,"%d %d",Sz,i);
				  sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_sz_%d_k_%d",Manager.GetString("peps-name"),NbrSites,Sz,i);
				}
			    }
			  else
			    {
			      if(SymmetryFlag)
				{
				  sprintf(TmpSzString,"%d %d %d",Sz, ZvalueBra, ZvalueKet);
				  sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_sz_%d_zbra_%d_zket_%d",Manager.GetString("peps-name"),NbrSites,Sz, ZvalueBra, ZvalueKet);
				}
			      else
				{
				  sprintf(TmpSzString,"%d",Sz);			  
				  sprintf(TmpEigenstateString,"TransfertMatrix_%s_l_%d_sz_%d",Manager.GetString("peps-name"),NbrSites,Sz);
				}
			    }
			}
		      Lanczos.SetComplexAlgorithms();
		      
		      
		      /*			  GenericComplexMainTask Task(&Manager, Space, &Lanczos, TransferMatrix, TmpSzString, SubspaceLegend, 0.0,  FullOutputFileName, FirstRunFlag, TmpEigenstateString);
			FirstRunFlag = false;
			MainTaskOperation TaskOperation (&Task);
			TaskOperation.ApplyOperation(Architecture.GetArchitecture());*/
		      
		      int NbrEigenvalues = Manager.GetInteger("nbr-eigen" );
		      GenericNonHermitianMainTask Task (&Manager,  TransferMatrix, NbrEigenvalues, Manager.GetBoolean("eigenstate"), LeftFlag, 1e-12, TmpSzString, SubspaceLegend,0.0,FirstRunFlag, FullOutputFileName,TmpEigenstateString);
		      FirstRunFlag = false;
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      
		      
		      cout << "------------------------------------" << endl;
		      //		}
		    }
		}
	      else
		{
		  RealVector TestVector(Space->GetHilbertSpaceDimension(),true);
		  RealVector DestinationVector(Space->GetHilbertSpaceDimension(),true);
		  for (int p = 0; p < Space->GetHilbertSpaceDimension(); p++)
		    {
		      TestVector[p] = (rand() - 32767) * 0.5;
		    }
		  TestVector /=  TestVector.Norm();
		  
		  TransferMatrix->Multiply(TestVector,DestinationVector);
		  if ( fabs(DestinationVector.Norm()) > 1e-14 )
		    { 
		      
		      sprintf(TmpSzString,"%d",Sz);
		      GenericRealMainTask Task(&Manager, Space, &Lanczos, TransferMatrix, TmpSzString, SubspaceLegend, 0.0,  FullOutputFileName, FirstRunFlag, TmpEigenstateString);
		      MainTaskOperation TaskOperation (&Task);
		      TaskOperation.ApplyOperation(Architecture.GetArchitecture());
		      FirstRunFlag = false;
		    }
		}
	      delete Space;
	      delete [] TmpEigenstateString;
	    }
	}
    }
  
  delete TransferMatrix;
  delete []  TmpSzString;
  delete [] SubspaceLegend;
  delete [] FullOutputFileName;
  return 0;
}
