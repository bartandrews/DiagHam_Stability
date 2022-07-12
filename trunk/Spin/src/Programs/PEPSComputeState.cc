#include "Options/Options.h"

#include "GeneralTools/MultiColumnASCIIFile.h"
#include "HilbertSpace/Spin0_1_2_ChainWithTranslations.h"


#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/ComplexMatrix.h"

#include "Operator/SpinS2Operator.h"



#include "MPSObjects/ComplexPEPSPBC.h"
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
  OptionManager Manager ("PEPSComputeState" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  
  Manager += SystemGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\n', "tensor-file", "name of the file containing the eigenstate to be displayed");
  (*SystemGroup) += new SingleStringOption  ('\n', "peps-name", "name of the peps used to form the output file name");
  (*SystemGroup) += new  SingleIntegerOption ('x',"Lx", "size of the lattice in x direction", 4);
  (*SystemGroup) += new  SingleIntegerOption ('y',"Ly", "size of the lattice in y direction (only 2 or 4)", 2);
  (*SystemGroup) += new  SingleIntegerOption ('sz',"sz", "sz value", -10000);
  (*SystemGroup) += new BooleanOption ('\n', "block-tensor", "compute the 2*2 blocked tensor (instead of computing the state)");

  (*SystemGroup) += new BooleanOption ('\n', "cylinder", "compute the state on a cylinder");
  (*SystemGroup) += new SingleStringOption  ('\n', "left-vector", "name of the file containing the left boundary vector (useful only in cylinder mode)");
  (*SystemGroup) += new SingleStringOption  ('\n', "right-vector", "name of the file containing the right boundary vector (useful only in cylinder mode)");

  (*SystemGroup) += new BooleanOption ('\n', "hstring", "add a horizontal string");
  (*SystemGroup) += new BooleanOption ('\n', "vstring", "add a vertical  string");
  (*SystemGroup) += new SingleDoubleOption ('\n',"theta", "angle in the U(1) inserted on the string", 0);

  (*SystemGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type PEPSComputeState -h" << endl;
      return -1;
    }

  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  int Lx = Manager.GetInteger("Lx");
  int Ly = Manager.GetInteger("Ly");
  int Sz = Manager.GetInteger("sz");
  bool BlockFlag = Manager.GetBoolean("block-tensor");
  bool CylinderFlag = Manager.GetBoolean("cylinder");

  bool HorinzontalStringFlag = Manager.GetBoolean("hstring");
  bool VerticalStringFlag = Manager.GetBoolean("vstring");
  
  MultiColumnASCIIFile TensorsElementsDefinition;
  if (TensorsElementsDefinition.Parse(Manager.GetString("tensor-file")) == false)
    {
      TensorsElementsDefinition.DumpErrors(cout) << endl;
      return -1;
    } 
  
  ComplexPEPSPBC PEPS(TensorsElementsDefinition);
  
  if (BlockFlag) 
    {
      PEPS.ComputeBlockTensor();
      return 0;
    }


  ComplexDiagonalMatrix MatrixOnString(PEPS.GetBondDimension(),true);
  if ( (HorinzontalStringFlag)|| (VerticalStringFlag ))
    {
      double Theta = Manager.GetDouble("theta");
      MatrixOnString.SetMatrixElement(0,0,Phase(2.0*M_PI*Theta));
      MatrixOnString.SetMatrixElement(1,1,Phase(2.0*M_PI*Theta));
      MatrixOnString.SetMatrixElement(2,2,Phase(-2.0*M_PI*3.0*Theta));
    }
    

  char * FullOutputFileName = new char [200];
  sprintf(FullOutputFileName,"PEPS_%s_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Lx,Ly); 

  if (CylinderFlag ==false)
    {
      if (Sz ==  -10000) 
	{
	  sprintf(FullOutputFileName,"PEPS_%s_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Lx,Ly); 
	  ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPS (Lx,Ly/2, MatrixOnString , MatrixOnString, HorinzontalStringFlag, VerticalStringFlag);
  cout <<State<<endl;

	}
      else
	{
	  sprintf(FullOutputFileName,"PEPS_%s_sz_%d_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Sz,Lx,Ly); 
	  ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPSSzConstraint (Lx,Ly/2,Sz,MatrixOnString , MatrixOnString, HorinzontalStringFlag, VerticalStringFlag);
	  State.WriteVector(FullOutputFileName);
  cout <<State<<endl;

	}
    }
  else
    {

      /*      Spin0_1_2_ChainWithTranslations Space  (Ly,100000,100000);
	      int NbrStates = Space.GetHilbertSpaceDimension();
      HermitianMatrix S2Matrix(NbrStates, true);
      SpinS2Operator TmpOperator(&Space, Ly);
      TmpOperator.GetOperator(S2Matrix) ;
      cout <<S2Matrix<<endl;
      RealDiagonalMatrix TmpS2Eigenvalues(NbrStates);
      ComplexMatrix TmpBasis (NbrStates, NbrStates);
      TmpBasis.SetToIdentity();
      S2Matrix.LapackDiagonalize(TmpS2Eigenvalues, TmpBasis);
      int NbrS0State = 0;
      for (int i = 0; i < NbrStates; ++i)
	{
	  double TmpS2 = TmpS2Eigenvalues[i];
	  cout << "<S^2>=" << TmpS2 << " <S>=" << (0.5 * (sqrt((4.0 * TmpS2) + 1.0) - 1.0)) << endl;
	  cout << "round(<2S>)=" <<  round(sqrt((4.0 * TmpS2) + 1.0) - 1.0) << endl; 
	  if ( fabs(TmpS2) < 1e-7  ) 
	    NbrS0State++;
	}
      cout << " Number of S = 0 boundary states : " <<  NbrS0State <<endl;

    for (int i = 0; i <  NbrS0State; i++)
	{
	  cout <<"i = "<< i <<endl;
	  cout << TmpBasis[i]<<endl;
	  
	  if (Sz ==  -10000) 
	    {
	      sprintf(FullOutputFileName,"PEPS_cylinder_%s_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Lx,Ly); 
	      ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPS (Lx,Ly/2, MatrixOnString , MatrixOnString,  TmpBasis[i],  TmpBasis[i], HorinzontalStringFlag, VerticalStringFlag);
	      State.WriteVector(FullOutputFileName);
	      cout <<State<<endl;
	    }
	  else
	    {
	      sprintf(FullOutputFileName,"PEPS_cylinder_%s_sz_%d_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Sz,Lx,Ly); 
	      ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPSSzConstraint (Lx, Ly/2,Sz, MatrixOnString , MatrixOnString,  TmpBasis[i],  TmpBasis[i], HorinzontalStringFlag, VerticalStringFlag);
	      State.WriteVector(FullOutputFileName);
	      cout <<State<<endl;
	    }
	}
    */
    
    ComplexVector LeftVector;
    if (LeftVector.ReadVector(Manager.GetString("left-vector")) == false)
      {
	cout << "error while reading " << Manager.GetString("left-vector") << endl;
	return -1;
      } 
    
    ComplexVector RightVector;
    if (RightVector.ReadVector(Manager.GetString("right-vector")) == false)
      {
	cout << "error while reading " << Manager.GetString("right-vector") << endl;
	return -1;
      } 
    
    if (Sz ==  -10000) 
      {
	sprintf(FullOutputFileName,"PEPS_cylinder_%s_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Lx,Ly); 
	ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPS (Lx,Ly/2, MatrixOnString , MatrixOnString, LeftVector, RightVector, HorinzontalStringFlag, VerticalStringFlag);
	State.WriteVector(FullOutputFileName);
	cout <<State<<endl;
      }
    else
      {
	sprintf(FullOutputFileName,"PEPS_cylinder_%s_sz_%d_lx_%d_ly_%d.vec",Manager.GetString("peps-name"),Sz,Lx,Ly); 
	ComplexVector State = PEPS.ComputeFockSpaceRepresentationOfAPEPSSzConstraint (Lx, Ly/2,Sz, MatrixOnString , MatrixOnString, LeftVector, RightVector, HorinzontalStringFlag, VerticalStringFlag);
	State.WriteVector(FullOutputFileName);
	cout <<State<<endl;
      }
    }
}
