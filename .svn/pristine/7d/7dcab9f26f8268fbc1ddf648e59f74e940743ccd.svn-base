#include "Vector/RealVector.h"

#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneBasisShort.h"
#include "HilbertSpace/BosonOnSphereHaldaneHugeBasisShort.h"
#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereHaldaneBasis.h"

#include "Options/Options.h"

#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"
#include "Tools/FQHEFiles/FQHESqueezedBasisTools.h"

#include "GeneralTools/MultiColumnASCIIFile.h"

#include "MathTools/WignerSmallDMatrix.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>
#include <string.h>


using std::cout;
using std::endl;
using std::ios;
using std::ofstream;


int main(int argc, char** argv)
{
  OptionManager Manager ("FQHESphereODLRORotation" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");
  OptionGroup* ToolsGroup  = new OptionGroup ("tools options");
  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += ToolsGroup;
  Manager += MiscGroup;
  (*SystemGroup) += new SingleStringOption  ('r', "right-state", "file describing right state decomposition into L^2 eigenstates");
  (*SystemGroup) += new SingleStringOption  ('l', "left-state", "file describing left state decomposition into L^2 eigenstates");
  (*SystemGroup) += new BooleanOption ('\n', "matrix-odlro", "compute matrix ODLRO");
  (*SystemGroup) += new SingleStringOption  ('m', "matrix-states", "file describing each state decompositions");
  (*SystemGroup) += new BooleanOption ('\n', "matrix-overlap", "indicates that an additional scalar product has to be taken into account when computing matrix ODLRO");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "matrix-eigenstates", "compute the n first eigenstates of the ODLRO matrix", 0);
  (*OutputGroup) += new SingleIntegerOption ('z', "lz-value", "twice the Lz value of the left state", 0);
  (*OutputGroup) += new SingleIntegerOption ('n', "nbr-points", "number of points that has to be computed", 100);

#ifdef __LAPACK__
  (*ToolsGroup) += new BooleanOption  ('\n', "use-lapack", "use LAPACK libraries instead of DiagHam libraries");
#endif

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereODLRORotation -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  int LzValue = Manager.GetInteger("lz-value");

  if (Manager.GetBoolean("matrix-odlro") == false)
    {
      MultiColumnASCIIFile LeftStateDecomposition;
      if (LeftStateDecomposition.Parse(Manager.GetString("left-state")) == false)
	{
	  LeftStateDecomposition.DumpErrors(cout);
	  return -1;
	}
      MultiColumnASCIIFile RightStateDecomposition;
      if (RightStateDecomposition.Parse(Manager.GetString("right-state")) == false)
	{
	  RightStateDecomposition.DumpErrors(cout);
	  return -1;
	}
      
      double* LeftStateComponents = LeftStateDecomposition.GetAsDoubleArray(1);
      double* RightStateComponents = RightStateDecomposition.GetAsDoubleArray(1);
      int* LeftStateLValues = LeftStateDecomposition.GetAsIntegerArray(2);
      int* RightStateLValues = RightStateDecomposition.GetAsIntegerArray(2);
      int NbrLValues = 0;
      for (int i = 0; i < LeftStateDecomposition.GetNbrLines(); ++i)
	{
	  int CurrentLValue = LeftStateLValues[i];	  
	  int j = 0;
	  for (; (j < RightStateDecomposition.GetNbrLines()) && (RightStateLValues[j] != CurrentLValue); ++j);
	  if (j < RightStateDecomposition.GetNbrLines())
	    ++NbrLValues;	
	}
      
      int* LValues = new int [NbrLValues];
      double* Coefficients = new double [NbrLValues];
      WignerSmallDMatrix** WignerCoefficients = new WignerSmallDMatrix*[NbrLValues];
      NbrLValues = 0;
      for (int i = 0; i < LeftStateDecomposition.GetNbrLines(); ++i)
	{
	  int CurrentLValue = LeftStateLValues[i];
	  int j = 0;
	  for (; (j < RightStateDecomposition.GetNbrLines()) && (RightStateLValues[j] != CurrentLValue); ++j);
	  if (j < RightStateDecomposition.GetNbrLines())
	    {
	      LValues[NbrLValues] = CurrentLValue;
	      WignerCoefficients[NbrLValues] = new WignerSmallDMatrix(CurrentLValue);
	      RealVector LeftVector;
	      if (LeftVector.ReadVector(LeftStateDecomposition(0, i)) == false)
		{
		  cout << "can't open vector file " << LeftStateDecomposition(0, i) << endl;
		  return -1;      	      
		}
	      RealVector RightVector;
	      if (RightVector.ReadVector(RightStateDecomposition(0, i)) == false)
		{
		  cout << "can't open vector file " << RightStateDecomposition(0, i) << endl;
		  return -1;      	      
		}
	      Coefficients[NbrLValues] = (LeftVector * RightVector) * (LeftStateComponents[i] * (RightStateComponents[i]));
	      ++NbrLValues;	
	    }
	}
      
      cout.precision(14);
      int NbrPoints = Manager.GetInteger("nbr-points");
      double Theta = 0.0;
      double ThetaInc = M_PI / ((double) NbrPoints);
      for (int i = 0; i <= NbrPoints; ++i)
	{
	  double Tmp = 0.0;
	  for (int j = 0; j < NbrLValues; ++j)
	    {
	      double Tmp2 = Coefficients[j] * (*(WignerCoefficients[j]))(LzValue, LzValue, Theta);// * (*(WignerCoefficients[j]))(LzValue, LzValue, Theta);
	      Tmp += Tmp2;
	    }      
	  cout << Theta << " " << Tmp << endl;
	  Theta += ThetaInc;
	}
    }
  else
    {
      MultiColumnASCIIFile States;
      if (States.Parse(Manager.GetString("matrix-states")) == false)
	{
	  States.DumpErrors(cout);
	  return -1;
	}
      int NbrStates = States.GetNbrLines();
      int NbrMaxLValues = 0;
      int TotalNbrComponents = 0;
      double** StateComponents = new double* [NbrStates];
      int** StateLValues = new int* [NbrStates];
      char*** StateFileNames = new char**[NbrStates];
      int* NbrStateComponents = new int [NbrStates];      
      int LShift = 0;
      RealSymmetricMatrix OverlapMatrix(NbrStates, true);
      if (Manager.GetBoolean("matrix-overlap") == true)
	{
	  if (States.GetNbrColumns() < 2)
	    {
	      cout << Manager.GetString("matrix-states") << " does not contain a second column with the original ODLRO states" << endl;
	      return -1;
	    }
	  RealVector TmpVector1;
	  if (TmpVector1.ReadVector(States(1, 0)) == false)
	    {
	      cout << "can't open " << States(1, 0) << endl;
	      return -1;
	    }
	  double Norm1 = TmpVector1.SqrNorm();
	  for (int i = 0; i < NbrStates; ++i)
	    {
	      if (TmpVector1.ReadVector(States(1, i)) == false)
		{
		  cout << "can't open " << States(1, i) << endl;
		  return -1;
		}
//	      double Norm1 = TmpVector1.Norm();
//	      OverlapMatrix(i, i) = 1.0;	      
	      OverlapMatrix(i, i) = TmpVector1.SqrNorm();// / Norm1;	      
	      for (int j = i + 1; j < NbrStates; ++j)
		{
		  RealVector TmpVector2;
		  if (TmpVector2.ReadVector(States(1, j)) == false)
		    {
		      cout << "can't open " << States(1, j) << endl;
		      return -1;
		    }
//		  double Norm2 = TmpVector2.Norm();
		  OverlapMatrix(i, j) = (TmpVector1 * TmpVector2);// / Norm1;	      
//		  OverlapMatrix(i, j) = (TmpVector1 * TmpVector2) / (Norm1 * Norm2);	      
		}
	    }
	}
      else
	{
	  for (int i = 0; i < NbrStates; ++i)
	    for (int j = i; j < NbrStates; ++j)
	      OverlapMatrix(i, j) = 1.0;
	}
//      cout << "overlap matrix : " << endl << OverlapMatrix << endl;
     
      for (int i = 0; i < NbrStates; ++i)
	{
	  MultiColumnASCIIFile TmpStateDecomposition;
	  if (TmpStateDecomposition.Parse(States(0, i)) == false)
	    {
	      TmpStateDecomposition.DumpErrors(cout);
	      return -1;
	    }
	  StateComponents[i] = TmpStateDecomposition.GetAsDoubleArray(1);
	  StateLValues[i] = TmpStateDecomposition.GetAsIntegerArray(2);
	  StateFileNames[i] = TmpStateDecomposition.GetAsStringArray(0);
	  NbrStateComponents[i] = TmpStateDecomposition.GetNbrLines();
	  if (NbrMaxLValues < StateLValues[i][NbrStateComponents[i] - 1])
	    NbrMaxLValues = StateLValues[i][NbrStateComponents[i] - 1];
	  LShift = StateLValues[i][NbrStateComponents[i] - 1] & 1;
	  TotalNbrComponents += NbrStateComponents[i];
	}

      NbrMaxLValues >>= 1;
      int* NbrStatesPerLValue = new int[NbrMaxLValues + 1];
      for (int i = 0; i <= NbrMaxLValues; ++i)
	NbrStatesPerLValue[i] = 0;
      
      for (int i = 0; i < NbrStates; ++i)
	for (int j = 0; j < NbrStateComponents[i]; ++j)
	  NbrStatesPerLValue[StateLValues[i][j] >> 1]++;

      WignerSmallDMatrix** WignerCoefficients = new WignerSmallDMatrix*[NbrMaxLValues + 1];
      RealSymmetricMatrix* CoefficientPerLSector = new RealSymmetricMatrix [NbrMaxLValues + 1];
      RealVector* TmpVectors = new RealVector[NbrStates];
      double* TmpComponents = new double [NbrStates];
      for (int i = 0; i <= NbrMaxLValues; ++i)
	if (NbrStatesPerLValue[i] > 0)
	  {	    
	    CoefficientPerLSector[i] = RealSymmetricMatrix (NbrStates, true);
	    WignerCoefficients[i] = new WignerSmallDMatrix((i << 1) + LShift);
	    int TmpNbrStates = NbrStatesPerLValue[i];	    
	    RealVector* TmpVectors = new RealVector[TmpNbrStates];
	    double* TmpComponents = new double [TmpNbrStates];
	    TmpNbrStates = 0;
	    for (int j = 0; j < NbrStates; ++j)
	      {
		TmpComponents[j] = 0.0;
		for (int k = 0; k < NbrStateComponents[j]; ++k)
		  if (StateLValues[j][k] == ((i << 1) + LShift))
		    {
		      if (TmpVectors[j].ReadVector(StateFileNames[j][k]) == false)
			{
			  cout << "can't open " << StateFileNames[j][k] << endl;
			  return -1;
			}
		      TmpComponents[j] = StateComponents[j][k];
		      ++TmpNbrStates;
		    }
	      }
	    for (int j = 0; j < NbrStates; ++j)
	      {
		if (TmpComponents[j] != 0.0)
		  {
		    for (int k = j; k < NbrStates; ++k)
		      {
			if (TmpComponents[j] != 0.0)
			  {
			    CoefficientPerLSector[i].SetMatrixElement(j, k, (TmpVectors[j] * TmpVectors[k]) * (TmpComponents[j] * TmpComponents[k]));
			  }
		      }
		  }
	      }
	  }
      delete[] TmpComponents;
      delete[] TmpVectors;
      
      cout.precision(14);
      int NbrPoints = Manager.GetInteger("nbr-points");
      int NbrEigenstate = Manager.GetInteger("matrix-eigenstates");
      if (NbrEigenstate > NbrStates)
	NbrEigenstate = NbrStates;
      RealMatrix PreviousEigenvector;
      double Theta = 0.0;
      double ThetaInc = M_PI / ((double) NbrPoints);
      RealSymmetricMatrix TmpMatrix(NbrStates, true);
      double Normalization = 0.0;
      bool NormalizationFlag = false;
      for (int i = 0; i <= NbrPoints; ++i)
	{
	  TmpMatrix.ClearMatrix();
	  for (int j = 0; j <= NbrMaxLValues; ++j)
	    if (NbrStatesPerLValue[j] > 0)
	      {
		TmpMatrix.MultiplyAndAdd((*(WignerCoefficients[j]))(LzValue, LzValue, Theta), CoefficientPerLSector[j]);
	      }
//	  cout << TmpMatrix << endl;

	  for (int j = 0; j < NbrStates; ++j)
	    for (int k = j; k < NbrStates; ++k)
	      TmpMatrix(j, k) *= OverlapMatrix(j, k);
//	  cout << TmpMatrix << endl;
	  if (NormalizationFlag == false)
	    {
	      NormalizationFlag = true;
	      for (int j = 0; j < NbrStates; ++j)
		Normalization += TmpMatrix(j, j);
	      Normalization = 1.0 / Normalization;
	    }
	  TmpMatrix *= Normalization;
	  cout << Theta;
	  if (NbrEigenstate == 0)
	    {
#ifdef __LAPACK__
	      if (Manager.GetBoolean("use-lapack") == true)
		{
		  RealDiagonalMatrix TmpDiag (NbrStates);
		  TmpMatrix.LapackDiagonalize(TmpDiag);
		  TmpDiag.SortMatrixDownOrder();
		  for (int j = 0; j < NbrStates; ++j)
		    cout << " " << TmpDiag[j];
		}
	      else
		{
#endif	    
		  RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrStates);
		  TmpMatrix.Householder(TmpTriDiag, 1e-7);
		  TmpTriDiag.Diagonalize();
		  TmpTriDiag.SortMatrixDownOrder();
		  for (int j = 0; j < NbrStates; ++j)
		    cout << " " << TmpTriDiag(j, j);
#ifdef __LAPACK__
		}
#endif	
	    }
	  else
	    {
	      RealMatrix TmpEigenvector (NbrStates, NbrStates, true);	      
	      for (int l = 0; l < NbrStates; ++l)
		TmpEigenvector(l, l) = 1.0;
#ifdef __LAPACK__
	      if (Manager.GetBoolean("use-lapack") == true)
		{
		  RealDiagonalMatrix TmpDiag (NbrStates);
		  TmpMatrix.LapackDiagonalize(TmpDiag, TmpEigenvector);
		  TmpDiag.SortMatrixDownOrder(TmpEigenvector);
		  for (int j = 0; j < NbrStates; ++j)
		    cout << " " << TmpDiag[j];
		}
	      else
		{
#endif	    
		  RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrStates);
		  TmpMatrix.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
		  TmpTriDiag.Diagonalize(TmpEigenvector);
		  TmpTriDiag.SortMatrixDownOrder(TmpEigenvector);
		  for (int j = 0; j < NbrStates; ++j)
		    cout << " " << TmpTriDiag(j, j);
#ifdef __LAPACK__
		}
#endif	
	      if (PreviousEigenvector.GetNbrRow() == 0)
		{
		  PreviousEigenvector = TmpEigenvector;
		}
	      double Overlap = (TmpEigenvector[0] * PreviousEigenvector[0]);
	      cout << " " << (Overlap * Overlap);
	    }

	  cout << endl;
	  Theta += ThetaInc;
	}

      if (NbrEigenstate > 0)
	{
	  if (strstr(States(0, 0), "_odlro_") == 0)
	    {
	      for (int i = 0; i < NbrStates; ++i)
		{
		  cout << States(0, i) << " : " << PreviousEigenvector[0][i] << endl;
		}
	    }
	  else
	    {
	      for (int i = 0; i < NbrStates; ++i)
		{
		  char* TmpStart = strstr(States(0, i), "_odlro_") + 7;
		  char* TmpEnd = strstr(TmpStart, "_n_");
		  char TmpChar = TmpEnd[0];
		  TmpEnd[0] = '\0';
		  cout << TmpStart << " : " << PreviousEigenvector[0][i] << endl;
		  TmpEnd[0] = TmpChar;
		}
	    }
	    
	}
      for (int i = 0; i < NbrStates; ++i)
	{
	  delete[] StateComponents[i];
	  delete[] StateLValues[i];
	  for (int j = 0; j < NbrStateComponents[i]; ++j)
	    delete[] StateFileNames[i][j] ;
	  delete[] StateFileNames[i];
	}
      delete[] StateComponents;
      delete[] StateLValues;
      delete[] StateFileNames;
      delete[] NbrStateComponents;
      delete[] CoefficientPerLSector;
    }
  return 0;
}


