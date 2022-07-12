#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("ExtractLinearlyIndependentVectors" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new  SingleStringOption ('b', "basis", "name of the file that contains the vector files used to describe the basis");
  (*SystemGroup) += new BooleanOption ('\n', "check-only", "check how many vectors are linearly independent without extracting the basis");
  (*SystemGroup) += new SingleDoubleOption ('\n', "error", "bound above which vectors are consider as linearly independent", 1e-10);
  (*SystemGroup) += new  SingleStringOption ('\n', "vector-prefix", "prefix to use for each vector of the basis", "vector");


  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ExtractLinearlyIndependentVectors -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  char* BasisDescription = ((SingleStringOption*) Manager["basis"])->GetString();  
  char* VectorPrefix = ((SingleStringOption*) Manager["vector-prefix"])->GetString();
  double Error = ((SingleDoubleOption*) Manager["error"])->GetDouble();

  ConfigurationParser ReducedBasis;
  if (ReducedBasis.Parse(BasisDescription) == false)
    {
      ReducedBasis.DumpErrors(cout) << endl;
      return -1;
    }
  int NbrVectors;
  char** VectorFileNames;
  if (ReducedBasis.GetAsStringArray("Basis", ' ', VectorFileNames, NbrVectors) == false)
    {
      cout << "Vectors are not defined or have a wrong value in " << BasisDescription << endl;
      return -1;
    }

  RealVector* Basis = new RealVector[NbrVectors];
  char* DirectoryName = ReducedBasis["Directory"];
  char* TmpName;
  for (int i = 0; i < NbrVectors; ++i)
    {
      TmpName = VectorFileNames[i];
      if (DirectoryName != 0)
	{
	  TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
	}
      cout << "reading vector " << TmpName << endl;
      if (Basis[i].ReadVector(TmpName) == false)
	{
	  cout << "error while reading " << TmpName << endl;
	  if (DirectoryName != 0)
	    delete[] DirectoryName;
	  for (int j = 0; j < NbrVectors; ++j)
	    delete[] VectorFileNames[j];
	  delete[] VectorFileNames;
	  return -1;
	}
      if (DirectoryName != 0)
	delete[] TmpName;
    }
 
  RealSymmetricMatrix HRep (NbrVectors);
  for (int i = 0; i < NbrVectors; ++i)
    for (int j = 0; j < NbrVectors; ++j)
      HRep(i ,j) = Basis[j] * Basis[i];

  RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrVectors);


  if (((BooleanOption*) Manager["check-only"])->GetBoolean() == true)	
    {
      HRep.Householder(TmpTriDiag, 1e-7);
      TmpTriDiag.Diagonalize();
      TmpTriDiag.SortMatrixUpOrder();
      int Count = 0;
      for (int i = 0; i < NbrVectors; ++i)
	{
	  cout << TmpTriDiag.DiagonalElement(i) << " ";
	  if (fabs( TmpTriDiag.DiagonalElement(i)) > Error)
	    Count++;
	}
      cout << endl;
      cout << Count << " linearly independent vectors" << endl;
    }
  else
    {
      RealMatrix TmpEigenvector (NbrVectors, NbrVectors, true);	      
      for (int l = 0; l < NbrVectors; ++l)
	TmpEigenvector(l, l) = 1.0;
      HRep.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
      TmpTriDiag.Diagonalize(TmpEigenvector);
      TmpTriDiag.SortMatrixUpOrder(TmpEigenvector);
      char* OutputVectorFileName = new char [strlen(VectorPrefix) + 32];
      int TmpDimension = Basis[0].GetVectorDimension();
      int Count = 0;
      for (int i = 0; i < NbrVectors; ++i)
	{
	  cout << TmpTriDiag.DiagonalElement(i) << " ";
	  if (fabs( TmpTriDiag.DiagonalElement(i)) > Error)
	    {
	      RealVector TmpVector (TmpDimension, true);
	      for (int j = 0; j < NbrVectors; ++j)
		for (int k = 0; k < TmpDimension; ++k)
		  {
		    TmpVector[k] += TmpEigenvector[i][j] * Basis[j][k];
		  }
	      TmpVector /= TmpVector.Norm();
	      sprintf (OutputVectorFileName, "%s_%d.vec", VectorPrefix, Count);
	      TmpVector.WriteVector(OutputVectorFileName);
	      Count++;
	    }
	}
      cout << endl;
      cout << Count << " linearly independent vectors" << endl;
      delete[] OutputVectorFileName;
    }

  if (DirectoryName != 0)
    delete[] DirectoryName;
  for (int j = 0; j < NbrVectors; ++j)
    delete[] VectorFileNames[j];
  delete[] VectorFileNames;
  return 0;
}
