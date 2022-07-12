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
#include "Options/MultipleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cstring>
#include <cstdlib>
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
  OptionManager Manager ("ReorthogonalizeVectorSet" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  Manager += SystemGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new  SingleStringOption ('b', "basis", "name of the file that contains the vector files of the orthonormalized input basis (InputBasis) and the files of the vectors to orthonalize upon (OrthogonalBasis)");
  (*SystemGroup) += new  MultipleStringOption ('o', "ortho-basis", "alternative input for orthogonal basis vectors(overrides basis file)");
  (*SystemGroup) += new  MultipleStringOption ('i', "input-basis", "alternative input for known input basis vectors");
  (*SystemGroup) += new SingleDoubleOption ('\n', "error", "bound above which vectors are considered as linearly independent", 1e-10);
  (*SystemGroup) += new  SingleStringOption ('\n', "vector-prefix", "prefix to use for each vector of the basis", "vector");


  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type ReorthogonalizeVectorSet -h" << endl;
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

  int NbrInputVectors;
  char** InputVectorFileNames;
  int NbrOrthogonalVectors;
  char** OrthogonalVectorFileNames;
  bool AddRandomVector=false;
  char* DirectoryName;

  if ((OrthogonalVectorFileNames=Manager.GetStrings("ortho-basis",NbrOrthogonalVectors))==NULL)
    {
      ConfigurationParser ReducedBasis;
      if (BasisDescription==NULL)
	{
	  cout << "Input of either --ortho-basis or --basis is required!"<<endl;
	  exit(1);
	}
      if (ReducedBasis.Parse(BasisDescription) == false)
	{
	  ReducedBasis.DumpErrors(cout) << endl;
	  return -1;
	}
      if (ReducedBasis.GetAsStringArray("OrthogonalBasis", ' ', OrthogonalVectorFileNames, NbrOrthogonalVectors) == false)
	{
	  cout << "Orthogonal basis vectors are not defined or have a wrong value in " << BasisDescription << endl;
	  return -1;
	}
      
      if (ReducedBasis.GetAsStringArray("InputBasis", ' ', InputVectorFileNames, NbrInputVectors) == false)
	{
	  cout << "Input basis vectors are not defined or have a wrong value in " << BasisDescription << endl;
	  return -1;
	}
      DirectoryName = ReducedBasis["Directory"];
    }
  else
    {
      if ((InputVectorFileNames=Manager.GetStrings("input-basis",NbrInputVectors))==NULL)
	{
	  NbrInputVectors=0;
	  AddRandomVector=true;
	}
      DirectoryName = new char[4];
      sprintf(DirectoryName,"./");
    }
      
  RealVector* InputBasis = new RealVector[NbrInputVectors+1];
  char* TmpName;
  long Dimension=0;
  for (int i = 0; i < NbrInputVectors; ++i)
    {
      TmpName = InputVectorFileNames[i];
      if (DirectoryName != 0)
	{
	  TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
	}
      cout << "reading vector " << TmpName << endl;
      if (InputBasis[i].ReadVector(TmpName) == false)
	{
	  cout << "error while reading " << TmpName << endl;
	  if (DirectoryName != 0)
	    delete[] DirectoryName;
	  for (int j = 0; j < NbrInputVectors; ++j)
	    delete[] InputVectorFileNames[j];
	  delete[] InputVectorFileNames;
	  return -1;
	}
      if (DirectoryName != 0)
	delete[] TmpName;
      long TmpDimension;
      if (InputBasis[i].IsLargeVector())
	TmpDimension = InputBasis[i].GetLargeVectorDimension();
      else
	TmpDimension = InputBasis[i].GetVectorDimension();
      if (Dimension==0)
	Dimension=TmpDimension;
      else
	if (TmpDimension!=Dimension)
	  {
	    cout << "Attention, vector "<<TmpName<<" has the wrong dimension"<<endl;
	    exit(1);
	  }
    }
  
  RealVector* OrthogonalBasis = new RealVector[NbrOrthogonalVectors];
  for (int i = 0; i < NbrOrthogonalVectors; ++i)
    {
      TmpName = OrthogonalVectorFileNames[i];
      if (DirectoryName != 0)
	{
	  TmpName = ConcatenatePathAndFileName(DirectoryName, TmpName);
	}
      cout << "reading vector " << TmpName << endl;
      if (OrthogonalBasis[i].ReadVector(TmpName) == false)
	{
	  cout << "error while reading " << TmpName << endl;
	  if (DirectoryName != 0)
	    delete[] DirectoryName;
	  for (int j = 0; j < NbrOrthogonalVectors; ++j)
	    delete[] OrthogonalVectorFileNames[j];
	  delete[] OrthogonalVectorFileNames;
	  return -1;
	}
      OrthogonalBasis[i]  /= OrthogonalBasis[i].Norm();

      if (DirectoryName != 0)
	delete[] TmpName;
      long TmpDimension;
      if (InputBasis[i].IsLargeVector())
	TmpDimension = OrthogonalBasis[i].GetLargeVectorDimension();
      else
	TmpDimension = OrthogonalBasis[i].GetVectorDimension();
      if (Dimension==0)
	Dimension=TmpDimension;
      else
	if (TmpDimension!=Dimension)
	  {
	    cout << "Attention, vector "<<TmpName<<" has the wrong dimension"<<endl;
	    exit(1);
	  }

    }
  if (Dimension==0)
    {
      cout << "Error: No input basis or orthogonal basis vectors found"<<endl;
      exit(1);
    }
  if (AddRandomVector)
    {
      InputBasis[NbrInputVectors].Resize(Dimension);
      for (int i = 0; i < Dimension; i++)
	InputBasis[NbrInputVectors][i] = (rand() - 32767) * 0.5;
      InputBasis[NbrInputVectors] /= InputBasis[NbrInputVectors].Norm();
      ++NbrInputVectors;
    }
  double* Coefficients = new double[NbrOrthogonalVectors];
  for (int i = 0; i < NbrInputVectors; ++i)
    {
     for (int j = 0; j < NbrOrthogonalVectors; ++j)
       Coefficients[j] = OrthogonalBasis[j] * InputBasis[i];     
     for (int j = 0; j < NbrOrthogonalVectors; ++j)
       InputBasis[i].AddLinearCombination(-Coefficients[j], OrthogonalBasis[j]);
     cout << i << " : " << InputBasis[i].Norm() << endl;
    }

  RealSymmetricMatrix HRep (NbrInputVectors);
  for (int i = 0; i < NbrInputVectors; ++i)
    for (int j = 0; j < NbrInputVectors; ++j)
      HRep(i ,j) = InputBasis[j] * InputBasis[i];

  RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrInputVectors);

  
  RealMatrix TmpEigenvector (NbrInputVectors, NbrInputVectors, true);	      
  for (int l = 0; l < NbrInputVectors; ++l)
    TmpEigenvector(l, l) = 1.0;
  HRep.Householder(TmpTriDiag, 1e-7, TmpEigenvector);
  TmpTriDiag.Diagonalize(TmpEigenvector);
  TmpTriDiag.SortMatrixUpOrder(TmpEigenvector);
  char* OutputVectorFileName = new char [strlen(VectorPrefix) + 32];
  int TmpDimension = InputBasis[0].GetVectorDimension();
  int Count = 0;
  for (int i = 0; i < NbrInputVectors; ++i)
    {
      cout << TmpTriDiag.DiagonalElement(i) << " ";
      if (fabs( TmpTriDiag.DiagonalElement(i)) > Error)
	{
	  RealVector TmpVector (TmpDimension, true);
	  for (int j = 0; j < NbrInputVectors; ++j)
	    for (int k = 0; k < TmpDimension; ++k)
	      {
		TmpVector[k] += TmpEigenvector[i][j] * InputBasis[j][k];
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

  if (DirectoryName != 0)
    delete[] DirectoryName;
//   int Offset = (AddRandomVector?1:0);
//   for (int j = 0; j < NbrInputVectors-Offset; ++j)
//     delete[] InputVectorFileNames[j];
//   delete[] InputVectorFileNames;
  return 0;
}
