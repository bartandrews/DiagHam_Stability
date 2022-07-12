#include "Vector/RealVector.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "HilbertSpace/FermionOnSphere.h"
#include "HilbertSpace/FermionOnSphereWithSpin.h"
#include "HilbertSpace/FermionOnSphereWithSpinLong.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleStringOption.h"
#include "Options/SingleDoubleOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <fstream>



using std::cout;
using std::endl;
using std::ios;
using std::ofstream;
using std::ifstream;



int main(int argc, char** argv)
{

  OptionManager Manager ("FQHESphereFermionsWithSpinEntanglementEntropy" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* OutputGroup = new OptionGroup ("output options");

  Manager += SystemGroup;
  Manager += OutputGroup;
  Manager += MiscGroup;

  int nbrFermions=0;
  int LzMax=0;
  int TotalLz=0;
  int Sz=0;
  bool Statistics = true;
  ParticleOnSphereWithSpin* Space = 0;
  
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "name of the file describing the system ground state");
  (*OutputGroup) += new SingleStringOption ('\n', "density-matrix", "store the eigenvalues of the reduced density matrices in the a given file");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsWithSpinEntanglementEntropy -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  char* DensityMatrixFileName = Manager.GetString("density-matrix");
  char* FileName = Manager.GetString("input-file");
  if (FileName == 0)
    {
      cout << " an input file has to be provided" << endl;
      return -1;
    }

  if (FQHEOnSphereWithSpinFindSystemInfoFromVectorFileName(FileName,nbrFermions,LzMax,TotalLz,Sz,Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << FileName << endl;
      return -1;
    }
  
#ifdef __64_BITS__
  if (LzMax <= 31)
#else
    if (LzMax <= 15)
#endif
      {
	Space=new FermionOnSphereWithSpin  (nbrFermions, TotalLz, LzMax, Sz);
      }
    else
      {
#ifdef __128_BIT_LONGLONG__
	if (LzMax <= 63)
#else
	  if (LzMax <= 31)
#endif
	    {
	      Space =new FermionOnSphereWithSpinLong (nbrFermions, TotalLz, LzMax, Sz);
	    }
	  else
	    {
	      cout << "States of this Hilbert space cannot be represented in a single word." << endl;
	      return 0;
	    }	
      }
  
  int nbrFermionsUp= (nbrFermions+Sz)>>1;
  
  RealVector GroundState;
  if (GroundState.ReadVector (FileName) == false)
    {
      cout << "can't open vector file " << FileName << endl;
      return -1;      
    }   
  
  if (Space->GetHilbertSpaceDimension() != GroundState.GetVectorDimension())
    {
      cout << "Number of rows of the vector is not equal to the Hilbert space dimension!";
      return -1;
    }

  if (DensityMatrixFileName != 0)
    {
      ofstream DensityMatrixFile;
      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out); 
      DensityMatrixFile << "# Lz    lambda" << endl;
      DensityMatrixFile.close();
    }

  double EntanglementEntropy = 0.0;
  double DensitySum =0.0;
  int MaxLzUp = (LzMax*nbrFermionsUp-nbrFermionsUp*(nbrFermionsUp-1)) + 1;
  for(int lzUp= (-LzMax*nbrFermionsUp+nbrFermionsUp*(nbrFermionsUp-1)); lzUp < MaxLzUp; lzUp+=2)
    {
      RealSymmetricMatrix PartialDensityMatrix;
      PartialDensityMatrix= Space->EvaluatePartialDensityMatrixSpinSeparation(lzUp,GroundState);

      if (PartialDensityMatrix.GetNbrRow() > 1)
	{
	  RealDiagonalMatrix TmpDiag(PartialDensityMatrix.GetNbrRow());
	  PartialDensityMatrix.Diagonalize(TmpDiag);
	  TmpDiag.SortMatrixDownOrder();
	  for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
	    {
	      if (TmpDiag[i] > 1e-14)
		{
		  EntanglementEntropy += TmpDiag[i] * log(TmpDiag[i]);
		  DensitySum += TmpDiag[i];
		}
	    }
	  if (DensityMatrixFileName != 0)
	    {
	      ofstream DensityMatrixFile;
	      DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
	      DensityMatrixFile.precision(14);
	      for (int i = 0; i < PartialDensityMatrix.GetNbrRow(); ++i)
		DensityMatrixFile << (0.5 * ((double) lzUp)) << " " << TmpDiag[i] << endl;
	      DensityMatrixFile.close();
	    }
	}
      else 
	{
	  if (PartialDensityMatrix(0,0) > 1e-14)
	    {
	      EntanglementEntropy += PartialDensityMatrix(0,0) * log(PartialDensityMatrix(0,0));
	      DensitySum += PartialDensityMatrix(0,0);
	      if (DensityMatrixFileName != 0)
		{
		  ofstream DensityMatrixFile;
		  DensityMatrixFile.open(DensityMatrixFileName, ios::binary | ios::out | ios::app); 
		  DensityMatrixFile.precision(14);
		  DensityMatrixFile << (0.5 * ((double) lzUp)) << " " << PartialDensityMatrix(0,0) << endl;
		  DensityMatrixFile.close();
	    }
	    }
	}
    }

  cout << (-EntanglementEntropy) << " " << DensitySum;
  return 0;
}
