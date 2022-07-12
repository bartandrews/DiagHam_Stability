#include "config.h"

#include "Vector/RealVector.h"

#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealMatrix.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/ArrayTools.h"

#include "Operator/ParticleOnSphereSquareTotalSpinOperator.h"
#include "Operator/ParticleOnSphereSquareTotalIsospinOperator.h"

#include "Tools/FQHESpectrum/QHEOnSphereLzSortedSpectrum.h"

#include "HilbertSpace/AbstractQHEParticle.h"
#include "HilbertSpace/FermionOnSphereWithSU4Spin.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>


using std::cout;
using std::endl;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereFermionsSU4Properties" , "0.01");
  OptionGroup* MainGroup = new OptionGroup ("main options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += MiscGroup;
  Manager += MainGroup;
 
  (*MainGroup) += new SingleStringOption  ('\n', "input-file", "name of the file which contains eigenstates description");
  (*MainGroup) += new SingleDoubleOption  ('\n', "ortho-error", "scalar product value below which two states are considered as orthogonal", 1e-12);
  (*MainGroup) += new SingleIntegerOption  ('\n', "output-precision", "numerical display precision", 14, true, 2, true, 14);
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereFermionsSU4Properties -h" << endl;
      return -1;
    }
  
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  double OrthogonalityError = ((SingleDoubleOption*) Manager["ortho-error"])->GetDouble();

  ConfigurationParser OverlapDefinition;
  if (OverlapDefinition.Parse(((SingleStringOption*) Manager["input-file"])->GetString()) == false)
    {
      OverlapDefinition.DumpErrors(cout) << endl;
      return -1;
    }

  unsigned long MemorySpace = (10ul) << 20;
  int LzMax = 0;
  int NbrFermions = 0;
  int TotalLz = 0;
  int SzTotal = 0;
  int IsoSzTotal = 0;
  int TotalEntanglement = 0;
  bool TotalEntanglementFlag = false;

  if ((OverlapDefinition.GetAsSingleInteger("LzMax", LzMax) == false) || (LzMax < 0))
    {
      cout << "LzMax is not defined or has a wrong value" << endl;
      return -1;
    }
  if ((OverlapDefinition.GetAsSingleInteger("NbrFermions", NbrFermions) == false) || (NbrFermions <= 0))
    {
      cout << "NbrFermions is not defined or has a wrong value" << endl;
      return -1;
    }
  if ((OverlapDefinition.GetAsSingleInteger("TotalLz", TotalLz) == false) || ((TotalLz & 1) != ((NbrFermions * LzMax) & 1)))
    {
      cout << "LzMax is not defined or has a wrong value" << endl;
      return -1;
    }
  if ((OverlapDefinition.GetAsSingleInteger("SzTotal", SzTotal) == false) || (SzTotal < (-NbrFermions)) || (SzTotal > NbrFermions))
    {
      cout << "SzTotal is not defined or has a wrong value" << endl;
      return -1;
    }
  if ((OverlapDefinition.GetAsSingleInteger("IsoSzTotal", IsoSzTotal) == false) || (IsoSzTotal < (-NbrFermions)) || (IsoSzTotal > NbrFermions))
    {
      cout << "IsoSzTotal is not defined or has a wrong value" << endl;
      return -1;
    }
  if (OverlapDefinition.GetAsSingleInteger("TotalEntanglement", TotalEntanglement) == true)
    {
      if ((TotalEntanglement< (-NbrFermions)) || (TotalEntanglement > NbrFermions))
	{
	  cout << "TotalEntanglement has a wrong value" << endl;
	  return -1;
	}
      TotalEntanglementFlag = true;
    }
  else
    TotalEntanglementFlag = false;
    
  int NbrUp = (NbrFermions + SzTotal) >> 1;
  int NbrDown = (NbrFermions - SzTotal) >> 1;
  if ((NbrUp < 0) || (NbrDown < 0))
    {
      cout << "This value of the spin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  int NbrPlus = (NbrFermions + IsoSzTotal) >> 1;
  int NbrMinus = (NbrFermions - IsoSzTotal) >> 1;
  if ((NbrPlus < 0) || (NbrMinus < 0))
    {
      cout << "This value of the isospin z projection cannot be achieved with this particle number!" << endl;
      return -1;
    }
  int NbrEntanglementPlus = (NbrFermions + TotalEntanglement) >> 1;
  int NbrEntanglementMinus = (NbrFermions - TotalEntanglement) >> 1;
  if ((TotalEntanglementFlag == true) && ((NbrEntanglementPlus < 0) || (NbrEntanglementMinus < 0)))
    {
      cout << "This value of the entanglement projection cannot be achieved with this particle number!" << endl;
      return -1;
    }

  char** VectorFileNames;
  int NbrVectors;
  if (OverlapDefinition.GetAsStringArray("States", ' ', VectorFileNames, NbrVectors) == false)
    {
      cout << "error while parsing States in " << ((SingleStringOption*) Manager["input-file"])->GetString() << endl;
      return -1;     
    }

  FermionOnSphereWithSU4Spin* Space;
  if (TotalEntanglementFlag == true)
    Space = new FermionOnSphereWithSU4Spin(NbrFermions, TotalLz, LzMax, SzTotal, IsoSzTotal, TotalEntanglement, MemorySpace);
  else
    Space = new FermionOnSphereWithSU4Spin(NbrFermions, TotalLz, LzMax, SzTotal, IsoSzTotal, MemorySpace);

  RealVector* InputVectors = new RealVector [NbrVectors];
  for (int i = 0; i < NbrVectors; ++i)
    if (InputVectors[i].ReadVector(VectorFileNames[i]) == false)
      {
	cout << "error while reading " << VectorFileNames[i] << endl;
	return -1;
      }
  
  ParticleOnSphereSquareTotalSpinOperator S2Oper (Space, LzMax, NbrFermions, SzTotal);
  RealSymmetricMatrix HRep (NbrVectors, NbrVectors);
  for (int k = 0; k < NbrVectors; ++k)
    for (int l = k; l < NbrVectors; ++l)
      HRep.SetMatrixElement(l, k,  S2Oper.MatrixElement(InputVectors[k], InputVectors[l]).Re);
  if (NbrVectors == 1)
    cout << "<S^2> = " << HRep(0, 0) << endl;
  else
    {
      RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrVectors);
      HRep.Householder(TmpTriDiag, 1e-7);
      TmpTriDiag.Diagonalize();
      TmpTriDiag.SortMatrixUpOrder();
      cout << "<S^2> =";
      for (int i = 0; i < NbrVectors; ++i)
	cout << " " << TmpTriDiag.DiagonalElement(i);
      cout << endl;
    }

  ParticleOnSphereSquareTotalIsospinOperator I2Oper (Space, LzMax, NbrFermions, IsoSzTotal);
  for (int k = 0; k < NbrVectors; ++k)
    for (int l = k; l < NbrVectors; ++l)
      HRep.SetMatrixElement(l, k, I2Oper.MatrixElement(InputVectors[k], InputVectors[l]).Re);
  if (NbrVectors == 1)
    cout << "<I^2> = " << HRep(0, 0) << endl;
  else
    {
      RealTriDiagonalSymmetricMatrix TmpTriDiag (NbrVectors);
      HRep.Householder(TmpTriDiag, 1e-7);
      TmpTriDiag.Diagonalize();
      TmpTriDiag.SortMatrixUpOrder();
      cout << "<I^2> =";
      for (int i = 0; i < NbrVectors; ++i)
	cout << " " << TmpTriDiag.DiagonalElement(i);
      cout << endl;
    }


  delete[] InputVectors;
  return 0;
}

