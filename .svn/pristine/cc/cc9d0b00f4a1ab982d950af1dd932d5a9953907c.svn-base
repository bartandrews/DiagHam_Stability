#include "Matrix/RealTriDiagonalSymmetricMatrix.h"
#include "Matrix/RealSymmetricMatrix.h"

#include "HilbertSpace/TrappedBosons.h"
#include "Hamiltonian/TrappedBosonHamiltonian.h"

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>


using std::cout;
using std::endl;
using std::ofstream;
using std::ios;

int main(int argc, char** argv)
{
  int NbrBosons = 4;
  int MMin = 1;
  int MMax = 1;
  if (argc >= 2)
    NbrBosons = atoi (argv[1]);
  if (argc >= 3)
    MMin = atoi (argv[2]);
   if (argc >= 4)
    MMax = atoi (argv[3]);
   char* OutputName = "test.dat";
   if (argc >= 5)
     OutputName = argv[4];
   ofstream File;
   File.open(OutputName, ios::binary | ios::out);
   for (int  L = MMin; L <= MMax; ++L)
     {
       TrappedBosons Space (NbrBosons, L);
       cout << "Nbr bosons = " << NbrBosons << "    L = " << L << "    Dimension = " << Space.GetHilbertSpaceDimension() << endl;
       TrappedBosonHamiltonian Hamiltonian(&Space, L);
       RealSymmetricMatrix HRep (Hamiltonian.GetHilbertSpaceDimension());
       Hamiltonian.GetHamiltonian(HRep);
       if (Hamiltonian.GetHilbertSpaceDimension() > 1)
	 {
	   RealTriDiagonalSymmetricMatrix TmpTriDiag (Hamiltonian.GetHilbertSpaceDimension());
	   HRep.Householder(TmpTriDiag, 1e-7);
	   TmpTriDiag.Diagonalize();
	   TmpTriDiag.SortMatrixUpOrder();
	   cout << "eigenvalues : " << endl;
	   for (int j = 0; j < Hamiltonian.GetHilbertSpaceDimension(); j++)
	     {
	       //	  cout << TmpTriDiag.DiagonalElement(j) << " ";
	       double TmpVal = (TmpTriDiag.DiagonalElement(j));
	       File << ((int) 10);//(int) L;
	       File << " ";
	       File << TmpVal;
	       File << endl;
	     }
	 }
       else
	 {
	   double TmpVal = HRep(0, 0);
	   File << L << " " << TmpVal << endl;
	 }
       cout << endl;
       cout << "//////////////////////////////////////////////////////" << endl;
     }
   File.close();
  return 0;
}

