#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"


#include "MPSObjects/MPOPeratorSixVertexModelTransferMatrixSquare.h"
#include "MPSObjects/RealMPSSite.h"
#include "MPSObjects/DMRGFiniteSizeRealOBCMainTask.h" 

#include "Options/Options.h"
#include "Matrix/RealMatrix.h"

#include "LanczosAlgorithm/LanczosManager.h"

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
  OptionManager Manager ("DMRGSixVertexModelTransferMatrixSquare" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;
  LanczosManager Lanczos(false);

  Architecture.AddOptionGroup(&Manager);
  Manager += SystemGroup;
  Manager += MiscGroup;
  Lanczos.AddOptionGroup(&Manager);

  (*SystemGroup) += new  SingleIntegerOption ('L', "length", "length of the spin chain", 4);
  (*SystemGroup) += new  SingleIntegerOption ('D', "bond-dimension", "bond dimension", 20);
  (*SystemGroup) += new  SingleIntegerOption ('s', "sweep", "number of sweep to be performed", 4);
  (*MiscGroup) += new  BooleanOption ('\n', "test-idmrg", "number of sweep to be performed", false);
  (*MiscGroup) += new  BooleanOption ('\n', "print-tensor", "number of sweep to be performed", false);

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
  
 
  if(Manager.GetBoolean("test-idmrg") == true)
  {
  
   return 0;
  }
 
  int NbrSites = Manager.GetInteger("length");
  int PhysicalDimension = 2;
  int MaxBondDimension = Manager.GetInteger("bond-dimension");;
  MPOPeratorSixVertexModelTransferMatrixSquare TransferMatrix(Architecture.GetArchitecture());
 
if(Manager.GetBoolean("print-tensor") == true)
  {
  TransferMatrix.PrintTensorElements();
   return 0;
  }

  RealMPSSite * Lattice = new RealMPSSite[NbrSites];
  Lattice[0] = RealMPSSite(PhysicalDimension, 0, &Lattice[1], MaxBondDimension,&TransferMatrix);
  for(int i = 1 ; i < NbrSites - 1 ; i++ )
    {
      Lattice[i] = RealMPSSite(PhysicalDimension, &Lattice[i-1], &Lattice[i+1], MaxBondDimension, &TransferMatrix);
    }
  Lattice[NbrSites-1] = RealMPSSite( PhysicalDimension, &Lattice[NbrSites-2], 0, MaxBondDimension,&TransferMatrix);
  int CurrentDimension = 1;
  int NextCurrentDimension = PhysicalDimension;

  for(int i = 0;  i < (NbrSites>>1) ; i++)
  {
    Lattice[i].SetBondDimension(CurrentDimension,NextCurrentDimension);
    Lattice[NbrSites - i - 1].SetBondDimension(NextCurrentDimension,CurrentDimension);

   CurrentDimension = NextCurrentDimension;
   NextCurrentDimension *=  PhysicalDimension;
   if(NextCurrentDimension > MaxBondDimension)
     { 
        NextCurrentDimension =   MaxBondDimension;
     }	
  }
  int NbrSweep = Manager.GetInteger("sweep");
  DMRGFiniteSizeRealOBCMainTask Algorithm (Lattice, &TransferMatrix, NbrSites, NbrSweep, MaxBondDimension, Architecture.GetArchitecture(), &Lanczos);
  Algorithm.RunAlgorithm();
  delete [] Lattice;
  return 0;
}
