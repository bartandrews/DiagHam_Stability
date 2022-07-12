#include "GeneralTools/ListIterator.h"
#include "MathTools/IntegerAlgebraTools.h"
#include "QuantumNumber/AbstractQuantumNumber.h"
#include "HilbertSpace/SubspaceSpaceConverter.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/MainTaskOperation.h"
#include "Architecture/ArchitectureOperation/VectorHamiltonianMultiplyOperation.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <fstream>

#define M1_12 0.08333333333333333

using std::cout;
using std::cin;
using std::endl;
using std::ofstream;
using std::ios;

double MadelungEnergy(double costheta, double aspect);


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("WignerEnergy" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  ArchitectureManager Architecture;

  Manager += SystemGroup; 
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles", 6);
  (*SystemGroup) += new SingleIntegerOption  ('l', "max-momentum", "maximum momentum for a single particle", 18);
  (*SystemGroup) += new SingleDoubleOption   ('r', "ratio", 
					      "ratio between lengths along the x and y directions (-1 if has to be taken equal to nbr-particles/4)", -1);
  (*SystemGroup) += new SingleDoubleOption   ('c', "costheta", 
					      "cosine of the angle between the vectors of the unit cell (0 for rectangular cell)", 0.0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type QHEBosonsDelta -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int NbrFermions = ((SingleIntegerOption*) Manager["nbr-particles"])->GetInteger();
  int MaxMomentum = ((SingleIntegerOption*) Manager["max-momentum"])->GetInteger();
  double XRatio = NbrFermions / 4.0;
  if (((SingleDoubleOption*) Manager["ratio"])->GetDouble() > 0)
    {
       XRatio = ((SingleDoubleOption*) Manager["ratio"])->GetDouble();
    }
  double CosTheta = ((SingleDoubleOption*) Manager["costheta"])->GetDouble();

  double WignerEnergy = (-NbrFermions * MadelungEnergy(CosTheta, XRatio)/(2.0 * sqrt(2.0 * M_PI * MaxMomentum)));
  cout<<"Calculating Wigner energy for N= "<<NbrFermions<<" Nphi= "<<MaxMomentum<<" Ratio= "<<XRatio<<" costheta= "<<CosTheta<<endl;
  cout << "Wigner Energy per electron: " << (WignerEnergy / NbrFermions) << endl;
  cout << "Total Wigner Energy : "<< WignerEnergy << endl;
  return 0;
}

double MadelungEnergy(double costheta, double aspect)
{
   double  rootpi = sqrt(M_PI);
   double sintheta = sqrt(1.0 - costheta * costheta); 
   double A = M_PI * aspect/sintheta;
   double B = M_PI/(aspect * sintheta);
   double C = 2.0 * M_PI * costheta/sintheta;

   int L = 0;
   double Energy = 4.0;
   double kmax, dcrmnt, q1, q2, qmax, E1, E2, E;
   int L1, L2;

   qmax = 9.25;

   int indicator = 0;
   while (indicator == 0)
    {
      L++;
      kmax = 2*L;
      dcrmnt = 0.0;
      for (int k = 1;  k <= kmax; ++k)
       {
         L1 = L;
         L2 = -L+k;
         q1 = A*L1*L1 + B*L2*L2 + C*L1*L2;
         q1 = sqrt(q1);

         L1 = L-k;
         L2 = L;
         q2 = A*L1*L1 + B*L2*L2 + C*L1*L2;
         q2 = sqrt(q2);

         E1 = 0.0;
         if (q1 <= qmax) 
           E1 = (erfc(q1))/q1;

         E2 = 0.0;
         if (q2 <= qmax) 
           E2 = (erfc(q2))/q2;

         E = 2.0 * rootpi * (E1+E2);
         E *= 2.0;
         dcrmnt += E;
       } 

      Energy = Energy - dcrmnt;

      if (L >= 3) indicator = 1;
      if (dcrmnt <= 1e-12) indicator = 1;
    }

  return Energy;

}
