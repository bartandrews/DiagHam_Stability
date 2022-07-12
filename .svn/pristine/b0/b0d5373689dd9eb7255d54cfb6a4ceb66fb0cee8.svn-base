#include "Options/Options.h"

#include "MathTools/ThreeJSymbol.h"

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
  OptionManager Manager ("PrintAllThreeJ" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");

  Manager += SystemGroup;
  Manager += MiscGroup;

  (*SystemGroup) += new SingleIntegerOption  ('j', "j1", "twice the first angular momentum to be coupled",0);
  (*SystemGroup) += new SingleIntegerOption  ('k', "j2", "twice the second angular momentum to be coupled",0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "j3", "twice the third angular momentum to be coupled",0);
  (*SystemGroup) += new SingleIntegerOption  ('m', "m3", "twice the z-component of the third angular momentum",0);

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);
  
  int J1=Manager.GetInteger("j1");
  int J2=Manager.GetInteger("j2");
  int J3=Manager.GetInteger("j3");
  int M3=Manager.GetInteger("m3");  

  if ((J3<abs(J1-J2))||(J3>J1+J2))
    {
      cout << "J1,J2,J3 do not satisfy the triangle relations"<<endl;
      exit(1);
    }
  if (abs(M3)>J3)
    {
      cout << "M3 is out of bounds (-J3<=M3<=J3)"<<endl;
      exit(1);
    }

  if ((M3^J3)&1)
    {
      cout << "M3 and J3 need to have the same parity"<<endl;
      exit(1);
    }

  ThreeJSymbol ThreeJ(J1,J2);
  
  int m2;
  for (int m1=-J1; m1<=J1; m1 += 2)
    {
      m2=-M3-m1;
      if((m2>=-J2)&&(m2<=J2))
	{
	  cout << "ThreeJSymbol[{";
	  if(J1&1)
	    cout<<J1<<"/2, "<<m1<<"/2}, {";
	  else
	    cout<<J1/2<<", "<<m1/2<<"}, {";
	  if(J2&1)
	    cout<<J2<<"/2, "<<m2<<"/2}, {";
	  else
	    cout<<J2/2<<", "<<m2/2<<"}, {";
	  if (J3&1)
	    cout<<J3<<"/2, "<<M3<<"/2}] = ";
	  else
	    cout<<J3/2<<", "<<M3/2<<"}] = ";
	  cout<<ThreeJ.GetCoefficient (m1, m2, J3)<<endl;
	}
    }

}
