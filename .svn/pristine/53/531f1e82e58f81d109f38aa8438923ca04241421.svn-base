#include "MCObservables/MCHistoryRecord.h"

#include "Options/Options.h"
#include <iostream>

using std::cout;
using std::endl;

using std::ofstream;
using std::ios;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("MCHistoryToText" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("miscellaneous options");

  Manager += SystemGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleStringOption  ('\0', "input-file", "history files to be combined");
  (*SystemGroup) += new SingleStringOption  ('o', "output-file", "file name that the combined history should be written to","history.txt");
  (*SystemGroup) += new SingleIntegerOption  ('p', "precision", "number of digits precision for output",10);

  

  (*MiscGroup) += new BooleanOption  ('v', "verbose", "be VERY verbose");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  Manager.StandardProceedings(argv, argc, cout);  

  int NbrPositions=0;
  char* InputFile = Manager.GetString("input-file");


  MCHistoryRecord Record(InputFile, NbrPositions /* could add additional observables here */);

  cout << "Converting History record with "<<NbrPositions<<" coordinate entries"<<endl;

  RealVector Positions(NbrPositions);
  int SampleCount;
  Complex Value;
  double SamplingAmplitude;

  ofstream Out;

  Out.open(Manager.GetString("output-file"),ios::out);

  // line 1: NbrPositions, line 2: two component output (complex)
  Out << NbrPositions <<"," << endl << "2," << endl;

  Out.precision(Manager.GetInteger("precision"));
  
  while ( Record.GetMonteCarloStep(SampleCount, SamplingAmplitude, &(Positions[0]), Value))
    {
      Out << Positions[0]<<",";
      for (int i=1; i<NbrPositions; ++i)
	Out << Positions[i]<<",";
      Out << endl << Value.Re << ","<<Value.Im<<","<<endl;
    }
  Out.close();
  cout << "Done." << endl;
}
