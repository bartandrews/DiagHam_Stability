#include "MCObservables/MultiFileMCHistoryRecord.h"

#include "Options/Options.h"
#include <iostream>

using std::cout;
using std::endl;

int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("MergeMCHistories" , "0.01");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  OptionGroup* MiscGroup = new OptionGroup ("miscellaneous options");

  Manager += SystemGroup;
  Manager += MiscGroup;
  
  (*SystemGroup) += new MultipleStringOption  ('\0', "input-files", "history files to be combined");
  (*SystemGroup) += new SingleStringOption  ('o', "output-file", "file name that the combined history should be written to","merged.samp");

  (*MiscGroup) += new BooleanOption  ('v', "verbose", "be VERY verbose");
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");


  Manager.StandardProceedings(argv, argc, cout);  

  int NbrFiles, NbrPositions=0;
  char** InputFiles = Manager.GetStrings("input-files",NbrFiles);  

  MultiFileMCHistoryRecord Merger(Manager.GetString("output-file"), NbrFiles, InputFiles, NbrPositions, NULL, Manager.GetBoolean("verbose"));
  
  cout << "Merged "<<NbrFiles<<" with a set of "<<NbrPositions<<" coordinate entries and a number of "<<Merger.GetProjectedSamples()<<" estimated entries" << endl;
}
