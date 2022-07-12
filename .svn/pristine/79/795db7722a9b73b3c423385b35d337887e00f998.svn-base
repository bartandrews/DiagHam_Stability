#include "HilbertSpace/ParticleOnSphereManager.h"
#include "HilbertSpace/BosonOnSphere.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasis.h"
#include "HilbertSpace/BosonOnSphereShort.h"
#include "HilbertSpace/BosonOnSphereLong.h"
#include "HilbertSpace/BosonOnSphereSymmetricBasisShort.h"
#include "HilbertSpace/BosonOnCP2.h"
#include "HilbertSpace/FermionOnCP2.h"
#include "HilbertSpace/BosonOnCP2TzSymmetry.h"
#include "HilbertSpace/BosonOnCP2TzZ3Symmetry.h"

#include "Matrix/ComplexLapackDeterminant.h"

#include "FunctionBasis/ParticleOnCP2FunctionBasis.h"

#include "Tools/FQHEWaveFunction/FQHECP2GeneralizedLaughlinWaveFunction.h"

#include "Tools/FQHEFiles/QHEOnSphereFileTools.h"

#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/RandomNumber/FileRandomNumberGenerator.h"
#include "GeneralTools/ArrayTools.h"
#include "GeneralTools/StringTools.h"
#include "GeneralTools/FilenameTools.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/MultiColumnASCIIFile.h"

#include "Architecture/ArchitectureManager.h"
#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticleWaveFunctionOperation.h"

#include "Options/Options.h"

#include "Options/OptionManager.h"
#include "Options/OptionGroup.h"
#include "Options/AbstractOption.h"
#include "Options/BooleanOption.h"
#include "Options/SingleIntegerOption.h"
#include "Options/SingleDoubleOption.h"
#include "Options/SingleStringOption.h"

#include "Vector/ComplexVector.h"
#include "GeneralTools/Endian.h"

#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <stdio.h>
#include <complex>


using std::ios;
using std::cout;
using std::endl;
using std::ofstream;


int main(int argc, char** argv)
{
  cout.precision(14);

  // some running options and help
  OptionManager Manager ("FQHESphereCP2CompareWaveFunctions" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");

  ArchitectureManager Architecture;

  Manager += SystemGroup;
  Architecture.AddOptionGroup(&Manager);
  Manager += MiscGroup;
  
  (*SystemGroup) += new SingleIntegerOption  ('p', "nbr-particles", "number of particles (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('l', "nbr-flux", "number of flux quanta (override autodetection from input file name if non zero)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('i', "nbr-iter", "number of Monte Carlo iterations", 10000);
  (*SystemGroup) += new SingleStringOption  ('\n', "use-exact", "file name of an exact state that has to be used as test wave function"); 
  (*SystemGroup) += new SingleStringOption  ('\n', "reference-state", "name of the file corresponding to the state calculated by exact diagonalization"); 
  (*SystemGroup) += new  SingleStringOption ('\n', "record-wavefunctions", "optional file where each wavefunctions will be tabulated and recorded");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "laughlin-exponent", "exponent of the generalized Laughlin wavefunction", 2);
  
  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");
  
  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type FQHESphereCP2CompareWaveFunctions -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }
  
  if (Manager.GetString("reference-state") == 0)
    {
      cout << "error, a reference state file should be provided. See man page for option syntax or type FQHESphereCP2CompareWaveFunctions -h" << endl;
      return -1;
    }
  if ((Manager.GetString("reference-state") != 0) && 
      (IsFile(Manager.GetString("reference-state")) == false))
    {
      cout << "can't open file " << Manager.GetString("reference-state") << endl;
      return -1;
    }
  char* RecordWaveFunctions = Manager.GetString("record-wavefunctions");
  if (RecordWaveFunctions != 0)
    {
      ofstream RecordFile;
      RecordFile.open(RecordWaveFunctions, ios::out | ios::binary);
      RecordFile.close();
    }
    
    int NbrBosons = Manager.GetInteger("nbr-particles");
    int NbrFluxQuanta = Manager.GetInteger("nbr-flux");
    char* ReferenceStateName = 0;
    bool Statistics = true;
    int TotalTz = 0;
    int TotalY = 0;
    int LaughlinExponent = Manager.GetInteger("laughlin-exponent");
    
    AbstractFunctionBasis* ReferenceBasis = 0;   
  
    ReferenceStateName = new char [strlen(Manager.GetString("reference-state")) + 1];
    strcpy (ReferenceStateName, Manager.GetString("reference-state"));    
    if (FQHEOnSphereFindSystemInfoFromFileName(ReferenceStateName, NbrBosons, NbrFluxQuanta, Statistics) == false)
    {
      cout << "error while retrieving system parameters from file name " << ReferenceStateName << endl;
      return -1; 
    }
    if ((NbrFluxQuanta % LaughlinExponent) == 0 && (NbrBosons != ((NbrFluxQuanta / LaughlinExponent + 1) * (NbrFluxQuanta / LaughlinExponent + 2) / 2)))
    {
      cout << " error : should provide a state with N = (Nphi/m + 1)(Nphi/m + 2)/2 " << endl;    
      return -1;
    }
    
    RealVector TmpPositions (NbrBosons * 4, true);
    double x1;
    double x2;
    double y1;
    double y2;
    double NormExact;
    double NormReference;
    double Argument;
    Complex NormalizationExact;
    RealVector* ReferenceState = new RealVector[1];
    ReferenceBasis = new ParticleOnCP2FunctionBasis(NbrFluxQuanta);
    
    if (ReferenceState[0].ReadVector (ReferenceStateName) == false)
      {
	cout << "can't open vector file " << ReferenceStateName << endl;
	return -1;      
      }
    
    FQHECP2GeneralizedLaughlinWaveFunction* BaseFunction = 0;
    ParticleOnSphere* ReferenceSpace = 0;
    ComplexLapackDeterminant* ExactDeterminant = 0;
    if ((NbrFluxQuanta % LaughlinExponent == 0))
       ExactDeterminant = new ComplexLapackDeterminant(NbrBosons);
    else
      ExactDeterminant = new ComplexLapackDeterminant(3);
    if (RecordWaveFunctions == 0)
      BaseFunction = new FQHECP2GeneralizedLaughlinWaveFunction(ExactDeterminant, NbrBosons, NbrFluxQuanta, LaughlinExponent);
    if (Statistics == false)
      {
	ReferenceSpace = new BosonOnCP2(NbrBosons, NbrFluxQuanta, TotalTz, TotalY);
      }
    else
      {
	ReferenceSpace = new FermionOnCP2(NbrBosons, NbrFluxQuanta, TotalTz, TotalY);
      }
    AbstractRandomNumberGenerator* RandomNumber = 0;
    RandomNumber = new StdlibRandomNumberGenerator (29457);
    Complex ValueExact;
    Complex ValueReference;
    
    
    for (int j = 0; j <= Manager.GetInteger("nbr-iter"); ++j)
      {
	for (int i = 0; i < NbrBosons; ++i)
	  {
	    x1 = acos (1.0 - (2.0 * RandomNumber->GetRealRandomNumber()));
	    x2 = acos (1.0 - (2.0 * RandomNumber->GetRealRandomNumber()));
	    y1 =  2.0 * M_PI * RandomNumber->GetRealRandomNumber();
	    y2 =  2.0 * M_PI * RandomNumber->GetRealRandomNumber();
	    TmpPositions[4*i] = x1;
	    TmpPositions[4*i + 1] = x2;
	    TmpPositions[4*i + 2] = y1;
	    TmpPositions[4*i + 3] = y2;
	  }
	QHEParticleWaveFunctionOperation Operation(ReferenceSpace, &(ReferenceState[0]), &TmpPositions, ReferenceBasis);
        Operation.ApplyOperation(Architecture.GetArchitecture());     
        ValueReference = Operation.GetScalar();
	if (RecordWaveFunctions != 0)
	  {
	    ofstream RecordFile;
	    RecordFile.precision(14);
	    RecordFile.open(RecordWaveFunctions, ios::out | ios::binary | ios::app);
	    for (int i = 0; i < NbrBosons; ++i)
	      RecordFile << TmpPositions[4*i] << " " << TmpPositions[4*i + 1] << " " << TmpPositions[4*i + 2] << " " << TmpPositions[4*i + 1] << " | ";
	    RecordFile << ValueReference;
	    RecordFile << endl;    
	    RecordFile.close();
	 }
	 else
	 {
	  ValueExact = BaseFunction->CalculateFromCoordinates(TmpPositions);
	  NormExact = (ValueExact.Re * ValueExact.Re) + (ValueExact.Im * ValueExact.Im);
	  NormReference = (ValueReference.Re * ValueReference.Re) + (ValueReference.Im * ValueReference.Im);
	  if (j == 0)
	    {
	      Argument = atan2(ValueReference.Im, ValueReference.Re) - atan2(ValueExact.Im, ValueExact.Re);
	      NormalizationExact.Re = cos(Argument);
	      NormalizationExact.Im = sin(Argument);
	      NormalizationExact *= sqrt(NormReference / NormExact);
	    }
	  ValueExact *= NormalizationExact;
	
	cout <<  (ValueReference.Re / ValueExact.Re) << " " << (ValueReference.Im / ValueExact.Im) << endl;
	 }
      }

//     x1 = M_PI / 4.0;
//     x2 =  M_PI / 4.0;
//     y1 =   0;
//     y2 =  0;
//     TmpPositions[0] = x1;
//     TmpPositions[0 + 1] = x2;
//     TmpPositions[0 + 2] = y1;
//     TmpPositions[0 + 3] = y2;
//     TmpPositions[4] = x1;
//     TmpPositions[5 ] = x2;
//     TmpPositions[6] = y1;
//     TmpPositions[7] = y2;
//     x1 = M_PI;
//     x2 = M_PI;
//     y1 = 0;
//     y2 = 0;
//     TmpPositions[4*2] = x1;
//     TmpPositions[4*2 + 1] = x2;
//     TmpPositions[4*2 + 2] = y1;
//     TmpPositions[4*2 + 3] = y2;
//     
//     QHEParticleWaveFunctionOperation Operation(ExactSpace, &(ReferenceState[0]), &TmpPositions, ReferenceBasis);
//     Operation.ApplyOperation(Architecture.GetArchitecture());     
//     ValueReference = Operation.GetScalar();
//     ValueExact = BaseFunction->CalculateFromCoordinates(TmpPositions);
//     NormExact = (ValueExact.Re * ValueExact.Re) + (ValueExact.Im * ValueExact.Im);
//     NormReference = (ValueReference.Re * ValueReference.Re) + (ValueReference.Im * ValueReference.Im);
//     NormalizationExact = NormReference / NormExact;
//     cout << ValueExact << " " << ValueReference << endl;
//     cout << "ratio between reference and test wave function : " << NormalizationExact << endl;

    return 0;
}
