#include "Options/Options.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
#include "Tools/FQHESpectrum/AbstractZDensityProfile.h"

#include "Vector/RealVector.h"
 
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "GeneralTools/FilenameTools.h"


#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <math.h>

using std::ofstream;
using std::ios;
using std::cout;
using std::endl;



int main(int argc, char** argv)
{
  OptionManager Manager ("CoulombPseudopotentials" , "0.01");
  OptionGroup* MiscGroup = new OptionGroup ("misc options");
  OptionGroup* SystemGroup = new OptionGroup ("system options");
  Manager += SystemGroup;
  Manager += MiscGroup;


  (*SystemGroup) += new SingleIntegerOption  ('l', "landau-level", "index of the Landau level (0 for the lowest Landau level)", 0, true, 0);
  (*SystemGroup) += new SingleIntegerOption  ('s', "nbr-flux", "number of flux quanta (i.e. twice the maximum momentum for a single particle)", 8);
  (*SystemGroup) += new BooleanOption ('\n', "add-impurities", "add two impurities (one at each pole)");
  (*SystemGroup) += new SingleDoubleOption ('\n', "north-potential", "potential assosciated to the impurity at the north pole", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "south-potential", "potential assosciated to the impurity at the south pole", 0.0);
  (*SystemGroup) += new  BooleanOption ('\n', "relativistic-fermions", "assume relativistic fermions");
  
  (*SystemGroup) += new  SingleDoubleOption ('d', "layer-separation", "assume finite layer separation",0.0);

  (*SystemGroup) += new  SingleDoubleOption ('t', "layer-thickness", "assume finite layer thickness",0.0);
  
  (*SystemGroup) += new  SingleStringOption ('\n', "profile-type", "type of density-profile (1=infiniteWell, 2=Fang-Howard, 3=infiniteWellExc)","1");

  (*SystemGroup) += new  SingleStringOption ('\n', "other-profile", "type of density-profile for second layer (1=infiniteWell, 2=Fang-Howard, 3=infiniteWellExc)",NULL);

  (*SystemGroup) += new  SingleIntegerOption ('\n', "profile-points", "number of points where profile is evaluated", 200);

  (*SystemGroup) += new  SingleDoubleOption ('\n', "profile-multiplier", "multiplier for number of integration points",1.0);

  (*SystemGroup) += new  BooleanOption ('n', "no-interpolation","do not use interpolation for finite width calculations");

  (*SystemGroup) += new  BooleanOption ('n', "nbody", "add n-body potentials");
  (*SystemGroup) += new  MultipleDoubleOption ('p', "nbody-potentials", "values of n-body potentials to be added (separated by ','",',');
  (*SystemGroup) += new SingleDoubleOption ('v', "add-v0", "add some amount to v_0",0.0);
  (*SystemGroup) += new SingleDoubleOption ('w', "add-v1", "add some amount to v_1",0.0);

  (*SystemGroup) += new BooleanOption  ('\n', "tab", "output a table with the results in columns (extension .tab)");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name (default is pseudopotential_coulomb_l_x_2s_y[_v_xxx_w_yyy].dat)");
  (*SystemGroup) += new BooleanOption ('\n', "std-output", "use standard output instead of an output file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type CoulombPseudopotentials -h" << endl;
      return -1;
    }
  if (((BooleanOption*) Manager["help"])->GetBoolean() == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }

  int LandauLevel = ((SingleIntegerOption*) Manager["landau-level"])->GetInteger();
  int NbrFlux = ((SingleIntegerOption*) Manager["nbr-flux"])->GetInteger();
  double layerSeparation = Manager.GetDouble("layer-separation");
  int MaxMomentum = NbrFlux + (LandauLevel << 1);
  
  char* OutputFile;
  ((SingleDoubleOption*)Manager["add-v0"])->SetStringFormat("%g");
  ((SingleDoubleOption*)Manager["add-v1"])->SetStringFormat("%g");
  ((SingleDoubleOption*)Manager["layer-thickness"])->SetStringFormat("%g");
  ((SingleDoubleOption*)Manager["layer-separation"])->SetStringFormat("%g");
  if (((SingleStringOption*) Manager["output"])->GetString() == 0l)
    {
      if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
	OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_relativistic_l_%landau-level%_2s_%nbr-flux%.dat");
      else
	{
	  if (Manager.GetDouble("layer-thickness")>0.0)
	    {
	      if (layerSeparation==0.0)
		OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_t_%layer-thickness%_%profile-type%.dat");
	      else
		OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_t_%layer-thickness%_d_%layer-separation%_%profile-type%.dat");
	    }
	  else
	    {	      
	      if ((Manager.GetDouble("add-v0")==0.0)&&(Manager.GetDouble("add-v1")==0.0))
		{
		  if (layerSeparation==0.0)
		    OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%.dat");
		  else
		    OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_d_%layer-separation%.dat");
		}
	      else
		{
		  if (Manager.GetDouble("add-v0")==0.0)
		    {
		      if (layerSeparation==0.0)
			OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_w_%add-v1%.dat");
		      else
			OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_d_%layer-separation%_w_%add-v1%.dat");
		    }
		  else if (Manager.GetDouble("add-v1")==0.0)
		    {
		      if (layerSeparation==0.0)
			OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_v_%add-v0%.dat");
		      else
			OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_d_%layer-separation%_v_%add-v0%.dat");
		    }
		  else
		    {
		      if (layerSeparation==0.0)
			OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_v_%add-v0%_w_%add-v1%.dat");
		      else
			OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_d_%layer-separation%_v_%add-v0%_w_%add-v1%.dat");
		    }
		}
	    }
	}
    }
  else
    {
      OutputFile = new char [strlen(((SingleStringOption*) Manager["output"])->GetString()) + 1];
      strcpy (OutputFile, ((SingleStringOption*) Manager["output"])->GetString());
    }

      double* Pseudopotentials;
  if (Manager.GetDouble("layer-thickness")>0.0)
    {      
      AbstractZDensityProfile *Profile = AbstractZDensityProfile::CreateZDensityProfile(Manager.GetString("profile-type"),Manager.GetDouble("layer-thickness"));
      AbstractZDensityProfile *OtherProfile=NULL;
      if (Manager.GetString("other-profile")!=NULL)
	OtherProfile = AbstractZDensityProfile::CreateZDensityProfile(Manager.GetString("other-profile"),Manager.GetDouble("layer-thickness"));
      if (Manager.GetBoolean("no-interpolation"))
	Pseudopotentials = EvaluateFiniteWidthPseudoPotentialsNoInterpolation(NbrFlux, LandauLevel, Profile,
							    /*, points=200 */ Manager.GetInteger("profile-points"),
							    /*, multiplier=5.0 */ Manager.GetDouble("profile-multiplier"),
							    Manager.GetDouble("layer-separation"), OtherProfile);
      else
	Pseudopotentials = EvaluateFiniteWidthPseudoPotentials(NbrFlux, LandauLevel, Profile,
							    /*, points=200 */ Manager.GetInteger("profile-points"),
							    /*, multiplier=5.0 */ Manager.GetDouble("profile-multiplier"),
							    Manager.GetDouble("layer-separation"), OtherProfile);
    }
  else
    {
      Pseudopotentials = EvaluatePseudopotentials(NbrFlux, LandauLevel, layerSeparation, true);
  
      if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
	{
	  double* PseudopotentialsNMinus1 = EvaluatePseudopotentials(NbrFlux, LandauLevel - 1, layerSeparation, true);
	  for (int i = 0; i <= MaxMomentum; ++i)
	    Pseudopotentials[i] = 0.5 * (Pseudopotentials[i] + PseudopotentialsNMinus1[i]);
	  delete[] PseudopotentialsNMinus1;
	}
    }
  double* OneBodyPotentials = 0;
  if (((BooleanOption*) Manager["add-impurities"])->GetBoolean() == true)
    OneBodyPotentials = EvaluateOneBodyPotentials(NbrFlux, LandauLevel, 
						  ((SingleDoubleOption*) Manager["north-potential"])->GetDouble(),
						  ((SingleDoubleOption*) Manager["south-potential"])->GetDouble());

  if (((BooleanOption*) Manager["std-output"])->GetBoolean() == false)
    {
      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      File << "# pseudopotentials on the sphere for coulomb interaction ";
      if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
	File << " for relativistic fermions";
      File << endl
	   << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
      if (layerSeparation != 0.0)
	File << "# with finite layer separation d=" << layerSeparation << endl;
      if (Manager.GetDouble("layer-thickness") > 0.0)
	{
	  File << "# with finite layer width W=" << Manager.GetDouble("layer-thickness") << endl;
	  File << "# and type: "<<AbstractZDensityProfile::DensityProfileName(Manager.GetString("profile-type"))<<endl;
	}
      if (OneBodyPotentials != 0)
	{
	  File << "# with two impurities (V_north = " << ((SingleDoubleOption*) Manager["north-potential"])->GetDouble() 
	       << " and V_south = " << ((SingleDoubleOption*) Manager["south-potential"])->GetDouble() << ")" << endl;
	}
      if (Manager.GetDouble("add-v0")!=0.0)
	File << "# with \\delta V_0="<<Manager.GetDouble("add-v0")<<endl;
      if (Manager.GetDouble("add-v1")!=0.0)
	File << "# with \\delta V_1="<<Manager.GetDouble("add-v1")<<endl;
      File << "#" << endl
	   << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
	   << "Pseudopotentials =";
      File << " " << Pseudopotentials[0]+Manager.GetDouble("add-v0");
      File << " " << Pseudopotentials[1]+Manager.GetDouble("add-v1");
      for (int i = 2; i <= MaxMomentum; ++i)
	File << " " << Pseudopotentials[i];
      File << endl;
      if (OneBodyPotentials != 0)
	{
	  File << endl << "Onebodypotentials =";
	  for (int i = 0; i <= MaxMomentum; ++i)
	    File << " " << OneBodyPotentials[i];
	  File << endl;
	}
      if (Manager.GetBoolean("nbody"))
	{
	  int Length;
	  double *potentials = Manager.GetDoubles("nbody-potentials",Length);
	  if (potentials != 0)
	    {
	      File << endl << "NbrNBody = "<<Length-1<<endl;
	      File << "Weights =";
	      for (int i=0; i<Length; ++i)
		File <<" "<<potentials[i];
	      File <<endl;
	    }
	  else
	    {
	      File << endl << "NbrNBody = 2"<<endl;
	      File << "Weights = 0 0 0";
	    }
	}
      
      File.close();
    }
  else
    {
      cout.precision(14);
      cout << "# pseudopotentials on the sphere for coulomb interaction ";
      if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
	cout << " for relativistic fermions";
      cout << endl
	   << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
      if (layerSeparation != 0.0)
	cout << "# with finite layer separation d=" << layerSeparation << endl;
      if (Manager.GetDouble("layer-thickness") > 0.0)
	{
	  cout << "# with finite layer width W=" << Manager.GetDouble("layer-thickness") << endl;
	  cout << "# and type: "<<AbstractZDensityProfile::DensityProfileName(Manager.GetString("profile-type"))<<endl;
	}
      if (OneBodyPotentials != 0)
	{
	  cout << "# with two impurities (V_north = " << ((SingleDoubleOption*) Manager["north-potential"])->GetDouble() 
	       << " and V_south = " << ((SingleDoubleOption*) Manager["south-potential"])->GetDouble() << ")" << endl;
	}
      cout << "#" << endl
	   << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
	   << "Pseudopotentials =";
      for (int i = 0; i <= MaxMomentum; ++i)
	cout << " " << Pseudopotentials[i];
      cout << endl;
      if (OneBodyPotentials != 0)
	{
	  cout << endl << "Onebodypotentials =";
	  for (int i = 0; i <= MaxMomentum; ++i)
	    cout << " " << OneBodyPotentials[i];
	  cout << endl;
	}
    }  

  if (Manager.GetBoolean("tab"))
    {
      ofstream File;
      char *TabFileName = ReplaceExtensionToFileName(OutputFile,"dat","tab");
      File.open(TabFileName, ios::binary | ios::out);
      File.precision(14);
      File << "# pseudopotentials on the sphere for coulomb interaction ";
      if (((BooleanOption*) Manager["relativistic-fermions"])->GetBoolean() == true)
	File << " for relativistic fermions";
      File << endl
	   << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
      if (layerSeparation != 0.0)
	File << "# with finite layer separation d=" << layerSeparation << endl;
      if (Manager.GetDouble("layer-thickness") > 0.0)
	{
	  File << "# with finite layer width W=" << Manager.GetDouble("layer-thickness") << endl;
	  File << "# and type: "<<AbstractZDensityProfile::DensityProfileName(Manager.GetString("profile-type"))<<endl;
	}      
      if (Manager.GetDouble("add-v0")!=0.0)
	File << "# with \\delta V_0="<<Manager.GetDouble("add-v0")<<endl;
      if (Manager.GetDouble("add-v1")!=0.0)
	File << "# with \\delta V_1="<<Manager.GetDouble("add-v1")<<endl;
      File << "#" << endl
	   << "# m\tV_m" << endl;      
      File << "0\t" << Pseudopotentials[0]+Manager.GetDouble("add-v0")<<endl;
      File << "1\t" << Pseudopotentials[1]+Manager.GetDouble("add-v1")<<endl;
      for (int i = 2; i <= MaxMomentum; ++i)
	File << i<<"\t" << Pseudopotentials[i]<<endl;
      File.close();
      delete [] TabFileName;
    }
  
  delete[] OutputFile;
  return 0;
}


