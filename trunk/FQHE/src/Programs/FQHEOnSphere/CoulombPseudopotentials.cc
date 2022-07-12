#include "Options/Options.h"

#include "Tools/FQHESpectrum/PseudoPotentials.h"
#include "Tools/FQHESpectrum/AbstractZDensityProfile.h"

#include "Vector/RealVector.h"
 
#include "FunctionBasis/ParticleOnSphereFunctionBasis.h"
#include "FunctionBasis/ParticleOnSphereGenericLLFunctionBasis.h"

#include "GeneralTools/FilenameTools.h"


#include <iostream>
#include <fstream>
#include <cstring>
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
  (*SystemGroup) += new SingleDoubleOption ('\n', "north-potential", "potential associated to the impurity at the north pole", 0.0);
  (*SystemGroup) += new SingleDoubleOption ('\n', "south-potential", "potential associated to the impurity at the south pole", 0.0);
  (*SystemGroup) += new MultipleDoubleOption ('\n', "north-charge", "potential associated to the a charge at the given distance above the north pole", ',');
  (*SystemGroup) += new MultipleDoubleOption ('\n', "south-charge", "potential associated to the a charge at the given distance above the south pole", ',');
  
  (*SystemGroup) += new  BooleanOption ('\n', "relativistic-fermions", "assume relativistic fermions");
  (*SystemGroup) += new  BooleanOption ('\n', "graphene-bilayer", "calculate pseudopotentials for bilayer graphene (need to specify the four Landau level indices)");
  (*SystemGroup) += new SingleIntegerOption  ('\n', "l1", "Landau level index for graphene bilayer pseudopotentials, can be 0 (LLL) or 1 (1st LL)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "l2", "Landau level index for graphene bilayer pseudopotentials, can be 0 (LLL) or 1 (1st LL)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "l3", "Landau level index for graphene bilayer pseudopotentials, can be 0 (LLL) or 1 (1st LL)", 0);
  (*SystemGroup) += new SingleIntegerOption  ('\n', "l4", "Landau level index for graphene bilayer pseudopotentials, can be 0 (LLL) or 1 (1st LL)", 0);


  (*SystemGroup) += new  SingleDoubleOption ('d', "layer-separation", "assume finite layer separation",0.0);
  
  (*SystemGroup) += new  SingleDoubleOption ('t', "layer-thickness", "assume finite layer thickness",0.0);
  (*SystemGroup) += new MultipleIntegerOption ('\n', "rescale", "rescale separation and thickness for a particular number of particles, filling factor and shift (4 arguments N,p,q,Sigma)",',');
  
  (*SystemGroup) += new SingleDoubleOption ('\n', "manual-scale", "scale all terms with a given prefactor",1.0);

  
  (*SystemGroup) += new  SingleStringOption ('\n', "profile-type", "type of density-profile (1=infiniteWell, 2=Fang-Howard, 3=infiniteWellExc)","1");

  (*SystemGroup) += new  SingleStringOption ('\n', "other-profile", "type of density-profile for second layer (1=infiniteWell, 2=Fang-Howard, 3=infiniteWellExc)",NULL);

  (*SystemGroup) += new  SingleIntegerOption ('\n', "profile-points", "number of points where profile is evaluated", 200);

  (*SystemGroup) += new  SingleDoubleOption ('\n', "profile-multiplier", "multiplier for number of integration points",1.0);
  (*SystemGroup) += new  SingleDoubleOption ('\n', "tolerance", "relative tolerance to be achieved in numerical integrations",1e-8);

  (*SystemGroup) += new  BooleanOption ('\n', "no-interpolation","do not use interpolation for finite width calculations");

  (*SystemGroup) += new  BooleanOption ('n', "nbody", "add n-body potentials");  (*SystemGroup) += new  MultipleDoubleOption ('p', "nbody-potentials", "values of n-body potentials to be added (separated by ','",',');
 
  (*SystemGroup) += new SingleDoubleOption ('v', "add-v0", "add some amount to v_0",0.0);
  (*SystemGroup) += new SingleDoubleOption ('w', "add-v1", "add some amount to v_1",0.0);
  (*SystemGroup) += new SingleDoubleOption ('k', "kappa", "add three-body terms for Landau-level mixing according to Bishara-Nayak",0.0);
  (*SystemGroup) += new BooleanOption ('\n', "BN-only", "include Landau-level mixing terms, only");

  (*SystemGroup) += new BooleanOption  ('\n', "tab", "output a table with the results in columns (extension .tab)");
  (*SystemGroup) += new SingleStringOption  ('o', "output", "output file name (default is pseudopotential_coulomb_l_x_2s_y[_v_xxx_w_yyy].dat)");
  (*SystemGroup) += new BooleanOption ('\n', "std-output", "use standard output instead of an output file");

  (*MiscGroup) += new BooleanOption  ('h', "help", "display this help");

  if (Manager.ProceedOptions(argv, argc, cout) == false)
    {
      cout << "see man page for option syntax or type CoulombPseudopotentials -h" << endl;
      return -1;
    }
  if (Manager.GetBoolean("help") == true)
    {
      Manager.DisplayHelp (cout);
      return 0;
    }


  bool Rescale=false;
  double NewScale=1.0;
  int NbrParticles=1, FillingP=1, FillingQ=1, ShiftSigma=0;
  {
    int TmpLength;
    int *ScaleParams=Manager.GetIntegers("rescale",TmpLength);
    if (TmpLength>0)
      {
	if (TmpLength!=4)
	  {
	    cout << "error: --rescale requires 4 parameters: NbrParticles,FillingP,FillingQ,ShiftSigma"<<endl;
	    exit(1);
	  }
	Rescale=true;
	NbrParticles=ScaleParams[0];
	FillingP=ScaleParams[1];
	FillingQ=ScaleParams[2];
	ShiftSigma=ScaleParams[3];
	NewScale=sqrt((double)NbrParticles/((double)NbrParticles-(double)FillingP/(double)FillingQ*ShiftSigma));
	cout << "Rescaling layer separation and thickness with scale factor: "<<NewScale<<endl;
      }
  }
  
  if (Manager.GetBoolean("graphene-bilayer") == false)
    {
      int LandauLevel = Manager.GetInteger("landau-level");
      int NbrFlux = Manager.GetInteger("nbr-flux");
      double layerSeparation = Manager.GetDouble("layer-separation");
      int MaxMomentum = NbrFlux + (LandauLevel << 1);
      
      char* OutputFile;
      ((SingleDoubleOption*)Manager["add-v0"])->SetStringFormat("%g");
      ((SingleDoubleOption*)Manager["add-v1"])->SetStringFormat("%g");
      ((SingleDoubleOption*)Manager["layer-thickness"])->SetStringFormat("%g");
      ((SingleDoubleOption*)Manager["layer-separation"])->SetStringFormat("%g");
      ((SingleDoubleOption*)Manager["kappa"])->SetStringFormat("%g");
      if (Manager.GetString("output") == 0l)
	{
	  if (Manager.GetBoolean("relativistic-fermions") == true)
	    OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_relativistic_l_%landau-level%_2s_%nbr-flux%.dat");
	  else
	    {
	      char* BaseString = new char[1024];
	      if (Manager.GetBoolean("BN-only"))
		{
		  OutputFile = Manager.GetFormattedString("pseudopotential_BN-mix_l_%landau-level%_2s_%nbr-flux%_k_%kappa%.dat");
		}
	      else
		{
		  if (Manager.GetDouble("kappa")>0.0)
		    BaseString = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%_k_%kappa%");	      
		  else
		    BaseString = Manager.GetFormattedString("pseudopotential_coulomb_l_%landau-level%_2s_%nbr-flux%");
		  if (Manager.HasDoubles("north-charge")||Manager.HasDoubles("south-charge"))
		    {
		      char * TmpStr=new char[512];
		      double Charge=0.0, D=0.0;
		      int L;
		      double *TmpD = Manager.GetDoubles("north-charge",L);
		      if (L>0) {Charge=1.0;D=TmpD[0];}
		      if (L>1) {Charge=TmpD[1];}
		      strcpy(TmpStr,BaseString);
		      if (L>0) sprintf(BaseString,"%s_QN_%g_%g",TmpStr,Charge,D);
		      TmpD = Manager.GetDoubles("south-charge",L);
		      if (L>0) {Charge=1.0;D=TmpD[0];}
		      if (L>1) {Charge=TmpD[1];}
		      strcpy(TmpStr,BaseString);
		      if (L>0) sprintf(BaseString,"%s_QS_%g_%g",TmpStr,Charge,D);
		      delete [] TmpStr;
		    }	      
		  char* FormatString=new char[1024];
		  if (Manager.GetDouble("layer-thickness")>0.0)
		    {
		      if (Rescale)
			{
			  if (layerSeparation==0.0)
			    sprintf(FormatString,"%s_t_%%layer-thickness%%_%%profile-type%%_scale_N_%d_nu_%d_%d_S_%d.dat", BaseString, NbrParticles, FillingP, FillingQ, ShiftSigma);
			  else
			    sprintf(FormatString,"%s_t_%%layer-thickness%%_%%profile-type%%_scale_N_%d_nu_%d_%d_S_%d.dat", BaseString, NbrParticles, FillingP, FillingQ, ShiftSigma);
			}
		      else
			{
			  if (layerSeparation==0.0)
			    sprintf(FormatString,"%s_t_%%layer-thickness%%_%%profile-type%%.dat",BaseString);
			  else
			    sprintf(FormatString,"%s_t_%%layer-thickness%%_d_%%layer-separation%%_%%profile-type%%.dat", BaseString);
			}
		    }
		  else
		    {	      
		      if ((Manager.GetDouble("add-v0")==0.0)&&(Manager.GetDouble("add-v1")==0.0))
			{
			  if (layerSeparation==0.0)
			    sprintf(FormatString,"%s.dat",BaseString);
			  else
			    {
			      if (Rescale)
				{
				  sprintf(FormatString, "%s_d_%%layer-separation%%_scale_N_%d_nu_%d_%d_S_%d.dat", BaseString, NbrParticles, FillingP, FillingQ, ShiftSigma);
				}
			      else
				sprintf(FormatString, "%s_d_%%layer-separation%%.dat", BaseString);
			    }
			}
		      else
			{
			  if (Manager.GetDouble("add-v0")==0.0)
			    {
			      if (layerSeparation==0.0)
				sprintf(FormatString, "%s_w_%%add-v1%%.dat", BaseString);
			      else
				{
				  if (Rescale)
				    {
				      sprintf(FormatString,"%s_d_%%layer-separation%%_w_%%add-v1%%_scale_N_%d_nu_%d_%d_S_%d.dat", BaseString, NbrParticles, FillingP, FillingQ, ShiftSigma);
				    }
				  else
				    sprintf(FormatString,"%s_d_%%layer-separation%%_w_%%add-v1%%.dat", BaseString);
				}
			    }
			  else if (Manager.GetDouble("add-v1")==0.0)
			    {
			      if (layerSeparation==0.0)
				sprintf(FormatString,"%s_v_%%add-v0%%.dat", BaseString);
			      else
				{
				  if (Rescale)
				    {
				      sprintf(FormatString,"%s_d_%%layer-separation%%_w_%%add-v0%%_scale_N_%d_nu_%d_%d_S_%d.dat", BaseString, NbrParticles, FillingP, FillingQ, ShiftSigma);
				    }
				  else
				    sprintf(FormatString,"%s_d_%%layer-separation%%_v_%%add-v0%%.dat", BaseString);
				}
			    }
			  else
			    {
			      if (layerSeparation==0.0)
				sprintf(FormatString,"%s_v_%%add-v0%%_w_%%add-v1%%.dat", BaseString);
			      else
				{
				  if (Rescale)
				    {
				      sprintf(FormatString,"%s_d_%%layer-separation%%_w_%%add-v0%%_w_%%add-v1%%_scale_N_%d_nu_%d_%d_S_%d.dat", BaseString, NbrParticles, FillingP, FillingQ, ShiftSigma);
				    }
				  else
				    sprintf(FormatString,"%s_d_%%layer-separation%%_v_%%add-v0%%_w_%%add-v1%%.dat",BaseString);
				}
			    }
			}
		    }
		  OutputFile = Manager.GetFormattedString(FormatString);
		  delete [] FormatString;
		  delete [] BaseString;
		}
	    }
	}
      else
	{
	  OutputFile = new char [strlen(Manager.GetString("output")) + 1];
	  strcpy (OutputFile, Manager.GetString("output"));
	}
    
      double* Pseudopotentials;
      double* OneBodyPotentials = 0;

      if (Manager.GetBoolean("BN-only"))
	{
	  Pseudopotentials = new double [MaxMomentum+1];
	  for (int i=0; i<=MaxMomentum; ++i)
	    Pseudopotentials[i]=0.0;
	}
      else
	{
	  if (Manager.GetDouble("layer-thickness")>0.0)
	    {      
	      AbstractZDensityProfile *Profile = AbstractZDensityProfile::CreateZDensityProfile(Manager.GetString("profile-type"),Manager.GetDouble("layer-thickness")*NewScale);
	      AbstractZDensityProfile *OtherProfile=NULL;
	      if (Manager.GetString("other-profile")!=NULL)
		OtherProfile = AbstractZDensityProfile::CreateZDensityProfile(Manager.GetString("other-profile"),Manager.GetDouble("layer-thickness")*NewScale);
	      if (Manager.GetBoolean("no-interpolation"))
		Pseudopotentials = EvaluateFiniteWidthPseudoPotentialsNoInterpolation(NbrFlux, LandauLevel, Profile,
										      /*, points=200 */ Manager.GetInteger("profile-points"),
										      /*, multiplier=5.0 */ Manager.GetDouble("profile-multiplier"),
										      Manager.GetDouble("layer-separation")*NewScale, OtherProfile,
										      Manager.GetDouble("tolerance"));
	      else
		Pseudopotentials = EvaluateFiniteWidthPseudoPotentials(NbrFlux, LandauLevel, Profile,
								       /*, points=200 */ Manager.GetInteger("profile-points"),
								       /*, multiplier=5.0 */ Manager.GetDouble("profile-multiplier"),
								       Manager.GetDouble("layer-separation")*NewScale, OtherProfile,
								       Manager.GetDouble("tolerance"));
	    }
	  else
	    {
	      Pseudopotentials = EvaluatePseudopotentials(NbrFlux, LandauLevel, layerSeparation*NewScale, true);
  
	      if (Manager.GetBoolean("relativistic-fermions") == true)
		{
		  double* PseudopotentialsNMinus1 = EvaluatePseudopotentials(NbrFlux, LandauLevel - 1, layerSeparation*NewScale, true);
		  for (int i = 0; i <= MaxMomentum; ++i)
		    Pseudopotentials[i] = 0.5 * (Pseudopotentials[i] + PseudopotentialsNMinus1[i]);
		  delete[] PseudopotentialsNMinus1;
		}
	    }
	  if (Manager.GetBoolean("add-impurities") == true)
	    OneBodyPotentials = EvaluateOneBodyPotentials(NbrFlux, LandauLevel, 
							  ((SingleDoubleOption*) Manager["north-potential"])->GetDouble(),
							  ((SingleDoubleOption*) Manager["south-potential"])->GetDouble());
	  if (Manager.HasDoubles("north-charge")||Manager.HasDoubles("south-charge"))
	    {
	      double ChargeN=0.0, ChargeS=0.0, DN=0.0, DS=0.0;
	      int L;
	      double *TmpD = Manager.GetDoubles("north-charge",L);
	      if (L>0) {ChargeN=1.0;DN=TmpD[0];}
	      if (L>1) {ChargeN=TmpD[1];}
	      TmpD = Manager.GetDoubles("south-charge",L);
	      if (L>0) {ChargeS=1.0;DS=TmpD[0];}
	      if (L>1) {ChargeS=TmpD[1];}
	      double *TmpPotentials = EvaluateOneBodyCoulombPotentials(NbrFlux, LandauLevel, ChargeN, DN, ChargeS, DS);
	      if (OneBodyPotentials!=NULL)
		{
		  for (int i=0; i<=NbrFlux+2*(LandauLevel-1);++i)
		    OneBodyPotentials[i]+=TmpPotentials[i];
		  delete [] TmpPotentials;
		}
	      else
		OneBodyPotentials=TmpPotentials;
	    }
	}
      if (fabs(Manager.GetDouble("kappa"))>0.0)
	{
	  double kappa=Manager.GetDouble("kappa");
	  if (LandauLevel==0)
	    {
	      Pseudopotentials[0]+= -0.1032*kappa;
	      Pseudopotentials[1]+= -0.0053*kappa;
	      Pseudopotentials[2]+= -0.0012*kappa;
	      Pseudopotentials[3]+= -0.0004*kappa;
	      Pseudopotentials[4]+= -0.0002*kappa;
	    }
	  else if (LandauLevel==1)
	    {
	      Pseudopotentials[0]+= -0.0325*kappa;
	      Pseudopotentials[1]+= 0.0042*kappa;
	      Pseudopotentials[2]+= -0.0174*kappa;
	      Pseudopotentials[3]+= 0.0023*kappa;
	      Pseudopotentials[4]+= -0.0034*kappa;
	    }
	}

      if (Manager.GetBoolean("std-output") == false)
	{
	  ofstream File;
	  File.open(OutputFile, ios::binary | ios::out);
	  File.precision(14);
	  File << "# pseudopotentials on the sphere for coulomb interaction ";
	  if (Manager.GetBoolean("relativistic-fermions") == true)
	    File << " for relativistic fermions";
	  File << endl
	       << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
	  if (layerSeparation != 0.0)
	    {
	      if (Rescale)
		File << "# with effective finite layer separation d=" << layerSeparation << " (rescaled by "<<NewScale<<")"<<endl;
	      else
		File << "# with finite layer separation d=" << layerSeparation << endl;
	    }
	  if (Manager.GetDouble("layer-thickness") > 0.0)
	    {
	      if (Rescale)
		File << "# with finite layer at effective width W=" << Manager.GetDouble("layer-thickness") << " (rescaled by "<<NewScale<<")"<<endl;
	      else
		File << "# with finite layer width W=" << Manager.GetDouble("layer-thickness") << endl;
	      if (Manager.GetString("other-profile")!=NULL)
		File << "# for type1 = "<<AbstractZDensityProfile::DensityProfileName(Manager.GetString("profile-type"))
		     << " and type2 = "<<AbstractZDensityProfile::DensityProfileName(Manager.GetString("other-profile"))<< endl;
	      else
		File << "# and type: "<<AbstractZDensityProfile::DensityProfileName(Manager.GetString("profile-type"))<<endl;	    
	    }
	  if (OneBodyPotentials != 0)
	    {
	      if (Manager.GetBoolean("add-impurities"))
		File << "# with two impurities (V_north = " << ((SingleDoubleOption*) Manager["north-potential"])->GetDouble() 
		     << " and V_south = " << ((SingleDoubleOption*) Manager["south-potential"])->GetDouble() << ")" << endl;
	      else
		{
		  int L;
		  double *TmpD = Manager.GetDoubles("north-charge",L);
		  if (L>0)
		    {
		      if (L>1)
			File << "# with a charge (Q_north = " << TmpD[1] << " at z=" << TmpD[0] << ")" << endl;
		      else
			File << "# with a charge (Q_north = +1 at z=" << TmpD[0] << ")" << endl;
		    }
		  TmpD = Manager.GetDoubles("south-charge",L);
		  if (L>0)
		    {
		      if (L>1)
			File << "# with a charge (Q_south = " << TmpD[1] << " at z=" << TmpD[0] << ")" << endl;
		      else
			File << "# with a charge (Q_south = +1 at z=" << TmpD[0] << ")" << endl;
		    }
		}

	    }
	  if (Manager.GetDouble("add-v0")!=0.0)
	    File << "# with \\delta V_0="<<Manager.GetDouble("add-v0")<<endl;
	  if (Manager.GetDouble("add-v1")!=0.0)
	    File << "# with \\delta V_1="<<Manager.GetDouble("add-v1")<<endl;
	  double ManualScale=Manager.GetDouble("manual-scale");
	  if (fabs(ManualScale-1.0)>1e-12)
	    {
	      cout << "Attention: overall scale factor "<<ManualScale<<" applied - but filename not changed in output"<<endl;
	      File << "# Additional overall manual scale factor: "<<ManualScale<<endl;
	    }

	  File << "#" << endl
	       << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
	       << "Pseudopotentials =";

	  File << " " << ManualScale*(Pseudopotentials[0]+Manager.GetDouble("add-v0"));
	  File << " " << ManualScale*(Pseudopotentials[1]+Manager.GetDouble("add-v1"));
	  
	  
	  for (int i = 2; i <= MaxMomentum; ++i)
	    File << " " << ManualScale*Pseudopotentials[i];
	  File << endl;
	  if (OneBodyPotentials != 0)
	    {
	      File << endl << "Onebodypotentials =";
	      for (int i = 0; i <= MaxMomentum; ++i)
		File << " " << ManualScale*OneBodyPotentials[i];
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
		    File <<" "<<ManualScale*potentials[i];
		  File <<endl;
		}
	      else
		{
		  File << endl << "NbrNBody = 2"<<endl;
		  File << "Weights = 0 0 0";
		}
	    }
	  if (fabs(Manager.GetDouble("kappa"))>0.0)
	    {
	      int LandauLevel = Manager.GetInteger("landau-level");
	      double kappa=Manager.GetDouble("kappa");
	      if (LandauLevel==0)
		{
		  File << "ThreebodyPseudopotentials = 0 0 0 "<<-0.0181*kappa<<" 0 "<<0.0033*kappa<<" "<<-0.0107*kappa<<" ";
		  File <<0.0059*kappa<<" "<<-0.0048*kappa<<endl;
		}
	      else if (LandauLevel==1)
		{
		  File << "ThreebodyPseudopotentials = 0 0 0 "<<-0.0147*kappa<<" 0 "<<-0.0054*kappa<<" "<<0.0099*kappa<<" ";
		  File <<0.0005*kappa<<" "<<-0.0009*kappa<<endl;
		}
	      else
		{
		  File << "#LL mixing not known for LandauLevel="<<LandauLevel<<endl;
		}
	      File << "NormalizeThreeBody = true"<<endl;
	    }

	  File.close();
	}
      else
	{
	  cout.precision(14);
	  cout << "# pseudopotentials on the sphere for coulomb interaction ";
	  if (Manager.GetBoolean("relativistic-fermions") == true)
	    cout << " for relativistic fermions";
	  cout << endl
	       << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
	  if (layerSeparation != 0.0)
	    {
	      if (Rescale)
		cout << "# with effective finite layer separation d=" << layerSeparation << " (rescaled by "<<NewScale<<")"<<endl;
	      else
		cout << "# with finite layer separation d=" << layerSeparation << endl;
	    }
	  if (Manager.GetDouble("layer-thickness") > 0.0)
	    {
	      if (Rescale)
		cout << "# with finite layer at effective width W=" << Manager.GetDouble("layer-thickness") << " (rescaled by "<<NewScale<<")"<<endl;
	      else
		cout << "# with finite layer width W=" << Manager.GetDouble("layer-thickness") << endl;
	      cout << "# and type: "<<AbstractZDensityProfile::DensityProfileName(Manager.GetString("profile-type"))<<endl;
	    }
	  if (OneBodyPotentials != 0)
	    {
	      cout << "# with two impurities (V_north = " << ((SingleDoubleOption*) Manager["north-potential"])->GetDouble() 
		   << " and V_south = " << ((SingleDoubleOption*) Manager["south-potential"])->GetDouble() << ")" << endl;
	    }
	  // xxx cout << "
	  double ManualScale=Manager.GetDouble("manual-scale");
	  if (fabs(ManualScale-1.0)>1e-12)
	    {
	      cout << "# Additional overall manual scale factor: "<<ManualScale<<endl;
	    }
	  cout << "#" << endl
	       << "# Pseudopotentials = V_0 V_1 ..." << endl << endl
	       << "Pseudopotentials =";
	  for (int i = 0; i <= MaxMomentum; ++i)
	    cout << " " << ManualScale*Pseudopotentials[i];
	  cout << endl;
	  if (OneBodyPotentials != 0)
	    {
	      cout << endl << "Onebodypotentials =";
	      for (int i = 0; i <= MaxMomentum; ++i)
		cout << " " << ManualScale*OneBodyPotentials[i];
	      cout << endl;
	    }
	  if (Manager.GetBoolean("nbody"))
	    {
	      int Length;
	      double *potentials = Manager.GetDoubles("nbody-potentials",Length);
	      if (potentials != 0)
		{
		  cout << endl << "NbrNBody = "<<Length-1<<endl;
		  cout << "Weights =";
		  for (int i=0; i<Length; ++i)
		    cout <<" "<<ManualScale*potentials[i];
		  cout <<endl;
		}
	      else
		{
		  cout << endl << "NbrNBody = 2"<<endl;
		  cout << "Weights = 0 0 0";
		}
	    }
	  if (fabs(Manager.GetDouble("kappa"))>0.0)
	    {
	      int LandauLevel = Manager.GetInteger("landau-level");
	      double kappa=Manager.GetDouble("kappa");
	      if (LandauLevel==0)
		{
		  cout << "ThreebodyPseudopotentials = 0 0 0 "<<-0.0181*kappa<<" 0 "<<0.0033*kappa<<" "<<-0.0107*kappa<<" ";
		  cout <<0.0059*kappa<<" "<<-0.0048*kappa<<endl;
		}
	      else if (LandauLevel==1)
		{
		  cout << "ThreebodyPseudopotentials = 0 0 0 "<<-0.0147*kappa<<" 0 "<<-0.0054*kappa<<" "<<0.0099*kappa<<" ";
		  cout <<0.0005*kappa<<" "<<-0.0009*kappa<<endl;
		}
	      else
		{
		  cout << "#LL mixing not known for LandauLevel="<<LandauLevel<<endl;
		}
	      cout << "NormalizeThreeBody = true"<<endl;
	    }
	}

      if (Manager.GetBoolean("tab"))
	{
	  ofstream File;
	  char *TabFileName = ReplaceExtensionToFileName(OutputFile,"dat","tab");
	  File.open(TabFileName, ios::binary | ios::out);
	  File.precision(14);
	  File << "# pseudopotentials on the sphere for coulomb interaction ";
	  if (Manager.GetBoolean("relativistic-fermions") == true)
	    File << " for relativistic fermions";
	  File << endl
	       << "# in the Landau level N=" << LandauLevel << " for 2S=" << NbrFlux << " flux quanta" << endl;
	  if (layerSeparation != 0.0)
	    {
	      if (Rescale)
		cout << "# with effective finite layer separation d=" << layerSeparation << " (rescaled by "<<NewScale<<")"<<endl;
	      else
		File << "# with finite layer separation d=" << layerSeparation << endl;
	    }
	  if (Manager.GetDouble("layer-thickness") > 0.0)
	    {
	      if (Rescale)
		File << "# with finite layer az effective width W=" << Manager.GetDouble("layer-thickness") << " (rescaled by "<<NewScale<<")"<<endl;
	      else
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
    }
  else //Calculate pseudopotentials for graphene bilayer
    {
      int NbrFlux = Manager.GetInteger("nbr-flux");
      int LLIndex1 = Manager.GetInteger("l1");
      int LLIndex2 = Manager.GetInteger("l2");
      int LLIndex3 = Manager.GetInteger("l3");
      int LLIndex4 = Manager.GetInteger("l4");

      if ((LLIndex1<0) || (LLIndex1>1) || (LLIndex2<0) || (LLIndex2>1) || (LLIndex3<0) || (LLIndex3>1) || (LLIndex4<0) || (LLIndex4>1) )
	{
	  cout<<"Invalid indices"<<endl;
	  exit(1);
	}
      cout<<"Calculating graphene bilayer pseudopotentials for Nflux= "<<NbrFlux<<" LL indices: "<<LLIndex1<<" "<<LLIndex2<<" "<<LLIndex3<<" "<<LLIndex4<<endl;

      //-----------------------------------------------------
      //---------  prepare output ---------------------------
      //-----------------------------------------------------

      char* OutputFile;

      if (Manager.GetString("output") == 0l)
	OutputFile = Manager.GetFormattedString("pseudopotential_coulomb_llmixing_2s_%nbr-flux%_l1_%l1_l2_%l2_l3_%l3_l4_%l4.dat");
      else
	{
	  OutputFile = new char [strlen(Manager.GetString("output")) + 1];
	  strcpy (OutputFile, Manager.GetString("output"));
	}


      double* Pseudopotentials;
      int NbrPseudopotentials = 0;
      Pseudopotentials = EvaluateGrapheneBilayerPseudopotentials(NbrFlux, NbrPseudopotentials, LLIndex1, LLIndex2, LLIndex3, LLIndex4);


      ofstream File;
      File.open(OutputFile, ios::binary | ios::out);
      File.precision(14);
      File << "Pseudopotentials=";        
      for (int i = 0; i <= NbrPseudopotentials - 1; ++i)
	File << " " << Pseudopotentials[i];
      File << endl;      
      File.close();

      cout.precision(14);
  
      delete[] OutputFile;
      delete[] Pseudopotentials;
    }

  return 0;
}
