////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2002 Nicolas Regnault                  //
//                                                                            //
//                                                                            //
//                     class of FQHE wave function manager                    //
//                                                                            //
//                        last modification : 18/01/2005                      //
//                                                                            //
//                                                                            //
//    This program is free software; you can redistribute it and/or modify    //
//    it under the terms of the GNU General Public License as published by    //
//    the Free Software Foundation; either version 2 of the License, or       //
//    (at your option) any later version.                                     //
//                                                                            //
//    This program is distributed in the hope that it will be useful,         //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of          //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           //
//    GNU General Public License for more details.                            //
//                                                                            //
//    You should have received a copy of the GNU General Public License       //
//    along with this program; if not, write to the Free Software             //
//    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////


#include "config.h"
#include "Options/Options.h"

#include "Tools/FQHEWaveFunction/QHEWaveFunctionManager.h"
#include "Tools/FQHEWaveFunction/JainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/ExtendedHalperinWavefunction.h"
#include "Tools/FQHEWaveFunction/SLBSWavefunction.h"
#include "Tools/FQHEWaveFunction/SLBSWavefunction2.h"
#include "Tools/FQHEWaveFunction/SLBSVariationalState.h"
#include "Tools/FQHEWaveFunction/SLBSWavefunctionUnprojected.h"
#include "Tools/FQHEWaveFunction/JainCFFilledLevelOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/ExplicitMooreReadOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/AdvancedReadRezayiOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnSphereTwoQuasiholeWaveFunction.h"
#include "Tools/FQHEWaveFunction/UnprojectedJainCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphere2QHWaveFunction.h"
#include "Tools/FQHEWaveFunction/LaughlinOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/MooreReadOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianOnDiskWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphereWithSpinWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSpherePermanentWaveFunction.h"
#include "Tools/FQHEWaveFunction/PairedCFOnSphereSpinSingletWaveFunction.h"
#include "Tools/FQHEWaveFunction/TwoThirdSingletState.h"
#include "Tools/FQHEWaveFunction/TwoThirdUnpolarizedCF.h"
#include "Tools/FQHEWaveFunction/HundRuleCFStates.h"
#include "Tools/FQHEWaveFunction/HundRuleBilayerSinglet.h"
#include "Tools/FQHEWaveFunction/CFOnSphereWithSpinPartonTunnellingWaveFunction.h"
#include "Tools/FQHEWaveFunction/PfaffianTimesPfaffianState.h"
#include "Tools/FQHEWaveFunction/FQHESphereSymmetrizedSU2ToU1WaveFunction.h"
#include "Tools/FQHEWaveFunction/SU4HalperinOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/NASSOnSphereWaveFunction.h"
#include "Tools/FQHEWaveFunction/TPfaffianCandidateOnSphereWaveFunction.h"
#include "MathTools/RandomNumber/StdlibRandomNumberGenerator.h"
#include "MathTools/NumericalAnalysis/AntisymmetrizedComplexFunction.h"

#include <iostream>
#include <cstring>

using std::endl;
using std::cout;

// constructor
//
// geometry = id of the geometry to use

QHEWaveFunctionManager::QHEWaveFunctionManager(int geometry, int geometryOptions)
{
  this->GeometryID = geometry;
  this->Options = 0;
  this->GeometryOptions = geometryOptions;
}

// destructor
//

QHEWaveFunctionManager::~QHEWaveFunctionManager()
{
}

// add an option group containing all options related to the wave functions
//
// manager = pointer to the option manager

void QHEWaveFunctionManager::AddOptionGroup(OptionManager* manager)
{
  this->Options = manager;
  OptionGroup* WaveFunctionGroup  = new OptionGroup ("analytical wave function options");
  (*(this->Options)) += WaveFunctionGroup;
  (*WaveFunctionGroup) += new SingleStringOption  ('\n', "test-wavefunction", "name of the test wave fuction",0);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "nbr-flux", "number of quantum flux attached to each particle (laughlin and *cf only)", 1, true, 0);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "cluster-size", "number of particles per cluster (read only)", 3, true, 2);
  (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "nbr-level", "number of pseudo Landau levels (filledcf only)", 1, true, 1);
  (*WaveFunctionGroup) += new SingleStringOption  ('\n', "cf-file", "name of the file describing the composite fermion state (genericcf only)");
  if (this->GeometryID & QHEWaveFunctionManager::SphereGeometry)
    {
      // PairedCFOptions:
      (*WaveFunctionGroup) += new SingleDoubleOption  ('\n', "MR-coeff", "coefficient for Moore-Read contribution (pairedcf only)",1.0);
      (*WaveFunctionGroup) += new MultipleDoubleOption  ('\n', "pair-coeff", "sequence of pairing coefficients (pairedcf only)",'+');
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "pair-compatibility", "adopt old conventions for normalisation (pairedcf only)");
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "fermion-state", "generate a fermionic wavefunction (RR)");
      (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "CF-levels", "number of CF levels to fill (><0) (SLBS)",-2);
      (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "Jz-Value", "Total angular momentum Jz (hund only)", 0);
      (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "hund-L2", "Total angular momentum J if less than maximal (hund only)", -1);
      (*WaveFunctionGroup) += new MultipleIntegerOption  ('\n', "JM-Values", "Angular momentum J and projection on z axis J,M (paired2QH only)",',');
      (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "QHMC-iter", "Number of MC steps for internal MC (paired2QH only)", 1000);
    }
  else if (this->GeometryID & QHEWaveFunctionManager::SphereWithSpinGeometry)
    {
      // PairedCF(CB)Options:
      (*WaveFunctionGroup) += new SingleDoubleOption  ('\n', "bosons", "coefficient for boson contribution (pairedcfcb only)",1.0);
      (*WaveFunctionGroup) += new SingleDoubleOption  ('\n', "MR-coeff", "coefficient for Moore-Read contribution (pairedcf only)",1.0);
      (*WaveFunctionGroup) += new MultipleDoubleOption  ('\n', "pair-coeff", "sequence of pairing coefficients (pairedcfp only)",'+');
      (*WaveFunctionGroup) += new SingleIntegerOption  ('\n', "pair-wave", "choose pairing channel s,+p,... (pairedcf* only)",1);
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "pair-compatibility", "adopt old conventions for normalisation (pairedcf* only)");
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "jastrow-inside", "move jastrow factors inside determinant (1s only)");
      (*WaveFunctionGroup) += new MultipleIntegerOption  ('\n', "XHC", "coefficients (k,l,p,q,r,s,t,u,v,b) of extended Halperin wavefunction  (XH only)",',' ,',', "2,2,-2,0");
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "antisymmetrize", "antisymmetrize wavefunction (XH only)");
      (*WaveFunctionGroup) += new MultipleIntegerOption  ('\n', "partonShells", "number of shells for parton & effective flux: lFup,lFdown,N_phi^eff",',');
      (*WaveFunctionGroup) += new BooleanOption  ('\n', "fermion-state", "generate a fermionic wavefunction (nass)");
    }
  else 
    if (this->GeometryID & QHEWaveFunctionManager::SphereWithSU3SpinGeometry)
      {	
      }
}

// get list of all available wave functions
// 
// str = reference on the output stream

ostream& QHEWaveFunctionManager::ShowAvalaibleWaveFunctions (ostream& str)
{
  str << "list of avalaible wave functions:" << endl;
  if (this->GeometryID == QHEWaveFunctionManager::SphereGeometry)
    {
      str << "  * laughlin : Laughlin wave function" << endl;
      str << "  * pfaffian : pfaffian state wave function" << endl;
      str << "  * pfaffian2qh : pfaffian state wave function with 2 quasiholes" << endl;
      str << "  * read : Read-Rezayi state wave function" << endl;
      str << "  * RR : Read-Rezayi wave function (using blocks)" << endl;
      str << "  * SLBS : Slingerland-Bonderson wave function" << endl;
      str << "  * SLBS-Unp : Slingerland-Bonderson wave function (unprojected)" << endl;
      str << "  * filledcf : composite fermions wave function (only with filled pseudo Landau levels)" << endl;
      str << "  * genericcf : generic composite fermions wave function" << endl;            
      str << "  * unprojectedcf : generic unprojected composite fermions wave function" << endl;
      str << "  * pairedcf : paired composite fermion wave function at flux 2N-3" << endl;
      str << "  * paired2QH : paired composite fermion wave function with two quasi-holes at flux 2N-2+PC" << endl;
      str << "  * hund : composite fermion state, if half filled highest shell using Hund's rule" << endl;
      str << "  * tpfaffian : a potential candidate for the t-Pfaffian state" << endl;
    }
  else
    if (this->GeometryID == QHEWaveFunctionManager::DiskGeometry)
      {
	str << "  * laughlin : Laughlin wave function" << endl;
	str << "  * pfaffian : pfaffian state wave function" << endl;
	str << "  * read : Read-Rezayi state wave function" << endl;	
      }
    else
      if (this->GeometryID == QHEWaveFunctionManager::SphereWithSpinGeometry)
	{	  
	  str << "  * 111 : 111-state" << endl;
	  str << "  * HR : Haldane-Rezayi state" << endl;
	  str << "  * XH : Extended Halperin wavefunctions (z-z)^k (w-w)^k (z-w)^l [det((z-w)^p)]^r [per((z-w)^q)]^s Pf[1/(r_i - r_j)]^t Pf[1/(z-z)]^u Pf[1/(w-w)]^v Pf[1/(z-w)]^b" << endl;
	  str << "  * 1s : Spin-singlet state at filling one" << endl;
	  str << "  * 2-3s : Spin-singlet state at filling two-thirds" << endl;
	  str << "  * 2-3up : Spin-unpolarized composite fermion state at two-third filling" << endl;
	  str << "  * nass : non-abelian spin singlet state" << endl;
	  str << "  * pairedcf : paired composite fermion wave function at flux 2N_1-1" << endl;	  
	  str << "  * pairedcfcb : paired composite fermion wave function at flux 2N_1-1 with CB component" << endl;
	  str << "  * pairedcfs : paired composite fermion spin singlet wave function at flux 2N-2" << endl;	  
	  str << "  * hund : singlet state with each layer formed according to Hund's rule" << endl;
	  str << "  * cfparton : bilayer state with parton tunnelling" << endl;
	}
    else
      if (this->GeometryID == QHEWaveFunctionManager::SphereWithSU3SpinGeometry)
	{
	  str << "  * halperin : SU(3)_generalized Halperin state" << endl;	  
	}
  return str;
}
  
// get the wave function corresponding to the option constraints
//
// return value = pointer to the wave function (null if an error occurs)

Abstract1DComplexFunction* QHEWaveFunctionManager::GetWaveFunction()
{
  if ((*(this->Options))["test-wavefunction"] == 0)
    {
      return 0;
    }
  if (this->Options->GetString("test-wavefunction") == 0)
    {
      return 0;
    }
  if (this->GeometryID == QHEWaveFunctionManager::SphereGeometry)
    {
      if (strcmp (this->Options->GetString("test-wavefunction"), "laughlin") == 0)
	{
	  return new LaughlinOnSphereWaveFunction(this->Options->GetInteger("nbr-particles"), 
						  this->Options->GetInteger("nbr-flux") + 1);
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "pfaffian") == 0)
	{
	  return new PfaffianOnSphereWaveFunction(this->Options->GetInteger("nbr-particles"), this->Options->GetBoolean("fermion-state"));
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "tpfaffian") == 0)
	{
	  return new TPfaffianCandidateOnSphereWaveFunction(this->Options->GetInteger("nbr-particles"));
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "pfaffian2qh") == 0)
	{
	  return new PfaffianOnSphereTwoQuasiholeWaveFunction(this->Options->GetInteger("nbr-particles"), 0.0, 0.0, M_PI, 0.0, this->Options->GetBoolean("fermion-state"));
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "read") == 0)
	{
	  return new MooreReadOnSphereWaveFunction(this->Options->GetInteger("nbr-particles"), 
						   this->Options->GetInteger("cluster-size"));
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "RR") == 0)
	{
	  int N=this->Options->GetInteger("nbr-particles");
	  int k=this->Options->GetInteger("cluster-size");
	  if (N%k != 0)
	    {
	      cout << "Number of particles has to be a multiple of the cluster-size!"<<endl;
	      exit(1);
	    }
	  return new AdvancedReadRezayiOnSphereWaveFunction(N/k, k, this->Options->GetBoolean("fermion-state"));
	}
      if ((strcmp (this->Options->GetString("test-wavefunction"), "SLBSV") == 0))
	{
	  int N= this->Options->GetInteger("nbr-particles");
	  int LL;
	  double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	  if (Coefficients==NULL)
	    {
	      Coefficients = new double[1];
	      Coefficients[0]=0.0;
	      LL=1;
	    }
	  int NbrLevels= this->Options->GetInteger("CF-levels");
	  double MR =this->Options->GetDouble("MR-coeff");
	  SLBSVariationalState* rst = new SLBSVariationalState(N, LL, MR, Coefficients, NbrLevels);
	  rst->AdaptAverageMCNorm();
	  delete [] Coefficients;
	  return rst;
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "SLBS2") == 0)
	{
	  int N=this->Options->GetInteger("nbr-particles");
	  SLBSWavefunction2 *rst = new SLBSWavefunction2(N);
	  rst->AdaptAverageMCNorm();
	  return rst;
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "SLBS-Unp") == 0)
	{
	  int N=this->Options->GetInteger("nbr-particles");
	  SLBSWavefunctionUnprojected *rst = new SLBSWavefunctionUnprojected(N);
	  rst->AdaptAverageMCNorm();
	  return rst;
	}      
      if (strcmp (this->Options->GetString("test-wavefunction"), "SLBS") == 0)
	{
	  int N=this->Options->GetInteger("nbr-particles");
	  int levels=this->Options->GetInteger("CF-levels");
	  SLBSWavefunction *rst = new SLBSWavefunction(N, levels);
	  rst->AdaptAverageMCNorm();
	  return rst;
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "BS38") == 0)
	{
	  int N=this->Options->GetInteger("nbr-particles");
	  SLBSWavefunction *rst = new SLBSWavefunction(N, -3);
	  rst->AdaptAverageMCNorm();
	  return rst;
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "ExplicitMR") == 0)
	{
	  int N=this->Options->GetInteger("nbr-particles");
	  int ClusterSize = this->Options->GetInteger("cluster-size");
	  if (N%ClusterSize != 0)
	    {
	      cout << "Number of particles has to be a multiple of the Cluster-size!"<<endl;
	      exit(1);
	    }
	  return new ExplicitMooreReadOnSphereWaveFunction(N/ClusterSize,ClusterSize,this->Options->GetBoolean("fermion-state"));
	}
      if (strcmp (this->Options->GetString("test-wavefunction"), "filledcf") == 0)
	{
	  return new JainCFFilledLevelOnSphereWaveFunction(this->Options->GetInteger("nbr-particles"), 
							   this->Options->GetInteger("nbr-level"),
							   this->Options->GetInteger("nbr-flux"));
	}
      if ((strcmp (this->Options->GetString("test-wavefunction"), "genericcf") == 0) && ((*(this->Options))["cf-file"] != 0))
	{	  
	  return new JainCFOnSphereWaveFunction(this->Options->GetString("cf-file"));
	}
      if ((strcmp (this->Options->GetString("test-wavefunction"), "unprojectedcf") == 0) && ((*(this->Options))["cf-file"] != 0))
	{	  
	  return new UnprojectedJainCFOnSphereWaveFunction(this->Options->GetString("cf-file"));
	}
      if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedcf") == 0))
	{
	  int N= this->Options->GetInteger("nbr-particles");
	  int LL;
	  double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	  if (Coefficients==NULL)
	    {
	      Coefficients = new double[1];
	      Coefficients[0]=0.0;
	      LL=1;
	    }
	  double MR =this->Options->GetDouble("MR-coeff");
	  bool conventions = this->Options->GetBoolean("pair-compatibility");
	  int attached = this->Options->GetInteger("nbr-flux");
	  if (attached % 2 != 0)
	    {
	      cout << "Attention, changing number of flux attached to 2 - please indicate the requested value with --nbr-flux"<<endl;
	      attached = 2;
	    }
	  PairedCFOnSphereWaveFunction* rst = new PairedCFOnSphereWaveFunction(N, LL, -1, MR,
									       Coefficients, conventions, attached);
	  rst->AdaptAverageMCNorm();
	  delete [] Coefficients;
	  return rst;
	}
      if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedQH") == 0))
	{
	  int N= this->Options->GetInteger("nbr-particles");
	  int J;
	  int M=0;
	  int i;
	  int *QuantumNumbers = this->Options->GetIntegers("JM-Values",i);
	  if (i>0) J=QuantumNumbers[0];
	  else J=N/2;
	  if (i>1) M=QuantumNumbers[1];
	  int LL;	  
	  double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	  if (Coefficients==NULL)
	    {
	      Coefficients = new double[1];
	      Coefficients[0]=0.0;
	      LL=1;
	    }
	  double MR =this->Options->GetDouble("MR-coeff");
	  int QHMCSteps = this->Options->GetInteger("QHMC-iter");
	  bool conventions = this->Options->GetBoolean("pair-compatibility");
	  PairedCFOnSphere2QHWaveFunction* rst = new PairedCFOnSphere2QHWaveFunction(N, LL, -1, J, M, MR, Coefficients, QHMCSteps, conventions, 2);
	  rst->AdaptAverageMCNorm();
	  delete [] Coefficients;
	  return rst;
	}
      if ((strcmp (this->Options->GetString("test-wavefunction"), "hund") == 0))
	{
	  int N = this->Options->GetInteger("nbr-particles");
	  int JastrowP = this->Options->GetInteger("nbr-flux")/2;
	  int RelevantLzMax = this->Options->GetInteger("lzmax");
	  if ((*this->Options)["product-state"]!=NULL)
	    if (this->Options->GetBoolean("product-state"))
	      {
		if (this->Options->GetInteger("second-lzmax")==0)
		  {
		    cout << "When using an analytic trial state for the second factor in a product state, please indicate --second-lzmax!"<<endl;
		    RelevantLzMax = this->Options->GetInteger("lzmax");
		    cout << "Defaulting to total lzmax = "<< RelevantLzMax << endl;
		  }
		else
		  RelevantLzMax = this->Options->GetInteger("second-lzmax");
	    }
	  int effectiveFlux = RelevantLzMax-2*JastrowP*(N-1);
	  int overrideL = this->Options->GetInteger("hund-L2");
	  if (JastrowP==0)
	    {
	      cout << "To obtain CF's, at least two flux need to be attached. Try:  --nbr-flux 2"<<endl;
	      exit(1);
	    }
	  cout << "Effective flux for CF's: "<<effectiveFlux<<endl;
	  HundRuleCFStates* rst = new HundRuleCFStates(N, effectiveFlux, JastrowP, overrideL);
	  rst->SelectMValue(this->Options->GetInteger("Jz-Value"));
	  rst->AdaptAverageMCNorm();
	  return rst;
	}
      return 0;
    }
  else
    if (this->GeometryID == QHEWaveFunctionManager::DiskGeometry)
      {
	if (strcmp (this->Options->GetString("test-wavefunction"), "laughlin") == 0)
	  {
	    return new LaughlinOnDiskWaveFunction(this->Options->GetInteger("nbr-particles"), 
						  this->Options->GetInteger("nbr-flux") + 1, 1.0, ((this->GeometryOptions & DiskWithBackground) != 0));
	  }
	if (strcmp (this->Options->GetString("test-wavefunction"), "pfaffian") == 0)
	  {
	    if (this->GeometryOptions & DiskWithBackground)
	      {
		cout << "Warning: exponential background not implemented for Pfaffian" << endl;
	      }
	    return new PfaffianOnDiskWaveFunction(this->Options->GetInteger("nbr-particles"));
	  }
	if (strcmp (this->Options->GetString("test-wavefunction"), "read") == 0)
	  {
	    if (this->GeometryOptions & DiskWithBackground)
	      {
		cout << "Warning: exponential background not implemented for Moore Read" << endl;
	      }
	    return new MooreReadOnDiskWaveFunction(this->Options->GetInteger("nbr-particles"), 
						   this->Options->GetInteger("cluster-size"));
	  }
	return 0;
      }
    else
      if (this->GeometryID == QHEWaveFunctionManager::SphereWithSpinGeometry)
	{
	  if ((strcmp ( this->Options->GetString("test-wavefunction"), "pairedcf") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "Paired CF bilayer states require equal population of both (pseudo-)spin species and even N!" << endl;
		  exit(1);
		}
	      int LL;
	      double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	      if (Coefficients==NULL)
		{
		  Coefficients = new double[1];
		  Coefficients[0]=0.0;
		  LL=1;
		}
	      bool conventions = Options->GetBoolean("pair-compatibility");
	      int pairWave = this->Options->GetInteger("pair-wave");
	      cout << "attention: (hardwired) testing conjugate down spin mode"<<endl;
	      PairedCFOnSphereWithSpinWaveFunction* rst = new PairedCFOnSphereWithSpinWaveFunction(N, LL, pairWave, false, 0.0, Coefficients, conventions, 2, true);
	      rst->AdaptAverageMCNorm();
	      delete [] Coefficients;
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedcfcb") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "Paired CF-CB bilayer states require equal population of both (pseudo-)spin species and even N!" << endl;
		  exit(1);
		}
	      int LL;
	      double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	      if (Coefficients==NULL)
		{
		  Coefficients = new double[1];
		  Coefficients[0]=0.0;
		  LL=1;
		}
	      double BC = this->Options->GetDouble("bosons");
	      bool conventions = this->Options->GetBoolean("pair-compatibility");
	      PairedCFOnSphereWithSpinWaveFunction* rst = new PairedCFOnSphereWithSpinWaveFunction(N, LL, 1, true, BC, Coefficients, conventions, 2);
	      rst->AdaptAverageMCNorm();
	      delete [] Coefficients;
	      return rst;
	    }
	  if ((strcmp ( this->Options->GetString("test-wavefunction"), "pairedcfs") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "Paired spin singlet CF states require equal population of both spin-species and even N!" << endl;
		  exit(1);
		}
	      int LL;
	      double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	      if (Coefficients==NULL)
		{
		  Coefficients = new double[1];
		  Coefficients[0]=0.0;
		  LL=1;
		}
	      int pairWave = this->Options->GetInteger("pair-wave");
	      PairedCFOnSphereSpinSingletWaveFunction* rst = new PairedCFOnSphereSpinSingletWaveFunction(N, LL, pairWave, Coefficients, 2);	      
	      rst->AdaptAverageMCNorm();
	      delete [] Coefficients;
	      return rst;
	    }
	  if ((strcmp ( this->Options->GetString("test-wavefunction"), "pairedcfp") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "Paired CF permanent states require equal population of both spin-species and even N!" << endl;
		  exit(1);
		}
	      int LL;
	      double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	      if (Coefficients==NULL)
		{
		  Coefficients = new double[1];
		  Coefficients[0]=0.0;
		  LL=1;
		}
	      double MR =this->Options->GetDouble("MR-coeff");
	      bool conventions = this->Options->GetBoolean("pair-compatibility");
	      int attached = this->Options->GetInteger("nbr-flux");
	      if (attached % 2 != 0)
		{
		  cout << "Attention, changing number of flux attached to 2 - please indicate the requested value with --nbr-flux"<<endl;
		  attached = 2;
		}
	      int pairWave = this->Options->GetInteger("pair-wave");
	      PairedCFOnSpherePermanentWaveFunction* rst = new PairedCFOnSpherePermanentWaveFunction(N, LL,
										 pairWave, MR,
										 Coefficients, conventions, attached);	      	      rst->AdaptAverageMCNorm();
	      delete [] Coefficients;
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "hund") == 0))
	    {
	      int N= Options->GetInteger("nbr-particles");
	      int n=N/2;
	      HundRuleBilayerSinglet* rst = new HundRuleBilayerSinglet(n);
	      rst->AdaptAverageMCNorm();
	      return rst;	      
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "cfparton") == 0))
	    {
	      int N= Options->GetInteger("nbr-particles");
	      int attached = this->Options->GetInteger("nbr-flux");
	      if (attached % 2 != 0)
		{
		  cout << "Attention, changing number of flux attached to 2 - please indicate the requested value with --nbr-flux"<<endl;
		  attached = 2;
		}
	      int Length;
	      int *PartonParameters = this->Options->GetIntegers("partonShells",Length);
	      int EffectiveFlux=1;
	      int NbrUp=-1;
	      if ((Length<2)||(Length>4))
		{
		  cout << "Need parameters for Parton wavefunction: --partonShells LFBonding,LFAntibonding[,EffectiveFlux][,NbrUp]"<<endl;
		  exit(-1);
		}
	      if (Length>2)
		EffectiveFlux = PartonParameters[2];
	      if (Length>3)
		NbrUp = PartonParameters[3];
	      CFOnSphereWithSpinPartonTunnellingWaveFunction* rst = new CFOnSphereWithSpinPartonTunnellingWaveFunction
		(N, PartonParameters[0], PartonParameters[1], EffectiveFlux, attached, NbrUp);
	      rst->AdaptAverageMCNorm();
	      return rst;	      
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "111-old") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      int Sz= Options->GetInteger("SzTotal");
	      if ((N&1)||(Sz!=0))
		{
		  cout << "For now, the implementation of the 111 state requires Sz=0 and even N!" << endl;
		  exit(1);
		}
	      double* Coefficients = new double[1];
	      Coefficients[0]=0.0;
	      bool conventions = this->Options->GetBoolean("pair-compatibility");
	      PairedCFOnSphereWithSpinWaveFunction* rst = new PairedCFOnSphereWithSpinWaveFunction(N, 1, 1, true, 1.0, Coefficients, conventions, 2);
	      rst->AdaptAverageMCNorm();
	      delete[] Coefficients;
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "111") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, 1, 1, 0);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "HR") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, 2, 2, -2);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "XH") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      int length;
	      int *Params = Options->GetIntegers("XHC",length);
	      cout <<"Read "<<length<<" coefficients for XHC"<<endl;
	      int K=1, L=1, P=0, Q=0, R=1, S=1, T=0, U=0, V=0, B=0;
	      if (length>0) K = Params[0];
	      if (length>1) L = Params[1];
	      if (length>2) P = Params[2];
	      if (length>3) Q = Params[3];
	      if (length>4) R = Params[4];
	      if (length>5) S = Params[5];
	      if (length>6) T = Params[6];
// warning: merge from testing_zr (revision 1159)
// old version is
//	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, K, L, P, Q, R, S, T);
// new version is 
	      if (length>7) U = Params[7];
	      if (length>8) V = Params[8];
	      if (length>9) B = Params[9];
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, K, L, P, Q, R, S, T, U, V, B);
	      delete [] Params;

// end warning
	      rst->AdaptAverageMCNorm();
	      if (Options->GetBoolean("antisymmetrize"))
		{
		  AntisymmetrizedComplexFunction* rst2 = new AntisymmetrizedComplexFunction(rst,N,2);		  
		  return rst2;
		}
	      else
		return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "1s") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      bool inside = Options->GetBoolean("jastrow-inside");
	      ExtendedHalperinWavefunction* rst = new ExtendedHalperinWavefunction(N, 0, 2, -2, 0, 1, 1, 0, 0, 0, 0, inside);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "2-3s") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      TwoThirdSingletState* rst = new TwoThirdSingletState(N);
	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "2-3up") == 0))
	    {
	      int N= this->Options->GetInteger("nbr-particles");
	      TwoThirdUnpolarizedCF* rst = new TwoThirdUnpolarizedCF(N);
	      rst->AdaptAverageMCNorm();
              return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "nass") == 0))
	    {
	      int N = this->Options->GetInteger("nbr-particles");
	      Abstract1DComplexFunctionOnSphere* BaseFunction = new SU4HalperinOnSphereWaveFunction(N >> 2, N >> 2, N >> 2, N >> 2, 
												    2, 2, 2, 2, 
												    1, 0, 0);
	      FQHESphereSymmetrizedSU2ToU1WaveFunction* rst = new FQHESphereSymmetrizedSU2ToU1WaveFunction (N, N >> 1, BaseFunction, true);
	      //	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "nass2") == 0))
	    {
	      int N = this->Options->GetInteger("nbr-particles");
	      bool Fermions = this->Options->GetBoolean("fermion-state");
	      if (!Fermions)
		cout << "Attention, using bosonic statistics for NASS state"<<endl;
	      Abstract1DComplexFunctionOnSphere* rst = new NASSOnSphereWaveFunction(N >> 2, Fermions);

	      //	      rst->AdaptAverageMCNorm();
	      return rst;
	    }
	  if ((strcmp (this->Options->GetString("test-wavefunction"), "pfaff2") == 0))
	    {
	      int N = this->Options->GetInteger("nbr-particles");
	      int LL;
	      bool conventions = Options->GetBoolean("pair-compatibility");
	      double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
	      double MR =this->Options->GetDouble("MR-coeff");
	      if (Coefficients==NULL)
		{
		  Coefficients = new double[1];
		  Coefficients[0]=0.0;
		  LL=1;
		}
	      PfaffianTimesPfaffianState* rst = new PfaffianTimesPfaffianState(N, LL, Coefficients, MR, conventions);
	      //	      rst->AdaptAverageMCNorm();
	      return rst;
	    }

	}
    else
      if (this->GeometryID == QHEWaveFunctionManager::SphereWithSU3SpinGeometry)
	{
	}
  return 0;
}

char* QHEWaveFunctionManager::GetDescription()
{
  if ((*(this->Options))["test-wavefunction"] == 0)
    {
      return 0;
    }
  if (this->Options->GetString("test-wavefunction") == 0)
    return 0;
  char * buffer = new char[1000];
  sprintf(buffer,"%s N=%ld",this->Options->GetString("test-wavefunction"), this->Options->GetInteger("nbr-particles"));
  if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedcf") == 0))
    {
      int LL;
      double *Coefficients = this->Options->GetDoubles("pair-coeff",LL);
      if (Coefficients==NULL)
	{
	  if (this->GeometryID & QHEWaveFunctionManager::SphereWithSpinGeometry)
	    sprintf(buffer,"%s, B: %g, c: 0",buffer, this->Options->GetDouble("bosons"));
	  else
	    sprintf(buffer,"%s, MR: %g, c: 0",buffer, this->Options->GetDouble("MR-coeff"));
	}
      else
	{
	  if (this->GeometryID & QHEWaveFunctionManager::SphereWithSpinGeometry)
	    sprintf(buffer,"%s, B: %g, c: %g",buffer, this->Options->GetDouble("bosons"), Coefficients[0]);
	  else
	    sprintf(buffer,"%s, MR: %g, c: %g",buffer, this->Options->GetDouble("MR-coeff"), Coefficients[0]);
	  for (int i=1; i<LL; ++i)
	    sprintf(buffer,"%s+%g",buffer, Coefficients[i]);
	  if(this->Options->GetBoolean("pair-compatibility"))
	    sprintf(buffer,"%s (c)",buffer);
	}
    }
  char *rst = new char[strlen(buffer)+1];
  strcpy(rst,buffer);
  delete [] buffer;
  return rst;
}


int QHEWaveFunctionManager::GetWaveFunctionType()
{
  if ((*(this->Options))["test-wavefunction"] == 0)
    return QHEWaveFunctionManager::InvalidWaveFunction;
  if (strcmp (this->Options->GetString("test-wavefunction"), "laughlin") == 0)
    return QHEWaveFunctionManager::Laughlin;
  if (strcmp (this->Options->GetString("test-wavefunction"), "pfaffian") == 0)
    return QHEWaveFunctionManager::Pfaffian;
  if (strcmp (this->Options->GetString("test-wavefunction"), "pfaffian2qh") == 0)
    return QHEWaveFunctionManager::Pfaffian2QH;
  if (strcmp (this->Options->GetString("test-wavefunction"), "read") == 0)
    return QHEWaveFunctionManager::ReadRezayi;
  if (strcmp (this->Options->GetString("test-wavefunction"), "filledcf") == 0)
    return QHEWaveFunctionManager::FilledCF;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "genericcf") == 0) && ((*(this->Options))["cf-file"] != 0))
    return QHEWaveFunctionManager::GenericCF;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "unprojectedcf") == 0) && ((*(this->Options))["cf-file"] != 0))
    return QHEWaveFunctionManager::UnprojectedCF;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedcf") == 0))
    return QHEWaveFunctionManager::PairedCF;  
  if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedcfcb") == 0))
    return QHEWaveFunctionManager::PairedCFCB;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "111") == 0))
    return QHEWaveFunctionManager::OneOneOne;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedcfs") == 0))
    return QHEWaveFunctionManager::PairedCFS;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "pairedcfp") == 0))
    return QHEWaveFunctionManager::PairedCFP;  
  if ((strcmp (this->Options->GetString("test-wavefunction"), "HR") == 0))
    return QHEWaveFunctionManager::HaldaneRezayi;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "XH") == 0))
    return QHEWaveFunctionManager::ExtendedHalperin;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "1s") == 0))
    return QHEWaveFunctionManager::OneS;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "2-3s") == 0))
    return QHEWaveFunctionManager::TwoThirdsS;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "2-3up") == 0))
    return QHEWaveFunctionManager::TwoThirdsUnpolarized;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "hund") == 0))
    return QHEWaveFunctionManager::HundRuleSinglet;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "halperin") == 0))
    return QHEWaveFunctionManager::Halperin;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "SLBSV") == 0))
    return QHEWaveFunctionManager::SLBSV;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "SLBS") == 0))
    return QHEWaveFunctionManager::SLBS;
  if ((strcmp (this->Options->GetString("test-wavefunction"), "pfaff2") == 0))
    return QHEWaveFunctionManager::Pfaff2;
  return QHEWaveFunctionManager::InvalidWaveFunction;
}
