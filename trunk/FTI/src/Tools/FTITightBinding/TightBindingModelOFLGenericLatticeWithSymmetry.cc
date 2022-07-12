////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2012 Nicolas Regnault                  //
//                                                                            //
//                        class author: Gunnar Möller                         //
//                                                                            //
//  class of tight binding model for the square lattice with homogeneous flux //
//                                                                            //
//                        last modification : 08/05/2012                      //
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
#include "Tools/FTITightBinding/TightBindingModelOFLGenericLatticeWithSymmetry.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cmath>
using std::cout;
using std::endl;


// default constructor
//
// nbrPoints1 = number of k-points in the 1-direction
// nbrPoints2 = number of k-points in the 2-direction
// gamma1 = boundary condition twisting angle along x
// gamma2 = boundary condition twisting angle along y
// architecture = pointer to the architecture
// cutOffMode = 0: Circular, 1: Square, 2: Periodic
// cutOffMomentum = maximum value of momenta considered (circular: norm)
// nMax1 = maximum unit cell index in symmetry enlarged k1 direction
// nMax2 = maximum unit cell index in symmetry enlarged k2 direction
// latticeDepth = parameter for the depth of the optical lattice in recoil energies
// nbrBandsToKeep = number of bands that should be calculated / kept for outside processing (negative = keep all)
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
TightBindingModelOFLGenericLatticeWithSymmetry::TightBindingModelOFLGenericLatticeWithSymmetry(int nbrPoints1, int nbrPoints2, double gamma1, double gamma2, AbstractArchitecture* architecture, CutOffModes cutOffMode, double cutOffMomentum, int nMax1, int nMax2, double latticeDepth, int nbrBandsToKeep, double symmetryThreshold, bool storeOneBodyMatrices)
{
  if (TightBindingModelOFLGenericLatticeWithSymmetry::Options==NULL)
    {
      cout << "Define the OptionManager, first, before creating any TightBindingModelOFLGenericLatticeWithSymmetry"<<endl;
      exit(1);
    }

  if (this->Options->GetString("lattice-definition") == NULL)
    {
      cout << "option --lattice-definition is required."<<endl;
      exit(1);
    }

  if (IsFile(this->Options->GetString("lattice-definition")) == false)
    {
      cout << "Could not open file: "<<this->Options->GetString("lattice-definition")<<"."<<endl;
      exit(1);
    }


  // default fields that are insensitive to the configuration file
  this->NbrSiteX= (this->NbrPoints1 = nbrPoints1);
  this->NbrSiteY= (this->NbrPoints2 = nbrPoints2);
  // set periodicity mode
  this->Mode = cutOffMode;
  this->CutOffMomentum = cutOffMomentum;
  this->LatticeDepth = latticeDepth;

  //this->DeltaK1 = 2.0*M_PI/(double)this->NbrPoints1;  // note unusual definition with respect to other models
  //this->DeltaK2 = 2.0*M_PI/(double)this->NbrPoints2;
  this->DeltaK1 = 1.0/(double)this->NbrPoints1;  // note unusual definition with respect to other models
  this->DeltaK2 = 1.0/(double)this->NbrPoints2;


  this->GammaX = gamma1;
  this->GammaY = gamma2;
  this->NbrStatePerBand = this->NbrPoints1 * this->NbrPoints2;

  this->Architecture = architecture;

  
  ConfigurationParser LatticeDefinition;
  if (LatticeDefinition.Parse(this->Options->GetString("lattice-definition")) == false)
    {
      LatticeDefinition.DumpErrors(cout) << endl;
      exit(-1);
    }
  if (LatticeDefinition["Descriptor"] == NULL)
    {
      cout << "Attention, 'Descriptor' is not defined, unnamed lattice geometry!" << endl;
      Descriptor = new char[10];
      sprintf(Descriptor,"unnamed");
    }
  else
    {
      this->Descriptor = new char[strlen(LatticeDefinition["Descriptor"])+1];
      strcpy(this->Descriptor, LatticeDefinition["Descriptor"]);
    }
  if ((LatticeDefinition.GetAsSingleInteger("NbrSubLattices", this->NbrSubLattices) == false) || (this->NbrSubLattices <= 0))
    {
      cout << "NbrSubLattices is not defined or has invalid value" << endl;
      exit(-1);
    }

  int NbrComponents;
  double *Components;
  int Dimension = 2;
  LatticeVector1.ResizeAndClean(Dimension);
  if (LatticeDefinition.GetAsDoubleArray("LatticeVector0", ',', Components, NbrComponents) == false)
    {
      cout << "error while parsing LatticeVector0 in " << this->Options->GetString("lattice-definition") << endl;
      exit(-1);
    }
  if (Dimension!=NbrComponents)
    {
      cout << "Lattice Vectors need to have dimension "<<Dimension<<"!"<<endl;
      exit(-1);
    }
  for (int j=0; j<Dimension; ++j)
    LatticeVector1[j]=Components[j];
  delete [] Components;

  LatticeVector2.ResizeAndClean(2);
  if (LatticeDefinition.GetAsDoubleArray("LatticeVector1", ',', Components, NbrComponents) == false)
    {
      cout << "error while parsing LatticeVector1 in " << this->Options->GetString("lattice-definition") << endl;
      exit(-1);
    }
  if (Dimension!=NbrComponents)
    {
      cout << "Lattice Vectors need to have dimension "<<Dimension<<"!"<<endl;
      exit(-1);
    }
  for (int j=0; j<Dimension; ++j)
    LatticeVector2[j]=Components[j];
  delete [] Components;
  
  
  char *FieldName = new char[255];

  if ((LatticeDefinition.GetAsSingleInteger("NbrSubLatticeFlavours", this->NbrSubLatticeFlavours) == false))
    {
      cout << "NbrSubLatticeFlavours is not defined: assuming all sublattices are distinct" << endl;
      this->NbrSubLatticeFlavours=-1;
    }
   
 if ((LatticeDefinition.GetAsSingleInteger("SymmetryMultiplier1", this->SymmetryMultiplier1) == false))
    {
      cout << "SymmetryMultiplier1 is not defined: no symmetry along k1" << endl;
      this->SymmetryMultiplier1=1;
    }

 if ((LatticeDefinition.GetAsSingleInteger("SymmetryMultiplier2", this->SymmetryMultiplier2) == false))
    {
      cout << "SymmetryMultiplier2 is not defined: no symmetry along k2" << endl;
      this->SymmetryMultiplier2=1;
    }
 
 if ((LatticeDefinition.GetAsSingleInteger("FlavourOffset1", this->FlavourOffset1) == false))
   {
     cout << "FlavourOffset1 is not defined: taking generic value NbrSubLattices" << endl;
     this->FlavourOffset1=NbrSubLattices;
   }

 if ((LatticeDefinition.GetAsSingleInteger("FlavourOffset2", this->FlavourOffset2) == false))
   {
     cout << "FlavourOffset2 is not defined: taking generic value SymmetryMultiplier1*NbrSubLattices" << endl;
     this->FlavourOffset2=SymmetryMultiplier1*NbrSubLattices;
   } 
 // assign array for sublattice flavours of model
 this->ExtNbrSubLattices = NbrSubLattices * SymmetryMultiplier1 * SymmetryMultiplier2;
 this->SubLatticeFlavours = new int[ExtNbrSubLattices];
 if (NbrSubLatticeFlavours == -1)
   NbrSubLatticeFlavours = ExtNbrSubLattices;

  
  // Parsing Reciprocal Sublattice Vectors
  SubLatticeVectors.Resize(Dimension, NbrSubLattices);
  for (int i=0; i<NbrSubLattices; ++i)
    {
      sprintf(FieldName,"SubLatticeVector%d",i);
      if (LatticeDefinition.GetAsDoubleArray(FieldName, ',', Components, NbrComponents) == false)
	{
	  cout << "error while parsing "<<FieldName<< " in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);     
	}      
      if (Dimension>NbrComponents)
	{
	  cout << "SubLattice Vectors need to have the dimension of the lattice!"<<endl;
	  exit(-1);
	}
      for (int j=0; j<Dimension; ++j)
	{
	  SubLatticeVectors[i][j]=Components[j];
	}
      if (Dimension < NbrComponents)
	{
          SubLatticeFlavours[i] = Components[Dimension];
	  cout << "read sublattice flavour for sublattice i"<<endl;
        }
      else
	{
          SubLatticeFlavours[i]=i % NbrSubLatticeFlavours;
	  cout << "Assigned SubLatticeFlavours["<<i<<"]="<<SubLatticeFlavours[i]<<endl;
        }
      delete [] Components;
    }

  // assign symmetry related sublattice indices
  for (int n1=0; n1<SymmetryMultiplier1; ++n1)
    for (int n2=0; n2<SymmetryMultiplier2; ++n2)
      for (int i=0; i<NbrSubLattices; ++i)
	{
	  if ( (n1>0) || (n2>0))
	    {
	      SubLatticeFlavours[LinearizedSublatticeIndex(i,n1,n2)] = (SubLatticeFlavours[i] + n1 * FlavourOffset1 + n2 * FlavourOffset2) % NbrSubLatticeFlavours;
	      cout << "Assigned SubLatticeFlavours["<<LinearizedSublatticeIndex(i,n1,n2)<<"]="<<SubLatticeFlavours[LinearizedSublatticeIndex(i,n1,n2)]<<endl;
	    }
        }

  // Parsing Lattice Parameters
  if (LatticeDefinition.GetAsSingleInteger ("NbrExtParameters", this->NbrExtParameters)==false)
    {
      this->NbrExtParameters=1;
      this->ExtParameters=new double[1];
      this->ExtParameters[0]=1.0;
    }
  else
    {
      this->ExtParameters=new double[NbrExtParameters];
      int TmpDimension = 0;
      double *TmpParameters = NULL;
      if ( (*TightBindingModelOFLGenericLatticeWithSymmetry::Options)["external-parameters"] != NULL)
	TmpParameters = TightBindingModelOFLGenericLatticeWithSymmetry::Options->GetDoubles("external-parameters",TmpDimension);
      for (int i=0; i<TmpDimension; ++i)
	this->ExtParameters[i]=TmpParameters[i];
      if (TmpDimension>0)
	{
	  cout << "Read system parameters: "<<this->ExtParameters[0];
	  for (int i=1; i<TmpDimension; ++i)
	    cout <<" "<<this->ExtParameters[i];
	  cout << endl;
	}
      if (TmpDimension<NbrExtParameters)
	{
	  if (TmpDimension>0)
	    cout << "Warning: reading only partial parameter list from command line"<<endl;
	  for (int i=TmpDimension; i<NbrExtParameters; ++i)
	    {
	      sprintf(FieldName,"ExtParameter%d",i);
	      if (LatticeDefinition.GetAsSingleDouble(FieldName, this->ExtParameters[i])==false)
		{
		  cout << "Require field "<<FieldName<<" to be defined, or given on command-line via --external-parameters"<<endl;
		  exit(1);
		}
	      else
		{
		  cout << "Read lattice parameter "<<FieldName<<"="<<this->ExtParameters[i]<<endl;
		}
	    }
	}
      if (TmpDimension>0)
	delete [] TmpParameters;

      // check for parameter labels in lattice definition
      char **ParameterLabels=NULL;
      int NbrLabels;
      if ((LatticeDefinition.GetAsStringArray ("ExtParametersLabels", ',', ParameterLabels, NbrLabels)==false) || (NbrLabels!=NbrExtParameters))
	{
	  // enhance lattice descriptor in default manner:
	  if (NbrLabels!=NbrExtParameters)
	    cout << "Incorrect number of parameter labels - ignoring ExtParametersLabels"<<endl;
	  char *NewDescriptor = new char[strlen(this->Descriptor)+4+15*NbrExtParameters];
	  sprintf(NewDescriptor,"%s_ep_%g",Descriptor,ExtParameters[0]);
	  for (int i=1; i<NbrExtParameters; ++i)
	    {
	      sprintf(NewDescriptor,"%s_%g",NewDescriptor,ExtParameters[i]);
	    }
	  delete [] this->Descriptor;
	  this->Descriptor=NewDescriptor;
	}
      else
	{
	  // customize parameter labels
	  int AddLength = 0;
	  for (int i=0; i<NbrExtParameters; ++i) AddLength += strlen(ParameterLabels[i]);
	  
	  char *NewDescriptor = new char[strlen(this->Descriptor)+AddLength+15*NbrExtParameters];
	  int PrefixOffset=0;
	  PrefixOffset+=sprintf(NewDescriptor, "%s_%s_%g", Descriptor, ParameterLabels[0], ExtParameters[0]);
	  delete [] ParameterLabels[0];
	  for (int i=1; i<NbrExtParameters; ++i)
	    {
	      PrefixOffset+=sprintf(NewDescriptor+PrefixOffset, "_%s_%g", ParameterLabels[i], ExtParameters[i]);
	      delete [] ParameterLabels[i];
	    }
	  delete [] ParameterLabels;
	  delete [] this->Descriptor;
	  this->Descriptor=NewDescriptor;
	}
    }
  
  // parsing jump parameters:
  bool TRSymmetrize=false;
  if ((LatticeDefinition["TimeReversalSymmetrize"]!=NULL)&&
      (    (strstr(LatticeDefinition["TimeReversalSymmetrize"],"yes")!=NULL)
	|| (strstr(LatticeDefinition["TimeReversalSymmetrize"],"YES")!=NULL)
	|| (strstr(LatticeDefinition["TimeReversalSymmetrize"],"true")!=NULL)
	|| (strstr(LatticeDefinition["TimeReversalSymmetrize"],"TRUE")!=NULL) ))
    {
      cout << "Generating time-reversal symmetrized hoppings"<<endl;
      TRSymmetrize=true;
    }
  
  if (LatticeDefinition.GetAsSingleInteger ("NbrJumpTerms", this->NbrJumpTerms)==false)
    {
      cout << "Error: NbrJumpTerms needs to be defined"<<endl;
      std::exit(1);
    }

  int MaxJumpTerms = NbrJumpTerms;
  if (TRSymmetrize) MaxJumpTerms*=2; // account for time-reversal inverse processes
  MaxJumpTerms *= SymmetryMultiplier1 * SymmetryMultiplier2; // jumps are replicated for each unit cell

  this->DeltaG1 = new int [MaxJumpTerms];
  this->DeltaG2 = new int [MaxJumpTerms];
  this->InitialSublattice = new int [MaxJumpTerms];
  this->FinalSublattice = new int [MaxJumpTerms];
  this->JumpAmplitudes  = new Complex [MaxJumpTerms];
  
  char **JumpString;
  int NbrFields;
  int dG1, dG2, sub1=0, sub2=0;
  double PhaseOffset, PhaseIncrementX, PhaseIncrementY, Amplitude;
  int CurrentJumpTerm=0;
  for (int i=0; i<NbrJumpTerms; ++i)
    {
      sprintf(FieldName,"JumpTerm%d",i);
      
      if (LatticeDefinition.GetAsStringArray (FieldName, ',', JumpString, NbrFields)==false)
	{
	  cout << "error while parsing "<<FieldName<<" in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);
	}

      if (NbrFields<4)
	{
	  cout << "error while decoding "<<FieldName<<" in " << this->Options->GetString("lattice-definition") <<endl ;
	  cout << "Indicate jump terms in the following format:\n"
	       << "JumpTermX = deltaG1,deltaG2,SublatticeI,SublatticeF[,phaseOffset=0.0,phaseIncrementX=0.0,phaseIncrementY=0.0,amplitudeParamaterID=0]\n"
	       << "Phases in units of pi (magnitudes) can either be indicated explitly or are assumed to be zero (one), otherwise."<<endl;
	  exit(-1);
	}
      dG1=strtol(JumpString[0], NULL,0);
      dG2=strtol(JumpString[1], NULL,0);
      if ((abs(dG1)>CutOffMomentum) || (abs(dG2)>CutOffMomentum))
	{
	  cout << "Attention: jump term with dG=("<<dG1<<", "<<dG2<<") is of same magnitude as cut-off momentum!"<<endl;
	}

      sub1 = strtol(JumpString[2], NULL,0);
      sub2 = strtol(JumpString[3], NULL,0);

      PhaseOffset=0.0;
      PhaseIncrementX=0.0;
      PhaseIncrementY=0.0;
      Amplitude=ExtParameters[0];
      // read phase
      if (NbrFields>4)
	{
          PhaseOffset = strtod(JumpString[4], NULL);
	}
      if (NbrFields>5)
	{
          PhaseIncrementX = strtod(JumpString[5], NULL);
	}
      if (NbrFields>6)
	{
          PhaseIncrementY = strtod(JumpString[6], NULL);
	}
      if (NbrFields>7)
	{
	  double AmplitudeNbr = strtod(JumpString[7], NULL);
	  int index = (int)nearbyint(AmplitudeNbr);
	  Amplitude = ExtParameters[index];
	}

      if (fabs(Amplitude)>1e-15)
	{
          for (int n1=0; n1<SymmetryMultiplier1; ++n1)
	    for (int n2=0; n2<SymmetryMultiplier2; ++n2)
	      {
                //G1, G2 here in terms of the symmetry extended BZ! deduce the effective values of the offsets and sublattices in this new frame
                int dG1Eff, dG2Eff, sub1Eff, sub2Eff;
		this->ExpressInExtendedCell(dG1, dG2, sub1, sub2, n1, n2, dG1Eff, dG2Eff, sub1Eff, sub2Eff);
                DeltaG1[CurrentJumpTerm]=dG1Eff;
	        DeltaG2[CurrentJumpTerm]=dG2Eff;
		// Sublattice flavours here in terms of extended sublattices!
	        InitialSublattice[CurrentJumpTerm]=sub1Eff;
	        FinalSublattice[CurrentJumpTerm]=sub2Eff;
		JumpAmplitudes[CurrentJumpTerm]=Polar(Amplitude,(PhaseOffset+n1*PhaseIncrementX+n2*PhaseIncrementY)*M_PI);
		//		cout <<"PhaseOffset="<<PhaseOffset<<", PhaseIncrementX="<<PhaseIncrementX<<", PhaseIncrementY="<<PhaseIncrementY<<", JumpAmplitudes["<<CurrentJumpTerm<<"]="<<JumpAmplitudes[CurrentJumpTerm]<<endl;
		++CurrentJumpTerm;
		if (TRSymmetrize)
		  {
		    this->ExpressInExtendedCell(-dG1, -dG2, sub2, sub1, n1, n2, dG1Eff, dG2Eff, sub1Eff, sub2Eff);   // changed argument from  (...sub1, sub2...) to  (...sub2, sub1...)
                    DeltaG1[CurrentJumpTerm]=dG1Eff;
		    DeltaG2[CurrentJumpTerm]=dG2Eff;
		    InitialSublattice[CurrentJumpTerm]=sub1Eff;
		    FinalSublattice[CurrentJumpTerm]=sub2Eff;
		    JumpAmplitudes[CurrentJumpTerm]=Polar(Amplitude,-(PhaseOffset+n1*PhaseIncrementX+n2*PhaseIncrementY)*M_PI);
		    //		    cout <<"PhaseOffset="<<PhaseOffset<<", PhaseIncrementX="<<PhaseIncrementX<<", PhaseIncrementY="<<PhaseIncrementY<<", JumpAmplitudes["<<CurrentJumpTerm<<"]="<<JumpAmplitudes[CurrentJumpTerm]<<endl;
		    ++CurrentJumpTerm;
                  }
              }
         }
    }
  if (CurrentJumpTerm < MaxJumpTerms)
    {
      int *TmpI1 = new int[CurrentJumpTerm];
      int *TmpI2 = new int[CurrentJumpTerm];
      int *TmpI3 = new int[CurrentJumpTerm];
      int *TmpI4 = new int[CurrentJumpTerm];
      for (int i=0; i<CurrentJumpTerm; ++i)
	{
	  TmpI1[i]=DeltaG1[i];
	  TmpI2[i]=DeltaG2[i];
	  TmpI3[i]=InitialSublattice[i];
	  TmpI4[i]=FinalSublattice[i];
	}
      delete [] DeltaG1; DeltaG1 = TmpI1;
      delete [] DeltaG2; DeltaG2 = TmpI2;
      delete [] InitialSublattice; InitialSublattice = TmpI3;
      delete [] FinalSublattice; FinalSublattice = TmpI4;

      Complex *TmpC = new Complex[CurrentJumpTerm];
      for (int i=0; i<CurrentJumpTerm; ++i)
	TmpC[i]=JumpAmplitudes[i];
      delete [] JumpAmplitudes; JumpAmplitudes=TmpC;
    }
  this->NbrJumpTerms = CurrentJumpTerm;
  // cleaning up things equal to numerical zero
  for (int i=0; i<this->NbrJumpTerms; ++i)
    {
      if (fabs(JumpAmplitudes[i].Re) < MACHINE_PRECISION)
	JumpAmplitudes[i].Re=0.0;
      if (fabs(JumpAmplitudes[i].Im) < MACHINE_PRECISION)
	JumpAmplitudes[i].Im=0.0;
    }
  
#ifdef DEBUG_OUTPUT
  cout << "Jump terms:\n" <<"Idx\tdG1\tdG2\ts_i\ts_f\tt_ij\n";
  for (int i=0; i<this->NbrJumpTerms; ++i)
    cout << i << "\t"<<DeltaG1[i]<< "\t"<< DeltaG2[i]<< "\t"<<InitialSublattice[i]<< "\t"<<FinalSublattice[i]<<"\t"<<JumpAmplitudes[i]<<endl;
#endif
	 
  // end of parsing parameters; infer remaining variables:

  // choose NMax values such that all k-points within circle of radius CutOffMomentum
  // (might revise this to carefully look at all sublattices)
  if( this->Mode == Circular)
    {
      this->NMax1 = std::ceil(CutOffMomentum / (SymmetryMultiplier1*this->LatticeVector1.Norm()));
      this->NMax2 = std::ceil(CutOffMomentum / (SymmetryMultiplier2*this->LatticeVector2.Norm()));
    }
  else
    {
      this->NMax1=nMax1;
      this->NMax2=nMax2;
    }

  if (ExtNbrSubLattices % NbrSubLatticeFlavours != 0)
    {
      cout << "Error: The extended number of sublattices must be a multiple of the number of flavours";
      exit(1);
    }
  
  this->NbrSubLatticesPerFlavour = ExtNbrSubLattices / NbrSubLatticeFlavours;

  this->NbrStatePerBand *= SymmetryMultiplier1*SymmetryMultiplier2;

  cout << "this->NbrStatePerBand = "<<this->NbrStatePerBand<<endl;


  this->FullNbrBands = (2*this->NMax1+1)*(2*this->NMax2+1)*ExtNbrSubLattices;
  if (nbrBandsToKeep<0)
    this->NbrBands = this->FullNbrBands;
  else
    {
      this->NbrBands = (FullNbrBands > nbrBandsToKeep ? nbrBandsToKeep : this->FullNbrBands);
      if (NbrBands<2 && FullNbrBands>1) ++NbrBands;
    }

  cout << "Total number of bands: "<<this->FullNbrBands<<", NMax1="<<this->NMax1<<", NMax2="<<this->NMax2<<", extended sublattices="<<ExtNbrSubLattices<<" ("<<NbrBands<<" kept for output)"<<endl;

  if (storeOneBodyMatrices == true)
    {
      this->OneBodyBasis = new ComplexMatrix [this->NbrStatePerBand];
    }
  else
    {
      this->OneBodyBasis = 0;
    }
  this->EnergyBandStructure = new double*[this->NbrBands];
  for (int i = 0; i < this->NbrBands; ++i)
    {
      this->EnergyBandStructure[i] = new double[this->NbrStatePerBand];
    }

  // increase number of k-points by symmetry multipliers for compatibility with Abstract class
  this->NbrSiteX *= SymmetryMultiplier1;
  this->NbrSiteY *= SymmetryMultiplier2;

  this->ComputeBandStructure();  

  // check if symmetries are realized sufficiently accurately:

  for (int k1 = 0; k1 < this->NbrPoints1; ++k1)
    for (int k2 = 0; k2 < this->NbrPoints2; ++k2)
      {
	int Index1 = this->GetLinearizedMomentumIndex(k1, k2);
	for (int s1=0; s1<SymmetryMultiplier1; ++s1)
	  for (int s2=0; s2<SymmetryMultiplier2; ++s2)
	    {
	      int Index2 = this->GetLinearizedMomentumIndex(k1+s1*this->NbrPoints1, k2+s2*this->NbrPoints2);
	      for (int i = 0; i < this->NbrBands; ++i)
		{		   
		  if ( fabs(this->EnergyBandStructure[i][Index1] - this->EnergyBandStructure[i][Index2]) > symmetryThreshold)
		    {
		      cout << "Error: inaccurate energy bandstructure - try to increase cut-off radius"<<endl;
		      if ( fabs(this->EnergyBandStructure[i][Index1] - this->EnergyBandStructure[i][Index2]) > 100.0*symmetryThreshold)
			{
			  cout << "Large inaccuracy: DeltaE[n="<<i<<", Index1="<<Index1<<", Index2="<<Index2<<"]="<<(fabs(this->EnergyBandStructure[i][Index1] - this->EnergyBandStructure[i][Index2])) <<"! Aborting calculation"<<endl;
			  exit(1);
			}
		    }
		}
	    }
      }	    
}

// destructor
//

TightBindingModelOFLGenericLatticeWithSymmetry::~TightBindingModelOFLGenericLatticeWithSymmetry()
{
}

OptionManager* TightBindingModelOFLGenericLatticeWithSymmetry::Options=NULL;



// add an option group containing all options related to the LatticeGeometry options
//
// manager = pointer to the option manager
void TightBindingModelOFLGenericLatticeWithSymmetry::AddOptionGroup(OptionManager* manager)
{
  TightBindingModelOFLGenericLatticeWithSymmetry::Options = manager;
  OptionGroup* LatticeGroup  = new OptionGroup ("flux lattice options");
  (*(TightBindingModelOFLGenericLatticeWithSymmetry::Options)) += LatticeGroup;
  (*LatticeGroup) += new SingleStringOption ('L', "lattice-definition", "File defining the geometry of the lattice");
  (*LatticeGroup) += new MultipleDoubleOption ('\n', "external-parameters", "values of external parameters for lattice definition", ',');
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelOFLGenericLatticeWithSymmetry::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  double K2;

  for (int k1 = 0; k1 < this->NbrPoints1*SymmetryMultiplier1; ++k1)
    {
      for (int k2 = 0; k2 < this->NbrPoints2*SymmetryMultiplier2; ++k2)
	{
	  int Index = this->GetLinearizedMomentumIndex(k1, k2);
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      // assign matrix to collect eigenvectors
	      ComplexMatrix TmpMatrix2(this->FullNbrBands, this->NbrBands, true);

	      K1 = this->DeltaK1*(((double) k1) + this->GammaX);
	      K2 = this->DeltaK2*(((double) k2) + this->GammaY);

	      // construct magnetic unit cell:
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->FullNbrBands, true);
	      for (int q1=0; q1< 2*NMax1+1 ; q1++)
		{
		  for (int q2=0; q2< 2*NMax2+1; q2++)
		    {
		      bool InitialIndexInBounds, FinalIndexInBounds;
		      double InitialNorm=0.0, FinalNorm=0.0;
		      int InitialIndex, FinalIndex;
		      // kinetic energy:
		      for (int i=0; i<ExtNbrSubLattices; ++i)
			{
			  InitialIndex=this->LinearizedReciprocalSpaceIndex(q1, q2, K1, K2, i, InitialNorm, InitialIndexInBounds);
			  TmpOneBodyHamiltonian.SetMatrixElement(InitialIndex, InitialIndex, InitialNorm*InitialNorm/(M_PI*M_PI));
			}
		      // hopping terms:
		      for (int j=0; j<this->NbrJumpTerms; ++j)				
			{
			  InitialIndex=this->LinearizedReciprocalSpaceIndex(q1, q2, K1, K2, InitialSublattice[j], InitialNorm, InitialIndexInBounds);
			  if (InitialIndexInBounds)
			    {
			      FinalIndex=this->LinearizedReciprocalSpaceIndex(q1+DeltaG1[j], q2+DeltaG2[j], K1, K2, FinalSublattice[j], FinalNorm, FinalIndexInBounds);
			      if (FinalIndexInBounds)
				{
				  TmpOneBodyHamiltonian.SetMatrixElement(InitialIndex, FinalIndex, this->LatticeDepth*this->JumpAmplitudes[j]);
#ifdef DEBUG_OUTPUT
				  cout << "Added jump term "<<j<<": (" << q1 <<", "<<q2<<"; "<<InitialIndex <<")->("<<q1+DeltaG1[j]<<", "<<q2+DeltaG2[j]<<";"<<FinalIndex <<") = "<<this->LatticeDepth*this->JumpAmplitudes[j]<<endl;
#endif
				}
			    }
			}
		    }
		}
#ifdef DEBUG_OUTPUT
	      cout << TmpOneBodyHamiltonian << endl;
#endif
	      if (this->OneBodyBasis != 0)
		{
		  ComplexMatrix TmpMatrix(this->FullNbrBands, this->FullNbrBands, true);
		  TmpMatrix.SetToIdentity();
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif

		  ComplexMatrix TmpMatrix2(this->FullNbrBands, this->NbrBands, true);
		  for (int i = 0; i < this->NbrBands; ++i)
		    TmpMatrix2[i] = TmpMatrix[i];
		  
		  this->OneBodyBasis[Index] = TmpMatrix2;		  		 

		  cout << "Sector "<< (Index+1) <<"/"<<this->NbrStatePerBand<<endl;		 
		  for (int i = 0; i < this->NbrBands; ++i)
		    {
		      this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		      cout << TmpDiag(i, i) << " ";			      
		    }
		  cout << endl;
		  
		  // check weight of wavefunction at boundaries of simulation cell
		  double Norm, Weight=0.0;
		  bool InBounds;
		  for (int i=0; i<ExtNbrSubLattices; ++i)
		    {
		      for (int q1=0; q1< 2*NMax1+1 ; q1++)
			{
			  int RSIndex=this->LinearizedReciprocalSpaceIndex(q1, 0, 0.0, 0.0, i, Norm, InBounds);
			  Weight += SqrNorm(TmpMatrix2[0][RSIndex]);
			  RSIndex=this->LinearizedReciprocalSpaceIndex(q1, 2*NMax2, 0.0, 0.0, i, Norm, InBounds);
			  Weight += SqrNorm(TmpMatrix2[0][RSIndex]);
			}
		      for (int q2=1; q2< 2*NMax2 ; q2++)
			{
			  int RSIndex=this->LinearizedReciprocalSpaceIndex(0, q2, 0.0, 0.0, i, Norm, InBounds);
			  Weight += SqrNorm(TmpMatrix2[0][RSIndex]);
			  RSIndex=this->LinearizedReciprocalSpaceIndex(2*NMax1, q2, 0.0, 0.0, i, Norm, InBounds);
			  Weight += SqrNorm(TmpMatrix2[0][RSIndex]);
			}
		    }
		  cout << "Weight in extremal unit-cells = "<<Weight<<" (Norm: "<< TmpMatrix2[0].Norm()<<")"<<endl;
		}
	      else
		{
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
		  cout << "Sector "<< (Index+1) <<"/"<<this->NbrStatePerBand<<endl;
		  for (int i = 0; i < this->NbrBands; ++i)
		    {
		      this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		      cout << TmpDiag(i, i) << " ";
		    }
		  cout << endl;
		}
	    }
	}
    }
}


// remap a vector for a translation of a given momentum
//
// vectorIn = input vector
// vectorOut = input vector
// dG1 = amount of translation along G1
// dG2 = amount of translation along G2
void TightBindingModelOFLGenericLatticeWithSymmetry::TranslateVector(ComplexVector &vectorIn, ComplexVector &vectorOut, int dG1, int dG2)
{
  int dG1Eff, dG2Eff, sub1Eff, sub2Eff;
  for (int i=0; i<NbrSubLattices; ++i)
    for (int s1= 0; s1 < SymmetryMultiplier1; ++s1)
      for (int s2= 0; s2 < SymmetryMultiplier2; ++s2)
	{
	  this->ExpressInExtendedCell(dG1, dG2, /*sub1*/ i , /*sub2*/ i, s1, s2, dG1Eff, dG2Eff, sub1Eff, sub2Eff);
	  
	  for (int q1=0; q1< 2*NMax1+1 ; q1++)	    
	    for (int q2=0; q2< 2*NMax2+1; q2++)
	      {
		int InitialIndex = PeriodicLinearizedReciprocalSpaceIndex(q1, q2, sub1Eff);
		int FinalIndex = PeriodicLinearizedReciprocalSpaceIndex(q1+dG1Eff, q2+dG2Eff, sub2Eff);
		vectorOut[FinalIndex] = vectorIn[InitialIndex];
	      }
	}
#ifdef DEBUG_OUTPUT
  cout << "NormIn = "<<vectorIn.Norm()<<", NormOut="<<vectorOut.Norm()<<endl;
#endif
}



// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool TightBindingModelOFLGenericLatticeWithSymmetry::WriteAsciiSpectrum(char* fileName)
{
  ofstream File;
  File.open(fileName);
  this->WriteASCIIHeader(File, '#');
  int LimitOut = 100;
  if (this->NbrBands < LimitOut)
    LimitOut = this->NbrBands;
  File << "# k1    k2";
  //for (int i = 0; i < LimitOut ; ++i)
  //  File <<  "    E_" << i;
  File << endl;
  //for (int k1 = 0; k1 < this->NbrPoints1*SymmetryMultiplier1; ++k1)
  for (int k1 = 0; k1 < this->NbrPoints1; ++k1)
    {
      //for (int k2 = 0; k2 < this->NbrPoints2*SymmetryMultiplier2; ++k2)
      for (int k2 = 0; k2 < this->NbrPoints2; ++k2)
	{
	  int LinearizedMomentumIndex = this->GetLinearizedMomentumIndex(k1, k2);

	  for (int i = 0; i < LimitOut; ++i)
	    File << k1 << " " << k2 << " " << this->EnergyBandStructure[i][LinearizedMomentumIndex]<<endl;
	  File << endl; 
	}
      File << endl;
    }
  File.close();
  return true;
}
