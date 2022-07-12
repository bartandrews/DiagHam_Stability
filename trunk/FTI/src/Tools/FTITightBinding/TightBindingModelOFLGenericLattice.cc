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
#include "Tools/FTITightBinding/TightBindingModelOFLGenericLattice.h"
#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"
#include "GeneralTools/ConfigurationParser.h"
#include "GeneralTools/FilenameTools.h"

#include <iostream>
#include <cmath>
using std::cout;
using std::endl;

//#define DEBUG_OUTPUT


// default constructor
//
// nbrPointsX = number of k-points in the 1-direction
// nbrPointsY = number of k-points in the 2-direction
// gamma1 = boundary condition twisting angle along x
// gamma2 = boundary condition twisting angle along y
// architecture = pointer to the architecture
// cutOffMomentum = maximum (absolute value) of momenta considered
// latticeDepth = parameter for the depth of the optical lattice in recoil energies
// storeOneBodyMatrices = flag to indicate if the one body transformation matrices have to be computed and stored
TightBindingModelOFLGenericLattice::TightBindingModelOFLGenericLattice(int nbrPoints1, int nbrPoints2, double gamma1, double gamma2, AbstractArchitecture* architecture, double cutOffMomentum, double latticeDepth, bool storeOneBodyMatrices)
{
  if (TightBindingModelOFLGenericLattice::Options==NULL)
    {
      cout << "Define the OptionManager, first, before creating any TightBindingModelOFLGenericLattice"<<endl;
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
  this->NbrSiteX= this->NbrPoints1 = nbrPoints1;
  this->NbrSiteY= this->NbrPoints2 = nbrPoints2;
  this->CutOffMomentum = cutOffMomentum;
  this->LatticeDepth = latticeDepth;

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
      if (Dimension!=NbrComponents)
	{
	  cout << "SubLattice Vectors need to have the dimension of the lattice!"<<endl;
	  exit(-1);
	}
      for (int j=0; j<Dimension; ++j)
	{
	  SubLatticeVectors[i][j]=Components[j];
	}
      delete [] Components;
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
      if ( (*TightBindingModelOFLGenericLattice::Options)["external-parameters"] != NULL)
	TmpParameters = TightBindingModelOFLGenericLattice::Options->GetDoubles("external-parameters",TmpDimension);
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
  if (NbrSubLattices==2) MaxJumpTerms*=2; // account for possibility of sigma matrices

  this->DeltaG1 = new int [MaxJumpTerms];
  this->DeltaG2 = new int [MaxJumpTerms];
  this->InitialSublattice = new int [MaxJumpTerms];
  this->FinalSublattice = new int [MaxJumpTerms];
  this->JumpAmplitudes  = new Complex [MaxJumpTerms];
  
  char **JumpString;
  int NbrFields;
  int dG1, dG2, sub1=0, sub2=0;
  double Phase, Amplitude;
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
	       << "JumpTermX = deltaG1,deltaG2,SublatticeI,SublatticeF[,phaseIF=0.0,amplitudeParamaterID=0]\n"
	       << "Phases and magnitudes can either be indicated explitly or are assumed to be one, otherwise."<<endl;
	  exit(-1);
	}
      dG1=strtol(JumpString[0], NULL,0);
      dG2=strtol(JumpString[1], NULL,0);
      if ((abs(dG1)>CutOffMomentum) || (abs(dG2)>CutOffMomentum))
	{
	  cout << "Attention: jump term with dG=("<<dG1<<", "<<dG2<<") is of same magnitude as cut-off momentum!"<<endl;
	}

      int HaveSigmaMatrix=-1;
      if (NbrSubLattices==2)
	{
	  if ((JumpString[3][0]=='i')||(JumpString[3][0]=='I'))
	    HaveSigmaMatrix=0;
	  if ((JumpString[3][0]=='x')||(JumpString[3][0]=='X'))
	    HaveSigmaMatrix=1;
	  if ((JumpString[3][0]=='y')||(JumpString[3][0]=='Y'))
	    HaveSigmaMatrix=2;
	  if ((JumpString[3][0]=='z')||(JumpString[3][0]=='Z'))
	    HaveSigmaMatrix=3;
	}
      if (HaveSigmaMatrix<0)
	{
	  sub1 = strtol(JumpString[2], NULL,0);
	  sub2 = strtol(JumpString[3], NULL,0);
	}

      Phase=0.0;
      Amplitude=ExtParameters[0];
      // read phase
      if ((NbrFields>4) && (HaveSigmaMatrix<0))
	{
	  Phase = strtod(JumpString[4], NULL);
	}
      if (NbrFields>5)
	{
	  double AmplitudeNbr = strtod(JumpString[5], NULL);
	  int index = (int)nearbyint(AmplitudeNbr);
	  Amplitude = ExtParameters[index];
	}

      if (fabs(Amplitude)>1e-15)
	{
	  if (HaveSigmaMatrix < 0)
	    {
	      DeltaG1[CurrentJumpTerm]=dG1;
	      DeltaG2[CurrentJumpTerm]=dG2;
	      InitialSublattice[CurrentJumpTerm]=sub1;
	      FinalSublattice[CurrentJumpTerm]=sub2;
	      JumpAmplitudes[CurrentJumpTerm]=Polar(Amplitude,Phase*M_PI);
	      ++CurrentJumpTerm;
	      if (TRSymmetrize)
		{
		  DeltaG1[CurrentJumpTerm]=-dG1;
		  DeltaG2[CurrentJumpTerm]=-dG2;
		  InitialSublattice[CurrentJumpTerm]=sub2;
		  FinalSublattice[CurrentJumpTerm]=sub1;
		  JumpAmplitudes[CurrentJumpTerm]=Polar(Amplitude,-Phase*M_PI);
		  ++CurrentJumpTerm;
		}
	    }
	  else
	    {
	      switch(HaveSigmaMatrix)
		{
		case 0: // identity
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=0;
		  FinalSublattice[CurrentJumpTerm]=0;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(Amplitude,0.0);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=1;
		  FinalSublattice[CurrentJumpTerm]=1;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(Amplitude,0.0);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  break;
		case 1: // sigma_x
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=0;
		  FinalSublattice[CurrentJumpTerm]=1;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(Amplitude,0.0);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=1;
		  FinalSublattice[CurrentJumpTerm]=0;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(Amplitude,0.0);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  break;
		case 2: // sigma_y
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=0;
		  FinalSublattice[CurrentJumpTerm]=1;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(0.0,-Amplitude);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=1;
		  FinalSublattice[CurrentJumpTerm]=0;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(0.0,Amplitude);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  break;
		case 3: // sigma_z
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=0;
		  FinalSublattice[CurrentJumpTerm]=0;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(Amplitude, 0.0);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  DeltaG1[CurrentJumpTerm]=dG1;
		  DeltaG2[CurrentJumpTerm]=dG2;
		  InitialSublattice[CurrentJumpTerm]=1;
		  FinalSublattice[CurrentJumpTerm]=1;
		  JumpAmplitudes[CurrentJumpTerm]=Complex(-Amplitude,0.0);
		  ++CurrentJumpTerm;
		  if (TRSymmetrize)
		    {
		      DeltaG1[CurrentJumpTerm]=-dG1;
		      DeltaG2[CurrentJumpTerm]=-dG2;
		      InitialSublattice[CurrentJumpTerm]=FinalSublattice[CurrentJumpTerm-1];
		      FinalSublattice[CurrentJumpTerm]=InitialSublattice[CurrentJumpTerm-1];
		      JumpAmplitudes[CurrentJumpTerm]=Conj(JumpAmplitudes[CurrentJumpTerm-1]);
		      ++CurrentJumpTerm;
		    }
		  break;
		default:
		  cout << "Error: sigma matrices have indices 0,...,3 only"<<endl;
		  exit(1);
		  break;
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
  this->NMax1 = CutOffMomentum / this->LatticeVector1.Norm();
  this->NMax2 = CutOffMomentum / this->LatticeVector2.Norm();

  this->TmpVector.ResizeAndClean(2);
  TmpVector.AddLinearCombination((double)NMax1,this->LatticeVector1);
  TmpVector.AddLinearCombination((double)NMax2,this->LatticeVector2);
  if (TmpVector.Norm()<CutOffMomentum)
    {
      NMax1 = std::ceil(TmpVector.Norm()/CutOffMomentum * NMax1);
      NMax2 = std::ceil(TmpVector.Norm()/CutOffMomentum * NMax2);
    }
  TmpVector.ClearVector();
  TmpVector.AddLinearCombination((double)NMax1,this->LatticeVector1);
  TmpVector.AddLinearCombination(-(double)NMax2,this->LatticeVector2);
  if (TmpVector.Norm()<CutOffMomentum)
    {
      NMax1 = std::ceil(TmpVector.Norm()/CutOffMomentum * NMax1);
      NMax2 = std::ceil(TmpVector.Norm()/CutOffMomentum * NMax2);
    }
  ++NMax1;
  ++NMax2;
  
  this->NbrBands = (2*this->NMax1+1)*(2*this->NMax2+1);

  cout << "Total number of bands: "<<this->NbrBands<<endl;

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

  this->ComputeBandStructure();  
}

// destructor
//

TightBindingModelOFLGenericLattice::~TightBindingModelOFLGenericLattice()
{
}

OptionManager* TightBindingModelOFLGenericLattice::Options=NULL;



// add an option group containing all options related to the LatticeGeometry options
//
// manager = pointer to the option manager
void TightBindingModelOFLGenericLattice::AddOptionGroup(OptionManager* manager)
{
  TightBindingModelOFLGenericLattice::Options = manager;
  OptionGroup* LatticeGroup  = new OptionGroup ("flux lattice options");
  (*(TightBindingModelOFLGenericLattice::Options)) += LatticeGroup;
  (*LatticeGroup) += new SingleStringOption ('L', "lattice-definition", "File defining the geometry of the lattice");
  (*LatticeGroup) += new MultipleDoubleOption ('\n', "external-parameters", "values of external parameters for lattice definition", ',');
}

// core part that compute the band structure
//
// minStateIndex = minimum index of the state to compute
// nbrStates = number of states to compute

void TightBindingModelOFLGenericLattice::CoreComputeBandStructure(long minStateIndex, long nbrStates)
{
  if (nbrStates == 0l)
    nbrStates = this->NbrStatePerBand;
  long MaxStateIndex = minStateIndex + nbrStates;
  double K1;
  double K2;
  for (int k1 = 0; k1 < this->NbrPoints1 ; ++k1)
    {
      for (int k2 = 0; k2 < this->NbrPoints2; ++k2)
	{
	  int Index = this->GetLinearizedMomentumIndex(k1, k2);
	  cout << MaxStateIndex<<endl;
	  cout << Index <<endl;
	  if ((Index >= minStateIndex) && (Index < MaxStateIndex))
	    {
	      K1 = this->DeltaK1*(((double) k1) + this->GammaX);
	      K2 = this->DeltaK2*(((double) k2) + this->GammaY);
	      
	      // construct magnetic unit cell:
	      
	      HermitianMatrix TmpOneBodyHamiltonian(this->NbrBands, true);
	      for (int q1=0; q1< 2*NMax1+1 ; q1++)
		{
		  for (int q2=0; q2< 2*NMax2+1; q2++)
		    {
		      bool InitialIndexInBounds, FinalIndexInBounds;
		      double InitialNorm=0.0, FinalNorm=0.0;
		      int InitialIndex, FinalIndex;
		      // kinetic energy:
		      for (int i=0; i<NbrSubLattices; ++i)
			{
			  InitialIndex=this->LinearizedReciprocalSpaceIndex(q1, q2, K1, K2, i, InitialNorm, InitialIndexInBounds);
			  TmpOneBodyHamiltonian.SetMatrixElement(InitialIndex, InitialIndex, InitialNorm*InitialNorm);
			}
		      // hopping terms:
		      for (int j=0; j<this->NbrJumpTerms; ++j)
			{
			  InitialIndex=this->LinearizedReciprocalSpaceIndex(q1, q2, K1, K2, InitialSublattice[j], InitialNorm, InitialIndexInBounds);
			  if (InitialIndexInBounds)
			    {
			      FinalIndex=this->LinearizedReciprocalSpaceIndex(q1+DeltaG1[j], q2+DeltaG2[j], K1, K2, FinalSublattice[j], FinalNorm, FinalIndexInBounds);
			      if (FinalIndexInBounds)
				TmpOneBodyHamiltonian.SetMatrixElement(InitialIndex, FinalIndex, this->LatticeDepth*this->JumpAmplitudes[j]);
			    }
			}
		    }
		}
#ifdef DEBUG_OUTPUT
	      cout << TmpOneBodyHamiltonian << endl;
#endif
	      if (this->OneBodyBasis != 0)
		{
		  ComplexMatrix TmpMatrix(this->NbrBands, this->NbrBands, true);
		  TmpMatrix.SetToIdentity();
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag, TmpMatrix);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag, TmpMatrix);
#endif
		  this->OneBodyBasis[Index] = TmpMatrix;
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	      else
		{
		  RealDiagonalMatrix TmpDiag;
#ifdef __LAPACK__
		  TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag);
#else
		  TmpOneBodyHamiltonian.Diagonalize(TmpDiag);
#endif
		  for (int i = 0; i < this->NbrBands; ++i)
		    this->EnergyBandStructure[i][Index] = TmpDiag(i, i);
		}
	    }
	}
    }
}


// write the energy spectrum in an ASCII file
//
// fileName = name of the ASCII file 
// return value = true if no error occured

bool TightBindingModelOFLGenericLattice::WriteAsciiSpectrum(char* fileName)
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
  for (int k1 = 0; k1 < this->NbrPoints1; ++k1)
    {
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
