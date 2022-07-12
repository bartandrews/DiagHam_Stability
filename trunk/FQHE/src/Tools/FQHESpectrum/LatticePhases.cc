////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2009 Nicolas Regnault                  //
//                         class author: Gunnar Möller                        //
//                                                                            //
//                                                                            //
//  class that defines sites and nearest neighbor relationships on a lattice  //
//                                                                            //
//                        last modification : 05/26/2009                      //
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

#include "LatticePhases.h"
#include "Options/Options.h"
#include "Matrix/RealSymmetricMatrix.h"
#include "Matrix/RealAntisymmetricMatrix.h"
#include "MathTools/RandomNumber/NumRecRandomGenerator.h"
#include "GeneralTools/ConfigurationParser.h"


#include <iostream>
using std::cout;
using std::endl;

#include <cstdlib>
#include <cstring>
#include <cmath>
using std::fabs;

// switch for debugging statements
#define DEBUG_LATTICE_PHASES

// generate the object using options from Option Manager
//
LatticePhases::LatticePhases()
{
  if (LatticePhases::Options==NULL)
    {
      cout << "Define the OptionManager, first, before creating any LatticePhases"<<endl;
      exit(1);
    }
  
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
  if ((LatticeDefinition.GetAsSingleInteger("NbrSites", NbrSitesPerCell) == false) || (NbrSitesPerCell <= 0))
    {
      cout << "NbrSites is not defined or has invalid value" << endl;
      exit(-1);
    }
  if ((LatticeDefinition.GetAsSingleInteger("Dimension", Dimension) == false) || (Dimension <= 0))
    {
      cout << "Dimension is not defined or has invalid value" << endl;
      exit(-1);
    }
  int TmpDimension;
  PeriodicRep = LatticePhases::Options->GetIntegers("cells-repeat",TmpDimension);
  if (TmpDimension==0)
    if ((LatticeDefinition.GetAsIntegerArray("PeriodicRepeat",',',PeriodicRep, TmpDimension) == false))
      {
	cout << "PeriodicRepeat is not defined or has invalid value, simulating a single unit cell" << endl;
	PeriodicRep = new int[Dimension];
	TmpDimension=Dimension;
	for (int i=0; i<Dimension; ++i) PeriodicRep[i]=1;
      }
  if (TmpDimension!=Dimension)
    {
      cout << "PeriodicRepeat does not have the right number of components (separator: ',')" << endl;
      exit(-1);
    }
  NbrSites=NbrSitesPerCell;
  for (int i=0; i<Dimension; ++i)
    NbrSites*=PeriodicRep[i];
  LatticeVectors.Resize(Dimension, Dimension);
  char *FieldName = new char[255];
  for (int i=0; i<Dimension; ++i)
    {
      sprintf(FieldName,"LatticeVector%d",i);
      int NbrComponents;
      double *Components;
      if (LatticeDefinition.GetAsDoubleArray(FieldName, ',', Components, NbrComponents) == false)
	{
	  cout << "error while parsing "<<FieldName<< " in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);     
	}
      if (Dimension!=NbrComponents)
	{
	  cout << "Lattice Vectors need to have the dimension of the lattice!"<<endl;
	  exit(-1);
	}
      for (int j=0; j<Dimension; ++j)
	{
	  LatticeVectors[i][j]=Components[j];
	}
      delete [] Components;
    }
  double Rescale=1.0;
  if (Options->GetBoolean("normalize-lattice"))
    {
      double Area = LatticeVectors.Determinant();
      Rescale = pow((double)NbrSitesPerCell/Area,1.0/(double)Dimension);
      for (int i=0; i<Dimension; ++i)
	LatticeVectors[i]*=Rescale;
      Area = LatticeVectors.Determinant();
    }
  
  
  SubLatticeVectors.Resize(Dimension, NbrSitesPerCell);
  for (int i=0; i<NbrSitesPerCell; ++i)
    {
      sprintf(FieldName,"SubLatticeVector%d",i);
      int NbrComponents;
      double *Components;
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
      if (fabs(Rescale-1.0)>1e-10)
	SubLatticeVectors[i]*=Rescale;
    }

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
      if ( (*LatticePhases::Options)["external-parameters"] != NULL)
	TmpParameters = LatticePhases::Options->GetDoubles("external-parameters",TmpDimension);
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

      // enhance lattice descriptor if non-trivial parameters
      char *NewDescriptor = new char[strlen(this->Descriptor)+4+15*NbrExtParameters];
      sprintf(NewDescriptor,"%s_ep_%g",Descriptor,ExtParameters[0]);
      for (int i=1; i<NbrExtParameters; ++i)
	{
	  sprintf(NewDescriptor,"%s_%g",NewDescriptor,ExtParameters[i]);
	}
      delete [] this->Descriptor;
      this->Descriptor=NewDescriptor;
    }

  int TrueNbrExtParameters = NbrExtParameters;
  // determine connectivity, and tunnelling phases
  // two alternate methods: indicate individual tunnelling phases, or give a gauge choice
  this->HaveGauge=false;
  if ((LatticeDefinition["UseGauge"]!=NULL)&&
      ( (strcmp(LatticeDefinition["UseGauge"],"yes")==0) || (strcmp(LatticeDefinition["UseGauge"],"YES")==0)
	|| (strcmp(LatticeDefinition["UseGauge"],"true")==0) || (strcmp(LatticeDefinition["UseGauge"],"TRUE")==0) ))
    {
      this->HaveGauge=true;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAxx", this->GaugeAxx)==false)
	this->GaugeAxx=0.0;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAxy", this->GaugeAxy)==false)
	this->GaugeAxy=0.0;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAyx", this->GaugeAyx)==false)
	this->GaugeAyx=0.0;
      if (LatticeDefinition.GetAsSingleDouble ("GaugeAyy", this->GaugeAyy)==false)
	this->GaugeAyy=0.0;
      this->AbsBField=GaugeAxy-GaugeAyx;
      if (fabs(AbsBField)<1e-14)
	{
	  cout << "If using gauge mode, the magnetic field B=GaugeAxy - GaugeAyx needs to be non-zero" << endl;
	  exit(-1);
	}       
      if (this->Options->GetBoolean("normalize-lattice"))
	{
	  this->AbsBField=1.0/this->AbsBField;
	  this->GaugeAxx*=this->AbsBField;
	  this->GaugeAxy*=this->AbsBField;
	  this->GaugeAyx*=this->AbsBField;
	  this->GaugeAyy*=this->AbsBField;
	  this->AbsBField=1.0;
	}
      cout << "Gauge used : A= ("<<GaugeAxx<<"*x +"<<GaugeAyx<<"*y) e_x"
	   << " + ("<<GaugeAxy << "*x +"<<GaugeAyy<<"*y) e_y, "
	   << "field strength B="<<this->AbsBField<<endl;
    }
#ifdef DEBUG_LATTICE_PHASES
  cout << "Attention, the code is currently not functional for gauges involving both Axy and Ayx!"<<endl;

  for (int i=0; i<Dimension; ++i)
    cout << "LatticeVector["<<i<<"]="<<endl<<LatticeVectors[i];
  for (int i=0; i<NbrSitesPerCell; ++i)
    cout << "SubLatticeVector["<<i<<"]="<<endl<<SubLatticeVectors[i];
#endif

  this->NbrCells = this->PeriodicRep[0];
  for (int d=1; d<Dimension; ++d)
    this->NbrCells *= this->PeriodicRep[d];

  
  // Flux densities defined?
  if (LatticeDefinition["NbrFlux"]!=NULL)
    {
      this->PredefinedFluxFlag=true;
      this->PredefinedFlux = atoi(LatticeDefinition["NbrFlux"])*NbrCells;
    }
  else
    {
      this->PredefinedFluxFlag=false;
    }
  if (LatticeDefinition["ContinuousPhases"]!=NULL)
    {
      if (( (strcmp(LatticeDefinition["ContinuousPhases"],"yes")==0) || (strcmp(LatticeDefinition["ContinuousPhases"],"YES")==0)
	    ||(strcmp(LatticeDefinition["ContinuousPhases"],"true")==0) || (strcmp(LatticeDefinition["ContinuousPhases"],"TRUE")==0) ))
	{
	  this->ContinuousPhases=true;
	}
      else this->ContinuousPhases=false;  
    }
  else this->ContinuousPhases=false;

  char ***NeighborString;
  int NbrPairs;
  int *NbrValues;
  RealSymmetricMatrix NeighborsInCellMatrix(NbrSitesPerCell, true);
  RealSymmetricMatrix NeighborsInCellAmplitudes(NbrSitesPerCell, true);
  RealAntisymmetricMatrix TunnellingPhaseMatrix(NbrSitesPerCell, true);
  if (LatticeDefinition["NeighborsInCell"]!=NULL)
    {
      if (LatticeDefinition.GetAsStringMultipleArray ("NeighborsInCell", '|', ',', NeighborString, NbrPairs, NbrValues)==false)
	{
	  cout << "error while parsing NeighborsInCell in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);
	}
      for (int p=0; p<NbrPairs; ++p)
	{
	  int s1, s2;
	  if (NbrValues[p]<2)
	    {
	      cout << "error while decoding NeighborsInCell in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << "NeighborsInCell = s1,s2[,phaseA12=0.0,amplitudeParamaterID=0] | s3, s4[,phaseA34=0.0,amplitudeParamaterID=0] | ..."
		   << "Phases can either be indicated explitly, or will be deduced from gauge if defined"
		   << " or assumed to be one, otherwise"<<endl;
	      
	      exit(-1);
	    }
	  s1=strtod(NeighborString[p][0], NULL);
	  s2=strtod(NeighborString[p][1], NULL);
	  if ((s1<0)||(s1>=NbrSitesPerCell))
	    {
	      cout << "Attention: pair index "<<s1<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	    }
	  else
	    {
	      if ((s2<0)||(s2>=NbrSitesPerCell))
		cout << "Attention: pair index "<<s2<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	      else
		{
		  bool IsZero=false;
		  if (NbrValues[p]>3)
		    {
		      double AmplitudeNbr = strtod(NeighborString[p][3], NULL);
		      int index = (int)nearbyint(AmplitudeNbr);
		      if ((index < NbrExtParameters) && (fabs(ExtParameters[index]) < 1e-15))
			IsZero=true;
		    }
		  if (IsZero)
		    {
#ifdef DEBUG_LATTICE_PHASES
		      cout << "Skipping zero matrix element for pair ("<<s1<<", "<<s2<<")."<<endl;
#endif
		    }
		  else
		    {
		      double Index;
		      NeighborsInCellMatrix.GetMatrixElement(s1,s2,Index);
		      if (Index<0.9)
			{
			  cout << "New matrix element "<<s1<<", "<<s2<<" Index = "<<Index << endl;
			  NeighborsInCellMatrix.SetMatrixElement(s1,s2,1.0);
			  NeighborsInCellMatrix.GetMatrixElement(s1,s2,Index);
			  cout << "set new index "<<s1<<", "<<s2<<" Index = "<<Index << endl;
			  // determine tunnelling phase
			  if (NbrValues[p]<3)
			    TunnellingPhaseMatrix.SetMatrixElement(s1,s2,0.0);
			  else             
			    TunnellingPhaseMatrix.SetMatrixElement(s1,s2,strtod(NeighborString[p][2], NULL));
			  // determine tunnelling amplitude
			  if (NbrValues[p]<4)
			    NeighborsInCellAmplitudes.SetMatrixElement(s1,s2,0.0);
			  else
			    {
			      NeighborsInCellAmplitudes.SetMatrixElement(s1,s2,strtod(NeighborString[p][3], NULL));
			      if (strtod(NeighborString[p][3],NULL)+1.0>(double)TrueNbrExtParameters+1e-10)
				{
				  cout << "Not enough external parameters defined for Lattice Defition"<<endl;
				  exit(1);
				}
			    }
			}
		      else
			{
			  if ((this->HavePredefinedFlux()==false)||(this->GetPredefinedFlux()!=this->GetNbrSites()))
			    {
			      cout << "Cannot have multiple matrix elements connecting the same sites, if FluxDensity!=1.0"<<endl;
			      cout << "this->GetPredefinedFlux()=" << this->GetPredefinedFlux() << " this->GetNbrSites()=" << this->GetNbrSites()<<endl;
			      exit(1);
			    }

			  double Phase;
			  TunnellingPhaseMatrix.GetMatrixElement(s1, s2, Phase);
			  double AmplitudeNbr;
			  NeighborsInCellAmplitudes.GetMatrixElement(s1, s2, AmplitudeNbr);
			  Complex OldMatrixElement = ExtParameters[(int)nearbyint(AmplitudeNbr)] * Polar(1.0,-2.0*M_PI*Phase);
			  Complex NewMatrixElement = 1.0;
			  // determine tunnelling phase
			  Phase = 0.0;
			  AmplitudeNbr=0;
			  if (NbrValues[p]>2)
			    {
			      Phase = strtod(NeighborString[p][2], NULL);
#ifdef DEBUG_LATTICE_PHASES
			      cout << "read tunnelling phase for "<<s1<<","<<s2<<"= "<<Phase<<endl;
#endif
			    }
			  // determine tunnelling amplitude
			  if (NbrValues[p]>3)
			    {
			      AmplitudeNbr = strtod(NeighborString[p][3], NULL);
			      if (AmplitudeNbr+1.0>(double)TrueNbrExtParameters+1e-10)
				{
				  cout << "Not enough external parameters defined for Lattice Defition"<<endl;
				  exit(1);
				}
			      NewMatrixElement = ExtParameters[(int)nearbyint(AmplitudeNbr)];
			    }
			  NewMatrixElement *= Polar(1.0,-2.0*M_PI*Phase);
			  NewMatrixElement += OldMatrixElement;
			  Phase = -Arg(NewMatrixElement)/(2.0*M_PI);
			  double Amplitude = Norm(NewMatrixElement);
			  Index=-1;
			  for (int i=0; i<NbrExtParameters; ++i)
			    if (fabs(Amplitude-ExtParameters[i])<1e-13)
			      {
				Index = i;
				break;
			      }
			  if (Index<0)
			    {
			      double *NewParameters = new double[NbrExtParameters+1];
			      for (int i=0; i<NbrExtParameters; ++i)
				NewParameters[i]=ExtParameters[i];
			      NewParameters[NbrExtParameters]=Amplitude;
			      delete [] this->ExtParameters;
			      this->ExtParameters = NewParameters;
			      Index = NbrExtParameters;
			      ++NbrExtParameters;
#ifdef DEBUG_LATTICE_PHASES
			      cout << "Created new entry for absolute value of matrix elements: "<<ExtParameters[NbrExtParameters-1]<<endl;
#endif
			    }
			  TunnellingPhaseMatrix.SetMatrixElement(s1,s2,Phase);
			  NeighborsInCellAmplitudes.SetMatrixElement(s1,s2,(double)Index);
#ifdef DEBUG_LATTICE_PHASES
			  cout << "Merged matrix element: "<<s1<<", "<<s2<<" Index: "<<Index <<endl;
#endif
			}
		    }
		}
	    }
	}
      for (int i=0; i<NbrPairs; ++i)
	{
	  for (int j=0; j<NbrValues[i]; ++j)
	    delete [] NeighborString[i][j];
	  delete [] NeighborString[i];
	}
      delete [] NeighborString;
      delete [] NbrValues;
    }
  
#ifdef DEBUG_LATTICE_PHASES
  cout << "NeighborsInCell="<<endl<<NeighborsInCellMatrix;
  cout << "NeighborsInCellAmplitudes="<<endl<<NeighborsInCellAmplitudes;
  cout << "TunnellingPhaseMatrix="<<endl<<TunnellingPhaseMatrix;
#endif
  if (LatticeDefinition.GetAsStringMultipleArray ("NeighborCells", '|', ',', NeighborString, NbrPairs, NbrValues)==false)
    {
      cout << "error while parsing NeighborCells in " << this->Options->GetString("lattice-definition") << endl;
      exit(-1);
    }
  this->NbrNeighborCells = NbrPairs;
  this->NeighborCells = new int*[NbrPairs];
  for (int p=0; p<NbrNeighborCells; ++p)
    {
      if (NbrValues[p]!=Dimension)
	{
	  cout << "error while decoding NeighborCells in " << this->Options->GetString("lattice-definition") << endl;
	  cout << "Indicate coordinates of neighboring cells as vectors of length Dimension separated by commas"<<endl
	       << "Separate multiple neighboring cells by bars: "
	       << "NeighborCells = v_11,...,v_1d | ... | v_k1, v_kd | ..."<<endl;
	  exit(-1);
	}
      this->NeighborCells[p] = new int[Dimension];
      for (int i=0; i<Dimension; ++i)
	this->NeighborCells[p][i] = strtod(NeighborString[p][i], NULL);
    }
  for (int i=0; i<NbrPairs; ++i)
    {
      for (int j=0; j<NbrValues[i]; ++j)
	delete [] NeighborString[i][j];
      delete [] NeighborString[i];
    }
  delete [] NeighborString;
  delete [] NbrValues;
  RealMatrix **NeighborsAcrossBoundary = new RealMatrix*[NbrNeighborCells];
  RealMatrix **AmplitudesAcrossBoundary = new RealMatrix*[NbrNeighborCells];
  RealMatrix **PhasesAcrossBoundary = new RealMatrix*[NbrNeighborCells];
  for (int d=0; d<NbrNeighborCells; ++d)
    {
      sprintf(FieldName,"NeighborsAcrossBoundary%d",NeighborCells[d][0]);
      for (int i=1; i<Dimension; ++i)
	sprintf(FieldName,"%s_%d", FieldName, NeighborCells[d][i]);
      if (LatticeDefinition.GetAsStringMultipleArray (FieldName, '|', ',', NeighborString, NbrPairs, NbrValues)==false)
	{
	  cout << "could not parse "<<FieldName<<" in " << this->Options->GetString("lattice-definition")
	       << ", no connections added."<< endl;
	}
      NeighborsAcrossBoundary[d] = new RealMatrix(NbrSitesPerCell, NbrSitesPerCell, true);
      AmplitudesAcrossBoundary[d] = new RealMatrix(NbrSitesPerCell, NbrSitesPerCell, true);
      PhasesAcrossBoundary[d] = new RealMatrix(NbrSitesPerCell, NbrSitesPerCell, true);
      for (int p=0; p<NbrPairs; ++p)
	{
	  int s1, s2;
	  if (NbrValues[p]<2)
	    {
	      cout << "error while decoding "<<FieldName<<" in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate paires of neighboring sites separated by commas and different pairs by bars: "
		   << FieldName <<" = s1,s2[,phaseA12] | s3, s4 [,phaseA34] | ..."
		   << "Phases can either be indicated explitly, or will be deduced from gauge if defined"
		   << "Amplitudes are given as indices 0,1,2,... referring to the respective external parameters;"
		   << " both are assumed to be one, otherwise"<<endl;
	      exit(-1);
	    }
	  s1=strtod(NeighborString[p][0], NULL);
	  s2=strtod(NeighborString[p][1], NULL);
	  if ((s1<0)||(s1>=NbrSitesPerCell))
	    {
	      cout << "Attention: pair index "<<s1<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	    }
	  else
	    {
	      if ((s2<0)||(s2>=NbrSitesPerCell))
		cout << "Attention: pair index "<<s2<<" out of range, ignoring pair ("<<s1<<", "<<s2<<")."<<endl;
	      else
		{
		  bool IsZero=false;
		  if (NbrValues[p]>3)
		    {
		      double AmplitudeNbr = strtod(NeighborString[p][3], NULL);
		      int index = (int)nearbyint(AmplitudeNbr);
		      if ((index < NbrExtParameters) && (fabs(ExtParameters[index]) < 1e-15))
			IsZero=true;
		    }
		  if (IsZero)
		    {
#ifdef DEBUG_LATTICE_PHASES
		      cout << "Skipping zero matrix element for pair ("<<s1<<", "<<s2<<")."<<endl;
#endif
		    }
		  else
		    {
		      double Index;
		      NeighborsAcrossBoundary[d]->GetMatrixElement(s1,s2,Index);
		      if (Index<0.9)
			{
			  NeighborsAcrossBoundary[d]->SetMatrixElement(s1, s2, 1.0);
			  // determine tunnelling phase
			  if (NbrValues[p]<3)
			    {
			      PhasesAcrossBoundary[d]->SetMatrixElement(s1, s2, 0.0);
			    }
			  else
			    {
			      PhasesAcrossBoundary[d]->SetMatrixElement(s1,s2,strtod(NeighborString[p][2], NULL));
#ifdef DEBUG_LATTICE_PHASES
			      {
				double tmp;
				PhasesAcrossBoundary[d]->GetMatrixElement(s1,s2,tmp);
				cout << "Phase between neighbor sites "<<s1<<", "<<s2<<": "<<tmp<<endl;
			      }
#endif
			    }
			  // determine tunnelling amplitude
			  if (NbrValues[p]<4)
			    {
			      AmplitudesAcrossBoundary[d]->SetMatrixElement(s1, s2, 0.0);
			    }
			  else
			    {
			      AmplitudesAcrossBoundary[d]->SetMatrixElement(s1,s2,strtod(NeighborString[p][3], NULL));
			      if (strtod(NeighborString[p][3],NULL)+1.0>(double)NbrExtParameters+1e-10)
				{
				  cout << "Not enough external parameters defined for Lattice Defition"<<endl;
				  exit(1);
				}
#ifdef DEBUG_LATTICE_PHASES
			      {
				double tmp;
				AmplitudesAcrossBoundary[d]->GetMatrixElement(s1,s2,tmp);
				cout << "Amplitude between neighbor sites "<<s1<<", "<<s2<<": "<<tmp<<endl;
			      }
#endif
			    }
			}
		      else
			{
			  if ((this->HavePredefinedFlux()==false)||(this->GetPredefinedFlux()!=this->GetNbrSites()))
			    {
			      cout << "Cannot have multiple matrix elements connecting the same sites, if FluxDensity!=1.0"<<endl;
			      cout << "this->GetPredefinedFlux()=" << this->GetPredefinedFlux() << " this->GetNbrSites()=" << this->GetNbrSites()<<endl;
			      exit(1);
			    }

			  double Phase;
			  PhasesAcrossBoundary[d]->GetMatrixElement(s1, s2, Phase);
			  double AmplitudeNbr;
			  NeighborsAcrossBoundary[d]->GetMatrixElement(s1, s2, AmplitudeNbr);
			  Complex OldMatrixElement = ExtParameters[(int)nearbyint(AmplitudeNbr)] * Polar(1.0,-2.0*M_PI*Phase);
			  Complex NewMatrixElement = 1.0;
			  // determine tunnelling phase
			  Phase = 0.0;
			  AmplitudeNbr=0;
			  if (NbrValues[p]>2)
			    {
			      Phase = strtod(NeighborString[p][2], NULL);
#ifdef DEBUG_LATTICE_PHASES
			      cout << "read tunnelling phase for "<<s1<<","<<s2<<"= "<<Phase<<endl;
#endif
			    }
			  // determine tunnelling amplitude
			  if (NbrValues[p]>3)
			    {
			      AmplitudeNbr = strtod(NeighborString[p][3], NULL);
			      if (AmplitudeNbr+1.0>(double)TrueNbrExtParameters+1e-10)
				{
				  cout << "Not enough external parameters defined for Lattice Defition"<<endl;
				  exit(1);
				}
			      NewMatrixElement = ExtParameters[(int)nearbyint(AmplitudeNbr)];
			    }
			  NewMatrixElement *= Polar(1.0,-2.0*M_PI*Phase);
			  NewMatrixElement += OldMatrixElement;
			  Phase = -Arg(NewMatrixElement)/(2.0*M_PI);
			  double Amplitude = Norm(NewMatrixElement);
			  Index=-1;
			  for (int i=0; i<NbrExtParameters; ++i)
			    if (fabs(Amplitude-ExtParameters[i])<1e-13)
			      {
				Index = i;
				break;
			      }
			  if (Index<0)
			    {
			      double *NewParameters = new double[NbrExtParameters+1];
			      for (int i=0; i<NbrExtParameters; ++i)
				NewParameters[i]=ExtParameters[i];
			      NewParameters[NbrExtParameters]=Amplitude;
			      delete [] this->ExtParameters;
			      this->ExtParameters = NewParameters;
			      Index = NbrExtParameters;
			      ++NbrExtParameters;
			      cout << "Created new entry for absolute value of matrix elements: "<<ExtParameters[NbrExtParameters-1]<<endl;
			    }
			  PhasesAcrossBoundary[d]->SetMatrixElement(s1,s2,Phase);
			  AmplitudesAcrossBoundary[d]->SetMatrixElement(s1,s2,(double)Index);
#ifdef DEBUG_LATTICE_PHASES
			  cout << "Merged matrix element: "<<s1<<", "<<s2<<" Index: "<<Index <<endl;
#endif
			}
		    }
		}
	    }
	}
      for (int i=0; i<NbrPairs; ++i)
	{
	  for (int j=0; j<NbrValues[i]; ++j)
	    delete [] NeighborString[i][j];
	  delete [] NeighborString[i];
	}
      delete [] NbrValues;
      delete [] NeighborString;
      
#ifdef DEBUG_LATTICE_PHASES
      cout << FieldName<<"="<<endl<<*NeighborsAcrossBoundary[d];
      cout << "Amplitudes("<<FieldName<<")="<<endl<<*AmplitudesAcrossBoundary[d];
      cout << "Phases("<<FieldName<<")="<<endl<<*(PhasesAcrossBoundary[d]);
#endif
    }
  
  this->Neighbors = new int*[NbrSites];
  this->NeighborShift = new int**[NbrSites];
  this->TunnellingPhases = new double*[NbrSites];
  this->TunnellingAmplitudes = new double*[NbrSites];
  this->NbrNeighbors = new int[NbrSites];
  for (int i=0; i<NbrSites; ++i)
    this->NbrNeighbors[i] = 0;
  int *TmpNeighbors = new int[NbrNeighborCells*NbrSites];
  int **TmpNeighborShift = new int*[NbrNeighborCells*NbrSites];
  for (int i=0; i<NbrNeighborCells*NbrSites; ++i)
    {
      TmpNeighborShift[i]=new int[Dimension];
      for (int j=0; j<Dimension; ++j)
	TmpNeighborShift[i][j]=0;
    }
  double *TmpPhases = new double[NbrNeighborCells*NbrSites];
  int *TmpAmplitudes = new int[NbrNeighborCells*NbrSites];
  
  

  int *CellCoordinates = new int[Dimension];
  int *CellCoordinates2 = new int[Dimension];
  int *Translation3 = new int[Dimension];
  
  for (int c=0; c<NbrCells; ++c)
    {
      this->GetCellCoordinates(c, CellCoordinates);
      
#ifdef DEBUG_LATTICE_PHASES
      cout << "Cell "<<c<<":"<< CellCoordinates[0]<<", "<<CellCoordinates[1]<<endl;
#endif
      int Site1, Site2, Site3;
      double TmpD;
      for (int i=0; i<NbrSitesPerCell; ++i)
	{
	  Site1 = this->GetSiteNumber(c, i);
	  
#ifdef DEBUG_LATTICE_PHASES
	  cout << "Site 1="<<Site1<<endl;
#endif
	  for (int j=0; j<NbrSitesPerCell; ++j)
	    {
	      Site2 = this->GetSiteNumber(c, j);
	      
	      if (NeighborsInCellMatrix(i,j)>0.0)
		{		  
		  TmpNeighbors[NbrNeighbors[Site1]]=Site2;
		  NeighborsInCellAmplitudes.GetMatrixElement(Site1%NbrSitesPerCell,Site2%NbrSitesPerCell,TmpD);
		  TmpAmplitudes[NbrNeighbors[Site1]]=(int)nearbyint(TmpD);
		  if (this->HaveGauge)
		    TmpPhases[NbrNeighbors[Site1]] = GetTunnellingPhaseFromGauge(Site1, Site2);
		  else
		    TunnellingPhaseMatrix.GetMatrixElement(Site1%NbrSitesPerCell,Site2%NbrSitesPerCell,TmpPhases[NbrNeighbors[Site1]]);
		  for (int r=0; r<Dimension; ++r)
		    TmpNeighborShift[NbrNeighbors[Site1]][r]=0;
#ifdef DEBUG_LATTICE_PHASES
		  cout << "Neighbors "<<Site1<<"->"<<Site2<<" with phase "<<TmpPhases[NbrNeighbors[Site1]];
		  if (NbrExtParameters>0)
		    cout <<", and with amplitudeID"<<TmpAmplitudes[NbrNeighbors[Site1]];
		  cout<<endl;
#endif
		  ++NbrNeighbors[Site1];
		}
	      for (int d=0; d<NbrNeighborCells; ++d)
		{
		  if ((*(NeighborsAcrossBoundary[d]))(i,j)!=0.0)
		    {
		      for (int k=0; k<Dimension; ++k)
			{
			  CellCoordinates2[k]=CellCoordinates[k]+NeighborCells[d][k];
#ifdef DEBUG_LATTICE_PHASES
			  cout << "CellCoordinates["<<k<<"]="<<CellCoordinates[k]<< ", "
			       << "NeighborCells["<<d<<", "<<k<<"]="<<NeighborCells[d][k]<< ", "
			       << "CellCoordinates2["<<k<<"]="<<CellCoordinates2[k]<<endl;
#endif
			}
		      Site3 = this->GetSiteNumber(CellCoordinates2, j, Translation3);
#ifdef DEBUG_LATTICE_PHASES
		      cout << "Translation3= ["<<Translation3[0]<<", "<<Translation3[1]<<"]"<<endl;
#endif
		      TmpNeighbors[NbrNeighbors[Site1]]=Site3;
		      for (int r=0; r<Dimension; ++r)
			TmpNeighborShift[NbrNeighbors[Site1]][r]=Translation3[r]/PeriodicRep[r];
		      AmplitudesAcrossBoundary[d]->GetMatrixElement(Site1%NbrSitesPerCell,Site3%NbrSitesPerCell,TmpD);
		      TmpAmplitudes[NbrNeighbors[Site1]]=(int)nearbyint(TmpD);
		      if (this->HaveGauge)
			{
#ifdef DEBUG_LATTICE_PHASES
			  cout << "Evaluating translation phase:"<<endl;
#endif
			  TmpPhases[NbrNeighbors[Site1]] = GetTunnellingPhaseFromGauge(Site1, Site3, Translation3);
			}
		      else
			PhasesAcrossBoundary[d]->GetMatrixElement(Site1%NbrSitesPerCell, Site3%NbrSitesPerCell, TmpPhases[NbrNeighbors[Site1]]);

#ifdef DEBUG_LATTICE_PHASES
		      cout << "additional neighbors "<<Site1<<"->"<<Site3<<" from NeigborCell "<<d<<" at "<<
			CellCoordinates2[0]<<", "<<CellCoordinates2[1]<<", "<<j<<" : Site 3="<<Site3
			   <<" with translation "<<Translation3[0];
		      for (int i=1; i<Dimension; ++i) cout << " "<<Translation3[i];
		      if (TmpAmplitudes[NbrNeighbors[Site1]]!=0)
			cout << ", amplitudeID "<<TmpAmplitudes[NbrNeighbors[Site1]];
		      cout << " and phase "<<TmpPhases[NbrNeighbors[Site1]]<<" Shift "<<TmpNeighborShift[NbrNeighbors[Site1]][0]
			   << ","<<TmpNeighborShift[NbrNeighbors[Site1]][1]<<endl;
#endif
		      ++NbrNeighbors[Site1];
		    }
		}	      
	    }
	  if (NbrNeighbors[Site1]>0)
	    {
	      Neighbors[Site1] = new int[NbrNeighbors[Site1]];
	      TunnellingPhases[Site1] = new double[NbrNeighbors[Site1]];
	      TunnellingAmplitudes[Site1] = new double[NbrNeighbors[Site1]];
	      NeighborShift[Site1] = new int*[NbrNeighbors[Site1]];
	      for (int k=0; k<NbrNeighbors[Site1]; ++k)
		NeighborShift[Site1][k] = new int[Dimension];
	    }
	  else Neighbors[Site1] = NULL;
	  for (int k=0; k<NbrNeighbors[Site1]; ++k)
	    {
	      Neighbors[Site1][k] = TmpNeighbors[k];
	      TunnellingPhases[Site1][k] = TmpPhases[k];
	      TunnellingAmplitudes[Site1][k] = ExtParameters[TmpAmplitudes[k]];
	      for (int r=0; r<Dimension; ++r)
		NeighborShift[Site1][k][r]=TmpNeighborShift[k][r];
	    }
	}
    }
  this->OneParticlePotentials=NULL;
  NbrOneParticlePotentials=0;
  if (LatticeDefinition["RandomPotentials"]!=NULL)
    {
      double randomStrength = strtod(LatticeDefinition["RandomPotentials"],NULL);
      this->OneParticlePotentials = new double[NbrSites];
      this->OneParticlePotentialPositions = new int[NbrSites];
      NumRecRandomGenerator G;
      G.UseTimeSeed();
      for (int i=0; i<NbrSites; ++i)
	{
	  OneParticlePotentials[i]=randomStrength*(-0.5+G.GetRealRandomNumber());
	  this->OneParticlePotentialPositions[i] = i;
	}
      NbrOneParticlePotentials = NbrSites;
    }
  else if (LatticeDefinition["LocalPotentials"]!=NULL)
    {
      cout << "Reading local potentials"<<endl;
      if (LatticeDefinition.GetAsStringMultipleArray ("LocalPotentials", '|', ',', NeighborString, NbrPairs, NbrValues)==false)
	{
	  cout << "error while parsing LocalPotentials in " << this->Options->GetString("lattice-definition") << endl;
	  exit(-1);
	}
      this->OneParticlePotentials = new double[NbrSites];
      this->OneParticlePotentialPositions = new int[NbrSites];
      for (int i=0; i<NbrSites; ++i) OneParticlePotentials[i]=0.0;
      for (int p=0; p<NbrPairs; ++p)
	{
	  if (NbrValues[p]!=2)
	    {
	      cout << "error while decoding LocalPotentials in " << this->Options->GetString("lattice-definition") << endl;
	      cout << "Indicate index of sites and potential at this site separated by commas"<<endl
		   << "Separate multiple entries by bars: "
		   << "NeighborCells = s_1,V_1 | ... | s_N, V_N | ..."<<endl;
	      exit(-1);
	    }
	  int SiteIndex = atoi(NeighborString[p][0]);
	  if ((SiteIndex<0)||(SiteIndex>=NbrSites))
	    {
	      cout << "SiteIndex "<<SiteIndex<<" out of range in LatticePhases - ignoring potential indicated for this site"<<endl;
	    }
	  else
	    {
	      this->OneParticlePotentials[NbrOneParticlePotentials] = strtod(NeighborString[p][1], NULL);
	      this->OneParticlePotentialPositions[NbrOneParticlePotentials] = SiteIndex;
	      ++NbrOneParticlePotentials;
	    }

	}
      for (int p=0; p<NbrPairs; ++p)
	{
	  for (int q=0; q<NbrValues[p]; ++q)
	    delete [] NeighborString[p][q];
	  delete [] NeighborString[p];
	}
      delete [] NbrValues;
      delete [] NeighborString;
    }
  
  cout << "LatticePhases created"<<endl;

  for (int d=0; d<NbrNeighborCells; ++d)
    {
      delete NeighborsAcrossBoundary[d];
      delete AmplitudesAcrossBoundary[d];
      delete PhasesAcrossBoundary[d];
    }
  delete [] NeighborsAcrossBoundary;
  delete [] AmplitudesAcrossBoundary;
  delete [] PhasesAcrossBoundary;
  delete [] CellCoordinates;
  delete [] CellCoordinates2;
  delete [] Translation3;
  delete [] FieldName;
  delete [] TmpNeighbors;
  for (int i=0; i<NbrNeighborCells*NbrSites; ++i)
    delete [] TmpNeighborShift[i];
  delete [] TmpNeighborShift;
  delete [] TmpPhases;
  delete [] TmpAmplitudes;
}

// destructor
//
LatticePhases::~LatticePhases()
{
  if (NbrSites!=0)
    {
      delete [] this->PeriodicRep;
      for (int i=0; i<NbrSites; ++i)
	{
	  if (this->Neighbors[i]!=NULL)
	    delete [] this->Neighbors[i];
	  if (this->Neighbors[i]!=NULL)
	    {
	      for (int j=0; j<NbrNeighbors[i]; ++j)
		delete [] NeighborShift[i][j];
	      delete [] NeighborShift[i];
	    }
	  if (this->TunnellingPhases[i]!=NULL)
	    delete [] this->TunnellingPhases[i];
	  if (this->TunnellingAmplitudes[i]!=NULL)
	    delete [] this->TunnellingAmplitudes[i];
	}      
      delete [] this->Neighbors;
      delete [] NeighborShift;
      delete [] this->TunnellingPhases;
      delete [] this->TunnellingAmplitudes;
      delete [] this->NbrNeighbors;
      for (int i=0; i<NbrNeighborCells; ++i)
	delete [] this->NeighborCells[i];
      delete [] this->NeighborCells;
      delete [] this->Descriptor;
      delete [] this->ExtParameters;
      if (NbrOneParticlePotentials>0)
	{
	  delete[] OneParticlePotentials;
	  delete[] OneParticlePotentialPositions;
	}
    }
}


// get cell coordinates given the number of the unit cell
// nbrCell = cell to be looked up
// cellCoordinates = resulting coordinates, has to be reserved prior to call
void LatticePhases::GetCellCoordinates(int nbrCell, int *cellCoordinates)
{
  int Divisor=1;
  while (nbrCell<0) nbrCell+=NbrCells;
  for (int i=0; i<Dimension; ++i)
    {
      cellCoordinates[i] = (nbrCell/Divisor)%this->PeriodicRep[i];
      Divisor*=this->PeriodicRep[i];
    }
}

// get cell coordinates given the number of the unit cell
// nbrSite = cell to be looked up
// cellCoordinates = resulting coordinates, has to be reserved prior to call
// sublattice = resulting sublattice
void LatticePhases::GetSiteCoordinates(int nbrSite, int *cellCoordinates, int &sublattice)
{
  while (nbrSite<0) nbrSite+=NbrSites;
  sublattice = nbrSite%NbrSitesPerCell;
  int Divisor=NbrSitesPerCell;
  for (int i=0; i<Dimension-1; ++i)
    {
      cellCoordinates[i] = (nbrSite/Divisor)%this->PeriodicRep[i];
      Divisor*=this->PeriodicRep[i];
    }
  cellCoordinates[Dimension-1] = (nbrSite/Divisor)%this->PeriodicRep[Dimension-1];
}

// retrieve the position of a given site
// cellCoordinates = resulting coordinates, has to be reserved prior to call
// sublattice = resulting sublattice  
RealVector LatticePhases::GetSitePosition(int *cellCoordinates, int sublattice)
{
  RealVector Position(this->Dimension,true);
  for (int i=0; i<Dimension; ++i)
    {
      Position.AddLinearCombination((double)cellCoordinates[i],LatticeVectors[i]);
    }
  Position.AddLinearCombination(1.0,SubLatticeVectors[sublattice]);
  // cout << "Position of site "<<cellCoordinates[0]<<", "<<cellCoordinates[1]<<", "<<sublattice<<endl<<Position;
  return Position;
}



// get number of a site in cell nbrCell
// cellCoordinates = coordinates of cell to be addressed
// sublattice = sublattice index
int LatticePhases::GetSiteNumber(int *cellCoordinates, int sublattice)
{
  int Result=this->Periodize(cellCoordinates[Dimension-1],Dimension-1);
  for (int i=Dimension-2; i>-1; --i)
    {
      Result*=this->PeriodicRep[i];
      Result+=this->Periodize(cellCoordinates[i],i);
    }
  Result*=NbrSitesPerCell;
  Result+=sublattice;
  return Result%NbrSites;
}

// get number of a site in cell nbrCell, and return translation vector back into the simulation cell
// cellCoordinates = coordinates of cell to be addressed
// sublattice = sublattice index
// translation = vector of tranlation back into simulation cell
int LatticePhases::GetSiteNumber(int *cellCoordinates, int sublattice, int *translation)
{
#ifdef DEBUG_LATTICE_PHASES
  cout << "Periodizing entry"<<Dimension-1<<endl;
#endif
  int Result=this->Periodize(cellCoordinates[Dimension-1], Dimension-1, translation[Dimension-1]);
  for (int i=Dimension-2; i>-1; --i)
    {
      Result*=this->PeriodicRep[i];
      Result+=this->Periodize(cellCoordinates[i], i, translation[i]);
    }
  Result*=NbrSitesPerCell;
  Result+=sublattice;
  return Result%NbrSites;
}


// request address of partners of site
// nbrSite = number of site whose partners to request
// nbrNeighbors = number of partners found
// Neighbors = array to partner sites
// phases = values of phase for tunnelling matrix element
// periodicTranslations = translations into the fundamental domain
// amplitudes = (optional) amplitudes of tunneling terms
void LatticePhases::GetNeighbors(int nbrSite, int &nbrNeighbors, int * &neighbors, double * &phases, int **&periodicTranslations, double *&amplitudes)
{
  if ((nbrSite>-1)&&(nbrSite<NbrSites))
    {
      neighbors = this->Neighbors[nbrSite];
      nbrNeighbors = this->NbrNeighbors[nbrSite];
      phases = this->TunnellingPhases[nbrSite];
      periodicTranslations = this->NeighborShift[nbrSite];
      if (NbrExtParameters>0)
	amplitudes = this->TunnellingAmplitudes[nbrSite];
      else
	amplitudes = NULL;
    }
  else
    {
      nbrNeighbors = 0;
      neighbors = NULL;
      phases = NULL;
      amplitudes = NULL;
    }
}

// get total number of hopping terms
int LatticePhases::GetNbrHoppingTerms()
{
  int sum=0;
  for (int i=0; i<NbrSites; ++i)
    sum += this->NbrNeighbors[i];
  return sum;
}

// get total number of local potential terms
int LatticePhases::GetNbrLocalPotentials()
{
  return NbrOneParticlePotentials;
}


// get a string describing the lattice geometry
// 
char *LatticePhases::GeometryString()
{
  char *rst = new char[strlen(this->Descriptor)+20];
  sprintf(rst,"%s_%d", this->Descriptor, this->PeriodicRep[0]);
  for (int i=1; i<Dimension; ++i)
    sprintf(rst,"%sx%d", rst, this->PeriodicRep[i]);
  return rst;
}


// request if single-particle potentials are defined
bool LatticePhases::HaveOneParticlePotentials()
{
  if (this->OneParticlePotentials!=NULL)
    return true;
  else
    return false;
}

// request single-particle potentials 
double* LatticePhases::GetOneParticlePotentials(int &nbrPotentials, int* &positions)
{
  nbrPotentials = NbrOneParticlePotentials;
  positions = OneParticlePotentialPositions;
  return OneParticlePotentials;
}

// get mapping of lattice sites under translations by multiples of the lattice vectors
// t = vector indicating translations in units of lattice vectors
// mappings = mapping of site numbers under this translation
// phases = eventual phases picked up by this translation operator
// solenoidFlux = solenoid fluxes in each period of lattice
void LatticePhases::GetTranslations(int *t, int* mappings, Complex *phases, double* solenoidFlux)
{
  int CellCoordinates[Dimension];
  int Translation[Dimension];
  int Sublattice;
  for (int i=0; i<NbrSites; ++i)
    {
      this->GetSiteCoordinates(i, CellCoordinates, Sublattice);
      for (int j=0; j<Dimension; ++j)
	CellCoordinates[j]+=t[j];
      mappings[i]=this->GetSiteNumber(CellCoordinates, Sublattice, Translation);
      if (this->HaveGauge)
	phases[i]=GetTranslationPhaseFromGauge(i, mappings[i], Translation);
      else
	{
	  phases[i]=1.0;
	  for (int j=0; j<Dimension; ++j)
	    {
	      //	      if (Translation[j]!=0)
	      //		cout <<"Translating ["<<i<<"] in "<<j<<": "<<Translation[j]<<" div: " << solenoidFlux[j]*(Translation[j]/PeriodicRep[j])<<endl;
	      if (abs(Translation[j])/PeriodicRep[j]>0)
		phases[i] *= Polar(1.0,solenoidFlux[j]*(Translation[j]/PeriodicRep[j]));
	    }
	}
    }
#ifdef DEBUG_LATTICE_PHASES
  cout << "Translating general lattice by "<<t[0]<<", "<<t[1]<<endl;
  for (int i=0; i<NbrSites; ++i)
    {
      cout << i<< " -> " << mappings[i];
      if (Norm(phases[i]-1.0)>1e-6)
	cout << " (phase: "<<Arg(phases[i])<<")";
      cout << endl;
    }
#endif
}




// add an option group containing all options related to the LatticeGeometry options
//
// manager = pointer to the option manager
void LatticePhases::AddOptionGroup(OptionManager* manager)
{
  LatticePhases::Options = manager;
  OptionGroup* LatticeGroup  = new OptionGroup ("lattice options");
  (*(LatticePhases::Options)) += LatticeGroup;

  (*LatticeGroup) += new SingleStringOption ('L', "lattice-definition", "File defining the geometry of the lattice");
  (*LatticeGroup) += new MultipleIntegerOption ('C', "cells-repeat", "number of times unit cell is repeated in the x-, y-,..., dim- directions of the lattice (overriding default given in definition)", ',');
  (*LatticeGroup) += new BooleanOption ('\n', "normalize-lattice", "normalize unit cell area to #of sites, and field strength to one");
  (*LatticeGroup) += new MultipleDoubleOption ('\n', "external-parameters", "values of external parameters for lattice definition", ',');
}



OptionManager* LatticePhases::Options=NULL;

int MyRound(double a) {
return int(a + 0.5);
}

// simple sort algorithm
// array = integer array to be sorted
// length = length of array
void LatticePhases::ArraySort(int* array, int length)
{
  int inc = MyRound(length/2.0);
  int tmpI;
  while (inc > 0)
    {
      for (int i = inc; i< length; ++i)
	{
	  tmpI = array[i];
	  int j = i;
	  while ((j>=inc) && ( array[j-inc] < tmpI) )
	    {
	      array[j] = array[j - inc];
	      j = j - inc;
	    }
	  array[j] = tmpI;
	}
      inc = round(inc / 2.2);
    }
}

// calculate the tunnelling phase between two given sites from the gauge
// s1 = start site
// s2 = end site
// cellTranslation = indicating whether translation across a boundary ocurred
// return = relative phase
double LatticePhases::GetTunnellingPhaseFromGauge(int s1, int s2, int *cellTranslation)
{
  if (this->HaveGauge)
    {
      double Result=0.0;
      // calculate site coordinates
      int *S1Coordinates = new int[this->Dimension];
      int *S2Coordinates = new int[this->Dimension];
      int S1Sublattice, S2Sublattice;
      this->GetSiteCoordinates(s1, S1Coordinates, S1Sublattice);
      this->GetSiteCoordinates(s2, S2Coordinates, S2Sublattice);
      RealVector Position1(this->Dimension,true);
      RealVector Position2(this->Dimension,true);      
      RealVector Translation(this->Dimension,true);
      for (int i=0; i<Dimension; ++i)
	{
	  Position1.AddLinearCombination((double)S1Coordinates[i],LatticeVectors[i]);
	  Position2.AddLinearCombination((double)S2Coordinates[i],LatticeVectors[i]);
	  if (cellTranslation!=NULL)
	    Translation.AddLinearCombination((double)cellTranslation[i],LatticeVectors[i]);
	}      
      Position2.AddLinearCombination(-1.0,Translation);
      RealVector CellPosition2(this->Dimension);
      CellPosition2.Copy(Position2);
      Position1.AddLinearCombination(1.0,SubLatticeVectors[S1Sublattice]);
      Position2.AddLinearCombination(1.0,SubLatticeVectors[S2Sublattice]);
      CellPosition2.AddLinearCombination(1.0,SubLatticeVectors[S2Sublattice]);
      
#ifdef DEBUG_LATTICE_PHASES
      cout << "Position1="<<endl<<Position1;
      cout << "Position2="<<endl<<Position2;
      cout << "Translation="<<endl<<Translation;
#endif
      if (this->Dimension==2)
	{
	  Result += 0.5*GaugeAxx*(Position2[0]*Position2[0]-Position1[0]*Position1[0]);
	  Result += 0.5*GaugeAyx*(Position2[0]-Position1[0])*(Position2[1]+Position1[1]);
	  Result += 0.5*GaugeAxy*(Position2[1]-Position1[1])*(Position2[0]+Position1[0]);
	  Result += 0.5*GaugeAyy*(Position2[1]*Position2[1]-Position1[1]*Position1[1]);
	  // xxx: check signs and prefactors here!
#ifdef DEBUG_LATTICE_PHASES
	  cout << "Raw phase="<<Result;
#endif
	  if (Translation.SqrNorm()>1e-14)
	    {
	      if (Translation[0]>0.0)
		{
		  // translate in x-direction first:
		  double MagneticTranslation =(GaugeAxy*Translation[0]+GaugeAyy*Translation[1])*CellPosition2[1]; // (CellPosition2[1]+0.5*Translation[1]);
		  // then translate in y-direction
		  MagneticTranslation+=(GaugeAxx*Translation[0]+GaugeAyx*Translation[1])*(CellPosition2[0]+Translation[0]); // (CellPosition2[0]+Translation[0]+0.5*Translation[0]);
		  
		  Result += MagneticTranslation;
#ifdef DEBUG_LATTICE_PHASES
		  cout << ", magnetic translation: "<<MagneticTranslation;
#endif
		}
	      else
		{
		  // translate in y-direction first:
		  double MagneticTranslation =(GaugeAxx*Translation[0]+GaugeAyx*Translation[1])*CellPosition2[0]; // (CellPosition2[0]+0.5*Translation[0]);
		  // translate in x-direction first:
		  MagneticTranslation+= (GaugeAxy*Translation[0]+GaugeAyy*Translation[1])*(CellPosition2[1]+Translation[1]); // (CellPosition2[1]+Translation[1]+0.5*Translation[1]);
		  Result += MagneticTranslation;
#ifdef DEBUG_LATTICE_PHASES
		  cout << ", magnetic translation: "<<MagneticTranslation;
#endif
		}
	    }
#ifdef DEBUG_LATTICE_PHASES
	  cout << ", after corrections"<<Result<<endl;
#endif
	}
      else
	{
	  cout << "Need to define LatticePhases::GetTunnellingPhaseFromGauge for dimension d>2"<<endl;
	}
      delete [] S1Coordinates;
      delete [] S2Coordinates;
      return Result;
    }
  else
    return 0.0;
}


// calculate the tunnelling phase between two given sites from the gauge
// s1 = start site
// s2 = end site
// cellTranslation = indicating whether translation across a boundary ocurred
// return = relative phase
double LatticePhases::GetTranslationPhaseFromGauge(int s1, int s2, int *cellTranslation)
{
  if (this->HaveGauge)
    {
      double Result=0.0;
      // calculate site coordinates
      int *S1Coordinates = new int[this->Dimension];
      int *S2Coordinates = new int[this->Dimension];
      int S1Sublattice, S2Sublattice;
      this->GetSiteCoordinates(s1, S1Coordinates, S1Sublattice);
      this->GetSiteCoordinates(s2, S2Coordinates, S2Sublattice);
      RealVector Position1(this->Dimension,true);
      RealVector Position2(this->Dimension,true);      
      RealVector Translation(this->Dimension,true);
      for (int i=0; i<Dimension; ++i)
	{
	  Position1.AddLinearCombination((double)S1Coordinates[i],LatticeVectors[i]);
	  Position2.AddLinearCombination((double)S2Coordinates[i],LatticeVectors[i]);
	  if (cellTranslation!=NULL)
	    Translation.AddLinearCombination((double)cellTranslation[i],LatticeVectors[i]);
	}      
      Position2.AddLinearCombination(-1.0,Translation);
      RealVector CellPosition2(this->Dimension);
      CellPosition2.Copy(Position2);
      Position1.AddLinearCombination(1.0,SubLatticeVectors[S1Sublattice]);
      Position2.AddLinearCombination(1.0,SubLatticeVectors[S2Sublattice]);
      CellPosition2.AddLinearCombination(1.0,SubLatticeVectors[S2Sublattice]);
      
      if (this->Dimension==2)
	{
	  if (Translation.SqrNorm()>1e-14)
	    {
	      if (Translation[0]>0.0)
		{
		  // translate in x-direction first:
		  double MagneticTranslation =(GaugeAxy*Translation[0]+GaugeAyy*Translation[1])*CellPosition2[1]; // (CellPosition2[1]+0.5*Translation[1]);
		  // then translate in y-direction
		  MagneticTranslation+=(GaugeAxx*Translation[0]+GaugeAyx*Translation[1])*(CellPosition2[0]+Translation[0]); // (CellPosition2[0]+Translation[0]+0.5*Translation[0]);
		  
		  Result += MagneticTranslation;
#ifdef DEBUG_LATTICE_PHASES
		  cout << ", magnetic translation: "<<MagneticTranslation;
#endif
		}
	      else
		{
		  // translate in y-direction first:
		  double MagneticTranslation =(GaugeAxx*Translation[0]+GaugeAyx*Translation[1])*CellPosition2[0]; // (CellPosition2[0]+0.5*Translation[0]);
		  // translate in x-direction first:
		  MagneticTranslation+= (GaugeAxy*Translation[0]+GaugeAyy*Translation[1])*(CellPosition2[1]+Translation[1]); // (CellPosition2[1]+Translation[1]+0.5*Translation[1]);
		  Result += MagneticTranslation;
		}
	    }
	}
      else
	{
	  cout << "Need to define LatticePhases::GetTunnellingPhaseFromGauge for dimension d>2"<<endl;
	}
      delete [] S1Coordinates;
      delete [] S2Coordinates;
      return Result;
    }
  else
    return 0.0;
}
