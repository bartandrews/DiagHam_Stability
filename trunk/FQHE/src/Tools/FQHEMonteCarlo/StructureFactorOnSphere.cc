////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2022 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//            class for a structure factor on the sphere geometry             //
//                                                                            //
//                        last modification : 29/04/2022                      //
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
#include "StructureFactorOnSphere.h"
#include <iostream>

using std::ios;
using std::ios_base;
using std::endl;


StructureFactorOnSphere::StructureFactorOnSphere()
{
  this->NbrStructureFactors=0;
}

// constructor
// resolution = total number of bins
// nbrStructureFactors = number of structure factors to calculate
StructureFactorOnSphere::StructureFactorOnSphere(int nbrFlux, int nbrStructureFactors)
{
  this->Type=AbstractObservable::VectorValued|AbstractObservable::RealObservableT;
  this->Measures = 0.0;
  
  this->LegendreBasis = new LegendrePolynomials(this->NbrStructureFactors);
  this->CurrentLegendrePolynomials = new double[this->NbrStructureFactors+1];
  this->StructureFactors = new WeightedRealVectorObservable(this->NbrStructureFactors, 16);
  this->NbrFlux=nbrFlux;
  this->NbrParticles=0;
  this->Radius = sqrt(0.5*(double)NbrFlux); // the radius is also the inverse magnetic length
  this->PrintFlag=true;
}
  
// destructor
StructureFactorOnSphere::~StructureFactorOnSphere()
{
  if (this->NbrStructureFactors!=0)
    {
      delete this->LegendreBasis;
      delete [] this->CurrentLegendrePolynomials; 
      delete this->StructureFactors;
    }
}

// call to make an observation
// weight = relative weight of this sample
void StructureFactorOnSphere::RecordValue(double weight)
{
  double r_ij;
  int index;
  this->Measures+=weight;
  const RealSymmetricMatrix &Distances = this->System->GetDistances();
  for( int i=1;i<this->NbrParticles;++i)
    for(int j=0; j<i;++j)
      {
	r_ij=Distances(i,j); // equals r_ij=sin(\theta_ij/2)
	double cosphi = 1.0-2.0*r_ij*r_ij;
	this->LegendreBasis->GetPolynomials(cosphi, CurrentLegendrePolynomials);
	this->StructureFactors->Observe(CurrentLegendrePolynomials,weight);
      }
}

// print legend to the given stream
// all = flag indicating whether to print all, or shortened information
void StructureFactorOnSphere::PrintLegend(std::ostream &output, bool all)
{
  output << "# n\tS_n"<<endl;
}

// print status to the given stream
// all = flag indicating whether to print all, or shortened information
void StructureFactorOnSphere::PrintStatus(std::ostream &output, bool all)
{
  // no action, for now
}

// request whether observable should be printed
//
bool StructureFactorOnSphere::IncludeInPrint()
{
  return this->PrintFlag;
}

// set print status
//
void StructureFactorOnSphere::IncludeInPrint(bool newStatus)
{
  this->PrintFlag=newStatus;
}


// print formatted data suitable for plotting
// ouput = the target stream
void StructureFactorOnSphere::WriteDataFile(std::ostream &output)
{
  if (output.flags() & ios_base::binary)
    {
      // write as binary file
      this->WriteBinaryData(output);
    }
  else
    {
      // write as textfile
      output << "# l\tS_l\terr(S_l)\n";
      for (int l=0;l<this->NbrStructureFactors; ++l)
	output << l <<"\t" <<  this->StructureFactors->Average(l) << this->StructureFactors->ErrorEstimate(l)<<endl;
    }
}


// write binary data 
// ouput = the target stream
void StructureFactorOnSphere::WriteBinaryData(std::ostream &output)
{
  std::cout << "Need to implement binary write"<<endl;
}

// set particle collection that the observable operates on
// system = particle collection
void StructureFactorOnSphere::SetParticleCollection(AbstractParticleCollection *system)
{
  this->System = (ParticleOnSphereCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  this->System->GetSpinorCoordinates(SpinorUCoordinates, SpinorVCoordinates);
}

// accessor function to return the legend and numerical values for legend
void StructureFactorOnSphere::GetVectorLegend(std::string &legendParameters, std::string &legendValue, RealVector &parameterValues)
{
  legendParameters = std::string("index 'l'");
  legendValue = std::string("Structure Factor 'S(l)'");
  parameterValues.Resize(this->NbrStructureFactors);
  for (int l=0;l<this->NbrStructureFactors; ++l)
    parameterValues[l] = l;
}

// accessor function for average and error for variables with real measurements
void StructureFactorOnSphere::GetRealVectorMeasurement(RealVector &values, RealVector &errors)
{
  values.Resize(this->NbrStructureFactors);
  errors.Resize(this->NbrStructureFactors);

  for (int l=0;l<this->NbrStructureFactors; ++l)
    {
      values[l] = this->StructureFactors->Average(l);
      errors[l] = this->StructureFactors->ErrorEstimate(l);
    }   
}
