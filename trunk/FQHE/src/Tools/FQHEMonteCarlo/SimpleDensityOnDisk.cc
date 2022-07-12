////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//          class for a binned density measurement on the disk geometry       //
//                                                                            //
//                        last modification : 19/07/2016                      //
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
#include "SimpleDensityOnDisk.h"
#include <iostream>

using std::ios;
using std::ios_base;
using std::endl;


SimpleDensityOnDisk::SimpleDensityOnDisk()
{
  this->Bins=0;
}

// constructor
// rMax = maximum radius to consider
// resolution = total number of bins
// highres = number of points in high resolution interval at small r
// range =  ranger over which high resolution is implemented
SimpleDensityOnDisk::SimpleDensityOnDisk(double rMax, int resolution, int highres, int range, double defectAngle)
{
  this->Type = RealObservableT | VectorValued;
  this->PrintFlag=true;
  this->Bins=resolution+highres-range+1;
  this->Resolution=resolution;
  this->Highres=highres;
  this->Range=range;
  if (this->Range > 0)
    this->HighResRatio = (double)this->Highres / (double) this->Range;
  else
    {
      this->Range = 0;
      this->HighResRatio = 1.0;
      this->Highres=0;
    }
  this->Measures = 0.0;
  this->Correlations=new double[Bins+1];
  for (int j=0;j<=Bins;++j)
    Correlations[j]=0.0;
  this->NbrParticles=0;
  this->MaxRadius = rMax; // the radius is also the inverse magnetic length
  this->DefectAngle = defectAngle;
  this->Gamma = 1.0-defectAngle;
  this->InvGamma = 1.0/this->Gamma;
}
  
// destructor
SimpleDensityOnDisk::~SimpleDensityOnDisk()
{
  if (Bins!=0)
    delete [] Correlations;
}

// call to make an observation
// weight = relative weight of this sample
void SimpleDensityOnDisk::RecordValue(double weight)
{
  double Ri,phi,tmp;
  int index;
  this->Measures+=weight;
  if (std::isnan(weight))
    {
      std::cout << "Error: irregular weight: "<<weight<<std::endl;
      exit(1);
    }
  
  for( int i=0;i<this->NbrParticles;++i)
    {
      //std::cout << CoordinatesZ[i] << " ";
      Ri=Norm(CoordinatesZ[i]);
      index=this->GetIndex(Ri);
      if (index>=Range)
	{
	  if (index>=this->Bins)
	    std::cout << "INdex = "<<index<<", Highres="<<Highres<<", Range="<<Range<<", Bins="<<this->Bins<<std::endl;
	  this->Correlations[index+Highres-Range]+=weight;  
	}
      else
	{
	  index = this->GetHighResIndex(Ri);
	  this->Correlations[index]+=weight;
	}
    }
}

// print legend to the given stream
// all = flag indicating whether to print all, or shortened information
void SimpleDensityOnDisk::PrintLegend(std::ostream &output, bool all)
{
  output << "# r\tg(r)"<<endl;
}

// print status to the given stream
// all = flag indicating whether to print all, or shortened information
void SimpleDensityOnDisk::PrintStatus(std::ostream &output, bool all)
{
  // no action, for now
}

// request whether observable should be printed
//
bool SimpleDensityOnDisk::IncludeInPrint()
{
  return this->PrintFlag;
}

// set print status
//
void SimpleDensityOnDisk::IncludeInPrint(bool newStatus)
{
  this->PrintFlag=newStatus;
}


// print formatted data suitable for plotting
// ouput = the target stream
void SimpleDensityOnDisk::WriteDataFile(std::ostream &output)
{
  if (output.flags() & ios_base::binary)
    {
      // write as binary file
      this->WriteBinaryData(output);
    }
  else
    {
      // write as textfile

      double Units;
      double Normalization=Measures/Resolution*NbrParticles;
      output << "# Rmax="<<this->MaxRadius<<", weight out of range = "<< this->Correlations[this->Bins]/Normalization<<"\n";
      output << "# r\tn(r)\n";
      Normalization=Measures/(Resolution*this->HighResRatio)*NbrParticles;
      for (int i=0;i<Highres;i++)
	output << this->GetHighResBinRadius(i+0.5)<<"\t"
	       << this->Correlations[i]/Normalization << endl;
      Normalization=Measures/Resolution*NbrParticles;
      for (int i=0;i<Resolution-Range; ++i)
	output << this->GetBinRadius(i+Range+0.5) <<"\t"
	       << this->Correlations[i+Highres]/Normalization << endl;
    }
}


// write binary data 
// ouput = the target stream
void SimpleDensityOnDisk::WriteBinaryData(std::ostream &output)
{
  std::cout << "Need to implement binary write"<<endl;
}

// set particle collection that the observable operates on
// system = particle collection
void SimpleDensityOnDisk::SetParticleCollection(AbstractParticleCollection *system)
{
  if (system->GetCollectionType() != AbstractParticleCollection::OnDiskCollection)
    {
      std::cerr << "Error: wrong type of particle collection!"<<endl;
      exit(1);
    }
  this->System = (ParticleOnDiskCollection*) system;
  this->NbrParticles = System->GetNbrParticles();
  this->System->GetCoordinates(CoordinatesZ);
  // std::cout << "Particle collection registered in TwoBody Correlations"<<endl;
}

// accessor function to return the legend and numerical values for legend
void SimpleDensityOnDisk::GetVectorLegend(std::string &legendParameters, std::string &legendValue, RealVector &parameterValues)
{
  legendParameters = std::string("Radius 'r'");
  legendValue = std::string("Density 'n'");
  parameterValues.Resize(this->Highres+this->Resolution-this->Range);
  int pos=0;
  for (int i=0;i<Highres;i++,pos++)
    parameterValues[pos] = this->GetHighResBinRadius(i+0.5);
  for (int i=0;i<Resolution-Range; ++i, ++pos)
    parameterValues[pos] = this->GetBinRadius(i+Range+0.5);
}

// accessor function for average and error for variables with real measurements
void SimpleDensityOnDisk::GetRealVectorMeasurement(RealVector &values, RealVector &errors)
{
  values.Resize(this->Highres+this->Resolution-this->Range);
  errors.Resize(0);
  double Normalization=Measures/Resolution*NbrParticles;
  Normalization=Measures/(Resolution*this->HighResRatio)*NbrParticles;
  int pos=0;
  for (int i=0;i<Highres;i++, pos++)
    values[pos] = this->Correlations[i]/Normalization;
  Normalization=Measures/Resolution*NbrParticles;
  for (int i=0;i<Resolution-Range; ++i, ++pos)
    values[pos] = this->Correlations[i+Highres]/Normalization;

}





