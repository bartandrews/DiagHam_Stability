////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2008 Gunnar Moeller                    //
//                                                                            //
//                                                                            //
//         class for a binned correlation function on the disk geometry       //
//                                                                            //
//                        last modification : 29/07/2016                      //
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
#include "SimpleTwoBodyCorrelatorOnDisk.h"
#include <iostream>

using std::ios;
using std::ios_base;
using std::endl;


SimpleTwoBodyCorrelatorOnDisk::SimpleTwoBodyCorrelatorOnDisk()
{
  this->Bins=0;
}

// constructor
// rMax = maximum radius to consider
// resolution = total number of bins
// highres = number of points in high resolution interval at small r
// range =  ranger over which high resolution is implemented
SimpleTwoBodyCorrelatorOnDisk::SimpleTwoBodyCorrelatorOnDisk(double rMax, int resolution, int highres, int range)
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
}
  
// destructor
SimpleTwoBodyCorrelatorOnDisk::~SimpleTwoBodyCorrelatorOnDisk()
{
  if (Bins!=0)
    delete [] Correlations;
}

// call to make an observation
// weight = relative weight of this sample
void SimpleTwoBodyCorrelatorOnDisk::RecordValue(double weight)
{
  double Rij,phi,tmp;
  int index;
  this->Measures+=weight;
  if (std::isnan(weight))
    {
      std::cout << "Error: irregular weight: "<<weight<<std::endl;
      exit(1);
    }
  
  for( int i=1;i<this->NbrParticles;++i)
    {
      //std::cout << CoordinatesZ[i] << " ";
      for(int j=0; j<i;++j)
	{
	  // set Rij=sin(\theta_ij/2)
	  Rij=Norm(CoordinatesZ[i]-CoordinatesZ[j]);
	  index=this->GetIndex(Rij);
	  if (index>=Range)
	    this->Correlations[index+Highres-Range]+=weight;  
	  else
	    {
	      index = this->GetHighResIndex(Rij);
	      this->Correlations[index]+=weight;
	    }
	}
    }
}

// print legend to the given stream
// all = flag indicating whether to print all, or shortened information
void SimpleTwoBodyCorrelatorOnDisk::PrintLegend(std::ostream &output, bool all)
{
  output << "# r\tg(r)"<<endl;
}

// print status to the given stream
// all = flag indicating whether to print all, or shortened information
void SimpleTwoBodyCorrelatorOnDisk::PrintStatus(std::ostream &output, bool all)
{
  // no action, for now
}

// request whether observable should be printed
//
bool SimpleTwoBodyCorrelatorOnDisk::IncludeInPrint()
{
  return this->PrintFlag;
}

// set print status
//
void SimpleTwoBodyCorrelatorOnDisk::IncludeInPrint(bool newStatus)
{
  this->PrintFlag=newStatus;
}


// print formatted data suitable for plotting
// ouput = the target stream
void SimpleTwoBodyCorrelatorOnDisk::WriteDataFile(std::ostream &output)
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
      double Normalization=Measures/Resolution*NbrParticles*NbrParticles;
      output << "# Rmax="<<this->MaxRadius<<", weight out of range = "<< this->Correlations[this->Bins]/Normalization<<"\n";
      output << "# r\tg(r)\n";
      Normalization=Measures/(Resolution*this->HighResRatio)*NbrParticles*NbrParticles;
      for (int i=0;i<Highres;i++)
	output << this->GetHighResBinRadius(i+0.5)<<"\t"
	       << this->Correlations[i]/Normalization << endl;
      Normalization=Measures/Resolution*NbrParticles*NbrParticles;
      for (int i=0;i<Resolution-Range; ++i)
	output << this->GetBinRadius(i+Range+0.5) <<"\t"
	       << this->Correlations[i+Highres]/Normalization << endl;
    }
}


// write binary data 
// ouput = the target stream
void SimpleTwoBodyCorrelatorOnDisk::WriteBinaryData(std::ostream &output)
{
  std::cout << "Need to implement binary write"<<endl;
}

// set particle collection that the observable operates on
// system = particle collection
void SimpleTwoBodyCorrelatorOnDisk::SetParticleCollection(AbstractParticleCollection *system)
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
void SimpleTwoBodyCorrelatorOnDisk::GetVectorLegend(std::string &legendParameters, std::string &legendValue, RealVector &parameterValues)
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
void SimpleTwoBodyCorrelatorOnDisk::GetRealVectorMeasurement(RealVector &values, RealVector &errors)
{
  values.Resize(this->Highres+this->Resolution-this->Range);
  errors.Resize(0);
  double Normalization=Measures/(Resolution*this->HighResRatio)*NbrParticles*NbrParticles;
  int pos=0;
  for (int i=0;i<this->Highres;i++, pos++)
    values[pos] = this->Correlations[i]/Normalization;
  Normalization=Measures/Resolution*NbrParticles*NbrParticles;
  for (int i=0;i<Resolution-Range; ++i, ++pos)
    values[pos] = this->Correlations[i+Highres]/Normalization;
}
