/***************************************************************************
 *            WeightedRealObservable.cc
 *
 *  Sat Dec  17 12:35:16 2005
 *  Copyright  2005  Gunnar Möller
 *  moller@lptms.u-psud.fr
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "WeightedRealObservable.h"
#include <cmath>

#include<iostream>
using std::cout;
using std::endl;

extern double dsqrarg;   // shift declaration to nrutil.cc to suppress warnings if unused in module that includes nrutil.h
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)


#define LARGE_BIN 256.0

WeightedRealObservable::WeightedRealObservable(int maxBins , double minObsPerBin)
{
  StartBins = maxBins;      
  StartObsPerBin = minObsPerBin; 
  BinSize = StartObsPerBin;   
  MaxBinNum = maxBins;   
  BinEntries=0;    
  BinWeight=0.0;
  PreviousWeight=0.0; 
  PreviousEntries=0;

  Values.clear();
  Values2.clear();

  Observations.clear();
  
  Weights.clear();

  Name = string("noname");
  HistoryAverages.clear();
  HistoryBinVariances.clear();
}

void WeightedRealObservable::operator<<(const double& x)
{
  this->Observe(x,1.0);
}

void WeightedRealObservable::Observe(const double& x, double w) // makes observation with different weight
{
  if (Values.empty())
    { 
    // start first bin
    Values.push_back(w*x);
    Values2.push_back(w*x*x);
    BinEntries = 1;
    BinWeight = w;
    Weights.push_back(w);
    Observations.push_back(1);
  }
  else
    {
      double AverageW= (PreviousWeight+BinWeight)/(PreviousEntries+BinEntries);
      if ( (BinWeight+0.5*w<BinSize*AverageW) || (BinEntries==0) )  // have space in current bin ?
	{
	  // continue filling it up!
	  int tmp;
	  Values[(tmp=Values.size()-1)] += w*x;
	  Values2[tmp] += w*x*x;
	  Weights[tmp] +=w;
	  ++BinEntries;
	  ++(Observations[tmp]);
	  BinWeight+=w;
	}
      else
	{
	  // new bin required!
	  if(Values.size()<MaxBinNum)
	    {
	      // start a new bin
	      PreviousWeight+=BinWeight;
	      PreviousEntries+=BinEntries;
	      Values.push_back(w*x);
	      Values2.push_back(w*x*x);
	      BinEntries = 1;
	      BinWeight = w;
	      Weights.push_back(w);
	      Observations.push_back(1);
	    }
	  else
	    {
	      // store average and bin variance, halve the bins and add
	      HistoryAverages.push_back(this->Average());
	      HistoryBinVariances.push_back(this->VarianceOfBins());
	      Rebin(2);
	      this->Observe(x,w); // and just call again
	      return;
	    }
	}
    }
}


double WeightedRealObservable::AverageBinSize()
{
  if (Values.size()==1)
    return ((double) BinEntries);
  else
    return PreviousEntries / (Values.size()-1);
}

double WeightedRealObservable::AverageBinWeight()
{
  if (Values.size()==1)
    return ((double) BinWeight);
  else
    return PreviousWeight / (Values.size()-1);
}


int WeightedRealObservable::BinNumber() const
{
  return Values.size();
}


int WeightedRealObservable::FilledBinNumber() const
{
  return Values.size()-1;
}


unsigned WeightedRealObservable::Measurements()
{
  return (PreviousEntries+BinEntries);
}

double WeightedRealObservable::TotalWeight()
{
  return (PreviousWeight+BinWeight);
}

double WeightedRealObservable::AverageWeight()
{
  return (PreviousWeight+BinWeight)/(PreviousEntries+BinEntries);
}


void WeightedRealObservable::Rescale(double factor)
{
  for (unsigned int i=0;i<Values.size();i++)
    {
      Values[i]*=factor;
      Values2[i]*=factor*factor;
    }
}


double WeightedRealObservable::Average()  // returns the average of the measurements
{
  if (Values.empty())
    return 0.0;
  double sum=Values[0];
  int binNum=Values.size();
  for (int i=1;i<binNum;i++)
    sum+=Values[i];
  return (sum/(PreviousWeight+BinWeight));
}

double WeightedRealObservable::PresentBinAverage()
{
  if (BinEntries!=0)
    return( 1.0/BinWeight * Values[Values.size()-1] );
  else
    return(  Values[Values.size()-2]/Weights[Values.size()-2] );
}

double WeightedRealObservable::PresentBinAverage2()
{
  if (BinEntries!=0)
    return( 1.0/BinWeight * Values2[Values.size()-1] );
  else
    return(  Values[Values2.size()-2]/Weights[Values.size()-2] );
}


void WeightedRealObservable::Rebin(unsigned int sequence)
{
  if ((Values.empty()) || (sequence<=1))
    return;

  double wbar = this->AverageWeight();

  int presentNewBin = 0;
  double storedInPresent = 0.0;
    // increase BinSize by argument of call or to average weight times that value...
  if (wbar>BinSize) BinSize = wbar*sequence;
  else BinSize *= sequence;

  // collect previous measurements into larger bins:
  for (unsigned int i=0; i<Values.size(); ++i)
  {
    if ((storedInPresent!=0.0)&&(storedInPresent + 0.5*Weights[i] > BinSize*wbar))
      {  // have no space -> start new bin!
	storedInPresent=0.0;
	++presentNewBin;
      }
	
    if (storedInPresent==0.0) 
      {       // starting new bin:
	Values[presentNewBin]=Values[i];
	Values2[presentNewBin]=Values2[i];
	Observations[presentNewBin]=Observations[i];
	Weights[presentNewBin]=Weights[i];
	storedInPresent=Weights[i];
      }
    else
      {
	// add to current bin:
	Values[presentNewBin]+=Values[i];
	Values2[presentNewBin]+=Values2[i];
	Observations[presentNewBin]+=Observations[i];
	Weights[presentNewBin]+=Weights[i];
	storedInPresent+=Weights[i];
      }

  } // end collecting

  PreviousWeight = PreviousWeight + BinWeight - Weights[presentNewBin];
  PreviousEntries = PreviousEntries + BinEntries - Observations[presentNewBin];
  BinEntries=Observations[presentNewBin];
  BinWeight=Weights[presentNewBin];

  Values.resize(presentNewBin+1);
  Values2.resize(presentNewBin+1);
  Observations.resize(presentNewBin+1);
  Weights.resize(presentNewBin+1);
}

double WeightedRealObservable::Variance(double typicalWgt)
{
  if ((Values.empty())||(Values.size()<2))
    return 0.0;
  if (typicalWgt==0.0) typicalWgt=AverageWeight();
  double sum2=0.0;
  double TotalWeight = this->TotalWeight();
  int binNum=Values.size();
  double Average = this->Average();
  for (int i=0;i<binNum;++i)
    sum2+= Values2[i];
  return ( (sum2 - Average*Average*TotalWeight) /(TotalWeight-typicalWgt));
}



double WeightedRealObservable::VarianceOfBins()  // treat as if all bins had the same size!
{
  if ((Values.empty())||(Values.size()<2))
    return 0.0;
  double sum2=0.0;
  int binNum=Values.size();
  double Average=this->Average();
  for (int i=0;i<binNum-1;i++)
    sum2+=DSQR(fabs(Values[i]/Weights[i]-Average));
  
  if (BinWeight / this->AverageWeight() > 0.5*this->AverageBinWeight() ) // take last bin into account, if big enough:
    return ( (sum2+ (Weights[binNum-1]/this->AverageBinWeight())*DSQR(Values[binNum-1]/Weights[binNum-1]-Average)) / (binNum-1.0) );
  else return (sum2/(binNum-2.0)); // else ignore it
}


double WeightedRealObservable::ErrorEstimate()
{
  return (std::sqrt(this->VarianceOfBins()/this->BinNumber()));
}


void WeightedRealObservable::SetName(string &newName)
{
  this->Name=newName;
}

void WeightedRealObservable::SetName(char *newName)
{
  this->Name = string(newName);
}

void WeightedRealObservable::PrintInternal()
{
  cout << "Average weight: "<<AverageWeight()<<endl;
  cout << "Total weight: "<<TotalWeight()<<endl;
  cout << "Size: "<<Values.size()<<endl;
  cout << "Average: "<<Average()<<endl;  
  cout << "Values2: [ " << Values2[0];
  for (unsigned int i=1; i<Values.size();++i)
    cout <<", "<<Values2[i];
  cout <<"]"<<endl;
}
  
ostream& operator << (ostream& str, WeightedRealObservable &O)
{
  str << O.Name << "= " << O.Average() << " +/- "  << O.ErrorEstimate() << " , has var= " << O.Variance() << std::endl;
  return str;
}
