////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                      Copyright (C) 2007 Gunnar Möller                      //
//                                                                            //
//                                                                            //
//           class implementing a paired CF wave function on the sphere       //
//                                                                            //
//                        last modification : 18/05/2007                      //
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


#include "RealObservable.h"

#include <cmath>
#include <iostream>

RealObservable::RealObservable(int Maxbins, int minObsPerBin)
{
  BinEntries=0;
  MaxBinNum=Maxbins;
  BinSize=minObsPerBin;
  StartObsPerBin=BinSize;
  StartBins=MaxBinNum;
  Values.clear();
  Values2.clear();
  name = string("noname");
}


void RealObservable::SetMinMaxbin(int newMaxbins, int newMinObsPerBin)
{
  if (Values.empty())
    {
      BinSize= newMinObsPerBin;
      MaxBinNum= newMaxbins;
      StartObsPerBin = BinSize;
      StartBins = MaxBinNum;
    }
}


void RealObservable::operator<<(const double& x)
{
  if (Values.empty())
    { 
    // start first bin
    Values.push_back(x);
    Values2.push_back(fabs(x*x));
    BinEntries= 1;
  }
  else if (BinEntries==BinSize) // have a full bin
  {
    if(Values.size()<MaxBinNum)
    {
      // start a new bin
      Values.push_back(x);
      Values2.push_back(fabs(x*x));
      BinEntries = 1;
    }
    else
    {
      // store Average and bin Variance, halve the bins and add
      HistoryAverages.push_back(this->Average());
      HistoryBinVariances.push_back(this->VarianceOfBins());
      Rebin(2);
      *this << x; // and just call again
      return;
    }
  }
  else
  {
    Values[Values.size()-1] += x;
    Values2[Values.size()-1] += fabs(x*x);
    ++BinEntries;
  }
}


int RealObservable::Measurements()
{
  return (BinSize*(Values.size()-1)+BinEntries);
}


void RealObservable::Rescale(double factor)
{
  for (unsigned int i=0;i<Values.size();i++)
    {
      Values[i]*=factor;
      Values2[i]*=factor*factor;
    }
}


void RealObservable::Reset()
{
  BinSize=StartObsPerBin;
  MaxBinNum=StartBins;
}




double RealObservable::Average()
{
  if (Values.empty())
    {
      cout << "Attention, average on empty observable!";
      return 0.0;
    }
  double sum;
  int BinNum=Values.size();
  sum=Values[0];
  for (int i=1;i<BinNum;i++)
    sum+=Values[i];
  return (1.0/(BinSize*(BinNum-1)+BinEntries)*sum);
}


double RealObservable::PresentBinAverage()
{
  if (BinEntries!=0)
    return( 1.0/BinEntries * Values[Values.size()-1] );
  else
    return( 1.0/BinSize * Values[Values.size()-2] );
}


double RealObservable::PresentBinAverage2()
{
  if (BinEntries!=0)
    return( 1.0/BinEntries * Values2[Values.size()-1] );
  else
    return( 1.0/BinSize * Values2[Values.size()-2] );
}


double RealObservable::Variance()
{
  if ((Values.empty())||(Values.size()<2))
    return 1e20; 
    //std::cerr<<"Variance of empty observable not defined!"<<std::endl;    
  double sum2=0.0;
  int BinNum=Values.size();
  int samples=(BinNum-1)*BinSize+BinEntries;
  double Average=this->Average();
  for (int i=0;i<BinNum;i++)
    sum2+=Values2[i];
  return ((sum2-samples*fabs(Average*Average))/(samples-1));
}




double RealObservable::VarianceOfBins()
{
  if ((Values.empty())||(Values.size()<2))
    return 1e20;
    //cerr<<"Variance of empty observable not defined!\n";    
  double sum2=0.0;
  int BinNum=Values.size();
  double Average=this->Average();
  double invbin=(double)1.0/BinSize;
  for (int i=0;i<BinNum-1;i++)
    sum2+=fabs(invbin*Values[i]-Average)*fabs(invbin*Values[i]-Average);
  if (BinEntries==BinSize)
    return ((sum2+fabs(invbin*Values[BinNum-1]-Average)*fabs(invbin*Values[BinNum-1]-Average))/(BinNum-1.0));
  else if (BinEntries > LARGE_N) // take partly filled bin into account
    return ((sum2+  (double)BinEntries/BinSize*
	     fabs(1.0/BinEntries*Values[BinNum-1]-Average)*fabs(1.0/BinEntries*Values[BinNum-1]-Average) )/(BinNum-1.0));
  else return (sum2/(BinNum-2.0)); // else ignore it
}


double RealObservable::ErrorEstimate()
{
  return (std::sqrt(this->VarianceOfBins()/this->BinNumber()));
}


void RealObservable::Rebin(unsigned int sequence)
{
  if (Values.empty() || sequence<=1)
    return;
    
  int newbins = (Values.size()+sequence-1)/sequence;
  
  // full bins
  for (unsigned int i=0;i<Values.size()/sequence;++i)
  {
    Values[i]=Values[sequence*i];
    Values2[i]=Values2[sequence*i];
    for (unsigned int j = 1 ; j<sequence;++j)
    {
      Values[i] += Values[sequence*i+j];
      Values2[i] += Values2[sequence*i+j];
    }
  }
  
  // last part of partly full bins
  Values[newbins-1]=Values[sequence*(newbins-1)];
  Values2[newbins-1]=Values2[sequence*(newbins-1)];
  for (unsigned int i=sequence*(newbins-1)+1;i<Values.size();++i){
    Values[newbins-1]+=Values[i];
    Values2[newbins-1]+=Values2[i];
  }
    
  // how many in last bin?
  BinEntries+=((Values.size()-1)%sequence)*BinSize;
  BinSize*=sequence;

  Values.resize(newbins);
  Values2.resize(newbins);
}


void RealObservable::SetBinNumber(unsigned binnum)
{  
  if(Values.size() > binnum)
    this->Rebin((Values.size()-1)/binnum+1);
  MaxBinNum=binnum;
}


int RealObservable::BinNumber() const 
{ 
  return Values.size();
}


int RealObservable::FilledBinNumber() const 
{ 
  if(Values.size()==0) return 0;
  else return Values.size() + (BinEntries ==BinSize ? 0 : -1);
}


void RealObservable::SetName(string &newName)
{
  this->name=newName;
}


void RealObservable::SetName(char *newName)
{
  this->name = string(newName);
}


ostream& operator << (ostream& str, RealObservable &B)
{
  str << B.name << "= " << B.Average() << " +/- "  << B.ErrorEstimate() << " , has var= " << B.Variance() << endl;
  return str;
}
