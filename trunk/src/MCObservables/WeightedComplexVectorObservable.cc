/***************************************************************************
 *            WeightedComplexVectorObservable.cc
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

#include "WeightedComplexVectorObservable.h"
#include <cmath>

#include<iostream>
using std::cout;
using std::endl;

#define LARGE_BIN 256.0

WeightedComplexVectorObservable::WeightedComplexVectorObservable(int length, int maxBins , WCVOIndividualWeightType MinObsPerBin)
{
  NumFields = length;
  StartBins = maxBins;      
  StartObsPerBin = MinObsPerBin;
  BinSize= new double[NumFields];       // number of measurements per bin
  for (int i=0; i<NumFields;++i) BinSize[i] = StartObsPerBin;
  MaxBinNum= maxBins;
  BinWeight= new double[NumFields];
  for (int i=0; i<NumFields;++i) BinWeight [i]=0.0;  
  BinEntries= new unsigned[NumFields];    // number of measurements in last bin
  for (int i=0; i<NumFields;++i) BinEntries [i]=0;
  PreviousWeight= new double[NumFields]; // weight of all observations to this point
  for (int i=0; i<NumFields;++i) PreviousWeight [i]=0.0;
  PreviousEntries= new unsigned[NumFields];
  for (int i=0; i<NumFields;++i) PreviousEntries [i]=0;

  Values = new std::vector<WCVOIndividualType>[NumFields];
  Values2 = new std::vector<double>[NumFields];
  Observations = new std::vector<int>[NumFields];
  Weights = new std::vector<WCVOIndividualWeightType>[NumFields];
  HistoryAverages = new std::vector<WCVOIndividualType>[NumFields];;
  HistoryBinVariances = new std::vector<double>[NumFields];
  for (int i=0; i<length; ++i)
    {
      Values[i].clear();
      Values2[i].clear();
      Observations[i].clear();
      Weights[i].clear();
      HistoryAverages[i].clear();
      HistoryBinVariances[i].clear();
    }

  Name = string("noname");
  Names = NULL;
}

WeightedComplexVectorObservable::~WeightedComplexVectorObservable()
{
  delete[]Values;
  delete[]Values2;
  delete[]Observations;
  delete[]Weights;
  delete[]HistoryAverages;
  delete[]HistoryBinVariances;
  delete []  BinSize;
  delete []  BinWeight;
  delete []  BinEntries;
  delete []  PreviousWeight;
  delete []  PreviousEntries;
  if (Names != 0) delete[]Names;
}

void WeightedComplexVectorObservable::operator<<(const WCVOType x)
{
  this->Observe(x,1.0);
}

// makes observation with default weight 1.0 on the same observable
void WeightedComplexVectorObservable::operator<<(const WCVOIndividualType x)
{
  this->Observe(x,NULL);
}
  

void WeightedComplexVectorObservable::Observe(const WCVOType allX, WCVOWeightType allW) // makes observation with different weight
{
  WCVOIndividualWeightType w;
  WCVOIndividualType x;
  for (int i=0; i<NumFields; ++i)
    {
    secondRun:
      if (allW!=NULL) w=allW[i]; else w=1.0;
      x=allX[i];
      if (Values[i].empty())
	{ 
	  // start first bin
	  Values[i].push_back(w*x);
	  Values2[i].push_back(w*SqrNorm(x));
	  BinEntries[i] = 1;
	  BinWeight[i] = w;
	  Weights[i].push_back(w);
	  Observations[i].push_back(1);
	}
      else
	{
	  double average_w= (PreviousWeight[i]+BinWeight[i])/(PreviousEntries[i]+BinEntries[i]);
	  if ( (BinWeight[i]+0.5*w<BinSize[i]*average_w) || (BinEntries[i]==0) )  // have space in current bin ?
	    {
	      // continue filling it up!
	      int tmp;
	      Values[i][(tmp=Values[i].size()-1)] += w*x;
	      Values2[i][tmp] += w*SqrNorm(x);
	      Weights[i][tmp] +=w;
	      ++BinEntries[i];
	      Observations[i][tmp]++;
	      BinWeight[i]+=w;
	    }
	  else
	    {
	      // new bin required!
	      if(Values[i].size()<MaxBinNum)
		{
		  // start a new bin
		  PreviousWeight[i]+=BinWeight[i];
		  PreviousEntries[i]+=BinEntries[i];
		  Values[i].push_back(w*x);
		  Values2[i].push_back(w*SqrNorm(x));
		  BinEntries[i] = 1;
		  BinWeight[i] = w;
		  Weights[i].push_back(w);
		  Observations[i].push_back(1);
		}
	      else
		{
		  // store average and bin variance, halve the bins and add
		  HistoryAverages[i].push_back(this->Average(i));
		  HistoryBinVariances[i].push_back(this->VarianceOfBins(i));
		  this->Rebin(i,2);
		  goto secondRun;
		  return;
		}
	    }
	}
    }
}

// observation with different weight on the same observable
void WeightedComplexVectorObservable::Observe(const WCVOIndividualType x, WCVOWeightType allW)
{
  WCVOIndividualWeightType w;
  for (int i=0; i<NumFields; ++i)
    {
    secondRun:
      if (allW!=NULL) w=allW[i]; else w=1.0;
      if (Values[i].empty())
	{ 
	  // start first bin
	  Values[i].push_back(w*x);
	  Values2[i].push_back(w*SqrNorm(x));
	  BinEntries[i] = 1;
	  BinWeight[i] = w;
	  Weights[i].push_back(w);
	  Observations[i].push_back(1);
	}
      else
	{
	  double average_w= (PreviousWeight[i]+BinWeight[i])/(PreviousEntries[i]+BinEntries[i]);
	  if ( (BinWeight[i]+0.5*w<BinSize[i]*average_w) || (BinEntries[i]==0) )  // have space in current bin ?
	    {
	      // continue filling it up!
	      int tmp;
	      Values[i][(tmp=Values[i].size()-1)] += w*x;
	      Values2[i][tmp] += w*SqrNorm(x);
	      Weights[i][tmp] +=w;
	      ++BinEntries[i];
	      Observations[i][tmp]++;
	      BinWeight[i]+=w;
	    }
	  else
	    {
	      // new bin required!
	      if(Values[i].size()<MaxBinNum)
		{
		  // start a new bin
		  PreviousWeight[i]+=BinWeight[i];
		  PreviousEntries[i]+=BinEntries[i];
		  Values[i].push_back(w*x);
		  Values2[i].push_back(w*SqrNorm(x));
		  BinEntries[i] = 1;
		  BinWeight[i] = w;
		  Weights[i].push_back(w);
		  Observations[i].push_back(1);
		}
	      else
		{
		  // store average and bin variance, halve the bins and add
		  HistoryAverages[i].push_back(this->Average(i));
		  HistoryBinVariances[i].push_back(this->VarianceOfBins(i));
		  this->Rebin(i,2);
		  goto secondRun;
		  return;
		}
	    }
	}
    }
}

// makes observation with same weight for different observables:
void WeightedComplexVectorObservable::Observe(const WCVOType allX, WCVOIndividualWeightType w) 
{
  WCVOIndividualType x;
  for (int i=0; i<NumFields; ++i)
    {
    secondRun:
      x=allX[i];
      if (Values[i].empty())
	{ 
	  // start first bin
	  Values[i].push_back(w*x);
	  Values2[i].push_back(w*SqrNorm(x));
	  BinEntries[i] = 1;
	  BinWeight[i] = w;
	  Weights[i].push_back(w);
	  Observations[i].push_back(1);
	}
      else
	{
	  double average_w= (PreviousWeight[i]+BinWeight[i])/(PreviousEntries[i]+BinEntries[i]);
	  if ( (BinWeight[i]+0.5*w<BinSize[i]*average_w) || (BinEntries[i]==0) )  // have space in current bin ?
	    {
	      // continue filling it up!
	      int tmp;
	      Values[i][(tmp=Values[i].size()-1)] += w*x;
	      Values2[i][tmp] += w*SqrNorm(x);
	      Weights[i][tmp] +=w;
	      ++BinEntries[i];
	      Observations[i][tmp]++;
	      BinWeight[i]+=w;
	    }
	  else
	    {
	      // new bin required!
	      if(Values[i].size()<MaxBinNum)
		{
		  // start a new bin
		  PreviousWeight[i]+=BinWeight[i];
		  PreviousEntries[i]+=BinEntries[i];
		  Values[i].push_back(w*x);
		  Values2[i].push_back(w*SqrNorm(x));
		  BinEntries[i] = 1;
		  BinWeight[i] = w;
		  Weights[i].push_back(w);
		  Observations[i].push_back(1);
		}
	      else
		{
		  // store average and bin variance, halve the bins and add
		  HistoryAverages[i].push_back(this->Average(i));
		  HistoryBinVariances[i].push_back(this->VarianceOfBins(i));
		  this->Rebin(i,2);
		  goto secondRun;
		  return;
		}
	    }
	}
    }
}



WCVOWeightType WeightedComplexVectorObservable::AverageBinSize()
{
  WCVOWeightType rst = new WCVOIndividualWeightType[NumFields];
  for (int i=0; i<NumFields; ++i)
  if (Values[i].size()==1)
    rst[i]= ((WCVOIndividualWeightType) BinEntries[i]);
  else
    rst[i]= PreviousEntries[i] / (Values[i].size()-1);
  return rst;
}

WCVOWeightType WeightedComplexVectorObservable::AverageBinWeight()
{
  WCVOWeightType rst = new WCVOIndividualWeightType[NumFields];
  for (int i=0; i<NumFields; ++i)
    if (Values[i].size()==1)
      rst[i]= ((WCVOIndividualWeightType) BinWeight[i]);
    else
      rst[i]= PreviousWeight[i] / (Values[i].size()-1);
  return rst;
}

WCVOIndividualWeightType WeightedComplexVectorObservable::AverageBinWeight(int field)
{
  if (Values[field].size()==1)
    return ((WCVOIndividualWeightType) BinWeight[field]);
  else
    return PreviousWeight[field] / (Values[field].size()-1);
}

int* WeightedComplexVectorObservable::BinNumber() const
{
  int *rst = new int[NumFields];
  for (int i=0; i<NumFields; ++i)
    rst[i] =  Values[i].size();
  return rst;
}

int WeightedComplexVectorObservable::SingleBinNumber(int i) const
{
  return  Values[i].size();
}

int* WeightedComplexVectorObservable::FilledBinNumber() const
{
  int* rst = new int[NumFields];
  for (int i=0; i<NumFields; ++i)
    rst[i] = Values[i].size()-1;
  return rst;
}


unsigned WeightedComplexVectorObservable::Measurements()
{
  return (PreviousEntries[0]+BinEntries[0]);
}

WCVOWeightType WeightedComplexVectorObservable::TotalWeight()
{
  WCVOWeightType rst = new WCVOIndividualWeightType[NumFields];
  for (int i=0; i<NumFields; ++i)
    rst[i] = (PreviousWeight[i]+BinWeight[i]);
  return rst;
}

WCVOIndividualWeightType WeightedComplexVectorObservable::TotalWeight(int field)
{
  return (PreviousWeight[field]+BinWeight[field]);
}


WCVOWeightType WeightedComplexVectorObservable::AverageWeight()
{
  WCVOWeightType rst = new WCVOIndividualWeightType[NumFields];
  for (int i=0; i<NumFields; ++i)
    rst[i] = (PreviousWeight[i]+BinWeight[i])/(PreviousEntries[i]+BinEntries[i]);
  return rst;
}

WCVOIndividualWeightType WeightedComplexVectorObservable::AverageWeight(int i)
{
  return (PreviousWeight[i]+BinWeight[i])/(PreviousEntries[i]+BinEntries[i]);
}


void WeightedComplexVectorObservable::Rescale(double factor)
{
  for (int i=0; i<NumFields;i++)
    {
      for (unsigned int c=0;c<Values[i].size();c++)
	{
	  Values[i][c]*=factor;
	  Values2[i][c]*=factor*factor;
	}
    }
}

void WeightedComplexVectorObservable::Rescale(int field, double factor)
{
  for (unsigned int c=0;c<Values[field].size();c++)
    {
      Values[field][c]*=factor;
      Values2[field][c]*=factor*factor;
    }
}


WCVOType WeightedComplexVectorObservable::Average()  // returns the average of the measurements
{
  WCVOType rst=new WCVOIndividualType[NumFields];
  WCVOIndividualType sum;
  int bin_num;
  for (int i=0; i<NumFields; ++i)
    {
      if (Values[i].empty())
	rst[i]= 0.0;
      else
	{
	  sum=Values[i][0];
	  bin_num=Values[i].size();
	  for (int j=1;j<bin_num;++j)
	    sum+=Values[i][j];
	  rst[i]=sum/(PreviousWeight[i]+BinWeight[i]);
	}
    }
  return rst;
}

WCVOIndividualType WeightedComplexVectorObservable::Average(int field)  // returns the average of the measurements
{
  WCVOIndividualType sum;
  int bin_num;
  if (Values[field].empty())
    return 0.0;
  else
    {
      sum=Values[field][0];
      bin_num=Values[field].size();
      for (int j=1;j<bin_num;++j)
	sum+=Values[field][j];
      return sum/(PreviousWeight[field]+BinWeight[field]);
    }
}


WCVOIndividualType WeightedComplexVectorObservable::PresentBinAverage(int i)
{
  if (BinEntries[i]!=0)
    return( 1.0/BinWeight[i] * Values[i][Values[i].size()-1] );
  else
    return(  Values[i][Values[i].size()-2]/Weights[i][Values[i].size()-2] );
}

WCVOIndividualType WeightedComplexVectorObservable::PresentBinAverage2(int i)
{
  if (BinEntries[i]!=0)
    return( 1.0/BinWeight[i] * Values2[i][Values[i].size()-1] );
  else
    return(  Values[i][Values2[i].size()-2]/Weights[i][Values[i].size()-2] );
}


void WeightedComplexVectorObservable::Rebin(int field, unsigned int sequence)
{
  if ((Values[field].empty()) || (sequence<=1))
    return;

  double wbar = this->AverageWeight(field);

  int present_new_bin = 0;
  WCVOIndividualWeightType stored_in_present = 0.0;
  // increase BinSize by argument of call:
  if (wbar>BinSize[field]) BinSize[field] = wbar*sequence;
  else BinSize[field] *= sequence;
  // collect previous measurements into larger bins:
  for (unsigned int i=0; i<Values[field].size(); ++i)
  {
    if ((stored_in_present!=0.0)&&(stored_in_present + 0.5*Weights[field][i] > BinSize[field]*wbar))
      {  // have no space -> start new bin!
	stored_in_present=0.0;
	++present_new_bin;
      }
	
    if (stored_in_present==0.0) 
      {       // starting new bin:
	Values[field][present_new_bin]=Values[field][i];
	Values2[field][present_new_bin]=Values2[field][i];
	Observations[field][present_new_bin]=Observations[field][i];
	Weights[field][present_new_bin]=Weights[field][i];
	stored_in_present=Weights[field][i];
      }
    else
      {
	// add to current bin:
	Values[field][present_new_bin]+=Values[field][i];
	Values2[field][present_new_bin]+=Values2[field][i];
	Observations[field][present_new_bin]+=Observations[field][i];
	Weights[field][present_new_bin]+=Weights[field][i];
	stored_in_present+=Weights[field][i];
      }

  } // end collecting

  PreviousWeight[field] = PreviousWeight[field] + BinWeight[field] - Weights[field][present_new_bin];
  PreviousEntries[field] = PreviousEntries[field] + BinEntries[field] - Observations[field][present_new_bin];
  BinEntries[field]=Observations[field][present_new_bin];
  BinWeight[field]=Weights[field][present_new_bin];

  Values[field].resize(present_new_bin+1);
  Values2[field].resize(present_new_bin+1);
  Observations[field].resize(present_new_bin+1);
  Weights[field].resize(present_new_bin+1);
}

double* WeightedComplexVectorObservable::Variance(WCVOWeightType typicalWgt)
{
  double * rst=new double[NumFields];
  double sum2;
  Complex average;
  double total_weight;
  int bin_num;
  for (int i=0; i< NumFields;++i)
    {
      if ((Values[i].empty())||(Values[i].size()<2))
	rst[i]= 0.0;
      if (typicalWgt[i]==0.0) typicalWgt[i]=AverageWeight(i);
      sum2 =0.0;
      total_weight = this->TotalWeight(i);
      bin_num=Values[i].size();
      average = this->Average(i);
      for (int j=0;j<bin_num;++j)
	sum2+= Values2[i][j];
      rst[i]=  ( (sum2 - SqrNorm(average)*total_weight) /(total_weight-typicalWgt[i]));
    }
  return rst;
}

double WeightedComplexVectorObservable::Variance(int field, WCVOIndividualWeightType typicalWgt)
{
  double sum2;
  Complex average;
  double total_weight;
  int bin_num;
  if ((Values[field].empty())||(Values[field].size()<2))
    return 0.0;
  if (typicalWgt==0.0) typicalWgt=AverageWeight(field);
  sum2 =0.0;
  total_weight = this->TotalWeight(field);
  bin_num=Values[field].size();
  average = this->Average(field);
  for (int j=0;j<bin_num;++j)
    sum2+= Values2[field][j];
  return ( (sum2 - SqrNorm(average)*total_weight) /(total_weight-typicalWgt));
}
 


double* WeightedComplexVectorObservable::VarianceOfBins()  // treat as if all bins had the same size!
{
  double *rst=new double[NumFields];
  double sum2;
  Complex average;
  int bin_num;
  for (int i=0; i< NumFields;++i)
    {
      if ((Values[i].empty())||(Values[i].size()<2))
	rst[i]=  0.0;
      sum2=0.0;
      bin_num=Values[i].size();
      average=this->Average(i);
      for (int j=0;j<bin_num-1;++j)
	sum2+=SqrNorm(Values[i][j]/Weights[i][j]-average);
      
      if (BinWeight[i] / this->AverageWeight(i) > 0.5*this->AverageBinWeight(i) ) // take last bin into account, if big enough:
	rst[i]=  ( (sum2+ (Weights[i][bin_num-1]/this->AverageBinWeight(i))*SqrNorm(Values[i][bin_num-1]/Weights[i][bin_num-1]-average)) / (bin_num-1.0) );
      else rst[i]=  (sum2/(bin_num-2.0)); // else ignore it
    }
  return rst;
}

double WeightedComplexVectorObservable::VarianceOfBins(int i)  // treat as if all bins had the same size!
{
  double sum2;
  Complex average;
  int bin_num;
  if ((Values[i].empty())||(Values[i].size()<2))
    return   0.0;
  sum2=0.0;
  bin_num=Values[i].size();
  average=this->Average(i);
  for (int j=0;j<bin_num-1;++j)
    sum2+=SqrNorm(Values[i][j]/Weights[i][j]-average);
  
  if (BinWeight[i] / this->AverageWeight(i) > 0.5*this->AverageBinWeight(i) ) // take last bin into account, if big enough:
    return ( (sum2+ (Weights[i][bin_num-1]/this->AverageBinWeight(i))*SqrNorm(Values[i][bin_num-1]/Weights[i][bin_num-1]-average)) / (bin_num-1.0) );
  else return  (sum2/(bin_num-2.0)); // else ignore it
}


double WeightedComplexVectorObservable::ErrorEstimate(int i)
{
  return  (sqrt(this->VarianceOfBins(i)/this->SingleBinNumber(i)));
}

double* WeightedComplexVectorObservable::ErrorEstimate()
{
  double *rst=this->VarianceOfBins();
  for (int i=0; i<NumFields;++i)
    rst[i]= sqrt(rst[i]/this->SingleBinNumber(i));
  return rst;
}


void WeightedComplexVectorObservable::SetName(string &newName)
{
  this->Name=newName;
}

void WeightedComplexVectorObservable::SetName(char *newName)
{
  this->Name = string(newName);
}

void WeightedComplexVectorObservable::SetFieldName(int field, string &newName)
{
  if ((field<0)||(field >=NumFields)) {std::cout<<"Field index out of range in SetFieldName!" << endl; return;}
  if (Names==0) Names = new string[NumFields];
  this->Names[field]=newName;
}

void WeightedComplexVectorObservable::SetFieldName(int field, char *newName)
{
  if ((field<0)||(field >=NumFields)) {std::cout<<"Field index out of range in SetFieldName!" << endl; return;}
  if (Names==0) Names = new string[NumFields];
  this->Names[field]=string(newName);
}



ostream& operator << (ostream& str, WeightedComplexVectorObservable &O)
{
  str << O.Name << endl;
  for (int i=0;  i<O.NumFields;++i)
    {
      if ((O.Names==NULL)||(O.Names[i].empty()))
	str << "X[" << i <<"]= ";
      else str << O.Names[i]<<"=";
      str << O.Average(i) << " +/- "  << O.ErrorEstimate(i) << " , has var= " << O.Variance(i) << std::endl;
    }
  return str;
}
