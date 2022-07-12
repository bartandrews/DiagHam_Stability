/***************************************************************************
 *            ComplexVectorObservable.cc
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

#include "ComplexVectorObservable.h"
#include <cmath>

#include<iostream>
using std::cout;
using std::endl;

#define LARGE_N 256  // defines the bin size from which a partly filled bin is used in Variance


ComplexVectorObservable::ComplexVectorObservable(int length, int maxBins , int minObsPerBin)
{
  NumFields = length;
  StartBins = maxBins;      
  StartObsPerBin = minObsPerBin;
  BinSize= minObsPerBin;
  MaxBinNum = maxBins;
  BinEntries = 0;    // number of measurements in last bin
  PreviousEntries=0;

  Values = new std::vector<WCVOIndividualType>[NumFields];
  Values2Re = new std::vector<double>[NumFields];
  Values2Im = new std::vector<double>[NumFields];
  HistoryAverages = new std::vector<WCVOIndividualType>[NumFields];
  HistoryBinVariances = new std::vector<double>[NumFields];
  for (int i=0; i<length; ++i)
    {
      Values[i].clear();
      Values2Re[i].clear();
      Values2Im[i].clear();
      HistoryAverages[i].clear();
      HistoryBinVariances[i].clear();
    }

  Name = string("noname");
  Names = NULL;
}

ComplexVectorObservable::~ComplexVectorObservable()
{
  delete[]Values;
  delete[]Values2Re;
  delete[]Values2Im;
  delete[]HistoryAverages;
  delete[]HistoryBinVariances;
  if (Names != 0) delete[]Names;
}  

void ComplexVectorObservable::Observe(const WCVOType allX) // makes observation with different weight
{
  if (Values[0].empty())
    {
      for (int i=0; i<NumFields; ++i)
	{ 
	  // start first bin
	  Values[i].push_back(allX[i]);
	  Values2Re[i].push_back(allX[i].Re*allX[i].Re);
	  Values2Im[i].push_back(allX[i].Im*allX[i].Im);
	}
      BinEntries = 1;
    }
  else if (BinEntries==BinSize) // have a full bin
    {
      if(Values[0].size()<MaxBinNum)
	{
	  for (int i=0; i<NumFields; ++i)
	    { 
	      // start a new bin
	      Values[i].push_back(allX[i]);
	      Values2Re[i].push_back(allX[i].Re*allX[i].Re);
	      Values2Im[i].push_back(allX[i].Im*allX[i].Im);
	    }
	  BinEntries = 1;
	}
      else
	{
	  // store Average and bin Variance, halve the bins and add
	  for (int i=0; i<NumFields; ++i)
	    { 
	      HistoryAverages[i].push_back(this->Average(i));
	      HistoryBinVariances[i].push_back(this->VarianceOfBins(i));
	    }
	  Rebin(2);
	  this->Observe(allX); // and just call again
	  return;
	}
    }
  else
    {
      int CurrentBin=Values[0].size()-1;
      for (int i=0; i<NumFields; ++i)
	{ 
	  Values[i][CurrentBin] += allX[i];
	  Values2Re[i][CurrentBin] += allX[i].Re*allX[i].Re;
	  Values2Im[i][CurrentBin] += allX[i].Im*allX[i].Im;
	}
      ++BinEntries;
    }
}


// observation with different weight on the same observable
void ComplexVectorObservable::Observe(const WCVOIndividualType x)
{
  if (Values[0].empty())
    {
      for (int i=0; i<NumFields; ++i)
	{ 
	  // start first bin
	  Values[i].push_back(x);
	  Values2Re[i].push_back(x.Re*x.Re);
	  Values2Im[i].push_back(x.Im*x.Im);
	}
      BinEntries = 1;
    }
  else if (BinEntries==BinSize) // have a full bin
  {
    if(Values[0].size()<MaxBinNum)
      {
	for (int i=0; i<NumFields; ++i)
	  { 
	    // start a new bin
	    Values[i].push_back(x);
	    Values2Re[i].push_back(x.Re*x.Re);
	    Values2Im[i].push_back(x.Im*x.Im);
	  }
	BinEntries = 1;
      }
    else
      {
	// store Average and bin Variance, halve the bins and add
	for (int i=0; i<NumFields; ++i)
	  { 
	    HistoryAverages[i].push_back(this->Average(i));
	    HistoryBinVariances[i].push_back(this->VarianceOfBins(i));
	  }
	Rebin(2);
	this->Observe(x); // and just call again
      return;
    }
  }
  else
  {
    int CurrentBin=Values[0].size()-1;
    for (int i=0; i<NumFields; ++i)
      { 
	Values[i][CurrentBin] += x;
	Values2Re[i][CurrentBin] += x.Re*x.Re;
	Values2Im[i][CurrentBin] += x.Im*x.Im;
      }
    ++BinEntries;
  }
}

int ComplexVectorObservable::BinNumber() const
{
  return Values[0].size();
}

int ComplexVectorObservable::FilledBinNumber() const
{
  return Values[0].size()-1;
}


unsigned long ComplexVectorObservable::Measurements()
{
  return (BinSize*(Values[0].size()-1)+BinEntries);
}


void ComplexVectorObservable::Rescale(double factor)
{
  for (int i=0; i<NumFields;i++)
    {
      for (unsigned int c=0;c<Values[i].size();c++)
	{
	  Values[i][c]*=factor;
	  Values2Re[i][c]*=factor*factor;
	  Values2Im[i][c]*=factor*factor;
	}
    }
}

void ComplexVectorObservable::Rescale(int field, double factor)
{
  for (unsigned int c=0;c<Values[field].size();c++)
    {
      Values[field][c]*=factor;
      Values2Re[field][c]*=factor*factor;
      Values2Im[field][c]*=factor*factor;
    }
}


WCVOType ComplexVectorObservable::Average()  // returns the average of the measurements
{
  WCVOType rst=new WCVOIndividualType[NumFields];
  WCVOIndividualType sum;
  int bin_num;
  double InvMeasurements = 1.0/(double)this->Measurements();
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
	  rst[i]=sum*InvMeasurements;
	}
    }
  return rst;
}


WCVOIndividualType ComplexVectorObservable::Average(int field)  // returns the average of the measurements
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
      return sum/(double)this->Measurements();
    }
}


WCVOIndividualType ComplexVectorObservable::PresentBinAverage(int i)
{
  if (BinEntries!=0)
    return( 1.0/(double)BinEntries * Values[i][Values[i].size()-1] );
  else
    if (Values[i].size()>1)
      return(  Values[i][Values[i].size()-2]/BinSize);
  else
    return Complex(0.0);
}

WCVOIndividualType ComplexVectorObservable::PresentBinAverage2(int i)
{
  if (BinEntries!=0)
    return( 1.0/(double)BinEntries * (Values2Re[i][Values[i].size()-1] + Values2Im[i][Values[i].size()-1]));
  else
    if (Values2Re[i].size()>1)
      return(  (Values2Re[i][Values[i].size()-2]+Values2Im[i][Values[i].size()-2])/BinSize);
  else
    return Complex(0.0);
}




void ComplexVectorObservable::Rebin(unsigned int sequence)
{
  if (Values[0].empty() || sequence<=1)
    return;
    
  int newbins = (Values[0].size()+sequence-1)/sequence;
  
  // full bins
  for (int field=0; field<NumFields; ++field)
    {
      for (unsigned int i=0;i<Values[field].size()/sequence;++i)
	{
	  Values[field][i]=Values[field][sequence*i];
	  Values2Re[field][i]=Values2Re[field][sequence*i];
	  Values2Im[field][i]=Values2Im[field][sequence*i];
	  for (unsigned int j = 1 ; j<sequence;++j)
	    {
	      Values[field][i] += Values[field][sequence*i+j];
	      Values2Re[field][i] += Values2Re[field][sequence*i+j];
	      Values2Im[field][i] += Values2Im[field][sequence*i+j];
	    }
	}
      
      // last part of partly full bins
      Values[field][newbins-1]=Values[field][sequence*(newbins-1)];
      Values2Re[field][newbins-1]=Values2Re[field][sequence*(newbins-1)];
      Values2Im[field][newbins-1]=Values2Im[field][sequence*(newbins-1)];
      for (unsigned int i=sequence*(newbins-1)+1;i<Values[field].size();++i){
	Values[field][newbins-1]+=Values[field][i];
	Values2Re[field][newbins-1]+=Values2Re[field][i];
	Values2Im[field][newbins-1]+=Values2Im[field][i];
      }
    }
  
  // how many in last bin?
  BinEntries+=((Values[0].size()-1)%sequence)*BinSize;
  BinSize*=sequence;
  
  for (int field=0; field<NumFields; ++field)
    {
      Values[field].resize(newbins);
      Values2Re[field].resize(newbins);
      Values2Im[field].resize(newbins);
    }
}

double* ComplexVectorObservable::Variance()
{
  double *rst=new double[NumFields];
  if ((Values[0].empty())||(Values[0].size()<2))
    return NULL; 
    //std::cerr<<"Variance of empty observable not defined!"<<std::endl;    
  int BinNum=Values[0].size();
  long samples=(BinNum-1)*BinSize+BinEntries;
  for (int f=0; f< NumFields;++f)
    {
      rst[f]=0.0;
      Complex Average=this->Average(f);
      for (int i=0;i<BinNum;i++)
	rst[f]+=Values2Re[f][i]+Values2Im[f][i];
      rst[f]= ((rst[f]-samples*SqrNorm(Average))/(samples-1));
    }
  return rst;
}

double ComplexVectorObservable::Variance(int field)
{
  double sum2;
  int bin_num;
  if ((Values[field].empty())||(Values[field].size()<2))
    return 0.0;
  sum2 =0.0;
  bin_num=Values[field].size();
  long samples=(bin_num-1)*BinSize+BinEntries;
  Complex average = this->Average(field);
  for (int j=0;j<bin_num;++j)
    sum2+= Values2Re[field][j] + Values2Im[field][j];
  return ( (sum2 - samples*SqrNorm(average))/(samples-1));
}
 
double* ComplexVectorObservable::VarianceOfBins()
{
  if ((Values[0].empty())||(Values[0].size()<2))
    return NULL;
    //cerr<<"Variance of empty observable not defined!\n";    
  double *rst=new double[NumFields];
  double sum2=0.0;
  int BinNum=Values[0].size();
  Complex Average;
  double invbin=(double)1.0/BinSize;

  for (int f=0; f< NumFields;++f)
    {
      Average=this->Average(f);
      for (int i=0;i<BinNum-1;i++)
	sum2+=SqrNorm(invbin*Values[f][i]-Average);
      if (BinEntries==BinSize)
	rst[f] = ((sum2+SqrNorm(invbin*Values[f][BinNum-1]-Average))/(BinNum-1.0));
      else if (BinEntries > LARGE_N) // take partly filled bin into account
	rst[f] = ((sum2+  (double)BinEntries/BinSize*
		 SqrNorm(1.0/BinEntries*Values[f][BinNum-1]-Average))
		/(BinNum-1.0));
      else rst[f] = (sum2/(BinNum-2.0)); // else ignore it
    }
  return rst;
}

double ComplexVectorObservable::VarianceOfBins(int f)
{
  if ((Values[0].empty())||(Values[0].size()<2))
    return 0.0;
    //cerr<<"Variance of empty observable not defined!\n";    
  double sum2=0.0;
  int BinNum=Values[0].size();
  double invbin=(double)1.0/BinSize;

  Complex Average=this->Average(f);
  for (int i=0;i<BinNum-1;i++)
    sum2+=SqrNorm(invbin*Values[f][i]-Average);
  if (BinEntries==BinSize)
    return ((sum2+SqrNorm(invbin*Values[f][BinNum-1]-Average))/(BinNum-1.0));
  else if (BinEntries > LARGE_N) // take partly filled bin into account
    return ((sum2+  (double)BinEntries/BinSize*
	       SqrNorm(1.0/BinEntries*Values[f][BinNum-1]-Average))
	      /(BinNum-1.0));
  else return (sum2/(BinNum-2.0)); // else ignore it
}



double ComplexVectorObservable::ErrorEstimate(int i)
{
  return  (sqrt(this->VarianceOfBins(i)/this->BinNumber()));
}

double* ComplexVectorObservable::ErrorEstimate()
{
  double *rst=this->VarianceOfBins();
  int BinNum=this->BinNumber();
  for (int i=0; i<NumFields;++i)
    rst[i]= sqrt(rst[i]/BinNum);
  return rst;
}


// variance and errors for real part:

double* ComplexVectorObservable::VarianceRe()
{
  double *rst=new double[NumFields];
  if ((Values[0].empty())||(Values[0].size()<2))
    return NULL; 
    //std::cerr<<"Variance of empty observable not defined!"<<std::endl;    
  int BinNum=Values[0].size();
  long samples=(BinNum-1)*BinSize+BinEntries;
  for (int f=0; f< NumFields;++f)
    {
      rst[f]=0.0;
      double Average=this->Average(f).Re;
      for (int i=0;i<BinNum;i++)
	rst[f]+=Values2Re[f][i];
      rst[f]= ((rst[f]-samples*Average*Average)/(samples-1));
    }
  return rst;
}

double ComplexVectorObservable::VarianceRe(int field)
{
  double sum2;
  int bin_num;
  if ((Values[field].empty())||(Values[field].size()<2))
    return 0.0;
  sum2 =0.0;
  bin_num=Values[field].size();
  long samples=(bin_num-1)*BinSize+BinEntries;
  double average = this->Average(field).Re;
  for (int j=0;j<bin_num;++j)
    sum2+= Values2Re[field][j];
  return ( (sum2 - samples*average*average)/(samples-1));
}
 
double* ComplexVectorObservable::VarianceOfBinsRe()
{
  if ((Values[0].empty())||(Values[0].size()<2))
    return NULL;
    //cerr<<"Variance of empty observable not defined!\n";    
  double *rst=new double[NumFields];
  double sum2=0.0;
  int BinNum=Values[0].size();
  double Average, tmp;
  double invbin=(double)1.0/BinSize;

  for (int f=0; f< NumFields;++f)
    {
      Average=this->Average(f).Re;
      for (int i=0;i<BinNum-1;i++)
	{
	  tmp = invbin*(Values[f][i].Re)-Average;
	  sum2+=tmp*tmp;
	}
      if (BinEntries==BinSize)
	{
	  tmp=invbin*Values[f][BinNum-1].Re-Average;
	  rst[f] = ((sum2+tmp*tmp)/(BinNum-1.0));
	}
      else if (BinEntries > LARGE_N) // take partly filled bin into account
	{
	  tmp = 1.0/BinEntries*Values[f][BinNum-1].Re-Average;
	  rst[f] = ((sum2+  (double)BinEntries/BinSize*tmp*tmp)
		    /(BinNum-1.0));
	}
      else rst[f] = (sum2/(BinNum-2.0)); // else ignore it
    }
  return rst;
}

double ComplexVectorObservable::VarianceOfBinsRe(int f)
{
  if ((Values[0].empty())||(Values[0].size()<2))
    return 0.0;
    //cerr<<"Variance of empty observable not defined!\n";    
  double sum2=0.0;
  int BinNum=Values[0].size();
  double tmp,invbin=(double)1.0/BinSize;

  double Average=this->Average(f).Re;
  for (int i=0;i<BinNum-1;i++)
    {
      tmp = invbin*(Values[f][i].Re)-Average;
      sum2+=tmp*tmp;
    }
  if (BinEntries==BinSize)
    {
      tmp=invbin*Values[f][BinNum-1].Re-Average;
      return ((sum2+tmp*tmp)/(BinNum-1.0));
    }
  else if (BinEntries > LARGE_N) // take partly filled bin into account
    {
      tmp = 1.0/BinEntries*Values[f][BinNum-1].Re-Average;
      return ((sum2+  (double)BinEntries/BinSize*tmp*tmp)
	      /(BinNum-1.0));
    }
  else return (sum2/(BinNum-2.0)); // else ignore it
}



double ComplexVectorObservable::ErrorEstimateRe(int i)
{
  return  (sqrt(this->VarianceOfBinsRe(i)/this->BinNumber()));
}

double* ComplexVectorObservable::ErrorEstimateRe()
{
  double *rst=this->VarianceOfBinsRe();
  int BinNum=this->BinNumber();
  for (int i=0; i<NumFields;++i)
    rst[i]= sqrt(rst[i]/BinNum);
  return rst;
}

// variance and error for imaginary part


double* ComplexVectorObservable::VarianceIm()
{
  double *rst=new double[NumFields];
  if ((Values[0].empty())||(Values[0].size()<2))
    return NULL; 
    //std::cerr<<"Variance of empty observable not defined!"<<std::endl;    
  int BinNum=Values[0].size();
  long samples=(BinNum-1)*BinSize+BinEntries;
  for (int f=0; f< NumFields;++f)
    {
      rst[f]=0.0;
      double Average=this->Average(f).Im;
      for (int i=0;i<BinNum;i++)
	rst[f]+=Values2Im[f][i];
      rst[f]= ((rst[f]-samples*Average*Average)/(samples-1));
    }
  return rst;
}

double ComplexVectorObservable::VarianceIm(int field)
{
  double sum2;
  int bin_num;
  if ((Values[field].empty())||(Values[field].size()<2))
    return 0.0;
  sum2 =0.0;
  bin_num=Values[field].size();
  long samples=(bin_num-1)*BinSize+BinEntries;
  double average = this->Average(field).Im;
  for (int j=0;j<bin_num;++j)
    sum2+= Values2Im[field][j];
  return ( (sum2 - samples*average*average)/(samples-1));
}
 
double* ComplexVectorObservable::VarianceOfBinsIm()
{
  if ((Values[0].empty())||(Values[0].size()<2))
    return NULL;
    //cerr<<"Variance of empty observable not defined!\n";    
  double *rst=new double[NumFields];
  double sum2=0.0;
  int BinNum=Values[0].size();
  double Average, tmp;
  double invbin=(double)1.0/BinSize;

  for (int f=0; f< NumFields;++f)
    {
      Average=this->Average(f).Im;
      for (int i=0;i<BinNum-1;i++)
	{
	  tmp = invbin*(Values[f][i].Im)-Average;
	  sum2+=tmp*tmp;
	}
      if (BinEntries==BinSize)
	{
	  tmp=invbin*Values[f][BinNum-1].Im-Average;
	  rst[f] = ((sum2+tmp*tmp)/(BinNum-1.0));
	}
      else if (BinEntries > LARGE_N) // take partly filled bin into account
	{
	  tmp = 1.0/BinEntries*Values[f][BinNum-1].Im-Average;
	  rst[f] = ((sum2+  (double)BinEntries/BinSize*tmp*tmp)
		    /(BinNum-1.0));
	}
      else rst[f] = (sum2/(BinNum-2.0)); // else ignore it
    }
  return rst;
}

double ComplexVectorObservable::VarianceOfBinsIm(int f)
{
  if ((Values[0].empty())||(Values[0].size()<2))
    return 0.0;
    //cerr<<"Variance of empty observable not defined!\n";    
  double sum2=0.0;
  int BinNum=Values[0].size();
  double tmp,invbin=(double)1.0/BinSize;

  double Average=this->Average(f).Im;
  for (int i=0;i<BinNum-1;i++)
    {
      tmp = invbin*(Values[f][i].Im)-Average;
      sum2+=tmp*tmp;
    }
  if (BinEntries==BinSize)
    {
      tmp=invbin*Values[f][BinNum-1].Im-Average;
      return ((sum2+tmp*tmp)/(BinNum-1.0));
    }
  else if (BinEntries > LARGE_N) // take partly filled bin into account
    {
      tmp = 1.0/BinEntries*Values[f][BinNum-1].Im-Average;
      return ((sum2+  (double)BinEntries/BinSize*tmp*tmp)
	      /(BinNum-1.0));
    }
  else return (sum2/(BinNum-2.0)); // else ignore it
}



double ComplexVectorObservable::ErrorEstimateIm(int i)
{
  return  (sqrt(this->VarianceOfBinsIm(i)/this->BinNumber()));
}

double* ComplexVectorObservable::ErrorEstimateIm()
{
  double *rst=this->VarianceOfBinsIm();
  int BinNum=this->BinNumber();
  for (int i=0; i<NumFields;++i)
    rst[i]= sqrt(rst[i]/BinNum);
  return rst;
}


void ComplexVectorObservable::SetName(string &newName)
{
  this->Name=newName;
}

void ComplexVectorObservable::SetName(char *newName)
{
  this->Name = string(newName);
}

void ComplexVectorObservable::SetFieldName(int field, string &newName)
{
  if ((field<0)||(field >=NumFields)) {std::cout<<"Field index out of range in SetFieldName!" << endl; return;}
  if (Names==0) Names = new string[NumFields];
  this->Names[field]=newName;
}

void ComplexVectorObservable::SetFieldName(int field, char *newName)
{
  if ((field<0)||(field >=NumFields)) {std::cout<<"Field index out of range in SetFieldName!" << endl; return;}
  if (Names==0) Names = new string[NumFields];
  this->Names[field]=string(newName);
}



ostream& operator << (ostream& str, ComplexVectorObservable &O)
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
