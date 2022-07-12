////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//                                                                            //
//                            DiagHam  version 0.01                           //
//                                                                            //
//                  Copyright (C) 2001-2007 Nicolas Regnault                  //
//                                                                            //
//                          class author: Ying-Hai Wu                         //
//                                                                            //
//                  class of Hofstadter model with any Chern number           //
//                            and two body interaction                        //
//                                                                            //
//                        last modification : 11/10/2013                      //
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

#include "Architecture/AbstractArchitecture.h"
#include "Architecture/ArchitectureOperation/QHEParticlePrecalculationOperation.h"

#include "Matrix/ComplexMatrix.h"
#include "Matrix/HermitianMatrix.h"
#include "Matrix/RealDiagonalMatrix.h"

#include "Hamiltonian/HofstadterWithAnyChernModelTwoBodyHamiltonian.h"

#include <iostream>
#include <sys/time.h>

using std::cout;
using std::endl;
using std::ostream;

// default constructor

HofstadterWithAnyChernModelTwoBodyHamiltonian::HofstadterWithAnyChernModelTwoBodyHamiltonian()
{

}

// constructor
//
// particles = Hilbert space associated to the system
// nbrParticles = number of particles
// nbrSiteX = number of sites in the x direction
// nbrSiteY = number of sites in the y direction
// uPotential = two body neareast neighbor strength
// vPotential = two body second neareast neighbor strength
// gammaX = boundary twist angle along x
// gammaY = boundary twist angle along y
// flatBandFlag = use flat band model
// architecture = architecture to use for precalculation
// memory = memory for fast multiplication (negative if there is no limit)

HofstadterWithAnyChernModelTwoBodyHamiltonian::HofstadterWithAnyChernModelTwoBodyHamiltonian(ParticleOnSphere* particles,int nbrParticles,int nbrSiteX,int nbrSiteY,double uPotential,double vPotential,double gammaX,double gammaY,int chern,int nbrBand,bool flatBandFlag,AbstractArchitecture* architecture,long memory)
{   
    this->Particles=particles;
    this->NbrParticles=nbrParticles;
    this->NbrSiteX=nbrSiteX;
    this->NbrSiteY=nbrSiteY;
    this->LzMax=nbrSiteX*nbrSiteY-1;
    this->KxFactor=2.0*M_PI/this->NbrSiteX;
    this->KyFactor=2.0*M_PI/this->NbrSiteY;
    this->UPotential=uPotential;
    this->VPotential=vPotential;
    this->GammaX=gammaX;
    this->GammaY=gammaY;
    this->Chern=chern;
    this->NbrBand=nbrBand;
    this->NbrCell=this->NbrBand/this->Chern;    
    
    this->Architecture=architecture;
    this->Memory=memory;
    this->HamiltonianShift=0.0;
    this->OneBodyInteractionFactors=0;
    
    this->FlatBand=flatBandFlag;
    this->FastMultiplicationFlag=false;
    
    long MinIndex,MaxIndex;
   
    this->Architecture->GetTypicalRange(MinIndex,MaxIndex);
    this->PrecalculationShift=(int)MinIndex;  
    this->EvaluateInteractionFactors();
    
    if(memory>0)
    {
        long TmpMemory=this->FastMultiplicationMemory(memory);
        
        if(TmpMemory<1024) cout<<"fast="<<TmpMemory<<"b"<<endl;
        else if(TmpMemory<(1<<20)) cout<<"fast="<<(TmpMemory>>10)<<"kb"<<endl;
	else if(TmpMemory<(1<<30)) cout<<"fast="<<(TmpMemory>>20)<<"Mb"<<endl;
	else
	{
	    cout<<"fast="<<(TmpMemory>>30)<<".";
	    TmpMemory-=((TmpMemory>>30)<<30);
	    TmpMemory*=100l;
	    TmpMemory>>=30;
	    
	    if(TmpMemory<10l) cout<<"0";
	    
	    cout<<TmpMemory<<"Gb"<<endl;
	}
        
        this->EnableFastMultiplication();
    }
}

// destructor

HofstadterWithAnyChernModelTwoBodyHamiltonian::~HofstadterWithAnyChernModelTwoBodyHamiltonian()
{

}

// evaluate all interaction factors

void HofstadterWithAnyChernModelTwoBodyHamiltonian::EvaluateInteractionFactors()
{
    long TotalNbrInteractionFactors=0;
    int TotNbrSite=this->NbrSiteX*this->NbrSiteY;
        
    if(this->FlatBand==false) this->OneBodyInteractionFactors=new double[TotNbrSite];
    
    ComplexMatrix* OneBodyBasis=new ComplexMatrix[TotNbrSite];
    this->ComputeOneBodyMatrix(OneBodyBasis);
 
    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {
        cout<<"====== only intra-color nearest neighbor and inter-color onsite terms are used for fermions ======"<<endl;       
    
        this->NbrSectorSums=this->NbrSiteX*this->NbrSiteY;
        this->NbrSectorIndicesPerSum=new int[this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i) this->NbrSectorIndicesPerSum[i]=0;      
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2) 
	{
	    int Index1=kx1*this->NbrSiteY+ky1;
            int Index2=kx2*this->NbrSiteY+ky2;
		
	    if(Index1<Index2) ++this->NbrSectorIndicesPerSum[((kx1+kx2)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY]; 
	}
         
        this->SectorIndicesPerSum=new int* [this->NbrSectorSums];
         
        for(int i=0;i<this->NbrSectorSums;++i)
	{
	    if(this->NbrSectorIndicesPerSum[i]>0)
	    {
	        this->SectorIndicesPerSum[i]=new int[2*this->NbrSectorIndicesPerSum[i]];
                this->NbrSectorIndicesPerSum[i]=0;
	    }
	}
         
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2) 
	{
	    int Index1=kx1*this->NbrSiteY+ky1;
	    int Index2=kx2*this->NbrSiteY+ky2;
            
            if(Index1<Index2)
            {
		int TmpSum=((kx1+kx2)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY;
		 
	    	this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum]<<1]=Index1;    
	    	this->SectorIndicesPerSum[TmpSum][1+(this->NbrSectorIndicesPerSum[TmpSum]<<1)]=Index2;
	    	
	    	++this->NbrSectorIndicesPerSum[TmpSum];    
            }
        }
        
        double FactorU=1.0/(this->NbrSiteX*this->NbrSiteY);          
        double FactorV=1.0*this->VPotential/(this->NbrSiteX*this->NbrSiteY);
        
        if(this->FlatBand==false) FactorU*=this->UPotential;
        
        this->InteractionFactors=new Complex* [this->NbrSectorSums];

        for(int i=0;i<this->NbrSectorSums;++i)
	{
	    this->InteractionFactors[i]=new Complex[this->NbrSectorIndicesPerSum[i]*this->NbrSectorIndicesPerSum[i]];
	    
	    int Index=0;
	    
	    for(int j1=0;j1<this->NbrSectorIndicesPerSum[i];++j1)
	    {
	        int Index1=this->SectorIndicesPerSum[i][j1<<1];
	        int Index2=this->SectorIndicesPerSum[i][(j1<<1)+1];
	        
	        for(int j2=0;j2<this->NbrSectorIndicesPerSum[i];++j2)
		{
		    int Index3=this->SectorIndicesPerSum[i][j2<<1];
   		    int Index4=this->SectorIndicesPerSum[i][(j2<<1)+1];
 		    		    
 		    this->InteractionFactors[i][Index]=0.0;

 		    this->InteractionFactors[i][Index]+=this->TwoBodyMatrixElement(OneBodyBasis,Index1,Index3,Index2,Index4);
 		    this->InteractionFactors[i][Index]-=this->TwoBodyMatrixElement(OneBodyBasis,Index2,Index3,Index1,Index4);
 		    this->InteractionFactors[i][Index]-=this->TwoBodyMatrixElement(OneBodyBasis,Index1,Index4,Index2,Index3);
 		    this->InteractionFactors[i][Index]+=this->TwoBodyMatrixElement(OneBodyBasis,Index2,Index4,Index1,Index3);

		    this->InteractionFactors[i][Index]*=-FactorU;
 		        
		    TotalNbrInteractionFactors++;
		    ++Index;
		}
	    }
	}
    }    
    else
    {
        cout<<"====== only intra-color onsite terms are used for bosons! ======"<<endl;
        
        this->NbrSectorSums=this->NbrSiteX*this->NbrSiteY;
        this->NbrSectorIndicesPerSum=new int[this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i) this->NbrSectorIndicesPerSum[i]=0;      
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2) 
	{
	    int Index1=kx1*this->NbrSiteY+ky1;
	    int Index2=kx2*this->NbrSiteY+ky2;
	    
	    if(Index1<=Index2) ++this->NbrSectorIndicesPerSum[((kx1+kx2)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY];
	
	}
        
        this->SectorIndicesPerSum=new int* [this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i)
	{
	    if(this->NbrSectorIndicesPerSum[i]>0)
	    {
	        this->SectorIndicesPerSum[i]=new int[2*this->NbrSectorIndicesPerSum[i]];      
	        this->NbrSectorIndicesPerSum[i]=0;
	    }
	}
        
        for(int kx1=0;kx1<this->NbrSiteX;++kx1) for(int kx2=0;kx2<this->NbrSiteX;++kx2)
	for(int ky1=0;ky1<this->NbrSiteY;++ky1) for(int ky2=0;ky2<this->NbrSiteY;++ky2) 
	{
	    int Index1=kx1*this->NbrSiteY+ky1;
            int Index2=kx2*this->NbrSiteY+ky2;
		
            if(Index1<=Index2)
            {
		int TmpSum=((kx1+kx2)%this->NbrSiteX)*this->NbrSiteY+(ky1+ky2)%this->NbrSiteY;
		
		this->SectorIndicesPerSum[TmpSum][this->NbrSectorIndicesPerSum[TmpSum]<<1]=Index1;		
		this->SectorIndicesPerSum[TmpSum][1+(this->NbrSectorIndicesPerSum[TmpSum]<<1)]=Index2;
		
	        ++this->NbrSectorIndicesPerSum[TmpSum];
            }
	}
        
        double FactorU=1.0/(this->NbrSiteX*this->NbrSiteY);
        double FactorV=1.0*this->VPotential/(this->NbrSiteX*this->NbrSiteY);
        
        if(this->FlatBand==false) FactorU*=this->UPotential;
                
        this->InteractionFactors=new Complex* [this->NbrSectorSums];
        
        for(int i=0;i<this->NbrSectorSums;++i)
	{
	    this->InteractionFactors[i]=new Complex[this->NbrSectorIndicesPerSum[i]*this->NbrSectorIndicesPerSum[i]];
	    
	    int Index=0;
	    
	    for(int j1=0;j1<this->NbrSectorIndicesPerSum[i];++j1)
	    {
	        int Index1=this->SectorIndicesPerSum[i][j1<<1];
	        int Index2=this->SectorIndicesPerSum[i][(j1<<1)+1];
	        
	        for(int j2=0;j2<this->NbrSectorIndicesPerSum[i];++j2)
		{
		    int Index3=this->SectorIndicesPerSum[i][j2<<1];
		    int Index4=this->SectorIndicesPerSum[i][(j2<<1)+1];

 		    this->InteractionFactors[i][Index]=0.0;
 		    
 		    this->InteractionFactors[i][Index]+=this->TwoBodyMatrixElement(OneBodyBasis,Index1,Index3,Index2,Index4);
 		    this->InteractionFactors[i][Index]+=this->TwoBodyMatrixElement(OneBodyBasis,Index2,Index3,Index1,Index4);
 		    this->InteractionFactors[i][Index]+=this->TwoBodyMatrixElement(OneBodyBasis,Index1,Index4,Index2,Index3);
 		    this->InteractionFactors[i][Index]+=this->TwoBodyMatrixElement(OneBodyBasis,Index2,Index4,Index1,Index3);
 
		    if(Index1==Index2) this->InteractionFactors[i][Index]*=0.5;
		    if(Index3==Index4) this->InteractionFactors[i][Index]*=0.5;
		  
		    this->InteractionFactors[i][Index]*=FactorU;

		    TotalNbrInteractionFactors++;
		    ++Index;
		}
	    }
	}
    }
    
    delete[] OneBodyBasis;
    
    cout<<"nbr two body interaction = "<<TotalNbrInteractionFactors<<endl;
    cout<<"===================================="<<endl;
    
}

// compute two body matrix
    
Complex HofstadterWithAnyChernModelTwoBodyHamiltonian::TwoBodyMatrixElement(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4)
{
    Complex Result;
    
    if(this->Chern==1) Result=TwoBodyMatrixElementCore1(oneBodyBasis,Index1,Index2,Index3,Index4);
    else if(this->Chern==2) Result=TwoBodyMatrixElementCore2(oneBodyBasis,Index1,Index2,Index3,Index4);
    else if(this->Chern==3) Result=TwoBodyMatrixElementCore3(oneBodyBasis,Index1,Index2,Index3,Index4);
    else if(this->Chern>=4) Result=TwoBodyMatrixElementCore4(oneBodyBasis,Index1,Index2,Index3,Index4);

    return Result;
} 
 
Complex HofstadterWithAnyChernModelTwoBodyHamiltonian::TwoBodyMatrixElementCore1(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4)
{
    Complex Result,Intra,Inter;
    double Xhase,Yhase;int bfirs,bsecd; 
    
    int ix1=Index1/this->NbrSiteY,iy1=Index1%this->NbrSiteY;
    int ix2=Index2/this->NbrSiteY,iy2=Index2%this->NbrSiteY;
    int ix3=Index3/this->NbrSiteY,iy3=Index3%this->NbrSiteY;
    int ix4=Index4/this->NbrSiteY,iy4=Index4%this->NbrSiteY;

    double kx1=ix1*this->KxFactor,ky1=iy1*this->KyFactor;
    double kx2=ix2*this->KxFactor,ky2=iy2*this->KyFactor;    
    double kx3=ix3*this->KxFactor,ky3=iy3*this->KyFactor;
    double kx4=ix4*this->KxFactor,ky4=iy4*this->KyFactor;
    
    Result=0.0;
    
    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {     
        for(int i=0;i<this->NbrCell;i++)
        {                  
            // intra-color nearest neighor term
                                   
            Intra=0.0;
            
            for(int j=0;j<this->Chern;j++) 
            {
                // +1+0  
                bfirs=i+j*this->NbrCell;bsecd=(i+1)%this->NbrCell+j*this->NbrCell;Xhase=((i+1)/this->NbrCell)*(kx3-kx4);Yhase=0;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
                                  
                // +0+1
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;Xhase=0;Yhase=ky3-ky4;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
            
            // no inter-color term
            
            Result+=Intra;
        }     
    }
    else
    {
        for(int i=0;i<this->NbrCell;i++)
        {
            // intra-color onsite term
         
            Intra=0.0;
            
            for(int j=0;j<this->Chern;j++)
            { 
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;                          
                Intra+=Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
                       
            // no inter-color term
            
            Result+=Intra;
        }
    }
        
    return Result;
} 
    
Complex HofstadterWithAnyChernModelTwoBodyHamiltonian::TwoBodyMatrixElementCore2(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4)
{   
    Complex Result,Intra,Inter;
    double Xhase,Yhase;int bfirs,bsecd; 
    
    int ix1=Index1/this->NbrSiteY,iy1=Index1%this->NbrSiteY;
    int ix2=Index2/this->NbrSiteY,iy2=Index2%this->NbrSiteY;
    int ix3=Index3/this->NbrSiteY,iy3=Index3%this->NbrSiteY;
    int ix4=Index4/this->NbrSiteY,iy4=Index4%this->NbrSiteY;

    double kx1=ix1*this->KxFactor,ky1=iy1*this->KyFactor;
    double kx2=ix2*this->KxFactor,ky2=iy2*this->KyFactor;    
    double kx3=ix3*this->KxFactor,ky3=iy3*this->KyFactor;
    double kx4=ix4*this->KxFactor,ky4=iy4*this->KyFactor;
    
    Result=0.0;
    
    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {     
        for(int i=0;i<this->NbrCell;i++)
        {               
            // intra-color nearest neighbor term
                        
            Intra=0.0;
            
            for(int j=0;j<this->Chern;j++) 
            {                
                // +1+0  
                bfirs=i+j*this->NbrCell;bsecd=(i+1)%this->NbrCell+j*this->NbrCell;Xhase=((i+1)/this->NbrCell)*(kx3-kx4);Yhase=0;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
                                  
                // +0+1
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;Xhase=0;Yhase=ky3-ky4;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
            
            // inter-color onsite term
                                                   
            Inter=0.0;
            
            for(int j=0;j<this->Chern;j++) for(int k=j+1;k<this->Chern;k++)
            {
                // +0+0                
                bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;Xhase=0;Yhase=0;
                Inter+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];        
            }

            // inter-color nearest neighbor term
                        
            for(int j=0;j<this->Chern;j++) for(int k=j+1;k<this->Chern;k++)
            {
                // +1+0
                bfirs=i+j*this->NbrCell;bsecd=(i+1)%this->NbrCell+k*this->NbrCell;Xhase=((i+1)/this->NbrCell)*(kx3-kx4);Yhase=0;
                // Inter+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
                
                // +0+1
                bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;Xhase=0;Yhase=ky3-ky4;
                // Inter+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];         
            }

            Result+=Intra+Inter;
        }     
    }
    else
    {
        for(int i=0;i<this->NbrCell;i++)
        {            
            // intra-color onsite term
         
            Intra=0.0;
            
            for(int j=0;j<this->Chern;j++)
            { 
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;                          
                Intra+=Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }   
                               
            // inter-color onsite term
            
            Inter=0.0;
                                    
            for(int j=0;j<this->Chern;j++) for(int k=j+1;k<this->Chern;k++)
            {
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;        
                Inter+=0.0*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
                                    
            Result+=Intra+Inter;
        }
    }
        
    return Result;
}
 
Complex HofstadterWithAnyChernModelTwoBodyHamiltonian::TwoBodyMatrixElementCore3(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4)
{
    Complex Result,Intra,Inter;
    double Xhase,Yhase;int bfirs,bsecd;
    
    int ix1=Index1/this->NbrSiteY,iy1=Index1%this->NbrSiteY;
    int ix2=Index2/this->NbrSiteY,iy2=Index2%this->NbrSiteY;
    int ix3=Index3/this->NbrSiteY,iy3=Index3%this->NbrSiteY;
    int ix4=Index4/this->NbrSiteY,iy4=Index4%this->NbrSiteY;

    double kx1=ix1*this->KxFactor,ky1=iy1*this->KyFactor;
    double kx2=ix2*this->KxFactor,ky2=iy2*this->KyFactor;    
    double kx3=ix3*this->KxFactor,ky3=iy3*this->KyFactor;
    double kx4=ix4*this->KxFactor,ky4=iy4*this->KyFactor;
        
    Result=0.0;

    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {      
        for(int i=0;i<this->NbrCell;i++)
        {              
            // intra-color nearest neighbor term
                     
            Intra=0.0;

            for(int j=0;j<this->Chern;j++) 
            {
                // +1+0  
                bfirs=i+j*this->NbrCell;bsecd=(i+1)%this->NbrCell+j*this->NbrCell;Xhase=((i+1)/this->NbrCell)*(kx3-kx4);Yhase=0;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
                                  
                // +0+1
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;Xhase=0;Yhase=ky3-ky4;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
              
            // inter-color onsite term

            Inter=0.0;
             
            for(int j=0;j<this->Chern;j++) for(int k=j+1;k<this->Chern;k++)
            {
                bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;Xhase=0;Yhase=0;
                Inter+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }

            Result+=Intra+Inter;
        }
    }   
    else
    {
        for(int i=0;i<this->NbrCell;i++)
        {   
            // intra-color onsite term

            Intra=0.0;
            
            for(int j=0;j<this->Chern;j++)
            { 
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;                          
                Intra+=Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }   
                          
            // inter-color onsite term

            Inter=0.0;
                        
            for(int j=0;j<this->Chern;j++) for(int k=j+1;k<this->Chern;k++)
            {
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;
                Inter+=0.0*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
                                    
            Result+=Intra+Inter;
        }
    }
    
    return Result;  
} 
 
Complex HofstadterWithAnyChernModelTwoBodyHamiltonian::TwoBodyMatrixElementCore4(ComplexMatrix* oneBodyBasis,int Index1,int Index2,int Index3,int Index4)
{
    Complex Result,Inter,Intra;
    double Xhase,Yhase;int bfirs,bsecd;
    
    int ix1=Index1/this->NbrSiteY,iy1=Index1%this->NbrSiteY;
    int ix2=Index2/this->NbrSiteY,iy2=Index2%this->NbrSiteY;
    int ix3=Index3/this->NbrSiteY,iy3=Index3%this->NbrSiteY;
    int ix4=Index4/this->NbrSiteY,iy4=Index4%this->NbrSiteY;

    double kx1=ix1*this->KxFactor,ky1=iy1*this->KyFactor;
    double kx2=ix2*this->KxFactor,ky2=iy2*this->KyFactor;    
    double kx3=ix3*this->KxFactor,ky3=iy3*this->KyFactor;
    double kx4=ix4*this->KxFactor,ky4=iy4*this->KyFactor;
        
    Result=0.0;

    if(this->Particles->GetParticleStatistic()==ParticleOnSphere::FermionicStatistic)
    {      
        for(int i=0;i<this->NbrCell;i++)
        {   
            // intra-color nearest neighbor
                                
            Intra=0.0;
            
            for(int j=0;j<this->Chern;j++) 
            {
                // +1+0  
                bfirs=i+j*this->NbrCell;bsecd=(i+1)%this->NbrCell+j*this->NbrCell;Xhase=((i+1)/this->NbrCell)*(kx3-kx4);Yhase=0;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
                                  
                // +0+1
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;Xhase=0;Yhase=ky3-ky4;
                Intra+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
            
            // inter-color onsite term

            Inter=0.0;
            
            for(int j=0;j<this->Chern;j++) for(int k=j+1;k<this->Chern;k++)
            {
                bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;Xhase=0;Yhase=0;
                Inter+=Phase(Xhase+Yhase)*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }

            Result+=Intra+Inter;
        }
    }   
    else
    {
        for(int i=0;i<this->NbrCell;i++)
        {         
            // intra-color onsite term
            
            Intra=0.0;
                     
            for(int j=0;j<this->Chern;j++)
            { 
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+j*this->NbrCell;                          
                Intra+=Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }   
                         
            // inter-color onsite term 

            Inter=0.0;
             
            for(int j=0;j<this->Chern;j++) for(int k=j+1;k<this->Chern;k++)
            {
                // +0+0
                bfirs=i+j*this->NbrCell;bsecd=i+k*this->NbrCell;        
                Inter+=0.0*Conj(oneBodyBasis[Index1][0][bfirs])*oneBodyBasis[Index2][0][bfirs]*Conj(oneBodyBasis[Index3][0][bsecd])*oneBodyBasis[Index4][0][bsecd];
            }
                                    
            Result+=Intra+Inter;
        }
    }
    
    return Result;  
} 

// compute the one body matrices and one body band stucture contribution
//
// oneBodyBasis = array of one body matrices

void HofstadterWithAnyChernModelTwoBodyHamiltonian::ComputeOneBodyMatrix(ComplexMatrix* oneBodyBasis)
{
    Complex **MatEle;
    double *phase,t10,t20,kx,ky;
    int Index,color1,color2;
    
    t10=1.0,t20=0.0;

    phase=new double[this->NbrBand];    
    for(int i=0;i<this->NbrBand;i++) phase[i]=2.0*M_PI*i/this->NbrBand;

    MatEle=new Complex *[this->NbrBand];
    for(int i=0;i<this->NbrBand;i++) MatEle[i]=new Complex[this->NbrBand];
    
    for(int ix=0;ix<this->NbrSiteX;ix++) for(int iy=0;iy<this->NbrSiteY;iy++) 
    {
        Index=ix*this->NbrSiteY+iy;
        
        kx=(ix+this->GammaX)*this->KxFactor;
        ky=(iy+this->GammaY)*this->KyFactor;
             
        for(int i=0;i<this->NbrBand;i++) for(int j=0;j<this->NbrBand;j++) MatEle[i][j]=0.0;

	for(int i=0;i<this->NbrBand;i++)
	{
	    MatEle[i][i]+=-2*t10*cos(ky+phase[i]);

            color1=i/this->NbrCell;color2=(i+1)/this->NbrCell;
            
            if(color1!=color2)
            {
	        MatEle[i][(i+1)%this->NbrBand]+=-t10*Phase(kx);
	        MatEle[(i+1)%this->NbrBand][i]+=-t10*Phase(-kx);
	    }    
	    else 
	    {
	        MatEle[i][(i+1)%this->NbrBand]+=-t10;	    
                MatEle[(i+1)%this->NbrBand][i]+=-t10;
	    }
	}
        
	HermitianMatrix TmpOneBodyHamiltonian(this->NbrBand,true);
        ComplexMatrix TmpMatrix(this->NbrBand,this->NbrBand,true);
        RealDiagonalMatrix TmpDiag;

        for(int i=0;i<this->NbrBand;i++) 
        {
            TmpMatrix[i][i]=1.0;
            for(int j=i;j<this->NbrBand;j++) TmpOneBodyHamiltonian.SetMatrixElement(i,j,MatEle[i][j]);
	}
	
#ifdef __LAPACK__
	TmpOneBodyHamiltonian.LapackDiagonalize(TmpDiag,TmpMatrix);
#else
	TmpOneBodyHamiltonian.Diagonalize(TmpDiag,TmpMatrix);
#endif  
	oneBodyBasis[Index]=TmpMatrix;

	if(this->FlatBand==false) this->OneBodyInteractionFactors[Index]=TmpDiag(0,0);
				
	cout<<TmpDiag(0,0)<<" "<<TmpDiag(1,1)<<" "<<TmpDiag(2,2)<<" "<<TmpDiag(3,3)<<endl;
    }
        
    for(int i=0;i<this->NbrBand;i++) delete[] MatEle[i];
    delete[] MatEle;

    delete[] phase;
}
