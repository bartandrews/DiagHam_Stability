#!/anaconda/bin/python
import sys
import math
import cmath
#sys.path.insert(1,'/anaconda/pkgs/numpy-1.7.1-py27_0/packages/lib/python2.7/site-packages/numpy')
import numpy as np
import scipy.linalg
#import matplotlib.pyplot as pp

class ComputeHofstadterBands:

    def __init__(self, numCellX, numCellY, unitCellX, unitCellY, numFlux, gauge):
        self.NumCellX = numCellX
        self.NumCellY = numCellY
        self.UnitCellX = unitCellX
        self.UnitCellY = unitCellY
        self.NumFlux = numFlux
        self.KXFactor = 2.0 * math.pi / (numCellX)
        self.KYFactor = 2.0 * math.pi / (numCellY)
        self.NumBand = unitCellX * unitCellY
        self.NumState = numCellX * numCellY
        self.FluxDensity = float(numFlux) / self.NumBand
        self.T1 = 1.0
        self.T2 = 0.0

        self.OneBodyBasis = [None]*(numCellX*numCellY)
        self.HaveBasis = False
        self.MixingAngle = 0
        self.MixingPhase = 0
        self.EnergyBandStructure = [[0 for i in range(self.NumBand)] for j in range(self.NumState)]
        #gauge choice. set to 0 for Landau gauge along y-axis. Set to 1 for symmetric gauge.
        self.Gauge = gauge
        if(gauge == 0):
            self.XTranslationPhase = cmath.rect(1.0, -2.0*math.pi*self.FluxDensity*unitCellX)
            self.YTranslationPhase = 1.0
        elif(gauge ==1):
            self.XTranslationPhase = cmath.rect(1.0,-math.pi*self.FluxDensity*unitCellX)
            self.YTranslationPhase = cmath.rect(1.0,math.pi*self.FluxDensity*unitCellY)
        else:
            print "Unrecognized gauge choice."
            sys.exit()


#returns a tuple (position index, translation phase)
    def GetIndexAndPhase(self, x, y, kx, ky):
        numXTrans = 0
        numYTrans = 0
        while(x < 0):
            x += self.UnitCellX
            numXTrans += 1
        while(x >= self.UnitCellX):
            x -= self.UnitCellX
            numXTrans -= 1  
        while(y < 0):
            y += self.UnitCellY
            numYTrans += 1 
        while(y >= self.UnitCellY):
            y -= self.UnitCellY
            numYTrans -= 1

        posIndex = x + self.UnitCellX*y
        tmpPhase1 = cmath.rect(1,0)
        tmpPhase2 = cmath.rect(0,0)
        phase = tmpPhase1

        if(numXTrans > 0):
            tmpPhase2 = self.XTranslationPhase
        else:
            tmpPhase2 = self.XTranslationPhase.conjugate()

        for i in range(0,int(math.fabs(numXTrans))):
            tmpPhase1 *= tmpPhase2

        for j in range(1,y+1):
            phase *= tmpPhase1

        phase *= cmath.rect(1,kx*numXTrans)
        phase = 1
        tmpPhase1 = cmath.rect(1,0)

        if(numYTrans > 0):
            tmpPhase2 = self.YTranslationPhase
        else:
            tmpPhase2 = self.YTranslationPhase.conjugate()

        for i in range(0,int(math.fabs(numYTrans))):
            tmpPhase1 *= tmpPhase2

        for j in range(1,x+1):
            phase *= tmpPhase1

        phase *= cmath.rect(1,ky*numYTrans)
        
        return (posIndex,phase)

    def MomentumIndex(self,kx,ky):
        return ky + kx*self.NumCellY

    def GetMomentumFromIndex(self, index):
        kx = int(index / self.NumCellY)
        ky = int(index % self.NumCellY)
        return (kx,ky)

    def GetPositionFromIndex(self,index):
        x = int(index / self.UnitCellX)
        y = int(index % self.UnitCellY)
        return (x,y)

    def ComputeOneBodyBasis(self):
        output_file = open("hofstadter_eigenvalues_N_" + str(self.NumBand) + ".dat",'w')
        output_file.write("index kx ky E\n")
        plot_indices =[]
        plot_evals = []
        a=1
        for kx in range(0,self.NumCellX):
            for ky in range(0,self.NumCellY):

                index = self.MomentumIndex(kx,ky)
                if(index >= 0 and index < self.NumState):

                    K1 = self.KXFactor * kx
                    K2 = self.KYFactor * ky
                    tmpHamiltonian = np.matrix(np.zeros(shape=(self.NumBand,self.NumBand)),dtype=complex)

                    # symmetric gauge
                    if(self.Gauge == 1):
                        for i in range(0,self.UnitCellX):
                            yPhaseFactor = cmath.rect(1.0,1.0*math.pi*self.FluxDensity*i)
                            for j in range(0,self.UnitCellY):
                                xPhaseFactor = cmath.rect(1.0,-1.0*math.pi*self.FluxDensity*j)
                                (initIndex,phase) = self.GetIndexAndPhase(i,j,K1,K2)

                                (finalIndex,phase) = self.GetIndexAndPhase(i+1,j,K1,K2)
                                if(initIndex >= finalIndex):
                                   tmpHamiltonian[initIndex,finalIndex] += phase* self.T1 * xPhaseFactor.conjugate()
                                   if(initIndex != finalIndex):
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate()* self.T1 * xPhaseFactor

                                (finalIndex,phase) = self.GetIndexAndPhase(i-1,j,K1,K2)

                                if(initIndex >= finalIndex):
                                   tmpHamiltonian[initIndex,finalIndex] += phase* self.T1 * xPhaseFactor
                                   if(initIndex != finalIndex):
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate()* self.T1 * xPhaseFactor.conjugate()

                                (finalIndex,phase) = self.GetIndexAndPhase(i,j+1,K1,K2)

                                if(initIndex >= finalIndex):
                                    tmpHamiltonian[initIndex,finalIndex] += phase * yPhaseFactor.conjugate()* self.T1
                                    if(initIndex != finalIndex):
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate() * yPhaseFactor* self.T1

                                (finalIndex,phase) = self.GetIndexAndPhase(i,j-1,K1,K2)

                                if(initIndex >= finalIndex):
                                    tmpHamiltonian[initIndex,finalIndex] += phase * yPhaseFactor* self.T1
                                    if(initIndex != finalIndex):
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate() * yPhaseFactor.conjugate() * self.T1
                            
                                  
                    # Landau gauge in y-axis
                    elif(self.Gauge ==0):
                        for i in range(0,self.UnitCellX):
                            yPhaseFactor = cmath.rect(1.0,2.0*math.pi*self.FluxDensity*(i+1))
                            offdiagonalphase = cmath.rect(1.0,K1)
                            #print "i=" + str(i)
                            #print "yphasefactor=" + str(yPhaseFactor)
                            for j in range(0,self.UnitCellY):
                                #print "j=" + str(j)
                                #xPhaseFactor = cmath.rect(1.0,-1.0*math.pi*self.FluxDensity*j)
                                (initIndex,phase) = self.GetIndexAndPhase(i,j,K1,K2)

                                if(initIndex == 0 or initIndex == (self.NumBand - 1)):
                                    offdiagonalphase = cmath.rect(1.0,-K1)
                                (finalIndex,phase) = self.GetIndexAndPhase(i+1,j,K1,K2)

                                if(initIndex >= finalIndex):
                                   tmpHamiltonian[initIndex,finalIndex] += phase* self.T1
                                   if(initIndex != finalIndex):
                                        tmpHamiltonian[initIndex,finalIndex] *= offdiagonalphase 
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate()* self.T1*offdiagonalphase.conjugate()

                                #print "[" + str(initIndex) + "," + str(finalIndex) + "]: " + str(tmpHamiltonian[initIndex,finalIndex])
                                (finalIndex,phase) = self.GetIndexAndPhase(i-1,j,K1,K2)

                                if(initIndex >= finalIndex):
                                   tmpHamiltonian[initIndex,finalIndex] += phase* self.T1 
                                   if(initIndex != finalIndex):
                                        tmpHamiltonian[initIndex,finalIndex] *= offdiagonalphase 
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate()* self.T1*offdiagonalphase.conjugate() 
                                #print "[" + str(initIndex) + "," + str(finalIndex) + "]: " + str(tmpHamiltonian[initIndex,finalIndex])
                                (finalIndex,phase) = self.GetIndexAndPhase(i,j+1,K1,K2)

                                if(initIndex >= finalIndex):
                                    tmpHamiltonian[initIndex,finalIndex] += phase * yPhaseFactor.conjugate()* self.T1
                                    if(initIndex != finalIndex):
                                        tmpHamiltonian[initIndex,finalIndex] *= offdiagonalphase 
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate()* self.T1*offdiagonalphase.conjugate()*yPhaseFactor
                                #print "[" + str(initIndex) + "," + str(finalIndex) + "]: " + str(tmpHamiltonian[initIndex,finalIndex])
                                (finalIndex,phase) = self.GetIndexAndPhase(i,j-1,K1,K2)

                                if(initIndex >= finalIndex):
                                    tmpHamiltonian[initIndex,finalIndex] += phase * yPhaseFactor* self.T1
                                    if(initIndex != finalIndex):
                                        tmpHamiltonian[initIndex,finalIndex] *= offdiagonalphase 
                                        tmpHamiltonian[finalIndex,initIndex] += phase.conjugate()* self.T1*offdiagonalphase.conjugate()*yPhaseFactor.conjugate()
                                #print "[" + str(initIndex) + "," + str(finalIndex) + "]: " + str(tmpHamiltonian[initIndex,finalIndex])
                                #(finalIndex,phase) = self.GetIndexAndPhase(i+2,j,K1,K2)

                                # modified Hofstadter with NNN hopping for simulating k^4 Hamiltonian:

                                # if(initIndex >= finalIndex):
                                #    tmpHamiltonian[initIndex,finalIndex] += phase * self.T2
                                #    if(initIndex != finalIndex):
                                #         tmpHamiltonian[finalIndex,initIndex] += phase.conjugate() * self.T2

                                # (finalIndex,phase) = self.GetIndexAndPhase(i-2,j,K1,K2)

                                # if(initIndex >= finalIndex):
                                #    tmpHamiltonian[initIndex,finalIndex] += phase * self.T2
                                #    if(initIndex != finalIndex):
                                #         tmpHamiltonian[finalIndex,initIndex] += phase.conjugate() * self.T2

                                # (finalIndex,phase) = self.GetIndexAndPhase(i,j+2,K1,K2)

                                # if(initIndex >= finalIndex):
                                #     tmpHamiltonian[initIndex,finalIndex] += phase * phaseFactor.conjugate()  * phaseFactor.conjugate() * self.T2
                                #     if(initIndex != finalIndex):
                                #         tmpHamiltonian[finalIndex,initIndex] += phase.conjugate() * phaseFactor * phaseFactor * self.T2

                                # (finalIndex,phase) = self.GetIndexAndPhase(i,j-2,K1,K2)

                                # if(initIndex >= finalIndex):
                                #     tmpHamiltonian[initIndex,finalIndex] += phase * phaseFactor * phaseFactor * self.T2
                                #     if(initIndex != finalIndex):
                                #         tmpHamiltonian[finalIndex,initIndex] += phase.conjugate() * phaseFactor.conjugate()  * phaseFactor.conjugate() * self.T2
                    

                
                eigenvalues,eigenvectors = scipy.linalg.eigh(tmpHamiltonian)
                self.OneBodyBasis[index] = eigenvectors
                #tmpDiag = np.array([])
                #for j in range(0,self.NumBand):
                #    x,y = self.GetPositionFromIndex(j)
                #    tmpDiag = np.append(tmpDiag,[cmath.rect(1,K1*kx*x + K2*ky*x)])
                #tmpU = np.diag(tmpDiag)
                #tmpBlochHamiltonian = np.dot(np.dot(tmpU,tmpHamiltonian),np.conjugate(tmpU))
                #print tmpHamiltonian
                #print tmpBlochHamiltonian
                #blochEigenvalues,blochEigenvectors = scipy.linalg.eigh(tmpBlochHamiltonian)
                
                #print "Index " + str(index)
                #print "kx ky"
                #print str(kx) + "  " +str(ky)
                #print "---------------------"
                #self.PrettyPrintMatrix(tmpHamiltonian)
                #print eigenvalues
                #print eigenvectors
                #print "==============="
                #self.PrettyPrintMatrix(tmpBlochHamiltonian)
                #print blochEigenvalues
                #print blochEigenvectors
                # #print eigenvalues
                # #print eigenvectors
                # for e in np.nditer(eigenvalues):
                # #    plot_indices.append(index)
                # #    plot_evals.append(e)
                #     output_file.write(str(index) + " " + str(kx) + " " + str(ky) + " " + str(e) + '\n')
                # #
                #for i in range (0,self.NumBand):
                #    print str(i) + "-th eigenvector"
                #    print self.OneBodyBasis[index][i]  

                self.HaveBasis = True

        #output_file.close()
        #pp.scatter(plot_indices,plot_evals)
        #pp.show()
    def GetBasisVector(self,index,energy):
        #self.OneBodyBasis[index][0] = self.OneBodyBasis[index][1]
        return self.OneBodyBasis[index][energy]

    def GetEnergy(self,index,energy):
        return self.EnergyBandStructure[energy][index]

    def PrettyPrintMatrix(self,matrix):
        for row in matrix:
            for e in np.nditer(row):
                (r,theta) = cmath.polar(e.item())
                if(r == 0.0 or abs(r) < 10**(-13)):
                    print 0,
                elif theta == 0.0 or abs(theta) < 10**(-14):
                    if(r != int(r)):
                        print "{0:.2f}".format(r),
                    else:
                        print int(r),
                else: 
                    if(r != int(r)):
                        print "{0:.2f}".format(r) + "exp(" + "{0:.2f}".format(theta) + "i) ",
                    else:
                        if (r==1):
                            print "exp(" + "{0:.2f}".format(theta) + "i) ",
                        else:
                            print str(int(r)) + "exp(" + "{0:.2f}".format(theta) + "i) ",
            print
        

#c = ComputeHofstadterBands(int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3]), int(sys.argv[4]), float(sys.argv[5]),int(sys.argv[6]))
#c.ComputeOneBodyBasis()

