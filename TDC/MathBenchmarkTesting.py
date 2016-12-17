from cmath import *
from random import *
from time import time
import FakeConfigurationFileMaker
import FakeLCAOFileMaker

LCAOmatrix=None
configmatrix=None

def main():
    seed(137)
    starttime=time()
    FakeConfigurationFileMaker.makefiles(nummolecul, numconfigu)
    FakeLCAOFileMaker.makefiles(nummolecul, numatomics)
    print("making files: "+str(time()-starttime));starttime=time()
    makeMatricies()
    print("making matricies: "+str(time()-starttime));starttime=time()
    for i in range(numtimes):
#        print("round "+str(i))
        coords=[randrange(-100,100), randrange(-100,100), randrange(-100,100)]
        sphercoords=[sqrt(coords[0]**2+coords[1]**2+coords[2]**2),
                     pi/2 if coords[0]==0 else atan(coords[1]/coords[0]),
                     acos(coords[2]/sqrt(coords[0]**2+coords[1]**2+coords[2]**2))]
        wav1=configurationwavefunc(0, sphercoords[0], sphercoords[1], sphercoords[2])       
    print(coords)
    print(sphercoords)
    endtime=time()
    print(endtime-starttime)

numtimes=1000
numatomics=14
nummolecul=2000
numconfigu=50

def makeMatricies():
    global LCAOmatrix, configmatrix
    LCAOfile=open("fakeLCAO.txt")
    configfile=open("fakeconfiguration.txt")
    LCAOmatrix=[[0]*numatomics for i in range(nummolecul)]
    configmatrix=[[0]*nummolecul for i in range(numconfigu)]
    for i in range(nummolecul):
        theline=LCAOfile.readline().split()
        LCAOmatrix[i]=[float(n) for n in theline]
    for i in range(numconfigu):
        theline=configfile.readline().split()
        configmatrix[i]=[float(n) for n in theline]
    
a0=1

def atomicwavefunc(n, l, ml, r, phi, theta):
    if (n,l,ml)==(1,0,0):   w=(1/sqrt(2*pi))*           (1/sqrt(2))*                        2/(a0**1.5)*exp(-r/(a0))
    elif (n,l,ml)==(2,0,0): w=(1/sqrt(2*pi))*           (1/sqrt(2))*                        (1/(2*a0)**1.5)*(2-r/a0)*exp(-r/(2*a0))
    elif (n,l,ml)==(2,1,0): w=(1/sqrt(2*pi))*           (sqrt(6)*cos(theta)/2)*             (3/(6*a0)**1.5)*(r/a0)*exp(-r/(2*a0))
    elif (n,l,ml)==(2,1,1): w=(exp(1j*phi)/sqrt(2*pi))* (sqrt(3)*sin(theta)/2)*             (3/(6*a0)**1.5)*(r/a0)*exp(-r/(2*a0))
    elif (n,l,ml)==(2,1,-1):w=(exp(-1j*phi)/sqrt(2*pi))*(sqrt(3)*sin(theta)/2)*             (3/(6*a0)**1.5)*(r/a0)*exp(-r/(2*a0))
    elif (n,l,ml)==(3,0,0): w=(1/sqrt(2*pi))*           (1/sqrt(2))*                        (2/(27*a0)**1.5)*(27-18*(r/a0)+2*(r/a0)**2)*exp(-r/(3*a0))
    elif (n,l,ml)==(3,1,0): w=(1/sqrt(2*pi))*           (sqrt(6)*cos(theta)/2)*             (8/(54*a0)**1.5)*(6-r/a0)*(r/a0)*exp(-r/(3*a0))
    elif (n,l,ml)==(3,1,1): w=(exp(1j*phi)/sqrt(2*pi))* (sqrt(3)*sin(theta)/2)*             (8/(54*a0)**1.5)*(6-r/a0)*(r/a0)*exp(-r/(3*a0))
    elif (n,l,ml)==(3,1,-1):w=(exp(-1j*phi)/sqrt(2*pi))*(sqrt(3)*sin(theta)/2)*             (8/(54*a0)**1.5)*(6-r/a0)*(r/a0)*exp(-r/(3*a0))
    elif (n,l,ml)==(3,2,0): w=(1/sqrt(2*pi))*           (sqrt(10)*(3*cos(theta)**2-1)/4)*   (40/(270*a0)**1.5)*(r/a0)**2*exp(-r/(3*a0))
    elif (n,l,ml)==(3,2,1): w=(exp(1j*phi)/sqrt(2*pi))* (sqrt(15)*sin(theta)*cos(theta)/2)* (40/(270*a0)**1.5)*(r/a0)**2*exp(-r/(3*a0))
    elif (n,l,ml)==(3,2,-1):w=(exp(-1j*phi)/sqrt(2*pi))*(sqrt(15)*sin(theta)*cos(theta)/2)* (40/(270*a0)**1.5)*(r/a0)**2*exp(-r/(3*a0))
    elif (n,l,ml)==(3,2,2): w=(exp(2j*phi)/sqrt(2*pi))* (sqrt(15)*sin(theta)**2/4)*         (40/(270*a0)**1.5)*(r/a0)**2*exp(-r/(3*a0))
    elif (n,l,ml)==(3,2,-2):w=(exp(-2j*phi)/sqrt(2*pi))*(sqrt(15)*sin(theta)**2/4)*         (40/(270*a0)**1.5)*(r/a0)**2*exp(-r/(3*a0))
    else: w=0
    return w



def atomicindextoquantumnumbers(ind):
    if ind==0:return (1,0,0)
    elif ind==1:return (2,0,0)
    elif ind==2:return (2,1,0)
    elif ind==3:return (2,1,1)
    elif ind==4:return (2,1,-1)
    elif ind==5:return (3,0,0)
    elif ind==6:return (3,1,0)
    elif ind==7:return (3,1,1)
    elif ind==8:return (3,1,-1)
    elif ind==9:return (3,2,0)
    elif ind==10:return(3,2,1)
    elif ind==11:return(3,2,-1)
    elif ind==12:return(3,2,2)
    elif ind==13:return(3,2,-2)
    else:return(0,0,0)

def molecularwavefunc(funcnum, r, phi, theta):
    result=0
    for i in range(numatomics):
        (n,l,ml)=atomicindextoquantumnumbers(i)
        result+=atomicwavefunc(n,l,ml,r,phi,theta)*LCAOmatrix[funcnum][i]
    return result

def configurationwavefunc(funcnum, r, phi, theta):
    result=0
    for i in range(nummolecul):
        result+=molecularwavefunc(i,r,phi,theta)*configmatrix[funcnum][i]
    return result

if __name__=="__main__":main()

