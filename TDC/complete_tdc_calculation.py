molecule_name="derpderp"
molecule_atoms=[[]]
molecule_symmetry="c1"

thesecoords="""
  C        0.00000        1.40272        0.00000
  H        0.00000        2.49029        0.00000
  C       -1.21479        0.70136        0.00000
  H       -2.15666        1.24515        0.00000
  C       -1.21479       -0.70136        0.00000
  H       -2.15666       -1.24515        0.00000
  C        0.00000       -1.40272        0.00000
  H        0.00000       -2.49029        0.00000
  C        1.21479       -0.70136        0.00000
  H        2.15666       -1.24515        0.00000
  C        1.21479        0.70136        0.00000
  H        2.15666        1.24515        0.00000
  """
import sys
from math import *
    
def check_os():
    if sys.platform.startswith('win32'):
         return ("Windows")
    elif sys.platform.startswith('linux2'):
         return("Linux")
    elif sys.platform.startswith('cygwin'):
         return("Windows/Cygwin")
    elif sys.platform.startswith('darwin'):
         return("Mac OS X")
    elif sys.platform.startswith('os2'):
         return("OS/2")
    elif sys.platform.startswith('os2emx'):
         return("OS/2 EMX")
    elif sys.platform.startswith('riscos'):
         return("RiscOS")
    elif sys.platform.startswith('atheos'):
        return("AtheOS")

is_mac=(check_os()=="Mac OS X")

print(check_os())

coordinateunits="au"

thispath="/scratch/owengray/"                                        #linux path
if is_mac:thispath="/Users/owengray/Desktop/scratch/owengray/"       #mac path

nwpath="/cm/shared/NWChem_old/bin/nwchem"
if is_mac:nwpath="/Users/owengray/nwchem-6.6/bin/MACX64/nwchem"

basis_set="STO-3G"

bounds_edge_size=2    # the transition density is calculated for a cube, with the edges
#   this number beyond the outermost nuclei in each of x, y, and z. Units are coordinateunits.

maxiter=1000    #this determines how many iterations are performed in the geometry optimization

#   nwchem places these files outside of the run scratch directory "thispath/thisname" for some reason,
#   files with extensions listed here get automatically moved back in. The TDC files are also placed outside,
#   but because there may be multiple thisname_root#.cube files, those are handled seperately
extensions_to_move=[".err",".drv.hess",".movecs",".civecs_singlet",".civecs_triplet",".db",".log"]

#   this list is used in extracting things from the nwchem output files. If your molecule contains an atom
#   that is not on this list, add it. Removing things from this list is never necessary
atomnames=("Oxygen","Hydrogen","Carbon","Nitrogen","Sulfur","Magnesium","Iron","Phosphorus","Nickel")
#   anything on the previous list should be on this list to. The number is the number of unique orbitals:
#   "P" is 9: 1s, 2s, 3s, 2px, 2py, 2pz, 3px, 3py, 3pz is nine different orbitals
uniqueorbitaldict={"H":1, "O":5, "C":5, "N":5, "P":9, "S":9, }


#   A list of which roots get TDCs printed for the transition groundstate:root#
whichroots=[1,2]

thisnum=0#    If this program is run multiple times on the same molecule, each run will be sequentially labelled
thisname=""#  The name of this run. Of the form <molecule>_TDC_<thisnum>
slurmnum=0#   The number of this run's slurm job. Used for tracking program completion

doNW=True

import subprocess
import time
import sys

def main_donw():
    make_nwfile()
    make_bashfile()
    run_nwchem()
    move_files()
    extract_orbitals()
    compute_TDC()
#    do_cleanup()    #this deletes the large files, like the TDC file. Comment this if those are desired
    print("DONE!!1!1")

def main():
    global molecule_name
#    try:molecule_atoms=convertcoordinatestolist(thesecoords.strip())
#    except:return
    ar=sys.argv
    if len(ar)==2:
        molecule_name=ar[1]
        make_filestructure(1)
        do_readxyz()
        main_donw()
    if len(ar)==3:
        if "xyzonly" in ar[2]:
            molecule_name=ar[1]
            do_readxyz()
            exit(0)
        if "readexisting" in ar[2]:
            molecule_name=ar[1]
            do_readxyz()
            make_nwfile()
            make_filestructure(0)
            extract_orbitals()
            print("TDM")
            print(transitiondensitymatricies)
            compute_TDC()

def do_readxyz():
    global molecule_atoms
    xyzfile=open(thispath+molecule_name+".xyz")
    xyzfile.readline()#number of atoms line
    xyzfile.readline()#molecule name line
    for line in xyzfile:
        molecule_atoms.append([trytomakenum(t) for t in line.strip().split()])
    molecule_atoms.pop(0)
    print(molecule_atoms)

def do_cleanup():
    #for whichroot in whichroots:subprocess.run(["rm","-rf",thispath+thisname+"/"+thisname+"_root"+str(whichroot)+".cube"])
    subprocess.run(["rm","-rf",thispath+thisname+"/"+thisname+".db"])


def move_files():
    for extens in extensions_to_move:
        subprocess.run(["mv",thispath+thisname+extens,thispath+thisname+"/"+thisname+extens],stderr=subprocess.DEVNULL)
    for rootnum in range(10):#this is bad, needs to be fixed.
        subprocess.run(["mv",thispath+thisname+"_root"+str(rootnum)+".cube",thispath+thisname+"/"+thisname+"_root"+str(rootnum)+".cube"],stderr=subprocess.DEVNULL)
    for i in range(maxiter+1):
        subprocess.run(["mv",thispath+"final-"+make_3digit(i)+".xyz",thispath+thisname+"/"+"final-"+str(i)+".xyz"],stderr=subprocess.DEVNULL)        

def make_3digit(num):
    num=str(num)
    return "0"*(3-len(num))+num
#def make_3digit(num):
#    return str(roundtodigits(num,3))

def roundtodigits(num,digits):
    factor=floor(-log(num)/log(10)+digits)
    return round(num*(10**factor))/(10**factor)

def wait_until_complete():
    sq=subprocess.run(["squeue"], stdout=subprocess.PIPE).stdout.decode("utf-8").split("\n")
    sq=[j.split()[0] for j in sq if len(j.split())>0]
    print("squeue returns:")
    print(sq)
    stillrunning=str(slurmnum) in sq
    print(stillrunning)
    t=0
    while stillrunning:
        time.sleep(10)
        t+=10
        sq=subprocess.run(["squeue"], stdout=subprocess.PIPE).stdout.decode("utf-8").split("\n")
        sq=[j.split()[0] for j in sq if len(j.split())>0]
        stillrunning=str(slurmnum) in sq
    print("finished within "+str(t)+" seconds")

def run_nwchem():
    global slurmnum
    if not is_mac:
        slurmnum=int(subprocess.run(["sbatch", thispath+thisname+"/"+thisname+".bash","-J","NwPY_"+thisname,"--shared"], stdout=subprocess.PIPE).stdout.decode("utf-8").split()[3])
        print(thisname+" has slurm job number "+str(slurmnum))
        wait_until_complete()
    else:
        subproc = subprocess.run(["bash",thispath+thisname+"/"+thisname+".bash"])
        t=time.time()
        subproc.communicate()#  waits until the process completes
        print("finished in "+str(time.time()-t)+" seconds")

def make_bashfile():

    bfile=open(thispath+thisname+"/"+thisname+".bash",'w')
    if not is_mac:
        bfile.write("#!/bin/bash\n")
        bfile.write("#SBATCH --partition=gpu\n#SBATCH -N 1\n#SBATCH --share\n#SBATCH --output="+thisname+".log\n#SBATCH --error="+thisname+".err\n")
        bfile.write("source /etc/profile.d/modules.sh\nmpirun -np 20 "+nwpath+" "+thispath+thisname+"/"+thisname+".nw\nexit 0\n")
    else:
        bfile.write(nwpath+" "+thispath+thisname+"/"+thisname+".nw\nexit 0\n")
    bfile.close()

cubecoords=[]

rootstotest=0

def make_nwfile():
    global cubecoords,rootstotest
    nwfile=open(thispath+thisname+"/"+thisname+".nw",'w')
    nwfile.write("echo\n\nstart "+thisname+"\ntitle \""+molecule_name+" TDC calculation number "+str(thisnum)+" in "+basis_set+" basis set\"\n")
    nwfile.write("scratch_dir "+thispath+thisname+"/\n")
    nwfile.write("permanent_dir "+thispath+thisname+"/\n")
    nwfile.write("geometry units "+coordinateunits+" noautoz nocenter\nsymmetry "+molecule_symmetry+"\n")
    for atom in molecule_atoms:nwfile.write("\t"+"\t".join([str(a) for a in atom])+"\n")
    #basis:
    nwfile.write("end\nbasis\n")
    atomtypes=list(set([at[0] for at in molecule_atoms]))
    for typ in atomtypes:nwfile.write("\t"+typ+" library "+basis_set+"\n")
    #driver:
    nwfile.write("end\ndriver\n\tclear\n\tmaxiter "+str(maxiter)+"\n\txyz final\n")
    #dft:
    nwfile.write("end\ncharge 0\ndft\n\txc b3lyp\n\titerations 100\n\tprint \"final vectors analysis\"\n\tvectors output "+thisname+".movecs\nend\ntask dft optimize\n")
    #tddft:
    rootstotest=10*ceil(2*len(whichroots)/10)
    nwfile.write("tddft\n\tnroots "+str(rootstotest)+"\n\tcis\n\tprint debug\n\tnotriplet\n\tcivecs\nend\ntask tddft energy\n")
    #dplot:
    bounds=[(roundtopoint2(min([at[i+1] for at in molecule_atoms])-bounds_edge_size),roundtopoint2(max([at[i+1] for at in molecule_atoms])+bounds_edge_size)) for i in range(3)]
    numpoints=[round((bounds[i][1]-bounds[i][0])/0.2)+1 for i in range(3)]
    cubecoords=bounds+numpoints
    for whichroot in whichroots:
        nwfile.write("dplot\n\tvectors "+thisname+".movecs\n\tcivecs "+thisname+".civecs_singlet\n\troot "+str(whichroot)+"\n\tlimitxyz\n")
        for i in range(3):nwfile.write("\t\t"+"\t".join([str(l) for l in bounds[i]])+"\t"+str(numpoints[i]-1)+"\n") #dplot wants # of SPACINGS not # of POINTS
        nwfile.write("\t\tgaussian\n\t\toutput "+thisname+"_root"+str(whichroot)+".cube\nend\ntask dplot\n")
    nwfile.close()

#if the argument is 1, creates a new folder for a new run. If the argument is 0, accesses the most recently created existing folder
def make_filestructure(donew):
    global thisnum,thisname
    ls=subprocess.run(["ls", thispath], stdout=subprocess.PIPE).stdout.decode("utf-8").split()
    matchnums=[trytomakenum(s.replace(molecule_name+"_TDC_","")) for s in ls if molecule_name+"_TDC" in s and "." not in s]
    matchnums=[a for a in matchnums if isinstance(a, int)]
    print(str(matchnums))
    thisnum=max(matchnums+[0])+int(bool(donew))
    thisname=molecule_name+"_TDC_"+str(thisnum)
    if donew:subprocess.run(["mkdir",thispath+thisname])

def convertcoordinatestolist(coords):
    outcoor=coords.split('\n')
    outcoor=[[trytomakenum(t) for t in oc.strip().split()] for oc in outcoor]
    return outcoor

def printmatrix(thismatrix):
    for line in thismatrix:
        for thing in line:
            print(makeshortifnum(thing),end='\t')
        print()
    print()

def hasanatomname(thisstring):
    for atomn in atomnames:
        if atomn in thisstring:return True
    return False

def roundtopoint2(num):
    return int(num*5)/5

def isanint(thisstring):
    try: 
        int(thisstring)
        return True
    except ValueError:
        return False

def makeshortifnum(thing):
    try:
        if int(thing)!=float(thing):raise ValueError
        return int(thing)
    except ValueError:pass
    try:
        return round(1000*float(thing))/1000
    except ValueError:pass
    return thing

def trytomakenum(thing):
    try:return int(thing)
    except ValueError:pass
    try:return float(thing)
    except ValueError:pass
    return thing


def findnumelectrons(molat):
    idents=[at[0] for at in molat]
    return sum([uniqueorbitaldict[a] for a in idents])

atomicorbitalparameters, molecularorbitalparameters, transitiondensitymatricies ,atomicorbitallist= [], [], [], []

def extract_orbitals():
    global atomicorbitalparameters, molecularorbitalparameters, transitiondensitymatricies, atomicorbitallist
    
    numatoms=len(molecule_atoms)
    numouterelectrons=findnumelectrons(molecule_atoms)

    atomicorbitalparameters={}

    atomicorbitallist=[""]*numouterelectrons

    print("number of a/m orbitals: "+str(numouterelectrons))

    molecularorbitalparameters=[[0 for i in range(numouterelectrons)] for j in range(numouterelectrons)]

    transitiondensitymatricies=[[[0 for i in range(numouterelectrons)] for j in range(numouterelectrons)] for k in range(rootstotest)]

    transitiondensitymatrixnodiagonal=[[0 for i in range(numouterelectrons)] for j in range(numouterelectrons)]

    excitedstatedensitymatrix=[[0 for i in range(numouterelectrons)] for j in range(numouterelectrons)]

    groundoccupation=[]

    excitedoccupation=[0]*max(numouterelectrons,rootstotest)

    cleanoutput=open(thispath+thisname+"/"+thisname+".cleanoutput",'w')
    with open(thispath+thisname+"/"+thisname+".log") as inF:
        for line in inF:
            if "Basis \"" in line:
                print("foundit: atomic orbitals")
                temp=""
                curratom=""
                whichpair=0;
                while("Summary of " not in temp):
                    temp=inF.readline().strip()
                    if hasanatomname(temp):curratom=temp[0]
                    if len(temp)>0 and isanint(temp[0]):
                        temp2=temp.split()
                        print(temp2)
                        atomicorbitalparameters[curratom+temp2[0]+temp2[1]+str(whichpair)]=(float(temp2[2]),float(temp2[3]))
                        whichpair=(whichpair+1)%3
 #               print(str(atomicorbitalparameters))
            if "Final Molecular Orbital Analysis" in line:
                print("foundit: molecular orbitals")
                temp=""
                currmo=""
                whichpair=0;
                while("center of mass" not in temp):
                    temp=inF.readline().strip()
                    if "Vector" in temp:
                        currmo=int(temp.split()[1])-1
 #                       molecularorbitalparameters[currmo]=[0]*len(atomicorbitalparameters)
                    if len(temp)>0 and isanint(temp[0]):
                        temp2=temp.split()
#                        print(temp2)
                        molecularorbitalparameters[currmo][int(temp2[0])-1]=float(temp2[1])
                        atomicorbitallist[int(temp2[0])-1]=str(temp2[2]+temp2[3]+temp2[4])
                        #print(int(temp2[0])-1)
                        #print(temp2[2]+temp2[3]+temp2[4])
                        if len(temp2)<=5:continue
 #                       print("len mop: "+str(len(molecularorbitalparameters))+" len mop[currmo] "+str(len(molecularorbitalparameters[currmo])))
                        molecularorbitalparameters[currmo][int(temp2[5])-1]=float(temp2[6])
                        atomicorbitallist[int(temp2[5])-1]=temp2[7]+temp2[8]+temp2[9]
##            if "global array: gdens1" in line:
##                continue # basically a comment
##                print("foundit: transition density matrix")
##                temp=""
##                rowwid=0
##                blockoffset=0
##                rowid=0
##                whichpair=0
##                while("global array: gdens+1" not in temp):
##                    temp=inF.readline().strip()
##                    if len(temp)<1:continue
##                    temp2=temp.split()
###                    print(temp2)
##                    if isanint(temp2[0]) and (len(temp2)==1  or isanint(temp2[1]) and int(temp2[1])-int(temp2[0])==1):
##                        rowwid=len(temp2)
##                        blockoffset=int(temp2[0])
## #                       print("rowwid="+str(rowwid))
###                        print("blockoffset="+str(blockoffset))
##                        continue
##                    if isanint(temp2[0]) and not isanint(temp2[1]):
##                        for i in range(len(temp2)):
##                            if i==0:rowid=int(temp2[i]);continue
##                            transitiondensitymatrix[rowid-1][i-1+blockoffset-1] = float(temp2[i])
##                line=temp
            if "global array: gdens+1" in line:
                print("foundit: excited state density matrix")
                temp=""
                rowwid=0
                blockoffset=0
                rowid=0
                whichpair=0
                while("grid_file" not in temp):
                    temp=inF.readline().strip()
                    if len(temp)<1:continue
                    temp2=temp.split()
#                    print(temp2)
                    if isanint(temp2[0]) and (len(temp2)==1  or isanint(temp2[1]) and int(temp2[1])-int(temp2[0])==1):
                        rowwid=len(temp2)
                        blockoffset=int(temp2[0])
#                        print("rowwid="+str(rowwid))
#                        print("blockoffset="+str(blockoffset))
                        continue
                    if isanint(temp2[0]) and not isanint(temp2[1]):
                        for i in range(len(temp2)):
                            if i==0:rowid=int(temp2[i]);continue
                            excitedstatedensitymatrix[rowid-1][i-1+blockoffset-1] = float(temp2[i])
            if "Root   " in line:
                print("foundit: an excited state occupation")
                lin=line.split()
                rootnum=int(lin[1])
                temp=""
                rowwid=0
                blockoffset=0
                rowid=0
                whichpair=0
                while("Occ." not in temp):
                    temp=inF.readline().strip()
                    if len(temp)<1:continue
                    if "global array: Transition" in temp:
                        print("foundit: transition density matrix")# these matrices are not symmetrized, to fix add the matrix to its transpose (doubles diag)
                        temp=""
                        rowwid=0
                        blockoffset=0
                        rowid=0
                        whichpair=0
                        while("--------------------------------------" not in temp):
                            temp=inF.readline().strip()
                            if len(temp)<1:continue
                            temp2=temp.split()
        #                    print(temp2)
                            if isanint(temp2[0]) and (len(temp2)==1  or isanint(temp2[1]) and int(temp2[1])-int(temp2[0])==1):
                                rowwid=len(temp2)
                                blockoffset=int(temp2[0])
         #                       print("rowwid="+str(rowwid))
        #                        print("blockoffset="+str(blockoffset))
                                continue
                            if isanint(temp2[0]) and not isanint(temp2[1]):
                                for i in range(len(temp2)):
                                    if i==0:rowid=int(temp2[i]);continue
                                    transitiondensitymatricies[rootnum][rowid-1][i-1+blockoffset-1] = float(temp2[i])
                        line=temp
                    temp2=temp.split()
#                    print(temp2)
                    if isanint(temp2[0]) and (len(temp2)==1  or isanint(temp2[1]) and int(temp2[1])-int(temp2[0])==1):
                        rowwid=len(temp2)
                        blockoffset=int(temp2[0])
#                        print("rowwid="+str(rowwid))
#                        print("blockoffset="+str(blockoffset))
                        continue
                    if isanint(temp2[0]) and not isanint(temp2[1]):
                        for i in range(len(temp2)):
                            if i==0:rowid=int(temp2[i]);continue
                            transitiondensitymatrixnodiagonal[rowid-1][i-1+blockoffset-1] = float(temp2[i])
                thisstate=[]
                while("Occ." in temp):
                    tem=temp.split()
                    print(tem)
                    thisstate+=[tem[1]+tem[2]+":"+tem[5]+tem[6]]
                    temp=inF.readline().strip()
                excitedoccupation[rootnum-1]=thisstate
    cleanoutput.write("atomicorbitalparameters\n")
    cleanoutput.write(str(atomicorbitalparameters))
    cleanoutput.write("\natomicorbitallist\n")
    cleanoutput.write(str(atomicorbitallist))
    cleanoutput.write("\nmolecularorbitalparameters\n")
    cleanoutput.write(str(molecularorbitalparameters))
    cleanoutput.write("\ntransitiondensitymatrixnodiagonal\n")
    cleanoutput.write(str(transitiondensitymatrixnodiagonal))
    cleanoutput.write("\nexcitedstatedensitymatrix\n")
    cleanoutput.write(str(excitedstatedensitymatrix))
    cleanoutput.write("\ngroundoccupation\n")
    cleanoutput.write(str(groundoccupation))
    cleanoutput.write("\nexcitedoccupation\n")
    cleanoutput.write(str(excitedoccupation))
#    computeMOlambdas(molecularorbitalparameters, atomicorbitalparameters, atomicorbitallist)
#    print(atomicorbitalparameters)
#    print(atomicorbitallist)
#    printmatrix(molecularorbitalparameters)
#    printmatrix(transitiondensitymatrix)
#    printmatrix(excitedstatedensitymatrix)
    cleanoutput.close();


def tdc_readin(filepath):
    tdcfile=open(filepath)
    tdc=[[[]]]
    tdcfile.readline()
    tdcfile.readline()
    tdcfile.readline()
    lin=[trytomakenum(n) for n in tdcfile.readline().strip().split()]
    xwid=lin[0]
    lin=[trytomakenum(n) for n in tdcfile.readline().strip().split()]
    ywid=lin[0]
    lin=[trytomakenum(n) for n in tdcfile.readline().strip().split()]
    zwid=lin[0]
    tdcfile.readline()
    tdcfile.readline()
    tdcfile.readline()
    buffer1=[]
    buffer2=[]
    for line in tdcfile:
        buffer2.extend([trytomakenum(n) for n in line.strip().split()])
        if len(buffer1)==xwid:buffer2.append(buffer1);buffer1=[]
        if len(buffer2)==ywid:tdc.append(buffer2);buffer2=[]
    return tdc

def mult_tdc():
    tdc1=tdc_readin(thispath+thisname+"/"+thisname+"_root"+whichroot+".cube")
    tdc2=tdc_readin(thispath+thisname+"/"+thisname+"_root"+whichroot+".cube")
    tdc_mult=[[[ tdc1[i][j][k]*tdc2[i][j][k] for k in range(len(tdc1[0][0]))] for j in range(len(tdc1[0]))] for i in range(len(tdc1))]
    return mult_tdc

STO3Gg={"1s":lambda a,x,r:(2*a/pi)**(3/4)*exp(-a*r*r),"2p":lambda a,x,r:(128*a**5/pi**3)**(1/4)*x*exp(-a*r*r)}

STO3Gd = {"K" : (0.44463454, 0.53532814, 0.15432897), 
   "L" : ((0.70011547, 0.39951283, -0.09996723), (0.39195739, 
      0.60768372, 0.15591627))};
STO3Ga = {"K" : (0.109818, 0.405771, 2.22766), 
   "L" : (0.0751386, 0.231031, 0.994203)};

StandardMolecularExponents = {"H" : {"K" : 1.24}, 
   "C" : {"K" : 5.67, "L" : 1.72}, 
   "N" : {"K" : 6.67, "L" : 1.95}, 
   "O" : {"K" : 7.66, "L" : 2.25}};

def sto3g(what, which, x, r):
#    print("sto3g call with args "+str(what)+","+str(which)+","+str(x)+","+str(r))
    if which=="1s":
        return sum([STO3Gd["K"][i]*STO3Gg["1s"](STO3Ga["K"][i]*StandardMolecularExponents[what]["K"]**2, x, r) for i in range(3)])
    elif which=="2s":
        return sum([STO3Gd["L"][0][i]*STO3Gg["1s"](STO3Ga["L"][i]*StandardMolecularExponents[what]["L"]**2, x, r) for i in range(3)])
    elif which=="2p":
        return sum([STO3Gd["L"][1][i]*STO3Gg["2p"](STO3Ga["L"][i]*StandardMolecularExponents[what]["L"]**2, x, r) for i in range(3)])
    else:print("DERPDERP sto3g fail with args "+str(what)+","+str(which)+","+str(x)+","+str(r))

print(sto3g("H","1s",1,1))
print(STO3Gg["1s"])

#this function turns x, y, z (from e.g. 5Npx1) into 1, 2, 3 and turns nothing (s orbital) into 0.
def gettheXnum(orbtype):
    thestr = orbtype[3]
    if thestr=='x':return 1
    if thestr=='y':return 2
    if thestr=='z':return 3
    return 0

def sublis4s(lis1, lis2):
    return [lis1[0]-lis2[0],lis1[1]-lis2[1],lis1[2]-lis2[2],lis1[3]-lis2[3]]

lowestpossible={"s":'1', "p":'2'}

import dis

#uses the MOs in terms of the AOs and the molecule_atoms array to compute a lambda function taking (x,y,z) to the value of the orbital at that point
def computeMOlambdas(MOinAO, aol):
#    print(MOinAO)
#    print(aop)
#    print(aol)
    highestoftype=[]
    result=[0]*len(MOinAO)
    thelambdas=[(0,)]*len(MOinAO)# a list of tuples. Each tuple is (atomtype, orbitaltype, [xoffset, yoffset, zoffset])
    for i in range(len(aol)):
        orbtype=""
        for j in range(len(highestoftype)):
            if aol[i] in highestoftype[j]:
                highestoftype[j]=aol[i]+str(int(highestoftype[j][-1:])+1)
                orbtype=highestoftype[j]
                break
        else:highestoftype.append(aol[i]+'1');orbtype=aol[i]+lowestpossible[aol[i][2]]
        if orbtype[2]=='s':
            print(orbtype)
            thelambdas[i]=(aol[i][1], orbtype, molecule_atoms[int(aol[i][0])-1][1:4])
        if orbtype[2:-1]=='px':
            print(orbtype)
            thelambdas[i]=(aol[i][1], orbtype, molecule_atoms[int(aol[i][0])-1][1:4])
        if orbtype[2:-1]=='py':
            print(orbtype)
            thelambdas[i]=(aol[i][1], orbtype, molecule_atoms[int(aol[i][0])-1][1:4])
        if orbtype[2:-1]=='pz':
            print(orbtype)
            thelambdas[i]=(aol[i][1], orbtype, molecule_atoms[int(aol[i][0])-1][1:4])
    for i in range(len(MOinAO)):
#        result[i]=lambda x,y,z,ii=i:sum([MOinAO[j][ii]*sto3g(thelambdas[j][0], thelambdas[j][1][-1]+thelambdas[j][1][2], sublis4s([0,x,y,z],[0]+(thelambdas[j][2]))[gettheXnum(thelambdas[j][1])], sqrt(x**2+y**2+z**2)) for j in range(len(MOinAO))])
        result[i]=lambda x,y,z,ii=i:
    return result

def computetransitionlambda(TDM, rootnum, MOinAO, aol):
    molamb=computeMOlambdas(MOinAO, aol)
    printmatrix(MOinAO)
#    return lambda x,y,z:sum([('''print(rootnum,i,TDM[rootnum][i],molamb[i](x,y,z))''',TDM[rootnum][i]*molamb[i](x,y,z))[1] for i in range(len(TDM[rootnum]))])/sqrt(x**2+y**2+z**2+10**(-30))
    return lambda x,y,z:sum([TDM[rootnum][i][j]*molamb[i](x,y,z)*molamb[j](x,y,z)/sqrt(x**2+y**2+z**2+10**(-30)) for i in range(len(TDM[rootnum])) for j in range(len(TDM[rootnum]))])

#import itertools

def compute_TDC():
    for whichroot in whichroots:
        outfile=open(thispath+thisname+"/"+thisname+"_manual_root"+str(whichroot)+".cube",'w')
        lambdas=computetransitionlambda(transitiondensitymatricies, whichroot, molecularorbitalparameters, atomicorbitallist)
        c=0
        print(cubecoords)
        print(transitiondensitymatricies[whichroot])
        outfile.write("Cube file generated by python program\n"+str(cubecoords[3]+1)+"\n"+str(cubecoords[4]+1)+"\n"+str(cubecoords[5]+1)+"\n")
        for i in floatrange(cubecoords[0][0], cubecoords[0][1], 0.2):
            for j in floatrange(cubecoords[1][0], cubecoords[1][1], 0.2):
                for k in floatrange(cubecoords[2][0], cubecoords[2][1], 0.2):
#                    print("writing cube datum "+'{: 06.5E}'.format(lambdas(i, j, k)) +" aka ("+str(lambdas(i, j, k))+ ") at point "+str(i)+", "+str(j)+", "+str(k))
#                    dis.dis(lambdas)
#                    exit(0)
                   # print(k)
                    outfile.write('{: 06.5E}'.format(lambdas(i, j, k))+' ')
                    c+=1
                    if c==6:outfile.write('\n');c=0
                outfile.write('\n');c=0
        outfile.close()

def floatrange(start, end, step):
    for i in range(round(start/step), round(end/step)+1, 1):
        yield round(i*step,1)

def matrix_mult(m1,m2):
    return [[sum([m1[i][k]*m2[k][j] for k in range(len(m1[0]))]) for j in range(len(m2[0]))] for i in range(len(m1))]

if __name__=="__main__":main()
