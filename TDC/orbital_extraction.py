
numatoms=3
numouterelectrons=7         #Not sure where this number comes from- 1/H * 2 + 4 or 6 O = 6 or 8, not 7

atomicorbitalparameters={}

atomicorbitallist=[""]*numouterelectrons

molecularorbitalparameters=[""]*numouterelectrons

transitiondensitymatrix=[[0 for i in range(numouterelectrons)] for j in range(numouterelectrons)]

transitiondensitymatrixnodiagonal=[[0 for i in range(numouterelectrons)] for j in range(numouterelectrons)]

excitedstatedensitymatrix=[[0 for i in range(numouterelectrons)] for j in range(numouterelectrons)]

groundoccupation=[]

excitedoccupation=[0]*numouterelectrons

atomnames=("Oxygen","Hydrogen","Carbon","Nitrogen")

def main():
    global atomicorbitalparameters, atomicorbitallist, molecularorbitalparameters, transitionmatrix
    #outfile=open("")
    with open("water_TDC_13.log") as inF:
        for line in inF:
            if 'Basis \"ao basis\" -> \"\" (cartesian)' in line:
                print("foundit: atomic orbitals")
                temp=""
                curratom=""
                whichpair=0;
                while("Summary of \"ao basis\" -> \"\" (cartesian)" not in temp):
                    temp=inF.readline().strip()
                    if hasanatomname(temp):curratom=temp[0]
                    if len(temp)>0 and isanint(temp[0]):
                        temp2=temp.split()
                        #print(temp2)
                        atomicorbitalparameters[curratom+temp2[0]+temp2[1]+str(whichpair)]=(float(temp2[2]),float(temp2[3]))
                        whichpair=(whichpair+1)%3
            if "Final Molecular Orbital Analysis" in line:
                print("foundit: molecular orbitals")
                temp=""
                currmo=""
                whichpair=0;
                while("center of mass" not in temp):
                    temp=inF.readline().strip()
                    if "Vector" in temp:
                        currmo=int(temp.split()[1])-1
                        molecularorbitalparameters[currmo]=[0]*len(atomicorbitalparameters)
                    if len(temp)>0 and isanint(temp[0]):
                        temp2=temp.split()
                        #print(temp2)
                        molecularorbitalparameters[currmo][int(temp2[0])-1]=temp2[1]
                        atomicorbitallist[int(temp2[0])-1]=str(temp2[2]+temp2[3]+temp2[4])
                        #print(int(temp2[0])-1)
                        #print(temp2[2]+temp2[3]+temp2[4])
                        if len(temp2)<=5:continue
                        molecularorbitalparameters[currmo][int(temp2[5])-1]=temp2[6]
                        atomicorbitallist[int(temp2[5])-1]=temp2[7]+temp2[8]+temp2[9]
            if "global array: gdens1" in line:
                continue # basically a comment
                print("foundit: transition density matrix")
                temp=""
                rowwid=0
                blockoffset=0
                rowid=0
                whichpair=0
                while("global array: gdens+1" not in temp):
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
                            transitiondensitymatrix[rowid-1][i-1+blockoffset-1] = float(temp2[i])
                line=temp
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
                tem=temp.split()
                print(tem)
                excitedoccupation[rootnum-1]=tem[1]+tem[2]+"->"+tem[5]+tem[6]
                
#    print(atomicorbitalparameters)
#    print(atomicorbitallist)
#    printmatrix(molecularorbitalparameters)
#    printmatrix(transitiondensitymatrix)
#    printmatrix(excitedstatedensitymatrix)
    #outfile.close();

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

if __name__ == "__main__": main()
