import sys

def main():
    ar=sys.argv
    if len(ar)!=3:print("something is rotten in the sys.argv! "+str(sys.argv))
    infile=ar[1]
    infile=open(infile)
    chemicalname=ar[2]
    outfile=open("/Users/owengray/Desktop/"+chemicalname+".xyz",'w')
    line=infile.readline();line=infile.readline()
    while len(line.strip())==0 or not isnum(line.strip().split()[0]):
        line=infile.readline()
    numatoms=trytomakenum(line.strip().split()[0])
    outfile.write(line.strip().split()[0]+"\n"+chemicalname+"\n")
    for i in range(numatoms):
        line=infile.readline()
        lin=[t for t in line.strip().split()]
        outfile.write(lin[3]+'\t'+lin[0]+'\t'+lin[1]+'\t'+lin[2]+"\n")
    outfile.close()

def isnum(thing):
    return not (thing == trytomakenum(thing))

def trytomakenum(thing):
    try:return int(thing)
    except ValueError:pass
    try:return float(thing)
    except ValueError:pass
    return thing

if __name__=="__main__":main()
