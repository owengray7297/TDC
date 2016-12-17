import sys

def main():
    ar=sys.argv
    if len(ar)!=3:print("something is rotten in the sys.argv! "+str(sys.argv))
    infile=ar[1]
    infile=open(infile)
    chemicalname=ar[2]
    outfile=open("/Users/owengray/Desktop/"+chemicalname+".xyz",'w')
    line=infile.readline()
    outfile.write("NaN\n"+chemicalname+"\n")
    while "_chem_comp_atom.pdbx_ordinal" not in line:line=infile.readline()
    line=infile.readline()
    while "#" not in line:
        lin=[t for t in line.strip().split()]
        outfile.write(lin[3]+'\t'+lin[9]+'\t'+lin[10]+'\t'+lin[11]+"\t\n")
        line=infile.readline()

def trytomakenum(thing):
    try:return int(thing)
    except ValueError:pass
    try:return float(thing)
    except ValueError:pass
    return thing

if __name__=="__main__":main()
