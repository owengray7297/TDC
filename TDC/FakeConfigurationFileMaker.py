from math import *
from random import *

# The n'th line consists of the coefficients of the molecular orbitals in the n'th configuration orbital.
def main():
    makefiles()

def makefiles(nummolecu=20, numconfig=50):
    seed(137)
    outfile = open("fakeconfiguration.txt",'w')
    for i in range(numconfig):
        rands=[0]*nummolecu
        for i in range(nummolecu):
            rands[i]=random()
        norm=sqrt(sum([r**2 for r in rands]))
        rands=[r/norm for r in rands]
        for r in rands:
            outfile.write(str(r)+' ')
        outfile.write('\n')
    outfile.close()
    print("make fake configuration coefficient file: done")

if __name__ =="__main__":main()
