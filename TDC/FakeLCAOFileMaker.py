from math import *
from random import *

# The n'th line consists of the coefficients of the atomic orbitals in the n'th molecular orbital.
# the atomic orbitals go 100, 200, 210, 211, 21m, 300, 310, 311, 31m, 320, 321, 321m, 322, 32n, etc. 
def main():
    makefiles()

def makefiles(nummolecu=20, numatomic=14):
    seed(137)
    outfile = open("fakeLCAO.txt",'w')
    for i in range(nummolecu):
        rands=[0]*numatomic
        for i in range(numatomic):
            rands[i]=random()
        norm=sqrt(sum([r**2 for r in rands]))
        rands=[r/norm for r in rands]
        for r in rands:
            outfile.write(str(r)+' ')
        outfile.write('\n')
    outfile.close()
    print("make fake LCAO coefficient file: done")

if __name__ =="__main__":main()
