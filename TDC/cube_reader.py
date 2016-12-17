def main():
#    tdcfile=open(thispath+thisname+"/"+thisname+"_root"+whichroot+".cube")
    tdcfile=open("/Users/owengray/Desktop/root1.cube")
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
        buffer1.extend([trytomakenum(n) for n in line.strip().split()])
        if len(buffer1)==xwid:buffer2.append(buffer1);buffer1=[]
        if len(buffer2)==ywid:tdc.append(buffer2);buffer2=[]
    printmatrix3(tdc)

def trytomakenum(thing):
    try:return int(thing)
    except ValueError:pass
    try:return float(thing)
    except ValueError:pass
    return thing

def makeshortifnum(thing):
    try:
        if int(thing)!=float(thing):raise ValueError
        return int(thing)
    except ValueError:pass
    try:
        return round(1000*float(thing))/1000
    except ValueError:pass
    return thing

def printmatrix(thismatrix):
    for line in thismatrix:
        for thing in line:
            print(makeshortifnum(thing),end='\t')
        print()
    print()

def printmatrix3(thismatrix):
    for line in thismatrix:
        for thing in line:
            for elem in thing:
                print(makeshortifnum(elem),end='\t')
            print('\t')
        print("------------------------------------------")
    print()

if __name__=="__main__":main()
