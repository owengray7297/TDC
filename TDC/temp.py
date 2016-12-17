def main():
    print(str(getlambdas('a')(2,3)))
    print(str(getlambdas('m')(2,3)))


thedict={"a":lambda x,y:x+y, "m":lambda x,y:x*y}

def getlambdas(which):
    return thedict[which]

if __name__=="__main__":main()
