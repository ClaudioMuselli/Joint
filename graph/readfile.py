#Fucntion to read from a c++ style file.dat

def read_from_file(name):
    """Function to read from a file called name

    it return a matrix of the various numbers"""
    f= open(name,'r')
    lines=f.readlines()
    lines2=[x.replace("\t"," ") for x in lines]
    lines3=[x.replace("\n","") for x in lines2]
    result=[[float(x) for x in res.split()] for res in lines3]
    return result


