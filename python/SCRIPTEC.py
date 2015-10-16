import os
import sys
import numpy as np
#   import scipy.io as sp
import hashlib
lineDensity=[]
linehar=[]
lineenergy=[]
linehash=[]
#C:\Users\user\Documents\skoltech\1year\innovation_workhsop\python_xuina\Coevo\4_TeTe
#POSCAR_dir= "C:/Users/user/Documents/skoltech/1year/innovation_workhsop/python_xuina/Coevo/4_TeTe/"


def SplitPoscars(POSCAR_dir):
    Denread = open(POSCAR_dir+"/Density","r")
    Harread = open(POSCAR_dir+"/Hardness","r")
    Enread = open(POSCAR_dir+"/Energies","r")
    fwriteHTML = open("res.txt", "w+")
    hash1 = hashlib.md5()
    count=0    
    for line in Denread:
        count+=1
        hash1.update(str(count))
        linehash.append(hash1.hexdigest())
        line=line[:len(line)-1]
        lineDensity.append(line)
    count=0       
    for line in Harread:
        if (count>=2):        
            line=line[:len(line)-1]
            linehar.append(line)
        count+=1
    for line in Enread:
        line=line[line.rfind(" ", 0,len(line)) : len(line)]
        lineenergy.append(line)
    i=0       
    for i in range(len(lineenergy)-1):
        fwriteHTML.write ("id:    "+str(linehash[i]) +"   density: " + lineDensity[i] + " hardness: "+ linehar[i]+" energy: "+lineenergy[i])   
    Denread.close()
    Harread.close()

    return

def splitPoscarDir(basedir):    
    withP=[]
    without=[]
    files=os.listdir(basedir)
    for name in files:
         basedir1=basedir+name
         if not os.path.isfile(basedir1):
            files1=os.listdir(basedir1)              
            for name1 in files1:
                if name1 == "POSCARS":
                    withP.append(name)
    #            else:
    #                count+=1
    #                if count== len(files1):
    #                    without.append(name)
         else:
             without.append(name)
    return withP


if __name__ == '__main__':  
 
#    basedir = sys.argv[2] ---Po xoroshemy wot tak
    
    basedir = "C:/Users/user/Documents/skoltech/1year/innovation_workhsop/python_xuina/Coevo/"
    basedir = "../res/"   
    poscarPath=splitPoscarDir(basedir)
    for name in poscarPath:     
        SplitPoscars(basedir+name)
    
    