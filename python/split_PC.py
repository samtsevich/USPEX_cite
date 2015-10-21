import os
import sys
import numpy as np
import scipy.io as sp
import hashlib


properties = ["Density", "Energies", "Hardness", "Individuals"]


def readEnergy(line):   # only last value
    line = line[line.rfind(" ", 0, len(line)) : len(line)]
    return line

def deleteEnd(line):
    tmpValue = float(line)
    return str(tmpValue)

def SplitPoscars(POSCAR_dir):
    fread = open(POSCAR_dir + "/POSCARS")
    count = 0
    count1 = 0
    hash1 = hashlib.md5()
    lineHash = str(hash1.hexdigest())

    # print(lineHash)

    fwritePOSCAR = open("./../out/" + lineHash, 'w')
    # fwritePOSCAR = open("./../out/{0}".format(hash1.hexdigest()), 'w')
    fwriteTable = open( "./../out/res.js", "w+")
    fwriteTable.write("var compounds = [\n")


    tmpLine = ""
    lineHardeness = ""
    lineDensity = ""
    lineEnergy = ""

    Denread = open(POSCAR_dir + "/Density", "r")
    Harread = open(POSCAR_dir + "/Hardness", "r")
    tmpLine = Harread.readline()
    tmpLine = Harread.readline()
    Enread = open(POSCAR_dir + "/Energies", "r")


    for line in fread:
        count1 += 1
        if "EA" in line:
            count1 = 0
            fwritePOSCAR.close()
            hash1.update(str(count))
            lineHash = str(hash1.hexdigest())
            fwritePOSCAR = open("./../out/" + lineHash, 'w')
            fwritePOSCAR.write(line)

            lineDensity = Denread.readline()
            # print(lineDensity)
            lineEnergy = readEnergy(Enread.readline())
            lineHardeness = Harread.readline()
 
            fwriteTable.write ( "{ id :    \"" + lineHash + "\" , density  :  " + lineDensity.rstrip() + ", hardness :  " + lineHardeness.rstrip() + ", energy: " + lineEnergy.rstrip())   
            # fwriteTable.write(lineHash + "\t\t" + str(count) + "\n")
            count += 1
        else:
            if count1 == 5:
                fwriteTable.write(" , atomTypes : " + str(list(set(line.split()))) + "},\n")
            fwritePOSCAR.write(line)



    fwritePOSCAR.close()
    fread.close()
    fwriteTable.write(" ]\n")
    fwriteTable.close()
    return

if __name__ == '__main__':  
    os.system('rm -rf ./../out/')
    os.system('mkdir ./../out')
    POSCARPath = sys.argv[1]
    SplitPoscars(POSCARPath)
    # os.system('rm ./../out/')
    
