import os
import sys
import numpy as np
import scipy.io as sp
import hashlib


properties = ["Density", "Energies", "Hardness", "Individuals"]

def SplitPoscars(POSCARPath):
    fread = open(POSCARPath)
    count = 0
    hash1 = hashlib.md5()
    fwritePOSCAR = open("./../out/{0}".format(hash1.hexdigest()), 'w')
    fwriteTable = open("./../out/res.txt", "w+")
    for line in fread:
        if "EA" in line:
            fwritePOSCAR.close()
            count += 1
            hash1.update(str(count))
            fwritePOSCAR = open("./out/{0}".format(hash1.hexdigest()), 'w')
            fwritePOSCAR.write(line)
            fwriteTable.write(str(hash1.hexdigest()) + "\t\t" + str(count) + "\n")
        else:
            fwritePOSCAR.write(line)
    fwritePOSCAR.close()
    fread.close()
    fwriteTable.close()
    return

if __name__ == '__main__':  
    POSCARPath = sys.argv[1]
    SplitPoscars(POSCARPath)
    os.system('rm ./out/0')
    
