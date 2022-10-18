# Program Name: main.py
# Author: Aravinthen Rajkumar
# Description: Main engine for MD results analysis of polymer chains

import numpy as np
import time
import os
import glob

class nanolysis:
    def __init__(self,):
        self.origin = False
        self.walks = {}

        
        self.header = False # this is used to build new files of only a few variables.

        self.sim_time = None # the time of the simulation
        self.sim_data = [] # the data of a simulation file, populated by read()

        
        self.hosted = [] # this list contains atom positions of interest.
                         # it is populated through a variety of methods.

    def verify(self,):
        if self.origin == False:
            raise SystemError("Origin file has not been set.")
        if self.header == False:
            raise SystemError("No files have been read in.")

    # READING =============================================================================
    def set_origin(self, filepath):
        # Set the origin file, used to encode the global bead numbers into their
        # original molecules.
        t0 = time.time()
        self.origin = filepath
        with open(self.origin) as f:
            lines = [j.split('\t') for j in ([i.strip() for i in f][1:])]
            
            for l in lines:
                if len(l) == 8:
                    # l[0] the global number of the bead
                    # l[1] the walk number that the bead belongs to.
                    ln = int(l[1])
                    if ln not in self.walks:
                        self.walks[int(ln)] = [int(l[0])]
                    else:
                        self.walks[int(ln)].append(int(l[0]))

        t1 = time.time()
        print(f"Origin file read in {t1-t0} seconds.")
        
    def read(self, filename):
        # Given an xmol file, read all atom data.
        self.sim_data = []
        with open(filename) as f:
            lines = [j.split(' ') for j in ([i.strip() for i in f])]            

        # build the header
        # can be used for reading.
        self.header = lines[0:9]
        self.sim_time = int(self.header[1][0])
        
        # read atom positions
        for l in lines:
            if len(l) == 5:
                self.sim_data.append([int(l[0]),
                                      int(l[1]),
                                      float(l[2]),
                                      float(l[3]),
                                      float(l[4])])
    
    def hostwalk(self, index):
        self.verify() # make sure all required values have been read in        
        walk_coords = self.walks[index]        
        for i in walk_coords:
            for j in self.sim_data:
                if i == j[0]:
                    self.hosted.append(j)

    def hosttype(self, typenum):
        self.verify() # make sure all required values have been read in
        for i in self.sim_data:
            if i[1] == typenum:
                self.hosted.append(i)
                    
    # CALCULATION =============================================================
    def get_data(self, walkID, index):
        self.verify() # make sure all required values have been read in

        # obtain global bead number
        glob = self.walks[walkID][index]
        gindex = 0
        for i in self.sim_data:
            if i[0] == glob:
                gindex = i
                break
            
        return i

    def bonds(self, walkID):
        bond_index = 0
        backb = np.array(self.get_data(walkID, 0)[2:5])

        
        for bead in range(1,len(self.walks[walkID])):
            frontb = np.array(self.get_data(walkID, bead)[2:5])            
            bond = frontb-backb
            
            bond_index+=1
            backb = frontb


    # VISUALISATION ===========================================================                
    def newfile(self, filename, flush=True):
        # build a structure file with only a few components of the original.
        # filename: the name of the new file produced.
        # the flush option wipes the hosted list
        
        self.verify() # make sure all required values have been read in
        
        if self.hosted == []:
            raise SystemError("Nothing hosted.")

        with open(filename, 'w') as f:
            self.header[3] = [str(len(self.hosted))]
            for i in self.header:
                f.write(" ".join(i) + "\n")
                
            for j in range(1,len(self.hosted)+1):
                strj = [str(i) for i in self.hosted[j-1][1::]]
                f.write(f"{str(j)} " + " ".join(strj) + "\n")
        
        if flush==True:
            self.hosted = []

test = nanolysis()
test.set_origin("structure-120.0.in")

current = os.getcwd()
# test.read("dump.settings-120.0.in.in_100000.in")

# test.bonds(30)
# test.hostwalk(30)
# test.hosttype(1)

current = os.getcwd()
path = current + "/simulation"

foldername = "cut3"
os.mkdir(foldername)

total = len(os.listdir(path))
count = 1
for i in os.listdir(path):
    t0 = time.time()
    test.read(path+"/"+i)
    test.hostwalk(1050)
    test.hostwalk(1000)
    test.newfile(f"test-{test.sim_time}.in", flush=True)
    t1 = time.time()
    os.rename(f"{current}/test-{test.sim_time}.in",
              f"{current}/{foldername}/test-{test.sim_time}.in")
    print(f"({count}/{total}) Predicted time: {int(((total - count)*(t1-t0))/60)} minutes.", end="\r")
    count+=1
print("Done!")

# test.newfile('test.in')
