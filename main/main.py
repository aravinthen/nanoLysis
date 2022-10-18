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
        self.xlen = None
        self.ylen = None
        self.zlen = None
        self.sim_data = [] # the data of a simulation file, populated by read()

        self.hosted = [] # this list contains atom positions of interest.
                         # it is populated through a variety of methods.

        self.bonds = []
        self.grid = None

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
        self.xlen = float(lines[5][1])
        self.ylen = float(lines[6][1])
        self.zlen = float(lines[7][1])
        
        # read atom positions
        for l in lines:
            if len(l) == 5:
                self.sim_data.append([int(l[0]),
                                      int(l[1]),
                                      float(l[2])*self.xlen,
                                      float(l[3])*self.ylen,
                                      float(l[4])*self.zlen])
    
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

                
    def get_data(self, walkID, index):
        # get the data of a single bead
        # 1: global num
        # 2: type
        # 2: x
        # 2: y
        # 2: z        
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
        backb = np.array(self.get_datapp(walkID, 0)[2:5])

        
        for bead in range(1,len(self.walks[walkID])):
            frontb = np.array(self.get_data(walkID, bead)[2:5])            
            bond = frontb-backb

            print(np.linalg.norm(bond))
            
            bond_index+=1
            backb = frontb

    def cell(self, xlen, ylen, zlen):
        # this function returns all beads within a particular cell
        # all arguments are tuples.
        self.verify()
        if len(xlen) != 2:
            raise SystemError("Must define a tuple of two range values in the x dimension.")
        if len(ylen) != 2:
            raise SystemError("Must define a tuple of two range values in the y dimension.")
        if len(zlen) != 2:
            raise SystemError("Must define a tuple of two range values in the z dimension.")

        condx = lambda x: ((x >= xlen[0]) and (x <= xlen[1]))
        condy = lambda y: ((y >= ylen[0]) and (y <= ylen[1]))
        condz = lambda z: ((z >= zlen[0]) and (z <= zlen[1]))

        cell = []
        for bead in self.sim_data:
            if (condx(bead[2]) and condy(bead[3]) and condz(bead[4])):
                cell.append(bead)
                
        return cell 

    def minl_maxl(self,):
        self.verify()
        xmax, xmin = self.sim_data[0][2], self.sim_data[0][2]
        ymax, ymin = self.sim_data[0][3], self.sim_data[0][3]
        zmax, zmin = self.sim_data[0][4], self.sim_data[0][4]
        for bead in self.sim_data:
            if bead[2] < xmin:
                xmin = bead[2]
            if bead[2] > xmax:
                xmax = bead[2]
                
            if bead[3] < ymin:
                ymin = bead[3]
            if bead[3] > ymax:
                ymax = bead[3]

            if bead[4] < zmin:
                zmin = bead[4]
            if bead[4] > zmax:
                zmax = bead[4]

        self.xmin, self.xmax = xmin, xmax
        self.ymin, self.ymax = ymin, ymax
        self.zmin, self.zmax = zmin, zmax

    # CALCULATION =============================================================
            
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
test.set_origin("../data/structure-120.0.in")

current = os.getcwd()
test.read("../data/simulation/dump.settings-120.0.in.in_100000.in")

test.minl_maxl()

# test.bonds(30)
# test.hostwalk(30)
# test.hosttype(1)

current = os.getcwd()
path = current + "/../data/simulation"
