# Program Name: spaghetti.py
# Author: Aravinthen Rajkumar
# Description: Testing out the nanoLysis suite capabilities

import sys
import os
import time
sys.path.insert(0, '../main/')

from main import nanolysis

test = nanolysis()
test.set_origin("../data/structure-120.0.in")

current = os.getcwd()
# test.read("dump.settings-120.0.in.in_100000.in")

# test.bonds(30)
# test.hostwalk(30)
# test.hosttype(1)

current = os.getcwd()
path = current + "/../data/simulation"

foldername = "cut3"
os.mkdir(foldername)

total = len(os.listdir(path))
count = 1
for i in os.listdir(path):
    t0 = time.time()
    test.read(path+"/"+i)
    test.hostwalk(1050)
    test.newfile(f"test-{test.sim_time}.in", flush=True)
    t1 = time.time()
    os.rename(f"{current}/test-{test.sim_time}.in",
              f"{current}/{foldername}/test-{test.sim_time}.in")
    print(f"({count}/{total}) Predicted time: {int(((total - count)*(t1-t0))/60)} minutes.", end="\r")
    count+=1
print("Done!")

# test.newfile('test.in')

