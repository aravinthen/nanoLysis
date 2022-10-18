# Program Name: test.py
# Author: Aravinthen Rajkumar
# Description: testing file

import sys
import os
import time
sys.path.insert(0, '../main/')

from main import nanolysis

test = nanolysis()
test.set_origin("../data/structure-120.0.in")

current = os.getcwd()
test.read("../data/simulation/dump.settings-120.0.in.in_100000.in")

# test.bonds(30)
# test.hostwalk(30)
# test.hosttype(1)
print(test.get_data(10, 1))

current = os.getcwd()
path = current + "/../data/simulation"
