import sys
import time
from math import *
import numpy as np
import os
#now set up the parameters to create the json file
Dn = 36.0
Chi = 0.0
Dc = 0.3 * Dn
#now create the json file
command = "./pFieldVisjson.sh " + str(Dn) + " " + str(Chi) + " " + str(Dc)
os.system(command)
#now run the main program
command = "./build/pFieldVis pFieldVis.json > output"
os.system(command)
print("script finished")
