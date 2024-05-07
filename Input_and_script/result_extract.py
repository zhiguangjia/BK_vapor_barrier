#!/usr/bin/python
import sys
import os
import os.path
import commands
import time
import math
import copy
import pprint
import collections
import numpy as np

try:
        infilename   = sys.argv[1]
        outfilename  = sys.argv[2]
        column_r     = int(sys.argv[3])
except:
     print "Usage:",sys.argv[0], "infile  ns per step  outfile"; sys.exit(1)
#########################################

ifile = open(infilename,'r') # open file for reading


print ' extracting column:  ' ,column_r 

ofile_1 = open(str(outfilename ),'w') # open file for writing


for line in ifile:
    columns = line.split()
    if (columns[0][0] != '#' and  columns[0][0] != '@' ):
        if len(columns) >= 5 :
            ofile_1.write( columns[column_r]  + " " + columns[0]  +"\n")












         
