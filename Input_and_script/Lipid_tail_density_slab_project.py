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
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
from  MDAnalysis import core
from  MDAnalysis.core import groups
from  MDAnalysis.core import flags

from collections import OrderedDict


try:

        structurefile = sys.argv[1]

        infilename1   = sys.argv[2]

        outfilename = sys.argv[3]

        thickness   = float(sys.argv[4])


        steps       = float(sys.argv[5])
        init_z      = float(sys.argv[6])
        end_z       = float(sys.argv[7])    #  this is cut in z axis, each step generate a file

        bins = int(sys.argv[8])
        binmin = int(sys.argv[9])
        binmax = int(sys.argv[10])          #  this is bins for histograme on the slab corrdinate 

        current_residue = int(sys.argv[11])  # the residue we cut a line


    #    start_value  = float (sys.argv[5])
    #    ns_per_step  = float (sys.argv[6])
        
except:
     print "Usage:",sys.argv[0], "infile  ns per step  outfile"; sys.exit(1)



#########################################

#ifile1 = open(infilename1,'r') # open file for reading


#total_lipid_select  = []
#lipid_frame_tracking = OrderedDict()
#lipid_time_record    = OrderedDict()


#                    cutting line
#                A  /
#                  /
#          D     F    B
#                  
#                C   
#
#

reference     = MDAnalysis.Universe(structurefile, structurefile)

VSDchain_1 = 'C'
VSDchain_2 = 'A'

anchor_residue = [288,current_residue,315,325]  # this is for later we need nome reference anchor to fit vmd to eps

anchor_1_x = OrderedDict()
anchor_1_y = OrderedDict()
anchor_1_z = OrderedDict()

anchor_2_x = OrderedDict()
anchor_2_y = OrderedDict()
anchor_2_z = OrderedDict()


for ts in reference.trajectory:
        # should be single frame 

    helix =reference.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(int(current_residue), VSDchain_1) )

    helix_ref_1_x = helix.center_of_mass()[0]
    helix_ref_1_y = helix.center_of_mass()[1]
    helix_ref_1_z = helix.center_of_mass()[2]

    helix =reference.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(int(current_residue), VSDchain_2) )
    
    helix_ref_2_x = helix.center_of_mass()[0]
    helix_ref_2_y = helix.center_of_mass()[1]
    helix_ref_2_z = helix.center_of_mass()[2]


    filter_gate =reference.select_atoms("  (resid 286 and name CA) or (resid 287 and name CA) or (resid 288 and name CA) ")
    filter_x = filter_gate.center_of_mass()[0]
    filter_y = filter_gate.center_of_mass()[1]
    filter_z = filter_gate.center_of_mass()[2]
     
    for temp_res in anchor_residue:

        temp_res_select =reference.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(int(temp_res), VSDchain_1) )

        anchor_1_x[temp_res] = temp_res_select.center_of_mass()[0]
        anchor_1_y[temp_res] = temp_res_select.center_of_mass()[1]
        anchor_1_z[temp_res] = temp_res_select.center_of_mass()[2]

        temp_res_select =reference.select_atoms(" (resid %d and name CA Ca ) and segid %s "%(int(temp_res), VSDchain_2) )

        anchor_2_x[temp_res] = temp_res_select.center_of_mass()[0]
        anchor_2_y[temp_res] = temp_res_select.center_of_mass()[1]
        anchor_2_z[temp_res] = temp_res_select.center_of_mass()[2]


##########################################
#
#  now we set up the slab and projection
#
####################   now, detrmine  the cutting slice line function
#
#   1  in current setting, filter should be center (  well, we will re center any way)
#
#    two point:  specific residue and filter    if y = ax + b,  
#
#     a  *  helix_ref_x    +  b  = helix_ref_y  ;  a* filter_average_x +b = filter_average_y
#
#     a= (filter_average_y - helix_ref_y)   / (filter_average_x- helix_ref_x)  ;  b =  ( helix_ref_y*filter_average_x - helix_ref_x*filter_average_y  ) / (filter_average_x- helix_ref_x)
#
#
#     tan (alpha) = (filter_average_y - helix_ref_y)   / (filter_average_x- helix_ref_x)
# 
#
#
#
#  2 No arbitry unit , direct prjection
#
#
#
#


#  y = slope_a * x + slope_b
#   pass filter_average and helix_ref

slope_a = (helix_ref_1_y - helix_ref_2_y)   / (helix_ref_1_x - helix_ref_2_x)
#slope_b = ( helix_ref_y*filter_average_x - helix_ref_x*filter_average_y  ) / (filter_average_x- helix_ref_x)  
slope_b = helix_ref_1_y - slope_a*helix_ref_1_x

#################################   now, we can decide the position of all anchor point in  artificial-slop  - z  plane  

def point_on_line(a, b, p):
    ap = p - a
    ab = b - a
    result = a + np.dot(ap, ab) / np.dot(ab, ab) * ab
    return result  # ## here we get is the xy of the projection point, we convert it to distance to center

####  we need this as the residue we check could not on the cutting plane


filter_project = point_on_line(np.array([helix_ref_1_x,helix_ref_1_y]), np.array([helix_ref_2_x,helix_ref_2_y]), np.array([filter_x,filter_y]))

for temp_res in anchor_residue:

   # temp_x =  - (filter_average_x -   anchor_residue_average_x[temp_res]  ) / numpy.cos(alpha)

    print '# current residue in check: ', temp_res , 'chain', VSDchain_1
    temp_project = point_on_line(np.array([helix_ref_1_x,helix_ref_1_y]), np.array([helix_ref_2_x,helix_ref_2_y]), np.array([anchor_1_x[temp_res],anchor_1_y[temp_res]]))

    #  now the problem is the distance and sign
    sign = abs(temp_project[0]-filter_project[0])/(temp_project[0]-filter_project[0])
    
    dist = sign * np.linalg.norm([ temp_project[0]-filter_project[0] , temp_project[1]-filter_project[1] ])

    print dist, -(anchor_1_z[temp_res] - filter_z )   # flip z axis

  
    print '# current residue in check: ', temp_res , 'chain', VSDchain_2
    temp_project = point_on_line(np.array([helix_ref_1_x,helix_ref_1_y]), np.array([helix_ref_2_x,helix_ref_2_y]), np.array([anchor_2_x[temp_res],anchor_2_y[temp_res]]))


    sign = abs(temp_project[0]-filter_project[0])/(temp_project[0]-filter_project[0])

    dist = sign * np.linalg.norm([ temp_project[0]-filter_project[0] , temp_project[1]-filter_project[1] ])

    print dist, temp_project-filter_project, -(anchor_2_z[temp_res] - filter_z)

  #  print '# for gnuplot anchor point '
  #  temp_gnuplot = str('set label "o" at ' + str(temp_x)[:6] + ',' + str(-anchor_residue_average_z[temp_res])[:6] + ' tc rgb "black" font ",30" front')
  #  print temp_gnuplot


###########
#         
#    current issue is   the data we read already centered on filter and z correct so   
#

current_z_ini = init_z
current_z_end = init_z - steps

while current_z_end >= end_z :

    ifile1 = open(infilename1,'r') # open file for reading

    ofile = open(str(outfilename + '_' + str(current_z_ini) + '_to_' + str(current_z_end) +  '.dat'),'w') # open file for writing

    projected_coordinate = []

    frame = 0
    for line in ifile1:

       if float(frame)%50000 == 0 :
          print  ' done data ', frame

       frame += 1

       columns = line.split()

       if len(columns) >= 3 and (columns[0][0] != '#' and  columns[0][0] != '@'   ):

          x = float(columns[0])
          y = float(columns[1])
          z = float(columns[2])

          # now project remember we already centered all the coordinate to filter and correct z 
  
          # !! also important is we select a large circle of lipid, now need use thickness to cared a rectangel centered at the slab

          # !!  and now the line cross the filter center so     #  y = slope_a * x      


        ## remember x y is  centerd and z has corrected direction

          if z <= current_z_ini and z >= current_z_end :
  
              # now find the projection of the lipid tail on the center

             temp_project = point_on_line(np.array([0.0,0.0]), np.array([1.0,slope_a]), np.array([x,y]))

             dist_to_slope = np.linalg.norm([ temp_project[0]-x , temp_project[1]-y ])
             
             if dist_to_slope <= thickness:

                 sign = abs(temp_project[0])/(temp_project[0])
 
                 dist = sign * np.linalg.norm([ temp_project[0] , temp_project[1] ])

             
                 projected_coordinate.append(dist)
 


    hist, bin_edges    = np.histogram(projected_coordinate,bins,(binmin, binmax))


    print len(hist), len(bin_edges)

    counter = 0
    sum_lipid = 0.0
    Max_lipid = 0.0

    for hist_temp in hist:
        sum_lipid += hist_temp  

        if hist_temp > Max_lipid:
           Max_lipid = hist_temp
        counter += 1

    counter = 0

    if sum_lipid == 0.0 :
           
       print 'between', current_z_ini , current_z_end , 'no lipid faound'       

       
    for hist_temp in hist:

        # hist_temp/sum  will be normalized data, but for plot, we scale it so later every thing can plat together 

        if Max_lipid != 0.0 :

            z_for_plot = current_z_end + (float(hist_temp)/float(Max_lipid))*steps*0.8  # *0.8 make it easier to plot 

            ofile.write(str(bin_edges[counter]) + "  " + str(z_for_plot )  +"\n")

        else:
            ofile.write(str(bin_edges[counter]) + "  " + str(current_z_end )  +"\n") 

        counter += 1

    ofile.close
    ifile1.close

    current_z_ini = current_z_ini - steps
    current_z_end = current_z_end - steps


#print 'average and std', np.average(lifetime_sum), np.std(lifetime_sum)


#       ofile.write(str(columns[0]) + "  " + str(np.average(value_1) ) + "  " + str(np.std(value_1) )  +"\n")





#for line in ifile2:
#    columns = line.split()
#    if (columns[0][0] != '#' and  columns[0][0] != '@' ):
#        value_2.append( float(columns[1]) )



#frame = 0
#for temp in value_1 :
#    ofile.write(str((start_value+frame)*ns_per_step) + "\t" + str( (value_1[frame] + value_2[frame] + value_3[frame]  )/3) +"\n")
#    frame += 1








         
