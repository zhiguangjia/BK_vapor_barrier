#!//usr/bin/python
import MDAnalysis
import numpy
import numpy as np
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
from collections import OrderedDict
from pprint import pprint
from numpy.linalg import norm



structurefile = sys.argv[1]
targettrajectory = sys.argv[2]

outfilename   = sys.argv[3]

residue_begin_1  = int(sys.argv[4])

residue_end_1  = int(sys.argv[5])


# discuss also see /home/zgjia/Project/TMEM16/TMEM16F/Paulino_6qp/Manuscript_2/Figures/TM_2D_plot/commond_tilt_and_rotate


try:
#    timestep     = float(sys.argv[4])
    initial_time  = int(sys.argv[6])
    final_time  = int(sys.argv[7])

except:
#    timestep     = 0.50   # suppose anton input  unit  ns
    initial_time  = 0
    final_time    = 999999999999


print  " "
print  " "
print  " "


ofile_tilt  = open(str(outfilename + '_tilt.xvg'),'w') # open file for writing
ofile_tilt_hist  = open(str(outfilename + '_tilt_hist.xvg'),'w') # open file for writing
#ofile_helix_rotation  = open(str(outfilename + '_helix_rotation.xvg'),'w') # open file for writing
#ofile_helix_rotation  = open(str(outfilename + '_helix_rotation_hist.xvg'),'w') # open file for writing
ofile_COM_rotation       = open(str(outfilename + '_COM_rotation.xvg'),'w') # open file for writing
ofile_COM_rotation_hist  = open(str(outfilename + '_COM_rotation_hist.xvg'),'w') # open file for writing
   #  based on the  z vector from filter center, is there a clocwise or anticlocwise rotation of helix ??
ofile_COM_move_xy           = open(str(outfilename + '_COM_move_xy.xvg'),'w') 
ofile_COM_move_xy_hist        = open(str(outfilename + '_COM_move_xy_hist.xvg'),'w') 
   #  based on the filter center, is there a move helix on xy  ?? ()
ofile_COM_move_z           = open(str(outfilename + '_COM_move_z.xvg'),'w')
ofile_COM_move_z_hist        = open(str(outfilename + '_COM_move_z_hist.xvg'),'w')

def rotation_angle_around_z(ini_vec, final_vec):

## clockwise is -
## counterclockwise is +

   ax = ini_vec[0]
   ay = ini_vec[1]
   bx = final_vec[0]
   by = final_vec[1]

   angle = numpy.rad2deg( math.atan2( ax*by-ay*bx, ax*bx+ay*by )  )

   return angle

def vecter_angle(A,B):
    theta = numpy.arccos(numpy.dot(A, B)/(norm(A)*norm(B)))
    current_angle = numpy.rad2deg(theta)

    ########################    special ,, our vector, as helix orientation issue , may up or but it will always less then 90 !!!!

    if current_angle > 90.0:

        current_angle = 180 - current_angle


    return current_angle


def point_projection(a, b, p):
    ##  calcualte projection of the pore_CB on helix_axis
    #           Cbeta (P)
    #             |
    #             | a2    ** a is  vector from helix COM to CB 
    #   .    .  _ | ____
    #   A    B    C   
    #           h helix_axis  ---->

    ap = p - a
    ab = b - a
    scaler = np.dot(ap, ab) / np.dot(ab, ab)
             #  |AP|*|AB|*cos(theta)/|AB|^2  then scale
    result = a + scaler * ab

    return result, numpy.sign(scaler)


######################################     reference  part   <==  1st frame 

reference     = MDAnalysis.Universe(structurefile, structurefile)

reference_helix_COM_shifted  = OrderedDict()
reference_helix_COM_distance = OrderedDict()
reference_helix_COM_distance_z = OrderedDict() 

RCKchain = ['B','C','D','A']   ## may be wrong ....
VSDchain = ['A','B','C','D']  

for ts in reference.trajectory:
    # should be single frame 

    filter_all = reference.select_atoms(" ( (resid 287:292 ) and name CA) ")
    reference_filter_COM = filter_all.center_of_mass()
    # this is used to translocate to origin point 0 0 0 

    for currect_VSDchain in VSDchain  :       
       helix_current = reference.select_atoms(" ( name CA and resid %d:%d ) and segid %s  "%(residue_begin_1, residue_end_1, currect_VSDchain) )

       # for toation and move, we only consider xy 
       reference_helix_COM_shifted[currect_VSDchain]  = helix_current.center_of_mass()[0:2] - reference_filter_COM[0:2]

       reference_helix_COM_distance[currect_VSDchain] = numpy.linalg.norm(reference_helix_COM_shifted[currect_VSDchain])

       reference_helix_COM_distance_z[currect_VSDchain] =   abs(helix_current.center_of_mass()[2] - reference_filter_COM[2]   )


#################################   reference  end 

targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)

totalframe = len(targetprotein.trajectory)

#filter_horizon_move  = open(str(outfilename + '_filter_2D.xvg'),'w') # open file for writing
#S1_horizon_move  = open(str(outfilename + '_S1_2D.xvg'),'w') # open file for writing

angle_tilt_list = []
angle_rotation_list = []
heix_move_list = []
heix_move_list_z = []
# initalize value

for ts in targetprotein.trajectory:

  # print ts.time/1000
 
   if initial_time != 0  and ( ts.time/1000 <= initial_time  ):
       pass
   elif ts.time/1000 >= final_time :

       break

   else :
      #(ts.frame == 0 or  float(ts.frame)%skip == 0 ) and (  ts.frame >= initial_time and ts.frame <= final_time ) :

      filter_all =targetprotein.select_atoms(" ( (resid 287:292 ) and name CA) ")
      filter_COM = filter_all.center_of_mass()
       # this is used to translocate to origin point 0 0 0 

      ##  a temperary list, later we first write average, then 4 chain seperate file
      temp_angle_tilt_list = []
      temp_angle_rotation_list = []
      temp_heix_move_list = []
      temp_heix_move_list_z = []
 
      BA = [0.0, 0.0, 1.0]

      for currect_VSDchain in VSDchain  :
         helix_current = targetprotein.select_atoms(" ( name CA and resid %d:%d ) and segid %s  "%(residue_begin_1, residue_end_1, currect_VSDchain) )
         helix_COM_shifted = helix_current.center_of_mass()- filter_COM

         #########################
         #  first, tilting angle 

         BC = helix_current.principal_axes()[2]   # see  /home/zgjia/test/test_mdanalysis_tilt_and_principal_axis
         current_angle =  vecter_angle(BA, BC)

         angle_tilt_list.append(float(current_angle))
         temp_angle_tilt_list.append(float(current_angle))

         #########################
         #   rotation around  cenrter axis clockwise (-) or anticlockwise (+)    ** reference is inital structure 

         # for toation and move, we only consider xy 

         helix_COM_shifted = helix_current.center_of_mass()[0:2] - filter_COM[0:2]
                             # make a vector start from 0 0  this also used in next xy move calculation
        # print helix_current.center_of_mass(),filter_COM,helix_current.center_of_mass()[0:2],filter_COM[0:2], helix_COM_shifted 
         current_angle = rotation_angle_around_z(reference_helix_COM_shifted[currect_VSDchain], helix_COM_shifted )

         angle_rotation_list.append(float(current_angle)) 
         temp_angle_rotation_list.append(float(current_angle))
        
         #########################
         #  move on xy ,  colser is -  away is +

         helix_COM_distance_changed =   numpy.linalg.norm(helix_COM_shifted) - reference_helix_COM_distance[currect_VSDchain]

         heix_move_list.append(float(helix_COM_distance_changed)) 
         temp_heix_move_list.append(float(helix_COM_distance_changed))

         #  move on z 

         helix_COM_distance_changed_z = abs( helix_current.center_of_mass()[2] - filter_COM[2] )  -reference_helix_COM_distance_z[currect_VSDchain]
         heix_move_list_z.append(float(helix_COM_distance_changed_z))
         temp_heix_move_list_z.append(float(helix_COM_distance_changed_z))
         
      ####################  

      ofile_tilt.write(         str(ts.time/1000) + "  " + str( numpy.average(temp_angle_tilt_list) )      + "  "  )
                                                              #  write average first
      ofile_COM_rotation.write( str(ts.time/1000) + "  " + str( numpy.average(temp_angle_rotation_list) )   + "  "  )
      ofile_COM_move_xy.write(     str(ts.time/1000) + "  " + str( numpy.average(temp_heix_move_list) )   + "  "  )
      ofile_COM_move_z.write(     str(ts.time/1000) + "  " + str( numpy.average(temp_heix_move_list_z) )   + "  "  )      

      ofile_tilt.write( "  " + str(temp_angle_tilt_list[0]) + "  " + str(temp_angle_tilt_list[1]) + "  " + str(temp_angle_tilt_list[2]) +"  " + str(temp_angle_tilt_list[3])  + "\n")
      ofile_COM_rotation.write( "  " + str(temp_angle_rotation_list[0]) + "  " + str(temp_angle_rotation_list[1]) + "  " + str(temp_angle_rotation_list[2]) + "  " + str(temp_angle_rotation_list[3])  + "\n")
      ofile_COM_move_xy.write( "  " + str(temp_heix_move_list[0]) + "  " + str(temp_heix_move_list[1]) + "  " + str(temp_heix_move_list[2]) + "  " + str(temp_heix_move_list[3]) + "\n")
      ofile_COM_move_z.write( "  " + str(temp_heix_move_list_z[0]) + "  " + str(temp_heix_move_list_z[1]) + "  " + str(temp_heix_move_list_z[2]) + "  " + str(temp_heix_move_list_z[3]) + "\n") 


      #  now  loop four  pair RCK-VSD    

      if float(ts.frame)%1000 == 0 :
          print  ' done frame ', ts.time/1000,'ns'


print 'average  tilt ', 'std' , numpy.average(angle_tilt_list), numpy.std(angle_tilt_list)
print 'average rotation around  cenrter axis clockwise (-) or anticlockwise (+)', 'std' , numpy.average(angle_rotation_list), numpy.std(angle_rotation_list)
print 'average  move on xy ,  colser is -  away is + ', 'std' , numpy.average(heix_move_list), numpy.std(heix_move_list)
print 'average  move on z ,  icrease is -  decrease is + ', 'std' , numpy.average(heix_move_list_z), numpy.std(heix_move_list_z)

quit()

binmin = -180
binmax = 180

bins = int((binmax-binmin)/0.5)

hist, bin_edges    = numpy.histogram(angle_list,bins,(binmin, binmax))
hist_sum = float( numpy.sum(hist) )# for  normalization

for x in range(bins-1):

   ofile_tilt_hist.write( str(bin_edges[x])   + " "  + str( float(hist[x]) /hist_sum )   +"\n")






























