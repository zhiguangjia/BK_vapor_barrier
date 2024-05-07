#!//usr/bin/python

import MDAnalysis
import numpy
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
import numpy.linalg
from numpy.linalg import norm


cutoff           = float(sys.argv[1])
structurefile    = sys.argv[2]
targettrajectory = sys.argv[3]

outfilename   = sys.argv[4]

#skip  =  int(sys.argv[4])

try:
    initial_frame = int(sys.argv[5])

    try:
       final_frame = int(sys.argv[6])

       try:
            skip  =  int(sys.argv[7])
       except:
            skip  =  1

    except:
       final_frame = 9999999999
       skip  =  1

except:
    initial_frame = 0
    final_frame = 9999999999
    skip  =  1

print '  '
print ' The idea of this script is   '
print '    1, define local COM canter by backbone of four  residues  e.g.315    (could not include side chain as flip may not symmetry)       '
print '    2, the angle between two vactors :  1  CA  of 315 to pore center;  2  CA of each 315 to COM of ring of PHE . This write to **********_***_orientation.xvg      '
print '                                                                                                                               inputname  315         '
print '  ' 
print ' note  :  '
print '    1, The PHE side chain should have atom name CG CD1 CD2 CE1 CE2 CZ  (true for charmm gmx, amber)  ' 
print '    2, for the angle, we not calculate 3D angle, we cal the project of angle on xy plane ......the 3D one will make unreal underestimate of the flip move ment '
print '                                                                                e.g.     CA----COM       |\ Z                 '
print '                                                                                           \             |    '
print '                                                                                            Sidechain    |    '
print '  '
print '    3 Suppose treating a 50 ps output trajectory, if not , modify line 122  time = ts.frame*50/1000     '
print '    4 for  GLY, we give a 0 angle ~l 204 '
print '  '
print '  '




ofile_water = open(str(outfilename + '_water_312-320_time.xvg'),'w') # open file for writing
ofile_pot = open(str(outfilename + '_pot_312-320_time.xvg'),'w') # open file for writing


# unit in A!


#######################################
#  relations :
#           RCK    VSD
#            A   -> C   
#            C   -> B
#            B   -> D 
#            D   -> A
VSD_chain_list = ['C', 'B', 'D', 'A']
RCK_chain_list = ['A', 'C', 'B', 'D']
  ############  whole calculation will be made ,VSD and RCK just means the contact method between two chain (who give VSD and who give RCK)

##########################

#  constants
# kB = 3.2976268E-24  cal/K
kB = 1.3806488E-23    # j/K
An = 6.02214179E23
T  = float(300)
R  = 1.987  # (cal/mol.degree)
#frame = int(0)
frame = float(0)
#  later it will be divided, so float may be better??

##########################

targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)

print  ' total frame number ',  targetprotein.trajectory.n_frames
#half_frame = int(targetprotein.trajectory.n_frames/2)

##  write file title  , set up average array 


frame_number  =  float(0)
for ts in targetprotein.trajectory:

   # if (ts.frame == 0 or  float(ts.frame)%skip == 0) and  ts.frame > intial_frame and  ts.frame  < last_frame:
   # if initial_frame != 0  and ( ts.frame < initial_frame or ts.frame > final_frame ):
    if ts.frame < initial_frame or float(ts.frame)%skip != 0 :
         pass
    elif ts.frame > final_frame :
         print ts.frame, 'exceed ', final_frame, 'fnishing and output final result' 
         break
    else :
         frame_number  +=  float(1)
         #time = float(ts.frame*50*skip) /1000  # suppose 50 ps step and final need ns unit
         # time = float(ts.frame*50) /1000  

         time = ts.time/1000  #  in ns unit  
         print 'frame', ts.frame, 'time (ns) ', time 

         ###################################  water
 
         gateresidue_CA   = targetprotein.select_atoms(" ( protein  and (resid 312 and name CA) ) ")
         x0 = gateresidue_CA.center_of_mass()[0]
         y0 = gateresidue_CA.center_of_mass()[1]
         current_radii = float(0)
         for temp_position       in  gateresidue_CA.positions:
              temp_x = float(temp_position[0])
              temp_y = float(temp_position[1])
              current_radii += math.sqrt( (temp_x-x0)*(temp_x-x0) + (temp_y-y0)*(temp_y-y0) )
         current_radii = current_radii/4
       
         min_radii        = current_radii            
         minz             = gateresidue_CA.center_of_mass()[2]

         gateresidue_CA   = targetprotein.select_atoms(" ( protein  and (resid 320 and name CA) ) ")
         x0 = gateresidue_CA.center_of_mass()[0]
         y0 = gateresidue_CA.center_of_mass()[1]
         current_radii = float(0)
         for temp_position       in  gateresidue_CA.positions:
              temp_x = float(temp_position[0])
              temp_y = float(temp_position[1])
              current_radii += math.sqrt( (temp_x-x0)*(temp_x-x0) + (temp_y-y0)*(temp_y-y0) )
         current_radii = current_radii/4

         max_radii        = current_radii
         maxz             = gateresidue_CA.center_of_mass()[2]

         # volumn = 3.1415926*(maxz-minz)*(min_radii*min_radii + max_radii*max_radii + max_radii*min_radii)/3
                #  V = pi*h*(R^2+r^2+R*r )/3
         #  switch to cut off
         volumn = 3.1415926*(maxz-minz)*(cutoff*cutoff)
                

         gateresidue = targetprotein.select_atoms("  (resid 286 and name CA) or (resid 287 and name CA) or (resid 288 and name CA) ")
         x0 = gateresidue.center_of_mass()[0]
         y0 = gateresidue.center_of_mass()[1]

         OWnear = targetprotein.select_atoms(" ( resname SOL TIP TIP3 WAT and name OW OH2 O ) and around 12.5 ( protein and resid 312:320 ) ")
                                                                                                  ### 25 is maxium pore size, should be enough
         current_OW_select= float(0)
         for OWatom_position     in  OWnear.positions :
              x = float(OWatom_position[0])
              y = float(OWatom_position[1])
              z = float(OWatom_position[2])
             #print x 
              if (x-x0)*(x-x0) + (y-y0)*(y-y0) <= cutoff*cutoff  and z >= minz  and z<= maxz :
                   #print x, y, z
                   current_OW_select += 1

         current_OW_density  =  current_OW_select/volumn
         #   reference water density :  /home/zgjia/Gromacs_test/TIP3P_box
                                      #  0.45 mM  :  7045 water 5.95086^3  nm box  = 0.03343
 
         #ofile_water.write( str(ts.time/1000) + "    "  + str((100*current_OW_density)/0.03343) + "    " + "\n")
                          #  ts time ps unit                                       
         print time, current_OW_select
         ofile_water.write( str(time) + "    "  + str(current_OW_select) + "    " + "\n")
                                                   # use absolute number ; for density, use current_OW_density

         current_POT_select= float(0)
         POTnear = targetprotein.select_atoms(" ( resname POT K ) and around 12.5 ( protein and resid 312:320 ) ")
         for POTatom_position  in  POTnear.positions :
              x = float(POTatom_position[0])
              y = float(POTatom_position[1])
              z = float(POTatom_position[2]) 
              if (x-x0)*(x-x0) + (y-y0)*(y-y0) <= cutoff*cutoff  and z >= minz  and z<= maxz :
                   #print x, y, z
                   current_POT_select += 1

         ofile_pot.write( str(time) + "    "  + str(current_POT_select) + "    " + "\n")




