#!/usr/bin/python
import MDAnalysis
import numpy
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis
from  MDAnalysis import core
from  MDAnalysis.core import groups
from  MDAnalysis.core import flags
import mdtraj as md
from collections import OrderedDict

structurefile    = sys.argv[1]

dxfile       = sys.argv[2]
ofile_name     = sys.argv[3]

e_file     = -float(sys.argv[4])   ##  ! correct here, we accumulate from low to up, see below

# 
#          cell/correct sim             our sim (BK)     1  because, z and E is reverse, when accumulate (from Z0 to Zmax), need  * -1
#       0                                                2  in situation 1, E  ????????????? 
#   Z   E                       Z   500 
#   |\  |\ ---------            |\  |  +++++++++
#   |   |  =========            |   |  =========  
#   |   |  =========            |   |  =========  
#   |   |  +++++++++            |   |/ ---------
#       500                         E  
#                                   0 
#                                                                
#


box_size_x      = float(sys.argv[5])
box_size_y      = float(sys.argv[6])
box_size_z      = float(sys.argv[7])

ofile         = open(str(ofile_name + '_converted.dx'),'w') # open file for writing
ofile_shifted = open(str(ofile_name + '_converted_shifted_by_e-fild.dx'),'w') # open file for writing

ofile_profile = open(ofile_name + '_profile.xvg','w')

try:

   range_sample = float(sys.argv[8])

except:

   range_sample = float(999)  

#########   Anton note:

# force.e_bias.E_applied       std::vector<double>  [0, 0, 0]
#    3-vector describing an applied electric field in units of kcal/mol/A/e
#
#   2   we need 500 mv, which  is 500 mv / 100A  =  5mv/A
#
#      1  kcal/mol  =   0.043 eV
#
#      1  kcal/mol/A/e   =  0.043 eV /A/e   = 0.043 V/A      = 43mV/A
#
#      our box ~ 100 A  so
#
#      x * 43 mv/A * 100A =  require voltage
#
#      500 mv / (43*100)   =  0.11627906976744186



#  output  PME potential (kT/e, T=300K)
# 
#  as Broux said :  V_total= V_reaction+ z*V_applied/L,
#


######################################################     reference point
#  

reference     = MDAnalysis.Universe(structurefile, structurefile)

for ts in reference.trajectory:
    # should be single frame 

 #   filter_gate =reference.select_atoms(" ( (resid 287:292 ) and name CA) ") 
 #   filter_gate =reference.select_atoms("  (resid 286 and name CA) or (resid 287 and name CA) or (resid 288 and name CA) ")

    filter_gate =reference.select_atoms(" ( resid 286 287 288 ) and  backbone  ")
    filter_gate_ref_x = filter_gate.center_of_mass()[0]
    filter_gate_ref_y = filter_gate.center_of_mass()[1] 
    filter_gate_ref_z = filter_gate.center_of_mass()[2]  

    membrane = reference.select_atoms(" ( resname OPC POPC OPE POPE OPS POPS POP) ")
    membrane_ref_z = membrane.center_of_mass()[2]

    print('filter position: ', filter_gate_ref_x,filter_gate_ref_y,filter_gate_ref_z)
    print('membrnae center  ', membrane_ref_z )

#####################################################   grid part 


total_data = 0
line_num   = 1
Nx_current = 1
Ny_current = 1
Nz_current = 1     #  tracker for profile

N_dx_shift = 1   # dx switch line every three data


e_fild_acummulated    = float(0)    

z_pozition         = OrderedDict()
potetial_e_filed   = OrderedDict()
potetial_pmepot_raw  = OrderedDict()   #  !  tthe content is array, average in end
potetial_pmepot    = OrderedDict()     #  !  the content is array, average in end
potetial_total     = OrderedDict()
temp_data_array    = []    # initialized every z slice
temp_potential_array = []    # initialized every z slice

for line in open(dxfile, 'rU'):
    columns = line.split()
    #print line_num
    if line_num == 1 :
        ofile.write(line)
        ofile_shifted.write(line)
        # line 1 always note
    elif line_num == 2 :
        Nx=int(columns[5])
        Ny=int(columns[6]) 
        Nz=int(columns[7]) 
        ofile.write(line)
        ofile_shifted.write(line)
        # 2nd line example 'object 1 class gridpositions counts 216 216 120'
    elif line_num == 3 :       
        # 3rd line example :   origin -21.5049 -20.0866 -10.3917
        x0 = float(columns[1]) 
        y0 = float(columns[2])
        z0 = float(columns[3])
        ofile.write(line) 
        ofile_shifted.write(line) 
        ## initalize the grid for later loop
        x  = x0
        y  = y0
        z  = z0 


    elif line_num == 4 :
        # 4th line example :  delta 0.930139 0 0
        dx = float(columns[1])
        ofile.write(line)
        ofile_shifted.write(line) 
    elif line_num == 5 :
        dy = float(columns[2])
        ofile.write(line)
        ofile_shifted.write(line)
    elif line_num == 6 :
        dz = float(columns[3])
        ofile.write(line)
        ofile_shifted.write(line)
    elif line_num == 7 :
        ofile.write(line)
        ofile_shifted.write(line)
    elif line_num == 8 : 
        total_data_expect = int(columns[9])
        ofile.write(line)
        ofile_shifted.write(line)
        # example object 3 class array type double rank 0 items 5598720 data follows
    elif line_num > 8 :

        ###  IMPORTANT VMD  dx output larger than the bnox, we need only average inbox region, also remove ? A from  edge (may not necessary, for safe)
        #      
        #    
        #    
        #                    
        #    
        #            ***************  inital test 
        #      
        #                  real box_size_x
        #              _____________________________________________________
        #        .               *                               . 
        #       x0            center = (x0 + Nx*dx/2)    edge of dx, x0 + Nx*dx     
        #



        #    now   use filter position as center
        #
        #

       #   2nd, use a small center region
    #    x_range_min = x0 + Nx*dx/2  - range_sample
    #    y_range_min = y0 + Ny*dy/2  - range_sample
        x_range_min = filter_gate_ref_x - range_sample
        y_range_min = filter_gate_ref_y - range_sample
        z_range_min = z0  - dz
    #    x_range_max = x0 + Nx*dx/2  + range_sample
    #    y_range_max = y0 + Ny*dy/2  + range_sample
        x_range_max = filter_gate_ref_x + range_sample
        y_range_max = filter_gate_ref_y + range_sample
        z_range_max = z0 + Nz*dz  + dz

                  #    we not really use  x_range_min/max x_range_min/max     <== this lead to ratangle region
                  #    later directlye radius as judgement



      #  print x0, Nx, dx, (x0 + Nx*dx/2)
      #  print x_range_min,y_range_min,x_range_max,y_range_max
      #  quit() 
        #
        #      end           

        if columns[0] == 'attribute'  or columns[0] == 'object' or columns[0] == 'component' :
            # example:  attribute "dep" string "positions"
            ofile.write("\n")
            ofile.write(line)
            ofile_shifted.write("\n")
            ofile_shifted.write(line)
            # this is end note part, not change    

        else: 
            for data in columns:
                total_data += 1

                ############################################################################################################
                #          
                #           voltage stage (electron included dx  |  data input     )  
                #
                ############################################################################################################

                e_fild_acummulated = e_file * (z-z0)   # unit  mv/A  * A     

                if Nz_current not in z_pozition:

                     z_pozition[Nz_current]          = z
                     potetial_e_filed[Nz_current]    = e_fild_acummulated
                     potetial_pmepot_raw[Nz_current] = []     #  !  the content is array, average in end
                     potetial_pmepot[Nz_current]     = []     #  !  the content is array, average in end


               ####################   grid out put specific, which require 3*? matrix

             #   dx_data = float(data)  # debug
                dx_data = float(data) * 25.7   #   + e_fild_acummulated
                                               ###  see   discuss in  /home/zgjia/test/test_gmx2019_electron_field_membrane/POPC_500mv
                                                                             ###  kbt/e = 25.7 mv     

                dx_data_shifted = dx_data + e_fild_acummulated

                dx_data =  "{:.5f}".format(dx_data)
                dx_data_shifted  =  "{:.5f}".format(dx_data_shifted)

                if N_dx_shift == 3 :  # dx switch line every three data
                    N_dx_shift =1
                    ofile.write( str(dx_data)  + " "  )
                    ofile.write( "\n" )
                    ofile_shifted.write( str(dx_data_shifted)  + " "  )
                    ofile_shifted.write("\n" )
                else:
                    N_dx_shift += 1
                    ofile.write( str(dx_data)  + " "  )
                    ofile_shifted.write( str(dx_data_shifted)  + " "  )

               ####################  profile data input 

              #  if ( x > x_range_min and x < x_range_max ) and ( y > y_range_min and y < y_range_max ) :  #   and (z  > z_range_min and z < z_range_max) :
              #                                                                                            #  z range is not so important  
                 ##  above lead to ratangle region
                 #    directlye radius as judgement
                if ( (x - filter_gate_ref_x)*(x - filter_gate_ref_x) + (y - filter_gate_ref_y)*(y - filter_gate_ref_y) ) < (range_sample*range_sample) : 

                     potetial_pmepot_raw[Nz_current].append(float(data))
                     potetial_pmepot[Nz_current].append(float(data)* 25.7)

                ####  continue on grid , update number 
                
                Nz_current += 1
                z          += dz
                    
               # print ' check z ',  z, ' NZ ', Nz_current,  ' e_fild_acummulated ', e_fild_acummulated

                if Nz_current == Nz +1 :   #####  
                    ##  re initlize:
                    Nz_current  = 1
                    z           = z0

                    ## increase y every time Z loop finished 
                    Ny_current += 1
                    y          += dy
                #    quick()

                if Ny_current == Ny +1 :

                 #  https://www.ics.uci.edu/~dock/manuals/apbs/html/user-guide/x2674.html
                 #   he data values, ordered with the z-index increasing most quickly, followed by the y-index, and then the x-index. 

                 # every time, y reach Ny +1 , Nx should add 1 , begin a new slice


                 #   although Z should already back to 1, we double initilize here, should no hurt

                   ##  re initlize:
                   Nx_current += 1
                   Ny_current = 1
                   Nz_current = 1

                   x += dx
                   y  = y0
                   z  = z0
    line_num += 1


print 'total read data: ', total_data, 'expect ', total_data_expect


                ############################################################################################################
                #          
                #        Profile stage   (profile)
                #
                ############################################################################################################


ofile_profile.write ( "# z_pozition             potetial_total                 potetial_pmepot                 potetial_e_filed (pure e feild )               energy (raw average)  "  + "\n" )

for Nz_current in z_pozition : 

    raw_average = numpy.mean(potetial_pmepot_raw[Nz_current])
    pot_average = numpy.mean(potetial_pmepot[Nz_current])
    
    total       = pot_average + potetial_e_filed[Nz_current]     

    print("current NZ ", Nz_current, " Z ", z_pozition[Nz_current]," data number", len(potetial_pmepot_raw[Nz_current])  ) 
    ofile_profile.write(str(   -(z_pozition[Nz_current] - membrane_ref_z )  ) + "               " + str(total) + "                " + str(pot_average) + "               " + str(potetial_e_filed[Nz_current]) +  "               "   + str(raw_average)   +"\n")
#    print('current z slice', Nz_current, 'at (A)', z_pozition[Nz_current] , ' average ',  potetial_pmepot_raw[Nz_current], ' to potential ', potetial_pmepot[Nz_current]  )

    #  *************   z first shifted with respect to membrane center, then reverse on z axis ....


