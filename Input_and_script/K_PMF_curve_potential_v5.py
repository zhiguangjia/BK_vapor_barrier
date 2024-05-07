#!//usr/bin/python

import MDAnalysis
import numpy
import sys
import os
import math
import MDAnalysis
import MDAnalysis.analysis.hbonds
from  MDAnalysis import analysis


#  fix a small bug in V4  L156 the selection should in outer circle


structurefile = sys.argv[1]
targettrajectory = sys.argv[2]


cutoff = float(sys.argv[3])

bins = int(sys.argv[4])
binmin = int(sys.argv[5])
binmax = int(sys.argv[6])

outfilename   = sys.argv[7]

try:
    intitial_time = float(sys.argv[8])

    try:
       final_time = float(sys.argv[9])
    except:
       final_time = 9999999.0

except:
    intitial_time = 0.0
    final_time    = 9999999.0



ofile1  = open(str(outfilename + '_K_PMF.dat'),'w') # open file for writing
ofilewater  = open(str(outfilename + '_water_PMF.dat'),'w') # open file for writing
ofile_k_number   = open(str(outfilename + '_K_density.dat'),'w') # open file for writing
ofile_w_number   = open(str(outfilename + '_water_density.dat'),'w') # open file for writing

# unit in A!


startpoint = 0

targetprotein = MDAnalysis.Universe(structurefile, targettrajectory)


#potassium = targetprotein.select_atoms(" ( resname K POT ) and not protein ")
 # the above selection generate a list in  mdanalysis form ,but how get index out ?
#OW        = targetprotein.select_atoms(" ( resname SOL TIP TIP3 WAT and name OW OH2 O ) and not protein ")

  ## put selection later, so we can narrow down only those near protein

#  filter:  T287  OG  O   V288 O G289   O 
#gateresidue = targetprotein.select_atoms(" (resid 287 and name OG1) or (resid 287 and name O) or (resid 288 and name O) or (resid 289 and name O) ")


##   use Ca
gateresidue = targetprotein.select_atoms("  (resid 286 and name CA) or (resid 287 and name CA) or (resid 288 and name CA) ")


#for t in gateresidue:
#   print t


slected_katom_z  = []
slected_OWatom_z = []
# this will store the z coordination of k falls in the x-y cutoff

##   for error analysis , we consider second half of the trajec
Fir_half_slected_katom_z  = []
Sec_half_slected_katom_z  = []
Fir_half_slected_OWatom_z = []
Sec_half_slected_OWatom_z = []


framenumber = float(0)
total_k_select  = 0
total_OW_select = 0


for ts in targetprotein.trajectory:

   current_time = ts.time/1000.0 # in ns

   if current_time < intitial_time :
      pass 

   elif current_time > final_time :
      break 

   #elif ts.frame >= intitial_frame and ts.frame <= final_frame ):
   else:
      framenumber += 1    
    #  print  ts.frame , current_time, framenumber  
      #position_array = potassium.coordinates() 
      #print      potassium.positions[:5]
      #print ts
      x0 = gateresidue.center_of_mass()[0]
      y0 = gateresidue.center_of_mass()[1]
      z0  = gateresidue.center_of_mass()[2]
      if framenumber == 1 or framenumber == 0:
           print  x0 , ' ', y0 ,  '  ',  z0,  ' at frame 0'
           print  'origin S6 position'
           s6_res = targetprotein.select_atoms("  (resid 286 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 286 ', z_temp,-( z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 287 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 287 ', z_temp, -( z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 288 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 288 ', z_temp, -(z_temp-z0)

           s6_res = targetprotein.select_atoms("  (resid 304 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 304 ', z_temp, -( z_temp-z0 )
           s6_res = targetprotein.select_atoms("  (resid 308 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 308 ', z_temp,-( z_temp-z0 )
           s6_res = targetprotein.select_atoms("  (resid 312 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 312 ', z_temp,-( z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 313 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 313 ', z_temp, -( z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 314 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 314 ', z_temp, -( z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 315 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 315 ', z_temp, -( z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 316 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 316 ', z_temp,-( z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 320 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 320 ', z_temp, -(z_temp-z0)
           s6_res = targetprotein.select_atoms("  (resid 324 and name CA)")
           z_temp = s6_res.center_of_mass()[2]
           print  ' 324 ', z_temp,-( z_temp-z0)

       #    s6_res = targetprotein.select_atoms("  (resid 682 and name CA)")
       #    z_temp = s6_res.center_of_mass()[2]
       #    print  '682, close to CTD bottom ', z_temp,-( z_temp-z0)


      potassium = targetprotein.select_atoms(" ( resname K POT ) and around 20  protein ")
      OW        = targetprotein.select_atoms(" ( resname SOL TIP TIP3 WAT and name OW OH2 O ) and around 20  protein ")

      for katom_position     in  potassium.positions :
           x = float(katom_position[0])
           y = float(katom_position[1])
           z = float(katom_position[2])
           #print x 
           if (x-x0)*(x-x0) + (y-y0)*(y-y0) <= cutoff*cutoff  and z >= (z0 + binmin)  and z<= (z0 + binmax) :
                                                                        # e.g.  + (-40)            e.g. + (+80)
                 slected_katom_z.append(z-z0)
                 total_k_select += 1


      for OWatom_position     in  OW.positions :
           x = float(OWatom_position[0])
           y = float(OWatom_position[1])
           z = float(OWatom_position[2])
           #print x 
           if (x-x0)*(x-x0) + (y-y0)*(y-y0) <= cutoff*cutoff  and z >= (z0 + binmin)  and z<= (z0 + binmax) :
                                                                        # e.g.  + (-40)            e.g. + (+80)
                 slected_OWatom_z.append(z-z0)
                 total_OW_select += 1

     # print 'frame ', framenumber
     # print 'totoal OW_selec', total_OW_select

print ' totoal frame : ', framenumber
print ' totoal k_selec : ', total_k_select
print ' totoal OW_selec : ', total_OW_select

hist, bin_edges    = numpy.histogram(slected_katom_z,bins,(binmin, binmax))
histOW,bin_edgesOW = numpy.histogram(slected_OWatom_z,bins,(binmin, binmax))


sum_for_normalize= sum(hist)


temp = bin_edges[:-1]


#t1 = temp[startpoint:]
#t2 = hist[startpoint:]
#time = numpy.array(t1)
#decay = numpy.array(t2)

#  if NOT want skip startpint 
zaxis    = numpy.array(temp)
decay   = numpy.array(hist)
decayOW = numpy.array(histOW)

print len(zaxis), len(decay)



# now change to energy

# kB = 3.2976268E-24  cal/K
kB = 1.3806488E-23    # j/K
An = 6.02214179E23
T = float(300)
R = 1.987  # (cal/mol.degree)
#frame = int(0)
frame = float(0)
#  later it will be divided, so float may be better??


###  bulk solvent K, is in   /home/zgjia/Project/BK_channel/Simulation_test_with_EM_result/Seperate_S0/K_reference_probability

##### Finding the maximum

Pmax = max(decay)/framenumber
#print 'Pmax ', max(decay)
#####

DG = numpy.zeros(bins+2)

DG_water = numpy.zeros(bins+2)


#####


for x in range(bins-1):
        #print x, ' ',  decay[x]
        #print x, ' ',  decayOW[x]
        if decay[x] == 0:
             DG[x] = 99

        else:
             #DG[x] = -0.001*R*T*(math.log( (float(decay[x])/framenumber)/(3.1415926*cutoff*cutoff)  ) )
             #DG[x] = -0.001*R*T*(math.log( (  ( float(decay[x])/float(decayOW[x]) )  /framenumber)/(3.1415926*cutoff*cutoff)  ) )
                                          #    water is used as a background
             DG[x] = -0.001*R*T*(math.log( (  ( float(decay[x]) )  /framenumber)/(3.1415926*cutoff*cutoff)  ) )
             # 0.001 change cal to kcal
 
        if decayOW[x] == 0:
             DG_water[x] = 99  

        else:
             DG_water[x]  = -0.001*R*T*(math.log( (  float(decayOW[x])    /framenumber)/(3.1415926*cutoff*cutoff)  ) )
             # 0.001 change cal to kcal   
   #     print x, ' ',  decayOW[x],DG_water[x]
 
tt = 0
while tt < bins-startpoint:

     #   please note the -1 in zaxis[tt] means CTD orientation correction (our sim is reversed)

     # FE 
     ofile1.write( str(DG[tt]) + "    " + str(-zaxis[tt]) + "    "  + "    "  +"\n")
     ofilewater.write( str(DG_water[tt]) + "    " + str(-zaxis[tt])  + "    "  +"\n")
     # density  
     ofile_k_number.write( str( (  float(decay[tt])  /framenumber )/(3.1415926*cutoff*cutoff)  ) + "    " + str(-zaxis[tt])  + "    "  +"\n")
     ofile_w_number.write( str( (  float(decayOW[tt]) /framenumber )/(3.1415926*cutoff*cutoff)  ) + "    " + str(-zaxis[tt])  + "    "  +"\n")   
     tt =tt+1
ofile1.close()



