
######################################################################################################

 note

 1) as our system is upside down, so most code haa an extra '-' sign hardcoded to correct the direction
 







####### pore  water number 

 theis calcualte pore water number in fig 1c

  ../Analysis_script/312-320_vapor_cutoff.py  5 Eqilibrated_system_with_chain.pdb  sim_1ns.xtc water_res-orientation/cutoff_5


  #  example usage:

  #  this script read a struc ture file Eqilibrated_system_with_chain.pdb ,  an trajectory file sim_1ns.xtc ,

  #   then it need parameter to define the pore region ( a cylinder), the axis is set at the center of filter entre (hard coded), 5 is pore radius, the center is located at the entrance pof the pore which sit on the middle on the membrane (hard coded)

  # number will be output in outp[ut name water_res-orientation/cutoff_5



#######  water density (fig 1d)

  the density profile in figure 1d is using script  , 
  example usage:

   ./K_PMF_curve_potential_v5.py Eqilibrated_system_with_chain.pdb Protein_0-1000ns.xtc \
  5     140     -120    120 PMF/Protein_ref-286-288_K-concentration_PMF_cutoff-5A_2A-bin_50-200ns 50 200
  
   #  this script read a struc ture file Eqilibrated_system_with_chain.pdb ,  an trajectory file Protein_0-1000ns.xtc ,
   #  then it need parameter to define the pore region ( a cylinder), the axis is set at the center of filter entre (hard coded), 5 is pore radius, the center is located at the entrance pof the pore which sit on the middle on the membrane (hard coded)
   #  then the sample range is 120 below (-120) and 120 above (120) with respect to membrane center, and the region is cutted to 140 bin/slices. 
   #  The output is PMF and density of K and water after the output name PMF/Protein_ref-286-288_K-concentration_PMF_cutoff-5A_2A-bin_50-200ns,  50 200 means only calculate 50 to 200 ns of the trajectory
   




#######  lipid density (fig S3)

 1st step, get a file crecord lipid tail position , run this in three replicas :

   ./Lipid_tail_density_v2.py Eqilibrated_system_with_chain_correctname.pdb sim_1ns.xtc analysis/lipid_density/50-200ns_25A 100 200

 2nd,  cat the output in 1st step to a tmeperory file  tt

   then, following script buuild the stab and slice it 

./Lipid_tail_density_slab_project.py Eqilibrated_system_with_chain.pdb tt analysis/lipid_density/50-200ns  12 5 0 -25   100 -50 50  310
                                                                                         input:  thickness, zstep, zmax, zmin, bin and range, residue to cut the line

    in the begging , it will also print the coordinate iof anchor residue, which should be put in a seperate anchor file for plot (so you know how to overlay this density to a structure file which the anchor residues are marked)

  
#######   time averaged electric field (fig S4)

  after get dx out put from vmd pmepot 

  use  -40 mv as an example: 

  psf for vmd pmepot has no voltage information, following script 1) select pore region 5A radius to calculate the profile;
                                                                  2) add the external electic field , our box is upside down, so -40 became 40, box z axis is 16.3 nm ; in mdp, unit is V/nm : 40/( 16.3   * 1000)  = 0.002453 ;  here unit is mv/A so * 100 ;

                                                                  3) out put will be a proile with raw dx, a profile with external -40 electric field and a dx file with external -40 electric field

   ./vmd_dx_to_potential_filter.py  step6.7_lipid_relaxtion_shifted.pdb  Protein_50-200ns_1A.dx Protein_50-200ns_1A.out_5Arange 0.2454    185.339 185.339 163.559    5

                                             # ref pdb                             ingrid       outname        unit mv/A, with direction     # box size    sample region, x y only
 
  now extract the data for plot :

   ./result_extract.py  Protein_50-200ns_1A.out_5Arange_profile.xvg   Protein_50-200ns_1A.out_5Arange.total.xvg   1


###########  helix tilting 

  ./helix_tilt_rotation.py Eqilibrated_system_with_chain.pdb  ~/Pikes_home/BK_channel/Simulation_hBK/6v3g_mutated/6v3g_PC_restrain_CTD_filter_from_syn_316I/Protein_0-1000ns.xtc helix_tilt_rotation/S6_up_50ps  313 324  50  200

    #  read a structure file, trjectory file, residuen range, time range (ns)   





















