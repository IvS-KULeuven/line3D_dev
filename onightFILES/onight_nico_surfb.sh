#CALCULATE SOME MODELS OVER NIGHT
function copy_to_directory () {
    #set indat_spec_nico.nml file
    dir_out_surfb=$1/surfb_opal
    
    mkdir $dir_out_surfb

    echo '########### copy to output directory ##################'
    echo 'output_dir: ' $dir_out_surfb
    echo

    cp output_spec.log $dir_out_surfb
    cp output_model.log $dir_out_surfb
    cp output_modelspec.log $dir_out_surfb
    cp modelspec.eo $dir_out_surfb
    cp spec.eo $dir_out_surfb
    cp model.eo $dir_out_surfb
    cp indat_sc3d_nico.nml $dir_out_surfb
    cp indat_modelspec_nico.nml $dir_out_surfb
    cp indat_spec_nico.nml $dir_out_surfb
    cp ./outputFILES/nico_wr3d/spec_surface* $dir_out_surfb
    cp ./outputFILES/nico_wr3d/photospheric_profile.dat $dir_out_surfb

    echo 'done'
    echo
}

function copy_fname () {


   echo '########### copying file from STER ####################'
   echo
    
   fname_base=$1
   isnap=$2
   
   if [ $isnap -lt 10 ]
      then
          fname=$file_base'_000'$isnap'.dat'
   elif [ $isnap -lt 100 ]
      then
         fname=$file_base'_00'$isnap'.dat'
   fi
   fname_in='/STER/levin/models_nico/models/nico_wr3d/'$fname
   fname_out='./models/nico_wr3d/'$fname

   echo $fname_in
   echo $fname_out
   cp $fname_in $fname_out

}

function delete_fname () {


   echo '########### deleting local model file #################'
   echo
   
   fname_base=$1
   isnap=$2
   
   if [ $isnap -lt 10 ]
      then
          fname=$file_base'_000'$isnap'.dat'
   elif [ $isnap -lt 100 ]
      then
         fname=$file_base'_00'$isnap'.dat'
   fi
   fname_out='./models/nico_wr3d/'$fname

   rm $fname_out

}
   
function set_spec () {
    #set indat_spec_nico.nml file
    vmicro=$1
    nalpha=$2
    ngamma=$3

    output_file='indat_spec_nico.nml'
    echo '########### setting indat file ########################'
    echo 'output_file: ' $output_file
    echo

   echo "&input_options" > $output_file
   echo "input_mod = 2" >> $output_file
   echo "input_file = './outputFILES/nico_wr3d/modspec_model00.h5'" >> $output_file
   echo "output_dir = './outputFILES/nico_wr3d'" >> $output_file
   echo "opt_photprof = 0" >> $output_file
   echo "opt_obsdir_read = t" >> $output_file
   echo "opt_surface = t" >> $output_file
   echo "opt_int2d = f" >> $output_file
   echo "opt_incl_gdark = f" >> $output_file
   echo "opt_incl_sdist = f" >> $output_file
   echo "nalpha =" $nalpha >> $output_file
   echo "ngamma =" $ngamma >> $output_file
   echo "/" >> $output_file   
   echo "" >> $output_file
#
   echo "&input_model" >> $output_file
   echo "vrot = 0.d0" >> $output_file
   echo "vth_fiducial = 1.d2" >> $output_file
   echo "vmicro =" $vmicro >> $output_file
   echo "rmin = 1.d0" >> $output_file
#   echo "rmax = 5.97d0" >> $output_file   
   echo "rmax = 10.97d0" >> $output_file
   echo "/" >> $output_file
   echo "" >> $output_file
#
   echo "&input_surface" >> $output_file
   echo "nsurfb = 1"  >> $output_file
   echo "alpha_surface = 1.570796d0" >> $output_file
   echo "gamma_surface = 0.d0" >> $output_file
   echo "xobs_surface = 0.d0" >> $output_file
   echo "/" >> $output_file   
   echo "" >> $output_file
#   
   echo "&dum" >> $output_file
   echo "/" >> $output_file
#
#
#    model.eo < $input_file > $output_file
    echo 'done'
    echo
}


function set_modelspec () {
    #set indat_modelspec_nico.nml file
    opt_opac=$1
    opt_scont=$2
    opt_sline=$3
    kcont=$4
    kline=$5
    opt_vlat=$6
    vmicro=$7

    output_file='indat_modelspec_nico.nml'
    echo '########### setting indat file ########################'
    echo 'output_file: ' $output_file
    echo

   echo "&input_options" > $output_file
   echo "input_file = './outputFILES/output_model00.h5'" >> $output_file
   echo "input_file2 = './inputFILES/model3d.h5'" >> $output_file
   echo "output_file = './outputFILES/nico_wr3d/modspec_model00.h5'" >> $output_file
   echo "input_mod = 19 /" >> $output_file
   echo "" >> $output_file
#
   echo "&input_model" >> $output_file
   echo "teff = 258390.7d0" >> $output_file
   echo "trad = 258390.7d0" >> $output_file
   echo "xlogg = 3.6d0" >> $output_file
   echo "rstar = 1.d0" >> $output_file
#   echo "rmax = 6.d0" >> $output_file
   echo "rmax = 11.d0" >> $output_file   
   echo "tmin = 1.d0" >> $output_file
   echo "xmloss = 1.d-6" >> $output_file
   echo "vmin = 10.d0" >> $output_file
   echo "vmax = 4.d3" >> $output_file
   echo "vmicro =" $vmicro >> $output_file
   echo "vth_fiducial=1.d2" >> $output_file
   echo "beta = 1.d0" >> $output_file
   echo "yhe = .98d0" >> $output_file
   echo "hei = 2.d0" >> $output_file
   echo "/" >> $output_file
   echo "" >> $output_file
#
   echo "&input_line" >> $output_file
   echo "iline = 11" >> $output_file
   echo "eps_line = 0.d0" >> $output_file
   echo "kline =" $kline >> $output_file
   echo "kappa0 = 1.d0" >> $output_file
   echo "alpha = 0.d0 /" >> $output_file
   echo "" >> $output_file
#
   echo "&input_usr" >> $output_file
   echo "opt_opac =" $opt_opac >> $output_file
   echo "opt_scont =" $opt_scont >> $output_file
   echo "opt_sline =" $opt_sline >> $output_file
   echo "kcont =" $kcont >> $output_file
   echo "opt_vlat =" $opt_vlat >> $output_file
   echo "tfloor = 20.d3" >> $output_file   
   echo "/" >> $output_file
   echo "" >> $output_file
#   
   echo "&dum" >> $output_file
   echo "/" >> $output_file
#
#
#    model.eo < $input_file > $output_file
    echo 'done'
    echo
}

function set_model () {
    #set indat_modelspec_nico.nml file
    file_base=$1
    isnap=$2

    output_file='indat_sc3d_nico.nml'
    echo '########### setting model indat file ##################'
    echo 'output_file: ' $output_file
    echo

   echo "&input_options" > $output_file
   echo "model_dir = 'inputFILES'" >> $output_file
   echo "output_file = 'output_model00.h5'" >> $output_file
   echo "input_mod = 17" >> $output_file
   echo "input_mod_dim = 3" >> $output_file
   echo "spatial_grid1d = 5" >> $output_file
   echo "spatial_grid3d = 4" >> $output_file
   echo "opt_opac = 0" >> $output_file
   echo "opt_opal = 0" >> $output_file
   echo "opt_angint_method = 9" >> $output_file
   echo "opt_method = 1" >> $output_file
   echo "opt_sol2d = f" >> $output_file
   echo "opt_ltec = 1" >> $output_file
   echo "opt_incl_cont = t" >> $output_file
   echo "opt_start_cont = t" >> $output_file
   echo "opt_ng_cont = t" >> $output_file
   echo "opt_ait_cont = f" >> $output_file
   echo "opt_incl_line = f" >> $output_file
   echo "opt_start_line = t" >> $output_file
   echo "opt_ng_line = t" >> $output_file
   echo "opt_ait_line = f" >> $output_file
   echo "opt_alo_cont = 3" >> $output_file
   echo "opt_alo_line = 3"  >> $output_file
   echo "opt_incl_gdark = f" >> $output_file
   echo "opt_incl_sdist = f /" >> $output_file   
   echo "" >> $output_file
   #
   echo "input_mod=0 if spherically symmetric model is used (beta-velocity-law)
               temperatures are read from JO's solution: RAD_JO.dat, TEMP_JO.dat" >> $output_file
   echo "input_mod=1 if Dylan Kee's 2d model is used (snapshot)" >> $output_file
   echo "input_mod=2 if Dylan Kee's 2d model is uised (wind ablation from initial
               conditions)" >> $output_file
   echo "input_mod=3 if ADM model is used" >> $output_file
   echo "input_mod=4 if ADM model is used (with spherical beta-velocities in wind)" >> $output_file
   echo "input_mod=5 if MHD snapshot is read in (from Asif's 3d atmosphere)" >> $output_file
   echo "" >> $output_file
   echo "input_mod_dim=[1,2,3] if 1d, 2d, 3d model shall be read in" >> $output_file
   echo "" >> $output_file
   echo "spatial_grid3d=0 if 3d grid is calculated from 1d grid with equidistant core
                 points " >> $output_file
   echo "spatial_grid3d=1 if 3d grid is calculated from a mean-value approach
                 (minimizing distance of subsequent coordinates from 1d-grid)" >> $output_file
   echo "spatial_grid3d=2 if 3d grid is calculated from a mean-value approach 
                 (minimizing distance of subsequent coordinates from original
                 input-grid)" >> $output_file
   echo "spatial_grid3d=3 if 3d grid is calculated completely equidistant" >> $output_file
   echo "spatial_grid3d=4 if 3d grid is calculated from 1d input model (with corresponding pdf)" >> $output_file
   echo "spatial_grid3d=5 if 3d grid is calculated from 3d input model (with corresponding pdf)" >> $output_file
   echo "" >> $output_file
   echo "spatial_grid1d=0 if equidistant radial grid is used (subroutine grid1d_r_equi)" >> $output_file
   echo "spatial_grid1d=1 if equidistant velocity grid is used (subroutine grid1d_vel_equi)" >> $output_file
   echo "spatial_grid1d=2 if equidistant tau_thomson grid is used (subroutine grid1d_tau_equi)" >> $output_file
   echo "spatial_grid1d=3 if equidistant log(tau_thomson) grid is used (subroutine grid1d_tau_log)" >> $output_file
   echo "spatial_grid1d=4 if combination is used (see subroutine grid1d_final for details)" >> $output_file
   echo "spatial_grid1d=5 if combination is used (see subroutine grid1d_final_2 for details)" >> $output_file
   echo "spatial_grid1d=6 if grid is calucalted equidistant in log-space (subroutine grid1d_r_log)" >> $output_file
   echo "" >> $output_file
   echo "&input_mod_1d" >> $output_file
   echo "teff = 258390.7d0" >> $output_file
   echo "trad = 258390.7d0" >> $output_file
   echo "xlogg = 3.6d0" >> $output_file
   echo "rstar = 1.d0" >> $output_file
   #   echo "rmax = 6.d0" >> $output_file
   echo "rmax = 11.d0" >> $output_file   
   echo "tmin = 1.d0" >> $output_file
   echo "xmloss = 1.d-6" >> $output_file
   echo "vmin = 10.d0" >> $output_file
   echo "vmax = 4000.d0" >> $output_file
   echo "vmicro = 1.d2" >> $output_file
   echo "vth_fiducial = 1.d2" >> $output_file
   echo "vrot = 0.d0" >> $output_file
   echo "beta = 1.d0" >> $output_file
   echo "yhe = .98d0" >> $output_file
   echo "hei = 2.d0 " >> $output_file
   echo "xnue0 = 1.93798d15" >> $output_file
   echo "na = 12 /" >> $output_file
   echo "" >> $output_file
   echo "" >> $output_file
   echo "TRANSITION-FREQUENCY (LAMBDA=1367 A):" >> $output_file
   echo "XNUE0=2.191983810829996D15" >> $output_file
   echo "FOR MODELS PULS DIPLOMA (LAMBDA=1548 A):" >> $output_file
   echo "XNUE0=1.93798D15" >> $output_file
   echo "H-ALPHA:" >> $output_file
   echo "XNUE0=4.5680294d14" >> $output_file
   echo "" >> $output_file
   echo "&input_infreg" >> $output_file
   echo "rmin = 1.d0" >> $output_file
   echo "rlim = 5.9d0 /" >> $output_file
   echo "" >> $output_file
   echo "&input_cont" >> $output_file
   echo "eps_cont = 0.d0" >> $output_file
   echo "kcont = 1.d0 /" >> $output_file
   echo "" >> $output_file
   echo "&input_line" >> $output_file
   echo "eps_line = 0.d0" >> $output_file
   echo "kline = 1.d0" >> $output_file
   echo "kappa0 = 5.d-1" >> $output_file
   echo "alpha = 0.d0 /" >> $output_file
   echo "" >> $output_file
   echo "&input_mod_adm" >> $output_file
   echo "ralfven = 2.7d0" >> $output_file
   echo "delta = 0.5d0" >> $output_file
   echo "chi_inf = 1.d-1" >> $output_file
   echo "obliquity = 0.d0 /" >> $output_file
   echo "obliquity in degree" >> $output_file
   echo "" >> $output_file
   echo "&input_mod_bc" >> $output_file
   echo "mstar = 52.5d0" >> $output_file
   echo "zeta = 4.18d0" >> $output_file
   echo "gamma = 0.35d0" >> $output_file
   echo "xi = -0.43d0" >> $output_file
   echo "alpha_cak = 0.66d0" >> $output_file
   echo "delta_cak = 0.07d0 /" >> $output_file
   echo "" >> $output_file
   echo "&input_mod_abl" >> $output_file
   echo "theta_d = 12.d0" >> $output_file
   echo "dtheta_abl = 4.d0" >> $output_file
   echo "beta_accr = 1.5d0" >> $output_file
   echo "tau_d = 1.d3 /" >> $output_file
   echo "" >> $output_file
   echo "&input_mod_be" >> $output_file
   echo "mdisc = 4.5d-8" >> $output_file
   echo "rdisc = 13.5d0" >> $output_file
   echo "tdisc = 30.d3" >> $output_file
   echo "dtheta = 45.d0" >> $output_file
   echo "slope = 1.5d0" >> $output_file
   echo "/" >> $output_file
   echo "" >> $output_file
   echo "&dimensions_1d" >> $output_file
   echo "n1d = 27" >> $output_file
   echo "n1d_t = 81" >> $output_file
   echo "n1d_r = 22" >> $output_file
   echo "delv = 0.3333333d0 /" >> $output_file
   echo "n1d=33" >> $output_file
   echo "" >> $output_file
   echo "!three different grids:" >> $output_file
   echo "n1d=60" >> $output_file
   echo "n1d=40" >> $output_file
   echo "n1d=17" >> $output_file
   echo "n1d_t = 81" >> $output_file
   echo "n1d_r = 22" >> $output_file
   echo "" >> $output_file
   echo "&dimensions_3d" >> $output_file
   echo "ncx=19" >> $output_file
   echo "ncy=19" >> $output_file
   echo "ncz=19" >> $output_file
   echo "delx_max=0.7d0" >> $output_file
   echo "dely_max=0.7d0 " >> $output_file
   echo "delz_max=0.7d0 /" >> $output_file
   echo "" >> $output_file
   echo "&dimensions_freq" >> $output_file
   echo "deltax = 0.3333333d0" >> $output_file
   echo "xcmf_max = 3.d0 /" >> $output_file
   echo "" >> $output_file
   echo "&dimensions_angles" >> $output_file
   echo "n_theta = 16 /" >> $output_file
   echo "16" >> $output_file
   echo "!low resolution:" >> $output_file
   echo "n_theta=6" >> $output_file
   echo "n_phi=5" >> $output_file
   echo "!high resolution:" >> $output_file
   echo "n_theta=16" >> $output_file
   echo "n_phi=5" >> $output_file
   echo "" >> $output_file
   echo "&benchmark" >> $output_file
   echo "benchmark_mod = 0" >> $output_file
   echo "im_source = 3" >> $output_file
   echo "im_opacity = 2" >> $output_file
   echo "im_vel = 0" >> $output_file
   echo "tau_min = 0.d0" >> $output_file
   echo "tau_max = 5.d0" >> $output_file
   echo "source_min = 0.1d0" >> $output_file
   echo "source_max = 1.d-6" >> $output_file
   echo "n_y = 0.d0" >> $output_file
   echo "n_z = 0.707107d0 /" >> $output_file
   echo "n_z = 0.707107d0 " >> $output_file
   echo "" >> $output_file
   echo "&input_usr" >> $output_file

   strdum='models/nico_wr3d/'$file_base
#   echo "fname_model='models/nico_wr3d/WR_3D_alpha_0.50'" >> $output_file
   echo "fname_model='$strdum'" >> $output_file
   echo "is_min="$isnap >> $output_file
   echo "is_max="$isnap"   !40" >> $output_file
   echo "unit_length = 1.d0" >> $output_file
   echo "unit_density = 2.5d-8" >> $output_file
   echo "unit_velocity = 1.d8" >> $output_file
   echo "unit_temperature = 1.d0" >> $output_file
   echo "opt_bvel=0" >> $output_file
   echo "max_refinement=4" >> $output_file
   echo "nd = 512, 64, 64" >> $output_file
   echo "/" >> $output_file
   echo "" >> $output_file
   echo "&test" >> $output_file
   echo "/" >> $output_file
   echo "test" >> $output_file
#
#  model.eo < $input_file > $output_file
    echo 'done'
    echo
}


function run_model () {
    #create model atmosphere
    input_file='in_sc3d'
    output_file='output_model.log'
    echo '########### calculating model atmosphere ##############'
    echo 'input_file: ' $input_file
    echo 'output_file: ' $output_file
    model.eo < $input_file > $output_file
    echo 'done'
    echo
}
#
function run_sc3d () {
    #calculate source functions
    input_file='in'
    output_file='output_sc3d.log'
    echo '######### calculating solution of source fct ##########'
    echo 'input_file: ' $input_file
    echo 'output_file: ' $output_file
    sc3d.eo < $input_file > $output_file
    echo 'done'
    echo    
}
#
function run_modelspec () {
    #calculate model-file for final formal solution
    input_file='in_modelspec'
    output_file='output_modelspec.log'
    echo '###### calculating model-file for formal solution #####'
    echo 'input_file: ' $input_file
    echo 'output_file: ' $output_file
    modelspec.eo < $input_file > $output_file
    echo 'done'
    echo    
}
#
function run_spec () {
    #calculate formal solution
    input_file='in_spec'
    output_file='output_spec.log'
    echo '############# calculating formal solution #############'
    echo 'input_file: ' $input_file
    echo 'output_file: ' $output_file
#    nohup spec.eo < $input_file > $output_file &
    spec.eo < $input_file > $output_file
    echo 'done'
    echo    
}
#
#
#
#
#run program
#
#-------------------------------------------------------------------
#
#opt_opac = 0 -> zero
#         = 1 -> opal opacities
#         = 2 -> thomson scattering
#opt_scont = 0 -> zero
#          = 1 -> planck function
#          = 2 -> dilution factor
#opt_sline = 0 -> zero
#          = 1 -> planck function
#          = 2 -> dilution factor
#opt_vlat = 0 -> only radial velocities
#         = 1 -> also theta,phi velocities

#set indat_modelspec_nico.nml file
vmicro='100.d0'
opt_opac='1'
opt_scont='1'
opt_sline='0'
kcont='1.d0'
kline='1.d0'
opt_vlat='0'
nalpha=1
ngamma=1
#
#--------------------------------------------
#
#--------------------------------------------
#
file_base='WR_3D_alpha_LTE_longbox'
#
#--------------------------------------------
#
isnap=21
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap21_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap21_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=22
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap22_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap22_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=23
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap23_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap23_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=24
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap24_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap24_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=25
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap25_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap25_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=26
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap26_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap26_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=27
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap27_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap27_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=28
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap28_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap28_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=29
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap29_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap29_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=30
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap30_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap30_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=31
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap31_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap31_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=32
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap32_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap32_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=33
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap33_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap33_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=34
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap34_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap34_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=35
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap35_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap35_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=36
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap36_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap36_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=37
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap37_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap37_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=38
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap38_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap38_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=39
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap39_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap39_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=40
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap40_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap40_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=41
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap41_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap41_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=42
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap42_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap42_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=43
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap43_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap43_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=44
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap44_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap44_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=45
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap45_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap45_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=46
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap46_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap46_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=47
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap47_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap47_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=48
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap48_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap48_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
#
#--------------------------------------------
#
isnap=49
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap49_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap49_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

#
#--------------------------------------------
#
isnap=40
copy_fname $file_base $isnap
set_model $file_base $isnap
run_model
delete_fname $file_base $isnap
#
#
vmicro='1.d2'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap50_kl1d13_vmicro1d2_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out

vmicro='1.d1'
dir_out='/lhome/levin/Postdoc/papers/paperIV/models/WR_3D_alpha_LTE_longbox/snap50_kl1d13_vmicro1d1_sresol'
set_modelspec $opt_opac $opt_scont $opt_sline $kcont $kline $opt_vlat $vmicro
set_spec $vmicro $nalpha $ngamma
run_modelspec
run_spec
copy_to_directory $dir_out
