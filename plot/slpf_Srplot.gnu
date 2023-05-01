#!/bin/sh
# This script will search the directory and generate plots
# for the SLPF executable and SLTC exec
# Written by: Terence Musho 2012
#

#Parse Tau*.dat extension
ls -la Tau*.dat | awk '{print($8)}' > fname.tmp

#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of Transmission
awk '{print "reset\nset term jpeg\nset output \""$1".jpg\"\nset grid\nset ylabel \"Scattering Rate 1/{/Symbol t} [1/s]\"\nset xlabel \"Energy [eV]\"\n#set xrange [-0.4:1]\nset yrange [1e5:1e18]\nset logscale y\nplot \""$1".dat\" using 1:2 lw 2 w lines title \"Phonon #1\",\""$1".dat\" using 1:3 lw 2 w lines title \"Phonon #2\",\""$1".dat\" using 1:4 lw 2 w lines title \"Phonon #3\",\""$1".dat\" using 1:5  lw 2 w lines title \"Phonon #4\",\""$1".dat\" using 1:6 lw 2 w lines title \"All Phonons\""}' ./fname1.tmp > makeplots.gp


#Parse Tr*.dat extension
ls -la Tr*.dat | awk '{print($8)}' > fname.tmp

#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of Transmission
awk '{print "reset\nset term jpeg\nset yrange [0:1]\nset output \""$1".jpg\"\nset grid\nset ylabel \"Transmission []\"\nset xlabel \"Energy [eV]\"\nplot \""$1".dat\" using 1:2 w lines title \"Subband #1\",\""$1".dat\" using 1:3 w lines title \"Subband #2\",\""$1".dat\" using 1:4 w lines title \"Subband #3\",\""$1".dat\" using 1:5 w lines title \"Subband #4\",\""$1".dat\" using 1:6 w lines title \"Subband #5\""}' ./fname1.tmp >> makeplots.gp

#Parse I*.dat extension
ls -la I_0*.dat | awk '{print($8)}' > fname.tmp

#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of I
awk '{print "reset\nset term jpeg\nset output \""$1".jpg\"\nset grid\nset xlabel \"Normalized Current [A/eV]\"\nset ylabel \"Energy [eV]\"\nplot \""$1".dat\" using 2:1 w lines title \"T21\",\""$1".dat\" using 3:1 w lines title \"T12\",\""$1".dat\" using 4:1 w lines title \"Scatter T13\""}' ./fname1.tmp >> makeplots.gp
awk '{print "reset\nset term jpeg\nset output \""$1"scat.jpg\"\nset grid\nset xlabel \"Normalized Current [A/eV]\"\nset ylabel \"Energy [eV]\"\nplot \""$1".dat\" using 4:1 w lines title \"Scatter T13\""}' ./fname1.tmp >> makeplots.gp

#Parse Rho*.dat extension
ls -la Rho_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of Rho
awk '{print "reset\nset term jpeg\nset output \""$1".jpg\"\nset grid\nset xlabel \"Device Length [m]\"\nset ylabel \"Electron Density [1/m^2]\"\nplot \""$1".dat\" using 1:2 w lines title \"Unbound Electron Density - Hori\",\""$1".dat\" using 1:3 w lines title \"Bound Electron Density - Hori\""}' ./fname1.tmp >> makeplots.gp

#Parse Scin*.dat extension
ls -la Scin_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of Sci
awk '{print "reset\nset term jpeg\nset output \""$1".jpg\"\nset grid\nset xlabel \"Device Length [m]\"\nset ylabel \"Scattered Electron Density [1/m^2]\"\nplot \""$1".dat\" using 1:2 w lines title \"In-Scatter\""}' ./fname1.tmp >> makeplots.gp
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Scattering States (LDOS) - Subband #1\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_p1.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:1e14]\nsplot \""$1".dat\" using 2:1:3"}' ./fname1.tmp >> makeplots.gp
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Scattering States (LDOS) - Subband #2\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_p2.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:1e14]\nsplot \""$1".dat\" using 2:1:4"}' ./fname1.tmp >> makeplots.gp
#awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Scattering States (LDOS) - Subband #3\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_p3.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:1e14]\nsplot \""$1".dat\" using 2:1:5"}' ./fname1.tmp >> makeplots.gp
#awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Scattering States (LDOS) - Subband #4\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_p4.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:1e14]\nsplot \""$1".dat\" using 2:1:6"}' ./fname1.tmp >> makeplots.gp
#awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Scattering States (LDOS) - Subband #4\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_all.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:1e14]\nsplot \""$1".dat\" using 2:1:7"}' ./fname1.tmp >> makeplots.gp

#Parse U*.dat extension
ls -la U_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of U
awk '{print "reset\nset title \"Field Emission - Free Electron Potential Plot - Horiz\"\nset term jpeg\nset output \""$1".jpg\"\nset grid\nset xlabel \"Device Length [m]\"\nset ylabel \"Device Potential [eV]\"\nplot \""$1".dat\" using 1:2 w lines title \"Subband #1\",\""$1".dat\" using 1:3 w lines title \"Subband #2\",\""$1".dat\" using 1:4 w lines title \"Subband #3\",\""$1".dat\" using 1:5 w lines title \"Subband #4\",\""$1".dat\" using 1:6 w lines title \"Subband #5\",\""$1".dat\" using 1:7 w lines title \"All Subbands\""}' ./fname1.tmp >> makeplots.gp

#Generate plot of Uo ionized donors
awk '{print "reset\nset title \"Field Emission - Ionized Donor Potential Plot - Horiz\"\nset term jpeg\nset output \""$1"_o.jpg\"\nset grid\nset xlabel \"Device Length [m]\"\nset ylabel \"Device Potential [eV]\"\nplot \""$1".dat\" using 1:8 w lines title \"Ionized Donors\""}' ./fname1.tmp >> makeplots.gp

#Generate plot of Ua band offset
awk '{print "reset\nset title \"Field Emission - Band-offset Potential Plot - Horiz\"\nset term jpeg\nset output \""$1"_a.jpg\"\nset grid\nset xlabel \"Device Length [m]\"\nset ylabel \"Device Potential [eV]\"\nplot \""$1".dat\" using 1:9 w lines title \"Ionized Donors\""}' ./fname1.tmp >> makeplots.gp

#Parse A*.dat extension
ls -la A_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of LDOS
awk '{print "reset\nset view map\nset samples 201, 201\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_raw.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:10]\nsplot \""$1".dat\" using 2:1:3"}' ./fname1.tmp >> makeplots.gp

#Parse A*.dat extension
ls -la Asub_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of LDOS
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of States (LDOS) - Subband #1\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_sub1.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:3]\nsplot \""$1".dat\" using 2:1:3"}' ./fname1.tmp >> makeplots.gp
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of States (LDOS) - Subband #2\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_sub2.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:3]\nsplot \""$1".dat\" using 2:1:4"}' ./fname1.tmp >> makeplots.gp
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of States (LDOS) - Subband #3\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_sub3.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:3]\nsplot \""$1".dat\" using 2:1:5"}' ./fname1.tmp >> makeplots.gp
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of States (LDOS) - Subband #4\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_sub4.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:3]\nsplot \""$1".dat\" using 2:1:6"}' ./fname1.tmp >> makeplots.gp

#Parse Gn*.dat extension
ls -la Gn_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of LDOS
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Filled States (LDOFS)\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_raw.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:10]\nsplot \""$1".dat\" using 2:1:3"}' ./fname1.tmp >> makeplots.gp

#Parse Gp*.dat extension
ls -la Gp_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of LDOS
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Unfilled States (LDOUS)\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_raw.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:10]\nsplot \""$1".dat\" using 2:1:3"}' ./fname1.tmp >> makeplots.gp

#Parse Gn*.dat extension
ls -la Gni_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of LDOS
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Filled States (LDOFS)\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_raw.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:3]\nsplot \""$1".dat\" using 2:1:3"}' ./fname1.tmp >> makeplots.gp

#Parse Gp*.dat extension
ls -la Gpi_*.dat | awk '{print($8)}' > fname.tmp
#Replace extension with nothing
sed "s/\./ /g" fname.tmp > fname1.tmp

#Generate plot of LDOS
awk '{print "reset\nset view map\nset title \"Superlattice - Local Density of Unfilled States (LDOUS)\"\nset samples 101, 101\nset isosamples 2, 2\nset style data pm3d\nset pm3d implicit at b\nset palette positive nops_allcF maxcolors 0 gamma 1.5 color model HSV \nset palette defined ( 0 0 1 0, 1 0 1 1, 6 0.8333 1 1, 7 0.8333 0 1)\nset term jpeg\nset output \""$1"_raw.jpg\"\nset xlabel \"Device Length [m]\"\nset ylabel \"Energy [eV]\"\nset zrange [0:-3]\nsplot \""$1".dat\" using 2:1:3"}' ./fname1.tmp >> makeplots.gp


echo "food" > fname.tmp

#Generate plot of IV
awk '{print "reset\nset term jpeg\nset output \"IV.jpg\"\nplot \"./IV.dat\" using 1:2"}' ./fname.tmp >> makeplots.gp

#Generate plot of H
awk '{print "reset\nset term jpeg\nset output \"H.jpg\"\nplot \"./H.dat\" using 1:3, \"./H.dat\" using 1:4, \"./H.dat\" using 1:5"}' ./fname.tmp >> makeplots.gp

#Generate plot of Eig
awk '{print "reset\nset title \"Superlattice Device - Eigenstates Plot\"\nset term jpeg\nset output \"Eig.jpg\"\nplot \"./Eig.dat\" using 1:2"}' ./fname.tmp >> makeplots.gp

#Plot *.prog
gnuplot ./makeplots.gp

#Delete temp files
#rm fname.tmp fname1.tmp

