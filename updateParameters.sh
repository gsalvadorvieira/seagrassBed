#!/bin/bash
filename='parameters.in'

#############

#EDIT GULF.IN
sed -i 4r<(sed '2,5!d' $filename) gulf.in
sed -i -e '1,4d' gulf.in

#EDIT MG.IN
sed -i 1r<(sed '9!d' $filename) mg.in
sed -i -e '1d' mg.in

# RUN PREPROC
gfortran preproc.f -o preproc.out
./preproc.out > batch

# READ NUMBERS TO BE USED
NI=$(awk 'NR==1 {print $9}' batch)
NJ=$(awk 'NR==1 {print $10}' batch)
NK=$(awk 'NR==1 {print $11}' batch)
NKg=$(awk 'NR==13 {print $1}' $filename)
drho=$(awk 'NR==17 {print $1}' $filename)

maxdim=$(awk 'NR==5 {print $3}' batch)
maxint=$(awk 'NR==5 {print $6}' batch)
ntint=$(awk 'NR==2 {print $3}' batch)

dampAmp=$(awk 'NR==21 {print $1}' $filename)
dampFrac=$(awk 'NR==21 {print $2}' $filename)
grassFrac=$(awk 'NR==21 {print $3}' $filename)

tol=$(awk 'NR==9 {print $7}' $filename)

dH=$(awk 'NR==24 {print $1}' $filename)

dz_0=$(awk 'NR==27 {print $1}' $filename)

# MODIFY DIMS.F
text="\ \ \ \ \ \ parameter(NI=$NI,NJ=$NJ,NK=$NK,NKg=$NKg,ntr=1)"
sed -i "6 c $text" dims.f

sed -i "7 c \ \ \ \ \ \ parameter(maxout=$maxdim)" dims.f

# MODIFY MGRID.F
text="\ \ \ \ \ \ parameter(ngrid=3,maxout=$maxdim,maxint=$maxint,int1=$ntint)"
sed -i "33 c $text" mgrid.f

# MODIFY grassshape_nobending.f
text="\ \ \ \ \ \ parameter(drho=$drho)"
sed -i "28 c $text" grassshape_nobending.f

# MODIFY advecn_openbc.f
text="\ \ \ \ \ \ parameter(dampAmp=$dampAmp, dampFrac=$dampFrac)"
sed -i "35 c $text" advecn_openbc.f

# MODIFY advecn_openbc.f
text="\ \ \ \ \ \ parameter(grassFrac=$grassFrac)"
sed -i "32 c $text" vdiffusion_damping_sponge.f

# MODIFY header.f
text="\ \ \ \ \ &     dHinit=$dH, dHfinal=$dH,"
sed -i "23 c $text" header.f

text="\ \ \ \ \ &     dz=$dz_0)"
sed -i "24 c $text" header.f

# MODIFY hsolve_openbc.f
text="\ \ \ \ \ \ tol=$tol"
sed -i "14 c $text" hsolve_openbc.f

#############
## Output on screen
#############

echo -e 'PARAMETERS MODIFIED ACCORDING TO PARAMETERS.IN:\n'
cat $filename
echo -e '\n'
