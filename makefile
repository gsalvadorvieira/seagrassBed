# This is a commmentary line in a makefile

# FFLAGS= -fast
# CFLAGS= -O3 -c
# FC=ifort

FLAGS=
CFLAGS= -O3 -c
FC=gfortran

### grass deforms (grasshape_nobending)
NHO  =  main.o init_grass.o momentum.o facediv.o intpol_3rdorder.o vcenter_grass.o vface_grass_openbc.o advecn_openbc.o cdiv.o cpcors.o cpfine_grass_openbc.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill_grass_openbc.o pbc_grass_openbc.o sor.o spline.o seval.o sigma_grass_openbc.o evalrho_none.o hsave.o linerelax_openbc.o dgtsl.o energy.o findzall_grass.o meanh.o vdiffusion_damping_sponge.o checks.o writeslice.o ran3.o tracerinit_grass.o mprove_openbc.o topog_grass.o vort_xz.o writeksurf.o writeyslice.o cfdiv.o steadystate_openbc.o pcorrect.o srcface_grass_openbc.o hsolve_openbc.o vhydro_openbc.o uvchy.o chfine_openbc.o hfill_openbc.o hbc_openbc.o calcpgrad.o grassshape_nobending.o openbc_fix.o updatedh.o outcdf.o 

### grass is fixed (grassshape_inflexible)
# NHO  =  main.o init_grass.o momentum.o facediv.o intpol_3rdorder.o vcenter_grass.o vface_grass_openbc.o advecn_openbc.o cdiv.o cpcors.o cpfine_grass_openbc.o efill.o mgrid.o prolong.o resid.o restrict.o mgpfill_grass_openbc.o pbc_grass_openbc.o sor.o spline.o seval.o sigma_grass_openbc.o evalrho_none.o hsave.o linerelax_openbc.o dgtsl.o energy.o findzall_grass.o meanh.o vdiffusion_damping_sponge.o checks.o writeslice.o ran3.o tracerinit_grass.o mprove_openbc.o topog_grass.o vort_xz.o writeksurf.o writeyslice.o cfdiv.o steadystate_openbc.o pcorrect.o srcface_grass_openbc.o hsolve_openbc.o vhydro_openbc.o uvchy.o chfine_openbc.o hfill_openbc.o hbc_openbc.o calcpgrad.o grassshape_inflexible.o openbc_fix.o updatedh.o outcdf.o

nhg: $(NHO)
	$(FC) $(FFLAGS) -o nhg $(NHO) -L/home/salvadorvieira.g/PSOM/netcdf-c-4.7.1_install/lib -lnetcdf -L/home/salvadorvieira.g/PSOM/netcdf-fortran-4.5.2_install/lib -lnetcdff
#	$(FC) $(FFLAGS) -o nhg $(NHO) -L/usr/local/lib -lnetcdf -lnetcdff

clean:
	mkdir -p output; mv -t output *.cdf *.out *.log *.err *.91; rm -f *.o nhg; rmdir --ignore-fail-on-non-empty output
