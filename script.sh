fname=run902g

for i in _half _double _quad
# for i in -0 -m3
    do
            cp -r run901g $fname$i
            cd $fname$i; make clean; rm -r output; cd ..
            # cp run900_model/advecn_openbc.f $fname$i
            echo $fname$i; sleep 2; vim $fname$i/parameters.in
            cd $fname$i; sbatch job.script; cd ..
            # vim -u $fname$i/output/*.log
    done


# fname=run901

# for i in a b c d e f g
# # for i in -0 -m3
#     do
#             # cp -r run900_model $fname$i
#             # cp run900_model/advecn_openbc.f $fname$i
#             # echo $fname$i; sleep 1; vim $fname$i/parameters.in
#             cd $fname$i; rm -r output; cd ..
#             # cd $fname$i; sbatch job.script; cd ..
#             # vim -u $fname$i/output/*.log
#     done

# fname=run81

# for i in 1-short 2-short 3-short
# # for i in -0 -m3
#     do
#             # cp -r ${fname}_model $fname$i
#             # cp run50/parameters.in $fname$i
#             echo $fname$i; sleep 1; vim $fname$i/parameters.in
#             cd $fname$i; sbatch job.script; cd ..
#             # vim -u $fname$i/output/*.log
#     done


# fname=run504

# for i in -0 -m1 -m2 -m3 -p1 -p2 -p3
# # for i in -0 -m3
#     do
#             cp -r ${fname}_model $fname$i
#             cp run501$i/parameters.in $fname$i
#             #  echo $fname$i; sleep 1; vim $fname$i/parameters.in
#             cd $fname$i; sbatch job.script; cd ..
#             # vim -u $fname$i/output/*.log
#     done

# fname=run211

# for i in a b c d e f; do cp -r run140_model $fname$i;
# for i in b c d e f; do cp -r run210_model $fname$i;
# cd $fname$i; nano parameters.in; sbatch job.script; cd ../; done

# for i in a b c d e; do mkdir $fname$i; cp run120_model/* $fname$i; done
# for i in a b c d e; do cd $fname$i; nano parameters.in; sbatch job.script; cd ../; done

# CLEAN FOLDER AND EDIT PARAMETERS
# for i in a b c d e; do cd $fname$i; make clean; rm -r output; nano parameters.in; cd ../; done

# SUBMIT JOBS
# for i in a b c d e; do cd $fname$i; sbatch job.script; cd ../; done

# COPY JOB SCRIPT
# for i in a b c d e f; do cp run120_model/job.script $fname$i; done


# #########
# fname1=run30
# fname2=run40

# N=1

# for i in a b c d e f g h i
# 	do
# 		cp -r ${fname1}$N$i ${fname2}$N$i
# 	done
# for i in a b c d e f g h i
# 	do
# 		cp run120_model/job.script ${fname2}$N$i
# 	done

#########
# fname1=run21
# fname2=run31

# N=1

# # for i in a b c d e f g h i
# for i in b c d e f
# 	do
# 		cp -r ${fname1}$N$i ${fname2}$N$i
# 		cd ${fname2}$N$i
# 		cp ../run200_model/updateParameters.sh .
# 		make clean
# 		rm -r output
# 		sbatch job.script
# 		cd ..
# 		echo ${fname2}$N$i done
# 	done


##########
#fname=run401

#for i in a b c d e f g h i
#	do
#		# cp -r run400_model $fname$i
#		# echo $fname$i; sleep 1; nano $fname$i/parameters.in
#		# cd $fname$i; sbatch job.script; cd ..
#		#nano $fname$i/output/*.log
#	done



