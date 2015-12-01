#!/bin/bash

SNPfiles_prefix=SNPfiles/Diplo_SNP-
matrix=matrix/av_matrix.matrix
k=100000        #number of iterations
p=4             #number of populations
env_file=env/Diplotaxodon_Morphometric_Data_raw.bayenv
out_dir=run_Bayenv      #prefix for results
job_prefix=Diplo
rep=10  #number of replicates
SNPs_per_run=10

bayenv_path=/home/c/ch380/src/Bayenv/bayenv2

##############################
#echo -e "find total number of SNPs"
SNPcount=$(ls -1 $SNPfiles_prefix* | grep "txt$" | wc -l)
echo -e "\nFound $SNPcount SNP files with prefix $SNPfiles_prefix\n"
echo -e "You specified $SNPs_per_run SNPs to be processed by each job"
#runs=$(($SNPcount/$SNPs_per_run)) #an integer will always be rounded down, the following line makes sure I get the round up value when appropriate
runs=$((($SNPcount+$SNPs_per_run/2)/$SNPs_per_run))
echo -e "and are therefore requesting a total of $runs parallel jobs\n"

#find number of environmental factors
number_env=$(cat $env_file | wc -l)
echo -e "Found $number_env environmental variables in $env_file\n"

#establish working directory
pwd=$(pwd)

for a in $(seq 1 $rep)
do
	ran=$RANDOM
	echo -e "Replicate $a (random seed: $ran)\n"
	cd $pwd
	mkdir $out_dir/replicate_$a
	cd $out_dir/replicate_$a
	echo -e "#!/bin/bash

#PBS -N bayenv-$job_prefix-$a
#PBS -l walltime=00:02:00
#PBS -l vmem=4gb
#PBS -l nodes=1:ppn=2
#PBS -m a
#PBS -M ch380@le.ac.uk
#PBS -t 1-$runs
#PBS -j oe

#LOAD MODULE
#
date
cd $pwd/$out_dir/replicate_$a

start=\$(((\$PBS_ARRAYID*$SNPs_per_run)-($SNPs_per_run-1)))
end=\$((\$start+$SNPs_per_run-1))

if [ \$end -gt $SNPcount ]
then 
	end=$SNPcount
fi

echo -e \"Processing SNPs \$start to \$end\\\n\"

for i in \$(seq \$start \$end)
do
	echo -e \"\\\n####\\\nSNP \$i\\\n\"
	date
	$bayenv_path -i $pwd/$SNPfiles_prefix\$i.txt -e $pwd/$env_file -m $pwd/$matrix -k $k -r -$ran -p $p -n $number_env -t -X -o $job_prefix-replicate-$a-\$PBS_ARRAYID
done
" > $pwd/$out_dir/replicate_$a/$job_prefix-replicate-$a.sh
	qsub $pwd/$out_dir/replicate_$a/$job_prefix-replicate-$a.sh
done
