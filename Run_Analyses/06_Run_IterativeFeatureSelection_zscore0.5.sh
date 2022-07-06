#!/bin/bash
filedir=""
zscore=""
usage() {
	echo
	echo ${0##/*}" Usage:"
	echo
	echo "${0##/} OPTIONS -f FILE"
	echo "[-f N]        inputfile dir"
	echo "[-z N]        zscore"
	exit
}


while getopts f:z: opt; do
	case ${opt} in
		f) filedir=${OPTARG};;
		z) zscore=${OPTARG};;
		*) echo "unrecognized option ${opt}"; usage;;
	esac
done

##create iterative directory and copy initial importance files
resultdir=${filedir}RandomForest_output/
cd $resultdir
scriptdir=/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/

for feature in class1 class2 class3 class4 
do
	mkdir iterative_zscore${zscore}_${featureclass}
	for i in ${featureclass}_*.txt
	do
		newfile=runiterative0_${i}
		cp ${i} iterative_zscore${zscore}_${featureclass}/${newfile}	
	done
	for i in {0..4..1}
	do	
		nice Rscript-4.0.3 ${scriptdir}IterativeFeatureSelection_zscore0.5.R --filedir ${resultdir}iterative_zscore${zscore}_${featureclass}/ --iterative ${i} --zscore ${zscore}
		i=$(($i+1))
		nice ionice -c 3 Rscript-4.0.3 ${scriptdir}RandomForest_iterative_zscore0.5.R --datadir ${filedir} --iterativefiledir ${filedir}iterative_zscore${zscore}_${featureclass}/ --iterative ${i} --outdir ${resultdir}iterative_zscore${zscore}_${featureclass}/ --zscore ${zscore}
	done
done



