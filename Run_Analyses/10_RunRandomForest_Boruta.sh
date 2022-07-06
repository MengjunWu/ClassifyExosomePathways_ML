#!/bin/bash
filedir=""
usage() {
	echo
	echo ${0##/*}" Usage:"
	echo
	echo "${0##/} OPTIONS -f FILE"
	echo "[-f N]        inputfile"
	exit
}


while getopts f: opt; do
	case ${opt} in
		f) filedir=${OPTARG};;
		*) echo "unrecognized option ${opt}"; usage;;
	esac
done

scriptdir=/binf-isilon/sandelin/people/mengjun//Exosome_ML/Run_Analyses/

for feature in class1 class2 class3 class4 TSS TES all
do
	echo ${feature}
	nice ionice -c 3 Rscript-4.0.3 ${scriptdir}RandomForest_Boruta.R --datadir ${filedir} --featureclass ${feature} 
done