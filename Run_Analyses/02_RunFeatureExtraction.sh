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

echo ${filedir}

############################################preparation for feature extraction###########################################

##get fasta file
genomefasta=/binf-isilon/sandelin/people/mengjun/genome.annotation/hg38/hg38.fa
bedfile1=${filedir}TSS.window1100.bed 
bedfile2=${filedir}TES.window1100.bed 

bedtools getfasta -fi ${genomefasta} -bed ${bedfile1} -name -s > ${filedir}TSS.window1100.fa
bedtools getfasta -fi ${genomefasta} -bed ${bedfile2} -name -s > ${filedir}TES.window1100.fa


##get all RBP overlapping files
RBPfile=/binf-isilon/sandelin/people/mengjun/Downloaded_data/RBP_POSTAR2/RBP_bindingSite_CLIP_sorted.bed
bedfileTES=${filedir}TES.upstream500.bed
bedfileTSS=${filedir}TSS.downstream500.bed

bedtools intersect -wa -wb -s -a ${bedfileTES} -b ${RBPfile} > ${filedir}RBP.TESupstream500.overlap.txt
bedtools intersect -wa -wb -s -a ${bedfileTSS} -b ${RBPfile} > ${filedir}RBP.TSSdownstream500.overlap.txt


###################################################Feature extraction####################################################

##output for all feature
mkdir ${filedir}features
featureout=${filedir}features/
scriptdir=/binf-isilon/sandelin/people/mengjun/Exosome_ML/ExtractFeatures/
fastafileTSS=${filedir}TSS.window1100.fa
fastafileTES=${filedir}TES.window1100.fa

##Chromatin enviroment
chredir=/binf-isilon/sandelin/people/mengjun/Downloaded_data/HistoneMarkers_hg38/
chromatinwindowsize=500
Rscript-3.6.1 ${scriptdir}Extract_Chromatin_Environment.R --bedfile ${bedfile1} --bigwigdir ${chredir} --windowsize ${chromatinwindowsize} --outdir ${featureout}TSS_
Rscript-3.6.1 ${scriptdir}Extract_Chromatin_Environment.R --bedfile ${bedfile2} --bigwigdir ${chredir} --windowsize ${chromatinwindowsize} --outdir ${featureout}TES_

##Transcription level
transdir=/binf-isilon/sandelin/people/mengjun/Downloaded_data/NET_seq_accessibility_hg38/
transwindowup=100
transwindowdn=500
Rscript-3.6.1 ${scriptdir}Extract_Transcription_Levels.R --bedfile ${bedfile1} --bigwigdir ${transdir} --windowup ${transwindowup} --windowdn ${transwindowdn} --outdir ${featureout}

##Calculating the CG content
slidingw=10
windowcg=500
Rscript-3.6.1 ${scriptdir}Promoter_SequenceContent_CGcontent.R --fastafiledir ${fastafileTSS} --slidingWindow ${slidingw} --windowcg ${windowcg} --outdir ${featureout}

##Calculate promoter TATA motifs
Rscript-3.6.1 ${scriptdir}Promoter_SequenceContent_TATA.R --fastafiledir ${fastafileTSS} --outdir ${featureout}

##Calculate U1 PAS RNAbinding motifs from TSS
window=500
Rscript-3.6.1 ${scriptdir}U1PAS_RBP_motifs_TSS.R --fastafiledir ${fastafileTSS} --onesidewindow ${window} --outdir ${featureout}

##Calculate U1 PAS RNAbinding motifs from TES
window1=500
window2=500
Rscript-3.6.1 ${scriptdir}U1PAS_RBP_motifs_TES.R --fastafiledir ${fastafileTES} --onesidewindow ${window1} --twosidewindow ${window2} --outdir ${featureout}


##Calculate CLIP RNA binding proteins
Rscript-3.6.1 ${scriptdir}RBP_CLIP_binding.R --datadir ${filedir} --overlapfile RBP.TSSdownstream500.overlap.txt --outdir ${featureout}TSS_overlap500_
Rscript-3.6.1 ${scriptdir}RBP_CLIP_binding.R --datadir ${filedir} --overlapfile RBP.TESupstream500.overlap.txt --outdir ${featureout}TES_upstream500_

#calculating PAS upstream TES
Rscript-3.6.1 ${scriptdir}PAS_upstream_TES.R --fastafiledir ${fastafileTES} --outdir ${featureout}
