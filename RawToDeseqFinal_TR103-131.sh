#!/bin/bash -x

SECONDS=0

#Run-specific user inputs.
Id=("TR103" "TR104" "TR105" "TR106" "TR107" "TR108" "TR109" "TR110" "TR111" "TR112" "TR113" "TR114" "TR115" "TR116" "TR117" "TR118" "TR119" "TR120" "TR121" "TR122" "TR123" "TR125" "TR126" "TR127" "TR128" "TR129" "TR130" "TR131") #Names of each sample in batch, respective to FastqList below (Id and FastqList MUST be the same length)
FastqList=("TR103_S28_L002" "TR104_S29_L002" "TR105_S30_L002" "TR106_S31_L002" "TR107_S32_L002" "TR108_S33_L002" "TR109_S34_L002" "TR110_S35_L002" "TR111_S36_L002" "TR112_S37_L002" "TR113_S38_L002" "TR114_S39_L002" "TR115_S40_L002" "TR116_S41_L002" "TR117_S42_L002" "TR118_S43_L002" "TR119_S44_L002" "TR120_S45_L002" "TR121_S46_L002" "TR122_S47_L002" "TR123_S48_L002" "TR125_S49_L002" "TR126_S50_L002" "TR127_S51_L002" "TR128_S52_L002" "TR129_S53_L002" "TR130_S54_L002" "TR131_S55_L002") #Names of raw .fastq.gz file names using the following format: <prefix>_R1_001 and <prefix>_R2_001 and located in loc directory defined below
loc=~/bin/OracleMount/TR103-131 #Full path to parent directory where all data will be stored and referenced. This must containing the raw fastq files within a subdirectory named "Raw" 
genome=~/bin/OracleMount/HumanGenome/grch38_tran/genome_tran #Full path to reference genome file(s)
annot=~/bin/OracleMount/HumanGenome/grch38_tran/Homo_sapiens.GRCh38.84.gtf #Full path to reference annotation file
mrgLstLoc=~/bin/OracleMount/TR103-131/merge_list.txt #location of "merge_list.txt" file containing <Id>.gtf (one per line) for each sample in batch

#More user options.
thr=20 #Number of threads to use in processing (based on system specs, default to DiazCrunch is 20)
fqc=true #should fastQC be run?
trm=true #should trimming be done?
mrg=false #should GTF be merged for novel transcipts?
cleanup=false #Do you want to remove all intermediate files that are not deleted by default? These include: entire ${loc}/Bam directory, entire ${loc}/GTF directory, entire ${loc}/TrimmedFastq directory, entire ${loc}/TableCounts directory. The only directories/files not removed by this are ${loc}/ForDESeq2, ${loc}/Raw, and ${loc}/FastQC. Do not enable this option if you are troubleshooting for interested in intermediate files. Enabling this will mean the entire process from raw fastq files needs to be restarted to make any changes. This will take a lot of compute time.

#Required tools and versions (Also required: FastQC, samtools, Python, and dependencies for all tools listed below)
TrimLoc=~/bin/Trimmomatic-0.38/trimmomatic-0.38.jar #Path to Trimmomatic .jar
HisatDir=~/bin/hisat2-2.1.0 #Path to HiSat2 directory
StringtieDir=~/bin/stringtie-1.3.4d.Linux_x86_64 #Path to Stringtie directory
prepLoc=~/bin/prepDE.py #Path to prepDE.py script provided by stringtie authors at: https://ccb.jhu.edu/software/stringtie/dl/prepDE.py


#
#Code below "should" no longer require user input.
#

#Basic sutomatic directory setup and variable definition
cd $loc
mkdir -p Bam Sam FastQC/Raw FastQC/Trimmed ForDESeq2 GTF TableCounts TrimmedFastq #make required directories within <loc> if they do not already exist
cd ~/bin/
mkdir -p TMP

#FastQC of raw fastq files if fqc is true. Output can be found in ${loc}/FastQC/Raw/
if [ "$fqc" = true ]
then
	fastq=()
	for i in ${FastqList[@]}
	do
		fastq+=(${loc}/Raw/${i}_R1_001.fastq.gz ${loc}/Raw/${i}_R2_001.fastq.gz)
	done
	cd ${loc}/Raw
	fastqc -o ${loc}/FastQC/Raw/ -t ${thr} ${fastq[@]}
fi

#Trimming of Raw fastq files if trm is true. Outputs resulting .fastq.gz files to ${loc}/TrimmedFastq/
if [ "$trm" = true ]
then
	inloc=${loc}/TrimmedFastq/
	in1=_P_R1_001.fastq.gz
	in2=_P_R2_001.fastq.gz
	#trimming loop
	cd ~/bin/	
	for i in ${FastqList[@]}
	do
		java -jar $TrimLoc PE -threads ${thr} ${loc}/Raw/${i}_R1_001.fastq.gz ${loc}/Raw/${i}_R2_001.fastq.gz ${inloc}${i}${in1} ${inloc}${i}_U_R1_001.fastq.gz ${inloc}${i}${in2} ${inloc}${i}_U_R2_001.fastq.gz ILLUMINACLIP:Trimmomatic-0.38/adapters/TruSeq3-PE-2.fa:2:30:10 SLIDINGWINDOW:4:15 LEADING:3 TRAILING:3 MINLEN:36
	done
else
	inloc=${loc}/Raw/
	in1=R1_001.fastq.gz
	in2=R2_001.fastq.gz
fi

#FastQC of trimmed fastq files. Output can be found in ${loc}/FastQC/Trimmed/
if [ "$fqc" = true ] && [ "$trm" = true ]
then
	fastq=()
	for i in ${FastqList[@]}
	do
		fastq+=(${inloc}${i}${in1} ${inloc}${i}${in2})
		#Since unpaired ("_U_") outputs are not used downstream, only paired ("_P_") outputs are analyzed by default. If you want to run FastQC on unpaired outputs as well then remove the hash from the code below this comment.
		#fastq+=(${inloc}${i}_U_R1_001.fastq.gz ${inloc}${i}_U_R2_001.fastq.gz)
	done
	cd ${loc}/TrimmedFastq
	fastqc -o ${loc}/FastQC/Trimmed/ -t ${thr} ${fastq[@]}
fi

#Sequence alignment, Sam>Bam conversion, and annotation of each Sample.
c=0
for i in ${FastqList[@]}
do
	cd $HisatDir
	hisat2 -p ${thr} --dta -x $genome -1 ${inloc}${i}${in1} -2 ${inloc}${i}${in2} -S ${loc}/Sam/${Id[c]}.sam
	samtools sort -@ ${thr} -o ${loc}/Bam/${Id[c]}.bam ${loc}/Sam/${Id[c]}.sam
	rm ${loc}/Sam/${Id[c]}.sam #removes sam files no longer needed to save space, hash this out if you want to keep sam files
	cd $StringtieDir
	./stringtie -p ${thr} -G $annot -o ${loc}/GTF/${Id[c]}.gtf -l ${Id[c]} ${loc}/Bam/${Id[c]}.bam
	cp ${loc}/GTF/${Id[c]}.gtf $StringtieDir
	((c++))
done

#Merge annotated .gtf files from all samples
cd $StringtieDir
if [ "$mrg" = true ]
then
	cp $mrgLstLoc $StringtieDir
	./stringtie --merge -p ${thr} -G $annot -o ${loc}/GTF/stringtie_merged.gtf merge_list.txt
	rm merge_list.txt
	annotver=${loc}/GTF/stringtie_merged.gtf
else
	annotver=$annot
fi

#Create table counts
for i in ${Id[@]}
do
	./stringtie -e -B -p ${thr} -G $annotver -o ${loc}/TableCounts/${i}.gtf -l ${i} ${loc}/Bam/${i}.bam
	rm ${i}.gtf
done

#Prepare files for DESeq2 using prepDE.py script
for i in ${Id[@]}
do
	cd ~/bin/TMP
	mkdir ${i}
	cp ${loc}/TableCounts/${i}.gtf ${i}
done
cd ~/bin/
python $prepLoc -i ~/bin/TMP
cd $loc
cp ~/bin/gene_count_matrix.csv ~/bin/transcript_count_matrix.csv ForDESeq2

#A little cleanup
rm ~/bin/gene_count_matrix.csv ~/bin/transcript_count_matrix.csv
cd ~/bin/TMP
for i in ${Id[@]}
do
	rm -r ${i}
done

#Optional thorough cleanup see comment in user input for details.
if [ "$cleanup" = true ]
then
	rm -r ${loc}/Bam ${loc}/GTF ${loc}/TrimmedFastq ${loc}/TableCounts ${loc}/Sam
fi

#Show total runtime
echo "Total Runtime: $(($SECONDS / 3600))hrs $((($SECONDS / 60) % 60))min $(($SECONDS % 60))sec"
echo "DONE!"
