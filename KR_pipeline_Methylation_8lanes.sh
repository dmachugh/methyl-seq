#!/bin/bash

###
# Genome previously indexed
###

# See KR_pipeline_Methylation_Epicentre-FINAL_and_bwameth.txt
# Index the genome 
#nohup bwameth.py index /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa > bwa-meth-index.out &

###
# Align raw reads using bwameth
###

mkdir /workspace/scratch/krue/Methylation/bwa_8lanes
cd /workspace/scratch/krue/Methylation/bwa_8lanes

ls /workspace/storage/krue/Methylation/Paired/*gz | wc -l
ls /workspace/storage/krue/Methylation/extra6lanes/20150209_D00238RB_BSSeq_PE/*gz | wc -l
ls /workspace/storage/krue/Methylation/extra6lanes/20150209_D00731A_BSSeq_PE/*gz | wc -l
ls /workspace/storage/krue/Methylation/extra6lanes/20150209_D00731B_BSSeq_PE/*gz | wc -l

# Alignment of raw unfiltered reads using bwameth.py
# Concatenate the four forward files, and four reverse files in a single alignment
# get the full file path to the forward read of the first lane
for fastq1_1 in `find /workspace/storage/krue/Methylation/Paired -name *R1_001.fastq.gz`
do
    # reverse read of first lane
    fastq1_2=`echo $fastq1_1 | perl -p -e 's/R1_001.fastq.gz$/R2_001.fastq.gz/'`
    # Full basename of fastq1_1
    baseFastq1=`basename $fastq1_1`
    # Sample name of fastq1_1
    sampleName=`echo $baseFastq1 | perl -p -e 's/_L00[12]_R1_001.fastq.gz$//'`
    sampleShortName=`echo $sampleName | perl -p -e 's/.{7}$//'`
    echo $sampleShortName
    echo $fastq1_1
    # Second run (2 lanes = 1 pool per lane = 1 pair of file for each sample)
    # find the name of the forward read
    fastq2_1=`find /workspace/storage/krue/Methylation/extra6lanes/20150209_D00238RB_BSSeq_PE -name "$sampleShortName"_[ATGC]*_L00[12]_R1_001.fastq.gz`
    echo $fastq2_1
    # reverse read of second lane
    fastq2_2=`echo $fastq2_1 | perl -p -e 's/R1_001.fastq.gz$/R2_001.fastq.gz/'`
    # Third run (2 lanes = same pool on both lanes = 2 pair of file for half the sample)
    # Fourth run (2 lanes = other pool on both lanes = 2 pair of file for half the sample)
    # find the name of the forward read
    fastq3_find=(`find /workspace/storage/krue/Methylation/extra6lanes/20150209_D00731A_BSSeq_PE -name "$sampleShortName"_[ATGC]*_R1_001.fastq.gz`)
    #echo 'fastq3_find '$fastq3_find
    #echo ${#fastq3_find[@]}
    #fastq4_find=(`find /workspace/storage/krue/Methylation/extra6lanes/20150209_D00731B_BSSeq_PE -name "$sampleShortName"*R1_001.fastq.gz`)
    fastq4_find=(`find /workspace/storage/krue/Methylation/extra6lanes/20150209_D00731B_BSSeq_PE -name "$sampleShortName"_[ATGC]*_R1_001.fastq.gz`)
    #echo 'fastq4_find '$fastq4_find
    #echo ${#fastq4_find[@]}
    if [ "${#fastq3_find[@]}" -eq "2" ]; then
        #echo 'Sample '$sampleShortName' found in run 3.'
        fastq3_1=${fastq3_find[0]}
        fastq4_1=${fastq3_find[1]}
    elif [ "${#fastq4_find[@]}" -eq "2" ]; then
        #echo 'Sample '$sampleShortName' found in run 4.'
        fastq3_1=${fastq4_find[0]}
        fastq4_1=${fastq4_find[1]}
    else
        echo 'Error. Sample '$sampleShortName' not found in either run 3 and 4.'
        break
    fi
    echo $fastq3_1
    echo $fastq4_1
    # reverse read of third lane
    fastq3_2=`echo $fastq3_1 | perl -p -e 's/R1_001.fastq.gz$/R2_001.fastq.gz/'` 
    # reverse read of fourth lane
    fastq4_2=`echo $fastq4_1 | perl -p -e 's/R1_001.fastq.gz$/R2_001.fastq.gz/'`
    echo 'bwameth.py -p '$sampleShortName' --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa '$fastq1_1','$fastq2_1','$fastq3_1','$fastq4_1' '$fastq1_2','$fastq2_2','$fastq3_2','$fastq4_2' -t 12 > '$sampleShortName'.bwameth.out' >> bwameth.sh
done

# Split the master script in subscripts of N lines
split -d -l 12 bwameth.sh bwameth.sh.

chmod +x bwameth.sh.*

for file in `ls bwameth.sh.*`
do
    nohup ./$file > ${file}.nohup &
done


###
# Count aligned reads
###

# Currently blocked here as Picard tools does not seem to detect adapter sequences

# See for samtools http://left.subtree.org/2012/04/13/counting-the-number-of-reads-in-a-bam-file/

# See for Picard tools: http://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics

# Adapter sequences
# Sequence of adapter ligated at 3' end
# GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT
# Reverse complement of adapter ligated at 3' end (expected in forward read)
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC      TruSeq_Adapter  Read1
# Sequence of adapter ligated at 5' end
# TACACTCTTTCCCTACACGACGCTCTTCCGATCT
# Reverse complement of adapter ligated at 5' end (expected in reverse read)
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA      Illumina_Single_End_PCR_Primer_1        Read2

# Proof that there is adapter (3489 read pairs / 100,000 = 3.4%)
samtools view C12.bam | head -n 100000 | grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC | wc -l
samtools view C12.bam | head -n 100 | grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# 10 times more reads, for assurance of stats: 35014 / 1,000,000 = 3.5 %
samtools view C12.bam | head -n 1000000 | grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC | wc -l


# Test with both reverse complement (%(1) = 0.000002 adapters detected)
java -jar /usr/local/src/picard-tools-1.128/picard.jar CollectAlignmentSummaryMetrics INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.test_picard_summary.txt ADAPTER_SEQUENCE=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ADAPTER_SEQUENCE=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA IS_BISULFITE_SEQUENCED=true STOP_AFTER=1000000

# Test with both adapter sequence (%(1) = 0.000002 adapters detected)
java -jar /usr/local/src/picard-tools-1.128/picard.jar CollectAlignmentSummaryMetrics INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.test_picard_summary.txt ADAPTER_SEQUENCE=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT ADAPTER_SEQUENCE=TACACTCTTTCCCTACACGACGCTCTTCCGATCT IS_BISULFITE_SEQUENCED=true STOP_AFTER=1000000

# Test removing default adapter sequences with both adapter sequence (%(1) = 0.000002 adapters detected)
java -jar /usr/local/src/picard-tools-1.128/picard.jar CollectAlignmentSummaryMetrics INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.test_picard_summary.txt ADAPTER_SEQUENCE=null ADAPTER_SEQUENCE=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT ADAPTER_SEQUENCE=TACACTCTTTCCCTACACGACGCTCTTCCGATCT IS_BISULFITE_SEQUENCED=true STOP_AFTER=100000

# Test removing default adapter sequences with both adapter sequence (%(1) = 0.000002 adapters detected)
java -jar /usr/local/src/picard-tools-1.128/picard.jar CollectAlignmentSummaryMetrics INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.test_picard_summary.txt ADAPTER_SEQUENCE=null ADAPTER_SEQUENCE=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ADAPTER_SEQUENCE=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA IS_BISULFITE_SEQUENCED=true STOP_AFTER=100000

# File seems valid
java -jar /usr/local/src/picard-tools-1.128/picard.jar ValidateSamFile INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam IS_BISULFITE_SEQUENCE=true

# Test removing default adapter sequences with both adapter and reverse complement (%(1) = 0. adapters detected)
java -jar /usr/local/src/picard-tools-1.128/picard.jar CollectAlignmentSummaryMetrics INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.test_picard_summary.txt ADAPTER_SEQUENCE=null ADAPTER_SEQUENCE=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT ADAPTER_SEQUENCE=TACACTCTTTCCCTACACGACGCTCTTCCGATCT ADAPTER_SEQUENCE=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ADAPTER_SEQUENCE=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA IS_BISULFITE_SEQUENCED=true STOP_AFTER=100000

# Test adding reference genome sequence - removing default adapter sequences with both adapter and reverse complement (%(1) = 0. adapters detected)
# Most statistics are OK now that genome is specified
# However, the adapters are still 0% instead of the expected 3%, approximately
java -jar /usr/local/src/picard-tools-1.128/picard.jar CollectAlignmentSummaryMetrics INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.test_picard_summary.txt ADAPTER_SEQUENCE=null ADAPTER_SEQUENCE=GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT ADAPTER_SEQUENCE=TACACTCTTTCCCTACACGACGCTCTTCCGATCT ADAPTER_SEQUENCE=AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC ADAPTER_SEQUENCE=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA IS_BISULFITE_SEQUENCED=true REFERENCE_SEQUENCE=/workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa STOP_AFTER=100000

mkdir /workspace/scratch/krue/Methylation/bwa_8lanes/tmp_picard

# Another Picard module to generate multiple stats from BAM file
# Generate tables and plots too !
# Although adapter sequences are not mentioned in the help as possible arguments
java -jar /usr/local/src/picard-tools-1.128/picard.jar CollectMultipleMetrics INPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa_8lanes/C12.test_picard_summary.txt REFERENCE_SEQUENCE=/workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa STOP_AFTER=100000 TMP_DIR=/workspace/scratch/krue/Methylation/bwa_8lanes/tmp_picard 


##
# Bias plot (bwameth)
##

# Paper on BioRXiv for bwameth
http://arxiv.org/pdf/1401.1129v2.pdf

# Extra grep to avoid a test file containing the word "top"
for bam in `find /workspace/scratch/krue/Methylation/bwa_8lanes -name '*\.bam' | grep -v 'top'`
do
    echo 'python /usr/local/src/bwa-meth-0.10/bias-plot.py '$bam' /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa' >> bias_plot.sh
done

chmod +x bias_plot.sh

./bias_plot.sh
# run the above command
# when the terminal connection crashes, the PNG images are not generated
# Then comment out the lines where the plots were generated, and re-run

###
# Bias plot (PileOMeth)
###

nohup PileOMeth mbias /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa C12.bam C12 > C12.PileOMeth.nohup &

# Extra grep to avoid a test file containing the word "top"
for bam in `find /workspace/scratch/krue/Methylation/bwa_8lanes -name '*\.bam' | grep -v 'top'`
do
    sample=`basename $bam | perl -p -e 's/.bam$//'`
    echo 'PileOMeth mbias /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa '$bam' '$sample >> bias_PileOMeth.sh
done

# Split the master script in subscripts of N lines
split -d -l 2 bias_PileOMeth.sh bias_PileOMeth.sh.

chmod +x bias_PileOMeth.sh.*

for file in `ls bias_PileOMeth.sh.*`
do
    nohup ./$file > ${file}.nohup &
done

###
# Summarise the suggested regions 
###

# Situation
# The suggested regions are outputted in the nohup log file without information about the corresponding sample
# while the corresponding sample can be found in the script submitted to nohup
# Solution:
# Give the list of scripts submitted to nohup
# assume that each line is a command
# fetch the sample name as the last word of the command
# assume that the corresponding nohup log file is named <script_name.nohup>
# assume that each line contains a suggested inclusion region, in the same order as the script
# fetch the suggested regions to include

# List all the PileOMeth script files in a file
ls bias_PileOMeth.sh.* | grep -v nohup > bias_PileOMeth_scripts.txt

# Summarise all the suggested regions in one simple table
python /home/krue/GitHub/NGS/PileOMeth_summary.py -l /workspace/scratch/krue/Methylation/bwa_8lanes/bias_PileOMeth_scripts.txt -o bias_PileOMeth_summary_test.txt -i NOT_BS
# Note that a bonus feature of the script would be to also output
# the max and min boundaries across all samples

# Look at the summarised statistics to decide the regions kept during the calling step
less bias_PileOMeth_summary_test.txt
# I took 1-2 bases more for some extra margin (based on the SVG graphs)
# OT: 10,0,0,135
# OB: 0,140,15,0


###
# CpG methylation calls
###

# Test line (not run)
#nohup PileOMeth extract /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa C12.bam --OT 10,0,0,135 --OB 0,140,15,0 -o C12 > C12.PileOMethExtract.nohup &

# Extra grep to avoid a test file containing the word "top"
# 10 trimmed from the start and 15 from the end
for bam in `find /workspace/scratch/krue/Methylation/bwa_8lanes -name '*\.bam' | grep -v 'top'`
do
    sample=`basename $bam | perl -p -e 's/.bam$//'`
    echo 'PileOMeth extract /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa '$bam' --OT 10,0,0,135 --OB 0,140,15,0 -o '$sample >> extract_PileOMeth.sh
done

# Split the master script in subscripts of N lines
split -d -l 2 extract_PileOMeth.sh extract_PileOMeth.sh.

chmod +x extract_PileOMeth.sh.*

for file in `ls extract_PileOMeth.sh.*`
do
    nohup ./$file > ${file}.nohup &
done


###
# The rest should be in R (see Methylation project)
###

# Results are fairly different from the preliminary data.
# Maybe we changed too many steps of the analysis.

###
# Tabulate methylation (bwameth)
###

#/usr/local/src/bwa-meth-0.09/bwameth.py tabulate âˆ’âˆ’trim 3 ,3 âˆ’âˆ’map-q 30 âˆ’âˆ’bissnp BisSNPâˆ’0.82.2.jar âˆ’âˆ’reference /path/to/ref.fasta input.bam # do for each file
