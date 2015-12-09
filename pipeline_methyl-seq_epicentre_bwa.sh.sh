# Adapters from Epicentre kit
# Expected in R1 (rc_adapter2 half before barcode): AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# Expected in R2 (rc_adapter1): AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# Folder for this analysis
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis

###
# Quality control of raw reads
###

# Quality control of raw files was done previously (/workspace/scratch/krue/Methylation/qc_MSU_data/Paired/)

###
# Filtering and trimming
###

cd /workspace/scratch/krue/Methylation/2014-08-27_analysis
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic

# Created a file describing the expected contaminant adapter sequences (rc)
# Search for rc of adapter 1 in mate1 and vice versa
nano adapters.fa
# >rc_adapter2/1
# AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# >rc_adapter1/2
# AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT

# test ona single sample
java -jar /usr/local/src/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -threads 8 -trimlog C11_CAGATC.log -baseout C11_CAGATC.fastq /workspace/storage/krue/Methylation/Paired/C11_CAGATC_L002_R1_001.fastq.gz /workspace/storage/krue/Methylation/Paired/C11_CAGATC_L002_R2_001.fastq.gz LEADING:30 HEADCROP:6 ILLUMINACLIP:/workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/adapters.fa:2:40:15 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:16 > C11_CAGATC.stdout
# TrimmomaticPE: Started with arguments: -phred33 -threads 8 -trimlog C11_CAGATC.log -baseout C11_CAGATC.fastq /workspace/storage/krue/Methylation/Paired/C11_CAGATC_L002_R1_001.fastq.gz /workspace/storage/krue/Methylation/Paired/C11_CAGATC_L002_R2_001.fastq.gz LEADING:30 HEADCROP:6 ILLUMINACLIP:/workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/adapters.fa:2:40:15 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:16
# Using Long Clipping Sequence: 'AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC'
# Using Long Clipping Sequence: 'AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT'
# ILLUMINACLIP: Using 0 prefix pairs, 0 forward/reverse sequences, 1 forward only sequences, 1 reverse only sequences
# Input Read Pairs: 11797190 Both Surviving: 11559470 (97.98%) Forward Only Surviving: 221326 (1.88%) Reverse Only Surviving: 7744 (0.07%) Dropped: 8650 (0.07%)
# TrimmomaticPE: Completed successfully

# Clenup
rm C11_CAGATC*


# Prepare script 
# adapted from the Epicentre guidebook (http://www.epibio.com/applications/bisulfite-sequencing/epignome-methyl-seq-kit?details)
# ILLUMINACLIP put before TRAILING as trimming sometimes left too short sequences
for file in `find /workspace/storage/krue/Methylation/Paired -name *R1_001.fastq.gz`; do file2=`echo $file | perl -p -e 's/R1(_001\.fastq.gz)$/R2$1/'`; baseout=`basename $file | perl -p -e 's/_L00[12]_R1_001\.fastq\.gz$//'`; echo "java -jar /usr/local/src/Trimmomatic-0.32/trimmomatic-0.32.jar PE -phred33 -threads 1 -trimlog $baseout.log -baseout $baseout.fastq $file $file2 LEADING:30 HEADCROP:6 ILLUMINACLIP:/workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/adapters.fa:2:40:15 TRAILING:30 SLIDINGWINDOW:4:15 MINLEN:16 > $baseout.stdout" >> trimmomatic.sh; done;


# Split and run all scripts on Stampede
split -d -l 1 trimmomatic.sh trimmomatic.sh.

# Make those subscripts executable
chmod +x trimmomatic.sh.*

# Run the 24 subscripts in parallel
for file in `ls trimmomatic.sh.*`
do
    nohup ./$file > ${file}.nohup &
done

# Check that all 24 pairs of files were fully processed (command successfully returns 24)
cat trimmomatic.sh.*.nohup | grep "TrimmomaticPE: Completed successfully" | c

# Compile all deconvolution reports into a single file
ls *.nohup > trimmomatic_reports.txt
python3 /home/krue/GitHub/NGS/Trimmomatic_summary.py --list-reports trimmomatic_reports.txt --output all_trimmomatics.txt --pattern '[MC]\d{1,2}.{0,7}?_[ATGC]{6}'

# Over 93% of read pairs are left

# Let's check that adapter sequences were removed (e.g. rc_adapter2 in mate1)
grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC C7_NOT_BS_TTAGGC_1P.fastq
# Make sure that the sequence was present before filtering
zcat /workspace/storage/krue/Methylation/Paired/C7_NOT_BS_TTAGGC_L002_R1_001.fastq.gz | head -n 100 | grep AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
# Note that this is the true adapter, as we can systematically see the barcode following it

# I forgot to create a subfolder for Paired data
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired
mv ./* /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired

# move fastq files of orphan reads to a subfolder and compress them(clean up)
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/Orphans
mv *1U.fastq *2U.fastq /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/Orphans
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/Orphans
for file in `ls *U*`; do echo gzip $file >> compress.sh ; done;
split -d -l 2 compress.sh compress.sh.
chmod +x compress.sh.*
for file in `ls compress.sh.*`
do
    nohup ./$file > ${file}.nohup &
done
rm compress.sh.*

###
# Quality control post filtering
###

# Create a directory for the fastqc output files and move into it
mkdir -p /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired

# Write commands to run fastqc on all 64 fastq.gz files in a master script
for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired -name *P.fastq`; do fileout=`basename $file | perl -p -e 's/(.*)P.fastq$/$1\.fastqc\.out/'`;echo "fastqc --noextract --nogroup -t 1 -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired $file > $fileout" >> fastqc.sh; done;

# Split the master script into 24 subscripts of 2 commands each
split -d -l 2 fastqc.sh fastqc.sh.

# Make those subscripts executable
chmod +x fastqc.sh.*

# Run the 22 subscripts in parallel
for file in `ls fastqc.sh.*`
do
    nohup ./$file > ${file}.nohup &
done

# Remove fastqc and nohup log files after checking that the fastqc log files all stated "analysis complete" 
cat *.out | grep "Analysis complete" | wc -l # 48
rm *.out
rm *.nohup

# Check all output from FastQC particularly for the level of reads duplication
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired/tmp
for file in `ls *_fastqc.zip`; do unzip $file -d /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired/tmp; done;
for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired/tmp -name summary.txt`; do cat $file | grep "FAIL" >> fastqc_MSU_data_fail.txt; done;
for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired/tmp -name summary.txt`; do cat $file | grep "WARN" >> fastqc_MSU_data_warning.txt; done;
rm -r /workspace/scratch/krue/Methylation/2014-08-27_analysis/qc_trimmomatic/Paired/tmp

# Count how many of each type of FAIL we had
cut -f2 fastqc_MSU_data_fail.txt | sort | uniq -c
# 32 Kmer Content
# 32 Per base sequence content

# Check which samples correspond to those (32 corresponds to number of BS-treated samples)
grep "Kmer Content" fastqc_MSU_data_fail.txt
grep "Per base sequence content" fastqc_MSU_data_fail.txt

# Count how many of each type of WARN(ING) we had
cut -f2 fastqc_MSU_data_warning.txt | sort | uniq -c
# 16 Kmer Content
# 1 Overrepresented sequences
# 23 Per base GC content
# 16 Per base sequence content
# 28 Per sequence GC content # only BS-treated samples: is expected
# 20 Sequence Duplication Levels 
# 48 Sequence Length Distribution # because different length
grep "Per sequence GC content" fastqc_MSU_data_warning.txt

###
# Bismark alignment
###

# Create a folder for the output of the bismark analysis
mkdir -p /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired/tmp

# default (= no) minimum fragment size
# do not print out unmapped and ambiguous pairs
# outputs files in a specific folder and uses another folder for temporary files
# 1 thread per bowtie2 alignment (2 alignments per sample). -p cannot be specified with value less than 2. Omitted.

for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired -path '*training*' -prune -o -name *1P.fastq | grep -v training`; do file2=`echo $file | perl -p -e 's/(.*)_1P\.fastq$/$1_2P\.fastq/'`; fileout=`basename $file | perl -p -e 's/(.*)P.fastq$/$1\.bismark\.out/'`;echo "bismark --bowtie2 --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/ -1 $file -2 $file2 > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired/$fileout" >> bismark.sh; done;

# Split the master script into 6 scripts of 4 commands
# (each command will run two bowtie2 jobs because it assumes directional libraries = strand specific)
split -d -l 2 bismark.sh bismark.sh.

# Make those subscripts executable
chmod +x bismark.sh.*

# Per manual instruction the current directory must contain the read files when launching bismark
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired

# Run the 6 subscripts in parallel
for file in `ls /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired/bismark.sh.*`
do
    nohup $file > ${file}.nohup &
done

# Very low alignment rate (30%). Could it be:
# too short insert max size (500 nts?)
# non-directional libraries

# Check the mapping efficiency report for each sample
grep efficiency *report*
# C11_CAGATC_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:       35.9%
# C11_NOT_BS_GATCAG_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        49.9%
# C12_NOT_BS_CTTGTA_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        51.9%
# C12_TAGCTT_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:       45.1%
# C1_ATCACG_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        28.4%
# C1_NOT_BS_TTAGGC_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency: 54.1%
# C4_NOT_BS_GCCAAT_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency: 49.3%
# C4_TGACCA_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        30.9%
# C5_CAGATC_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        29.4%
# C5_NOT_BS_GATCAG_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency: 52.4%
# C6_NOT_BS_CTTGTA_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency: 53.5%
# C6_TAGCTT_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        36.5%
# C7_ATCACG_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        41.4%
# C7_NOT_BS_TTAGGC_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency: 52.5%
# C8_NOT_BS_GCCAAT_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency: 55.5%
# C8_TGACCA_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        29.4%
# M11_ACTTGA_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:       41.2%
# M12_GGCTAC_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:       28.4%
# M1_CGATGT_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        41.4%
# M4_ACAGTG_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        35.3%
# M5_ACTTGA_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        30.7%
# M6_GGCTAC_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        44.6%
# M7_CGATGT_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        34.0%
# M8_ACAGTG_1P.fastq_bismark_bt2_PE_report.txt:Mapping efficiency:        40.3%

grep "Number of paired-end alignments with a unique best hit:" *report*



###
# Test zone
###

# Takes samples C1 BS-treated
# Align it in non-directional mode
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/non_directional
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/non_directional
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/non_directional/tmp
# Per manual instruction the current directory must contain the read files when launching bismark
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired
nohup bismark --bowtie2 --non_directional --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/non_directional/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/non_directional /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/ -1 /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_1P.fastq -2 /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_2P.fastq > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/non_directional/C1_ATCACG.out &

# Align it using max insert size of 1,000 nucleotides using 4 cores per alignment (8 total)
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/insert_1000
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/insert_1000
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/insert_1000/tmp
# Per manual instruction the current directory must contain the read files when launching bismark
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired
nohup bismark --bowtie2 --maxins 1000 --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/insert_1000/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/insert_1000 /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/ -1 /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_1P.fastq -2 /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_2P.fastq > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/insert_1000/C1_ATCACG.out &


# Not a massive improvement
# We will continue with the initial results, even if low for some samples

###
# Remove PCR duplicates (like we needed to lose more data...)
###

mkdir -p /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates

# Converts sam file to bam 
mkdir -p /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/bam_input
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/bam_input

for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired -name *sam`; do fileout=`basename $file | perl -p -e 's/(.*)sam$/$1bam/'`; report=`basename $file | perl -p -e 's/(.*)sam$/$1out/'`; echo "samtools view -S $file -b -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/bam_input/$fileout > $report" >> sam2bam.sh; done;

split -d -l 1 sam2bam.sh sam2bam.sh.
chmod +x sam2bam.sh.*
for file in `ls /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/bam_input/sam2bam.sh.*`
do
    nohup $file > ${file}.nohup &
done

# Removes PCR duplicates
mkdir -p /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/duplicate_removed
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/duplicate_removed

for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/bam_input -name *bam`; do fileout=`basename $file | perl -p -e 's/(.*)bam$/noDup_$1bam/'`; report=`basename $file | perl -p -e 's/(.*)bam$/$1out/'`; echo "samtools rmdup -S $file /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/duplicate_removed/$fileout > $report" >> remove_duplicates.sh; done;

split -d -l 1 remove_duplicates.sh remove_duplicates.sh.
chmod +x remove_duplicates.sh.*
for file in `ls /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/duplicate_removed/remove_duplicates.sh.*`
do
    nohup $file > ${file}.nohup &
done

# Converts bam file back to sam file (can be used in next step, methylation calls) 
mkdir -p /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/sam_output
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/sam_output

for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/duplicate_removed -name *bam`; do fileout=`basename $file | perl -p -e 's/(.*)bam$/$1sam/'`; report=`basename $file | perl -p -e 's/(.*)bam$/$1out/'`; echo "samtools view $file -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/sam_output/$fileout > $report" >> bam2sam.sh; done;

split -d -l 1 bam2sam.sh bam2sam.sh.
chmod +x bam2sam.sh.*
for file in `ls /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/sam_output/bam2sam.sh.*`
do
    nohup $file > ${file}.nohup &
done


###
# Methylation extractor
###

mkdir -p /workspace/scratch/krue/Methylation/2014-08-27_analysis/methyl_extracted
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/methyl_extracted

for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/PCRduplicates/sam_output -name '*sam' | grep -v NOT_BS`; do echo "bismark_methylation_extractor -s --comprehensive --scaffolds --bedGraph  --counts  --buffer_size  10G  --cytosine_report  --genome_folder  /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file $file" >> methylation.sh; done;
# Epicentre uses the single-end option.
# I tried using the paired-end, but it crashes because mates are not correctly sorted anymore it seems
# This was likely introduced by the remove duplicate step
# Anyway, I suppose that paired-end was most useful to correctly align samples

# Split the master script into 6 scripts of 4 commands
split -d -l 1 methylation.sh methylation.sh.

# Make those subscripts executable
chmod +x methylation.sh.*

# Run the 16 subscripts in parallel
for file in `ls methylation.sh.*`
do
    nohup $file > ${file}.nohup &
done

###
# Repeat alignment with truncated reads aligned in single-end mode
###

mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75/tmp

# truncate the forward reads of the sample with worst alignment rate to max 75 nts.
cat /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_1P.fastq | cut -c 1-75 > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75/C1_ATCACG_1P_truncated75.fastq

# Per manual instruction the current directory must contain the read files when launching bismark
nohup bismark --bowtie2 --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75 /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75/C1_ATCACG_1P_truncated75.fastq > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75/C1_ATCACG.out &

###
# Repeat alignment with 1 million full-length reads
###

mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/1million
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/1million
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/1million/tmp

# Per manual instruction the current directory must contain the read files when launching bismark
nohup bismark --bowtie2 -u 1000000 --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/1million/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/1million /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_1P.fastq > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/1million/C1_ATCACG.out &

###
# Repeat alignment with 1 million 75bp paired-end reads
###

mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/tmp

# truncate the forward reads of the sample with worst alignment rate to max 75 nts.
cat /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_2P.fastq | cut -c 1-75 > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/C1_ATCACG_2P_truncated75.fastq

# Per manual instruction the current directory must contain the read files when launching bismark
nohup bismark --bowtie2 -u 1000000 --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file -1 /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75/C1_ATCACG_1P_truncated75.fastq -2 /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/C1_ATCACG_2P_truncated75.fastq > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/C1_ATCACG.out &

###
# Repeat alignment with full-length forward read
###

rm -r /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/1million
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/single_end
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/single_end
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/single_end/tmp

# Per manual instruction the current directory must contain the read files when launching bismark
nohup bismark --bowtie2 --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/single_end/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/single_end /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_1P.fastq > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/single_end/C1_ATCACG.out &

###
# Repeat alignment with 75bp paired-end reads
###

rm -r /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated
mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/tmp

# Per manual instruction the current directory must contain the read files when launching bismark
nohup bismark --bowtie2 --temp_dir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/tmp -o /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file -1 /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75/C1_ATCACG_1P_truncated75.fastq -2 /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/C1_ATCACG_2P_truncated75.fastq > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated/C1_ATCACG.out &

###
# Alignment using bwa-meth
###

mkdir /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/bwa_meth
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/bwa_meth

# Index the genome 
nohup bwameth.py index /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa > bwa-meth-index.out &

# repeated  the following steps after installing the most recent version of bwa
# using the old version (0.7.2) caused crash in downstream step due to CIGAR string starting with deletion
# hopefully 0.7.10 solves that issue
rm ex* C* bwa-meth.bam*

# Align my favorite test sample
nohup bwameth.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_1P.fastq /workspace/scratch/krue/Methylation/2014-08-27_analysis/Trimmomatic/Paired/C1_ATCACG_2P.fastq -t 12 -p C1_ATCACG > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/bwa_meth/C1_ATCACG.out &

samtools flagstat C1_ATCACG.bam

python /home/krue/GitHub/bwa-meth/bias-plot.py C1_ATCACG.bam /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa

# Script by the author of bwa-meth using Bis-SNP
nohup bwameth.py tabulate \
            # removed trimming
             --map-q 60 --bissnp /usr/local/bin/BisSNP-0.82.2.jar \
             --prefix C1_ATCACG \
             -t 12 \
             --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa \
             C1_ATCACG.bam > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/bwa_meth/C1_ATCACG.tabulate.out &


# 17570633 + 0 in total (QC-passed reads + QC-failed reads)
# 0 + 0 duplicates
# 17493812 + 0 mapped (99.56%:-nan%)
# 17570633 + 0 paired in sequencing
# 8787708 + 0 read1
# 8782925 + 0 read2
# 16950242 + 0 properly paired (96.47%:-nan%)
# 17470980 + 0 with itself and mate mapped
# 22832 + 0 singletons (0.13%:-nan%)
# 209437 + 0 with mate mapped to a different chr
# 66124 + 0 with mate mapped to a different chr (mapQ>=5

# Python script by the author of bwa-meth
nohup python /home/krue/GitHub/bwa-meth/scripts/tabulate-methylation.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa --read-length 144 C1_ATCACG.bam > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/bwa_meth/C1_ATCACG.tabulate-py.out &

## NOT RUN yet ###
nohup bwameth.py tabulate \
             --trim 6,6 \
             --map-q 60 --bissnp /usr/local/bin/BisSNP-0.82.2.jar \
             --prefix C1_ATCACG \
             -t 12 \
             --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa \
             C1_ATCACG.bam > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/bwa_meth/C1_ATCACG.tabulate.out &

## NOT RUN yet ###
nohup python /home/krue/GitHub/bwa-meth/scripts/tabulate-methylation.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa --read-length 144 --skip-left 6 --skip-right 6 C1_ATCACG.bam > /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/bwa_meth/C1_ATCACG.tabulate-py.out &


###
# Simpler-is-the-killer solution? Raw reads using bwa-meth
###

mkdir /workspace/scratch/krue/Methylation/bwa-meth
cd /workspace/scratch/krue/Methylation/bwa-meth

nohup bwameth.py -p C1_ATCACG --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa /workspace/storage/krue/Methylation/Paired/C1_ATCACG_L001_R1_001.fastq.gz /workspace/storage/krue/Methylation/Paired/C1_ATCACG_L001_R2_001.fastq.gz -t 12 > C1_ATCACG.bwameth.out &

nohup samtools flagstat C1_ATCACG.bam > C1_ATCACG.flagstat.out &
nohup samtools view -f2 -F0x200 C1_ATCACG.bam | awk '$5 > 15' | wc -l >> C1_ATCACG.flagstat.out &

nohup python /home/krue/GitHub/bwa-meth/bias-plot.py C1_ATCACG.bam /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa > C1_ATCACG.bias-plot.out &

# Given the bias plot, I think we should trim 20 bases from each end
# because the first 20 of read1 and the last 20 of read 2 have a %GC methylation drifting from the trend.
nohup python /home/krue/GitHub/bwa-meth/scripts/tabulate-methylation.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa --read-length 144 --skip-left 20 --skip-right 20 C1_ATCACG.bam > C1_ATCACG.tabulate.out &

###
# Cleanup
###

# Compress all the Trimmomatic output fastq files in the bismark pipeline
cd /home/krue/scratch/Methylation/2014-08-27_analysis/Trimmomatic/Paired
for file in `find /home/krue/scratch/Methylation/2014-08-27_analysis/Trimmomatic/Paired -name '*.fastq'`; do echo "gzip $file >> gzip.out" >> gzip_fastq.sh; done;
chmod +x gzip_fastq.sh
./gzip_fastq.sh &

# Remove the SAM files obtained after removing PCR duplicates
# Code is available in this script to regenerate those files in one step
rm -r /home/krue/scratch/Methylation/2014-08-27_analysis/PCRduplicates/sam_output

# Convert test SAM files to BAM files
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/insert_1000/
samtools view -S C1_ATCACG_1P.fastq_bismark_bt2_pe.sam -b -o C1_ATCACG_1P.fastq_bismark_bt2_pe.bam > sam2bam.out &
# When the above command is done:
rm C1_ATCACG_1P.fastq_bismark_bt2_pe.sam
# More:
cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired_truncated
samtools view -S C1_ATCACG_1P_truncated75.fastq_bismark_bt2_pe.sam -b -o C1_ATCACG_1P_truncated75.fastq_bismark_bt2_pe.bam > sam2bam.out &
rm C1_ATCACG_1P_truncated75.fastq_bismark_bt2_pe.sam

cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/single_end
samtools view -S C1_ATCACG_1P.fastq_bismark_bt2.sam -b -o C1_ATCACG_1P.fastq_bismark_bt2.bam > sam2bam.out &
rm C1_ATCACG_1P.fastq_bismark_bt2.sam

cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/truncated_75
samtools view -S C1_ATCACG_1P_truncated75.fastq_bismark_bt2.sam -b -o C1_ATCACG_1P_truncated75.fastq_bismark_bt2.bam > sam2bam.out &
rm C1_ATCACG_1P_truncated75.fastq_bismark_bt2.sam

cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/non_directional
samtools view -S C1_ATCACG_1P.fastq_bismark_bt2_pe.sam -b -o C1_ATCACG_1P.fastq_bismark_bt2_pe.bam > sam2bam.out &
rm C1_ATCACG_1P.fastq_bismark_bt2_pe.sam

cd /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired
for file in `find /workspace/scratch/krue/Methylation/2014-08-27_analysis/bismark/paired -name '*.sam'`; do fileOut=`echo $file | perl -p -e 's/sam$/bam/'`; echo "samtools view -S $file -b -o $fileOut >> sam2bam.out; rm $file" >> sam2bam.sh; done;
chmod +x sam2bam.sh
nohup ./sam2bam.sh &



###
# Repeat the above strategy for all samples
###

mkdir /workspace/scratch/krue/Methylation/bwa-meth-all
cd /workspace/scratch/krue/Methylation/bwa-meth-all

for file in `find /workspace/storage/krue/Methylation/Paired -name *R1_001.fastq.gz | grep -v NOT_BS`; do file2=`echo $file | perl -p -e 's/R1_001.fastq.gz$/R2_001.fastq.gz/'`; sample=`basename $file | perl -p -e 's/_[ATGC]{6}.*$//'`; folderOut=/workspace/scratch/krue/Methylation/bwa-meth-all/$sample; echo "mkdir $folderOut; cd $folderOut; bwameth.py -p $sample --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa $file $file2 -t 12 > $sample.bwameth.out" >> bwameth.sh; done;


split -d -l 8 bwameth.sh bwameth.sh.
chmod +x bwameth.sh.*

for file in `ls bwameth.sh.*`
do
    nohup ./$file > ${file}.nohup &
done

###
# Check the number of read pairs aligned with quality score > 15
###

for file in `find /workspace/scratch/krue/Methylation/bwa-meth-all -name *.bam`; do sample=`basename $file | perl -p -e 's/.bam$//'`; echo $sample; reads=`samtools view -f2 -F0x200 $file | awk '$5 > 15' | wc -l`; echo $reads; read_pairs=`expr $reads / 2`; echo $read_pairs; echo -e "$sample\t$read_pairs" >> countaligned_q15.txt; done &

###
# Check the bias-plot of all bwa-meth alignments in previous section
###

cd /workspace/scratch/krue/Methylation/bwa-meth-all

for file in `find /workspace/scratch/krue/Methylation/bwa-meth-all -name *.bam `; do sample=`basename $file | perl -p -e 's/.bam$//'`; fileIn=`basename $file`; folder=`dirname $file`; echo "cd $folder; python /home/krue/GitHub/bwa-meth/bias-plot.py $fileIn /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa > $sample.bias-plot.out" >> biasplot.sh; done;


split -d -l 1 biasplot.sh biasplot.sh.
chmod +x biasplot.sh.*

# without nohup, this requires the terminal not to crash for about 10min (10 million paired-read)
for file in `ls biasplot.sh.*`
do
    ./$file > ${file}.sdout & # nohup causes crash prior to the PNG output
done

###
# Call methylation levels, while ignoring the first 20 and last 20 nucleotides of each read (pair?)
###
cd /workspace/scratch/krue/Methylation/bwa-meth-all

for file in `find /workspace/scratch/krue/Methylation/bwa-meth-all -name *.bam `; do sample=`basename $file | perl -p -e 's/.bam$//'`; fileIn=`basename $file`; folder=`dirname $file`; echo "cd $folder; python /home/krue/GitHub/bwa-meth/scripts/tabulate-methylation.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa --read-length 144 --skip-left 20 --skip-right 20 $fileIn > $sample.tabulate.out" >> tabulate.sh; done;
nohup python /home/krue/GitHub/bwa-meth/scripts/tabulate-methylation.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa --read-length 144 --skip-left 20 --skip-right 20 C1_ATCACG.bam > C1_ATCACG.tabulate.out &

split -d -l 1 tabulate.sh tabulate.sh.
chmod +x tabulate.sh.*

for file in `ls tabulate.sh.*`
do
    nohup ./$file > ${file}.sdout &
done


###
# Merging BAM files by groups of 4
###

# The low coverage of cytosines in individual samples is an issue for the statistical significance of methylation calls.
# Merging the reads of multiple samples into "master-samples" may allow better calling of methylation levels in each master-sample
# while reducing the accuracy of the biological variability estimate

# A randomisation script in R assigned the following 2 groups of 4 animals:
# Group     Element
# 1     1 (5, 1, 7, 11)
# 2     2 (12, 8, 6, 4)

cd /workspace/scratch/krue/Methylation/bwa-meth-all

# Using picard-tools
# I will merge the BAM files of the animals in each group, for control and infected samples separately
mkdir /workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped
nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C1/C1.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C5/C5.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C7/C7.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C11/C11.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/C_1-5-7-11.bam > C_1-5-7-11.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C4/C4.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C6/C6.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C8/C8.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C12/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/C_4-6-8-12.bam > C_4-6-8-12.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M1/M1.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M5/M5.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M7/M7.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M11/M11.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/M_1-5-7-11.bam > M_1-5-7-11.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M4/M4.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M6/M6.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M8/M8.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M12/M12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/M_4-6-8-12.bam > M_4-6-8-12.merge.nohup &

###
# Call methylation levels of 4-sample master samples
###

cd /workspace/scratch/krue/Methylation/bwa-meth-all

for file in `find /workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped -name *-*-*-*.bam `; do sample=`basename $file | perl -p -e 's/.bam$//'`; fileIn=`basename $file`; folder=`dirname $file`; echo "cd $folder; python /home/krue/GitHub/bwa-meth/scripts/tabulate-methylation.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa --read-length 144 --skip-left 20 --skip-right 20 $fileIn > $sample.tabulate.out" >> tabulate_groupsOf4.sh; done;

split -d -l 1 tabulate_groupsOf4.sh tabulate_groupsOf4.sh.
chmod +x tabulate_groupsOf4.sh.*

for file in `ls tabulate_groupsOf4.sh.*`
do
    nohup ./$file > ${file}.sdout &
done


###
# Merging BAM files by groups of 2
###

cd /workspace/scratch/krue/Methylation/bwa-meth-all

# A randomisation script in R assigned the following 4 pairs of samples:
# (5, 8) (11, 12) (4, 6) (1, 7)
# control samples
nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C5/C5.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C8/C8.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/C_5-8.bam > C_5-8.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C11/C11.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C12/C12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/C_11-12.bam > C_11-12.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C4/C4.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C6/C6.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/C_4-6.bam > C_4-6.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C1/C1.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/C7/C7.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/C_17.bam > C_1-7.merge.nohup &

# infected samples
nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M5/M5.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M8/M8.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/M_5-8.bam > M_5-8.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M11/M11.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M12/M12.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/M_11-12.bam > M_11-12.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M4/M4.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M6/M6.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/M_4-6.bam > M_4-6.merge.nohup &

nohup java -jar /usr/local/src/picard-tools-1.119/MergeSamFiles.jar INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M1/M1.bam INPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/M7/M7.bam OUTPUT=/workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped/M_17.bam > M_1-7.merge.nohup &


###
# Call methylation levels of 2-sample master samples (TODO later maybe)
###

cd /workspace/scratch/krue/Methylation/bwa-meth-all

for file in `find /workspace/scratch/krue/Methylation/bwa-meth-all/BAM_grouped -name [CM]_*-*.bam | grep -v .*-.*-.*-.*-.*.bam`; do sample=`basename $file | perl -p -e 's/.bam$//'`; fileIn=`basename $file`; folder=`dirname $file`; echo "cd $folder; python /home/krue/GitHub/bwa-meth/scripts/tabulate-methylation.py --reference /workspace/storage/genomes/bostaurus/UMD3.1.75/source_file/Bos_taurus.UMD3.1.75.dna.toplevel.fa --read-length 144 --skip-left 20 --skip-right 20 $fileIn > $sample.tabulate.out" >> tabulate_groupsOf2.sh; done;

split -d -l 1 tabulate_groupsOf2.sh tabulate_groupsOf2.sh.
chmod +x tabulate_groupsOf2.sh.*

for file in `ls tabulate_groupsOf2.sh.*`
do
    nohup ./$file > ${file}.sdout &
done

