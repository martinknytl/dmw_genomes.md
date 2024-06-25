This is a repo with Martin Knytl for the X. largeni and X. pygmaeus genomics project.

The goal of this project is to identify sex-linked regions in these species based on genomic data. 


### 1) dowload sequencing files to /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data**


```
mkdir 2024_larg_pygm

cp ~/projects/rrg-ben/for_martin/2024_larg_pygm* .
```
```.``` /home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm

### 2) Run trimmomatic

It was done by Ben


### 3) FastQC 

go to folder where you want to create folder with fastqc.html files in computecanada and in googledisk, make a directory 'mkdir fastqc/' 

computecanada: /home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm/fastqc
googledisk: /Users/martinknytl/Library/CloudStorage/GoogleDrive-knytlma@natur.cuni.cz/My Drive/pracovni_slozka/vyzkum/xenopus/clivii_largeni_pygmaeus/fastqc_trimmed_2024_larg_pygm

FasQC script of trimmed sequences:

```
#!/bin/sh
#SBATCH --job-name=fastQC
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=8gb
#SBATCH --output=fastQC.%J.out
#SBATCH --error=fastQC.%J.err
#SBATCH --account=def-ben

# run by passing a directory argument like this
# sbatch sbatch /home/knedlo/projects/rrg-ben/knedlo/martin_scripts/2024_fastQC.sh

module load fastqc/0.12.1
module load StdEnv/2023

for i in *fq.gz ; do fastqc $i
done
```

```
sbatch /home/knedlo/projects/rrg-ben/knedlo/martin_scripts/2024_fastQC.sh
```

dowload files to my googledisk account using command:

```
scp knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm/fastqc/*fastqc.html .
```


### 4) map trimmed sequences to X. laevis genome

use R1 and R2 trimmed reads for bwa mapping

make a folder for each individual sequences: e.g., `mkdir fem_larg_BJE1581/`, move the fem_larg_BJE1581_S8_trim_R1.fq.gz and fem_larg_BJE1581_S8_trim_R2.fq.gz files in the fem_larg_BJE1581 folder using `mv` command

`mkdir fem_larg_BJE1582/` and move the fem_larg_BJE1582_S9_trim_R1.fq.gz and fem_larg_BJE1582_S9_trim_R2.fq.gz files there etc.

use the '2024_align_paired_fq_to_ref.sh' script ('martin_csripts' folder because modules in the script were updated) for each individual separately 

```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=168:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2018_Austin_XB_genome/Austin_genome/Xbo.v1.fa.gz pathtofqfilez

module load bwa/0.7.17
# module load samtools/1.10
module load StdEnv/2023  gcc/12.3 samtools/1.20

for file in ${2}/*_trim_R2.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
	bwa mem ${1} ${file::-14}_trim_R1.fq.gz ${file::-14}_trim_R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${file::-14}_sorted.bam
	samtools index ${file::-14}_sorted.bam
  fi
done
```

script was executed:

```
sbatch ../../ben_scripts/2020_align_paired_fq_to_ref.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa .
```

### 5) read groups using picard

a script for read groups was modified with modules. It is not possible to change only one module because they require certain modules which can be found using 'spider' as it is recommended in the error file

```
#!/bin/sh
#SBATCH --job-name=readgroups
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --mem=8gb
#SBATCH --output=readgroups.%J.out
#SBATCH --error=readgroups.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this
# sbatch 2022_picard_add_read_groups_and_index.sh bamfile_prefix

# module load picard/2.23.3
module load picard/3.1.0

    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${1}.bam O=${1}_rg.bam RGID=4 RGLB=${1} RGPL=ILLUMINA RGPU=${1} RGSM=${1}

# module load StdEnv/2020 samtools/1.12
module load StdEnv/2023 samtools/1.20
samtools index ${1}_rg.bam
```

script was executed:
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2022_picard_add_read_groups_and_index.sh fem_pygm_ELI2372_S3_sorted
```
