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


### 4) map trimmed sequences to X. laevis genome (bwa alignment)

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

after bwa mapping is done, trimed files moved to the "trimmed_trimmomatic" folder

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

after picard script is done, all files moved to the 'bwa_alignment' file. For example:

```
mv fem_larg_BJE1581/*sorted.bam fem_larg_BJE1581/*sorted.bam.bai bwa_alignment/
```
The 'sorted_rg.bam' and 'sorted_rg.bam.bai' files will be used for haplotype caller in next step.

### 6) call haplotype aka haplotype caller

```
#!/bin/sh
#SBATCH --job-name=HaplotypeCaller
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=30gb
#SBATCH --output=HaplotypeCaller.%J.out
#SBATCH --error=HaplotypeCaller.%J.err
#SBATCH --account=def-ben


# This script will read in the *_sorted.bam file names in a directory, and 
# make and execute the GATK command "RealignerTargetCreator" on these files. 

# execute like this:
# sbatch 2021_HaplotypeCaller.sh ref dir_of_bam chr
# sbatch 2021_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa bamfiledirectory chr

# module load nixpkgs/16.09 gatk/4.1.0.0
module load gatk/4.4.0.0 StdEnv/2023

for file in ${2}*_rg.bam
do
    gatk --java-options -Xmx24G HaplotypeCaller  -I ${file} -R ${1} -L ${3} -O ${file}_${3}.g.vcf -ERC GVCF
done
```

commands for the execution of scripts for all chromosomes in sample 'fem_pygm_ELI2081':

```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2024_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa /home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm/fem_pygm_ELI2081/ Chr$x\L; done
```
```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2024_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa /home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm/fem_pygm_ELI2081/ Chr$x\S; done
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2024_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa /home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm/fem_pygm_ELI2081/ Chr9_10L
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2024_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa /home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm/fem_pygm_ELI2081/ Chr9_10S
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2024_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa /home/knedlo/projects/rrg-ben/knedlo/2024_larg_pygm/fem_pygm_ELI2081/ Scaffolds
```

### 7) CombineGVCFs

```
#!/bin/sh
#SBATCH --job-name=CombineGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=12gb
#SBATCH --output=CombineGVCFs.%J.out
#SBATCH --error=CombineGVCFs.%J.err
#SBATCH --account=rrg-ben

# for graham.computecanada change account def-ben to rrg-ben 
# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute job for Chr9_10L:
# sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_CombineGVCFs.sh /home/knedlo/projects/def-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ Chr9_10L 

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx10G CombineGVCFs -R ${1}"
for file in ${2}*${3}*g.vcf
do
    commandline+=" -V ${file}"
done

commandline+=" -L ${3} -O ${2}allsites_${3}.g.vcf.gz"

${commandline}
```

using commands:

```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_CombineGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ Chr$x\L; done
```
```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_CombineGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ Chr$x\S; done
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_CombineGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ Chr9_10L
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_CombineGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ Chr9_10S
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_CombineGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ Scaffolds
```

### 8) GenotypeGVCFs

```
#!/bin/sh
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=12gb
#SBATCH --output=GenotypeGVCFs.%J.out
#SBATCH --error=GenotypeGVCFs.%J.err
#SBATCH --account=rrg-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:

# sbatch ../../ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr1L.g.vcf.gz Chr1L
# or loop for chromosomes 1-8L: for x in {1..8}; do sbatch ../../ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr$x\L.gvcf.gz Chr$x\L; done

# ./ = /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/bams_combined/combined_GVCFs

module load StdEnv/2023 gatk/4.4.0.0

commandline="gatk --java-options -Xmx10G GenotypeGVCFs -R ${1} -V ${2}${3} -L ${4} -O ${2}${3}_${4}.vcf.gz"

${commandline}
```

```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr9_10L.g.vcf.gz Chr9_10L
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr9_10S.g.vcf.gz Chr9_10S
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Scaffolds.g.vcf.gz Scaffolds
```
```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr$x\L.g.vcf.gz Chr$x\L; done
```
```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr$x\S.g.vcf.gz Chr$x\S; done
```

### 9) VariantFiltration

```
#!/bin/sh
#SBATCH --job-name=VariantFiltration
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=10gb
#SBATCH --output=VariantFiltration.%J.out
#SBATCH --error=VariantFiltration.%J.err
#SBATCH --account=def-ben


# This script will execute the GATK command "VariantFiltration" on a gvcf file

# execute like this:
# sbatch 2021_VariantFiltration.sh path_and_file

module load nixpkgs/16.09 gatk/4.1.0.0

gatk --java-options -Xmx8G VariantFiltration -V ${1}\
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 20.0" --filter-name "QUAL20" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -O ${1}_filtered.vcf.gz
```

```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_VariantFiltration.sh ./allsites_Chr$x\L.g.vcf.gz_Chr$x\L.vcf.gz; done
```
```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_VariantFiltration.sh ./allsites_Chr$x\S.g.vcf.gz_Chr$x\S.vcf.gz; done
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_VariantFiltration.sh ./allsites_Chr9_10L.g.vcf.gz_Chr9_10L.vcf.gz
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_VariantFiltration.sh ./allsites_Chr9_10S.g.vcf.gz_Chr9_10S.vcf.gz
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_VariantFiltration.sh ./allsites_Scaffolds.g.vcf.gz_Scaffolds.vcf.gz
```

### 10) SelectVariants

```
#!/bin/sh
#SBATCH --job-name=SelectVariants
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=10gb
#SBATCH --output=SelectVariants.%J.out
#SBATCH --error=SelectVariants.%J.err
#SBATCH --account=rrg-ben


# This script will execute the GATK command "SelectVariants" on a file

# execute like this:
# sbatch 2021_SelectVariants.sh pathandfile

# example of loop: for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_SelectVariants.sh/2021_SelectVariants.sh ./allsites_Chr$x\L.g.vcf.gz_Chr$x\L.vcf.gz_filtered.vcf.gz; done
# no loop: sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_SelectVariants.sh ./allsites_Chr9_10L.g.vcf.gz_Chr9_10L.vcf.gz_filtered.vcf.gz


module load nixpkgs/16.09 gatk/4.1.0.0

gatk --java-options -Xmx8G SelectVariants \
        --exclude-filtered \
        -V ${1} \
        -O ${1}_filtered_removed.vcf
```

from the directory "variantFiltration"
```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_SelectVariants.sh ./allsites_Chr$x\L.g.vcf.gz_Chr$x\L.vcf.gz_filtered.vcf.gz; done
```
```
for x in {1..8}; do sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_SelectVariants.sh ./allsites_Chr$x\S.g.vcf.gz_Chr$x\S.vcf.gz_filtered.vcf.gz; done
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_SelectVariants.sh ./allsites_Chr9_10L.g.vcf.gz_Chr9_10L.vcf.gz_filtered.vcf.gz
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_SelectVariants.sh ./allsites_Chr9_10S.g.vcf.gz_Chr9_10S.vcf.gz_filtered.vcf.gz
```
```
sbatch /home/knedlo/projects/rrg-ben/knedlo/ben_scripts/2021_SelectVariants.sh ./allsites_Scaffolds.g.vcf.gz_Scaffolds.vcf.gz_filtered.vcf.gz
```

### 11) convert vcf to tab

all `*filtered_removed.vcf` and *filtered_removed.vcf.idx` were moved to the `selectVariants` directory

Script for the L chromosomes:

```
#!/bin/sh
#SBATCH --job-name=vcf-to-tab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=10gb
#SBATCH --output=vcf-to-tab.%J.out
#SBATCH --error=vcf-to-tab.%J.err
#SBATCH --account=rrg-ben


# This script will execute zcat file.vcf.gz | vcf-to-tab > out.tab
# sbatch /home/knedlo/projects/rrg-ben/knedlo/martin_scripts/2024_vcf-to-tab_for_chrL.sh

module load StdEnv/2020 vcftools/0.1.16

# loop:

for x in {1..8}; do zcat allsites_Chr$x\L.g.vcf.gz_Chr$x\L.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chr$x\L.tab; done
allsites_Chr1L.g.vcf.gz_Chr1L.vcf.gz_filtered.vcf.gz_filtered_removed.vcf
# no loop: 

zcat allsites_Chr9_10L.g.vcf.gz_Chr9_10L.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chr9_10L.tab
```

Script for the S chromosomes:

```
#!/bin/sh
#SBATCH --job-name=vcf-to-tab
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=10gb
#SBATCH --output=vcf-to-tab.%J.out
#SBATCH --error=vcf-to-tab.%J.err
#SBATCH --account=def-ben


# This script will execute zcat file.vcf.gz | vcf-to-tab > out.tab
# sbatch 2024_vcf-to-tab_for_chrS.sh 

module load StdEnv/2020 vcftools/0.1.16

# loop:

for x in {1..8}; do zcat allsites_Chr$x\S.g.vcf.gz_Chr$x\S.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chr$x\S.tab; done

# no loop: 

zcat allsites_Chr9_10S.g.vcf.gz_Chr9_10S.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chr9_10S.tab
```

