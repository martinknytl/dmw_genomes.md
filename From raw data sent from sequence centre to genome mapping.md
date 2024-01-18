This is a repo with Martin Knytl for the X. clivii, X. largeni, X. pygmaeus genomics project.

The goal of this project is to identify sex-linked regions in these species based on genomic data. 


### 1) dowload sequencing files to /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data**

each sample has R1 and R2 raw reads and is present three times = 28 samples *2 = 56 files *3 = 168 files in total

check it using `ls -1 | wc -l` in the folder `for_martin`

```
mkdir 2023_clivii_largeni_pygmaeus
mkdir raw_data
cp ~/projects/rrg-ben/for_martin/2023_clivii_largeni_pygmaeus/raw_data/* .
```
```.``` /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data

### 2)FastQC 

```
for i in *fastq.gz ; do fastqc $i; done
```
fastqc files with .html prefix will be created in same folder as raw data

go to folder where you want to create folder with fastqc.html files, make a directory 'mkdir fastqc_raw' ; dowload files to my googledisk account using command:

```
scp knedlo@graham.computecanada.ca:/home/knedlo... .. .. .. /*fastqc.html .
```

### 3) Run trimmomatic
```
cp -r ~/projects/rrg-ben/for_martin/2023_clivii_largeni_pygmaeus/ben_scripts/ .
```
```.``` /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/ben_scripts

**2020_trimmomatic.sh script**
```
#!/bin/sh
#SBATCH --job-name=trimmomatic
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=96:00:00
#SBATCH --mem=8gb
#SBATCH --output=trimmomatic.%J.out
#SBATCH --error=trimmomatic.%J.err
#SBATCH --account=def-ben

# run by passing a directory argument like this
# sbatch ./2020_trimmomatic.sh ../raw_data/plate1


module load StdEnv/2020
module load trimmomatic/0.39
# _R1.fq.gz
#v=1
#  Always use for-loop, prefix glob, check if exists file.
for file in $1/*_R1.fastq.gz; do         # Use ./* ... NEVER bare *
  if [ -e "$file" ] ; then   # Check whether file exists.
      java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.39.jar PE ${file::-12}_R1.fastq.gz ${file::-12}_R2.fastq.gz ${file::-12}_trim_R1.fq.gz ${file::-12}_trim_single_R1.fq.gz ${file::-12}_trim_R2.fq.gz ${file::-12}_trim_single_R2.fq.gz ILLUMINACLIP:TruSeq2_and_3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  fi
done
```
Trimmomatic can be executed by:
```
sbatch ../ben_scripts/2020_trimmomatic.sh .
```
```
cp ~/projects/rrg-ben/for_martin/2023_clivii_largeni_pygmaeus/ben_scripts/TruSeq2_and_3-PE-2.fa .
```

### 4) map trimmed sequences to X. laevis genome

each sample has (in the `raw_data` folder) R1 and R2 raw reads, R1 and R2 trimmed reads, and R1 and R2 single reads

first 10 R1 and R2 trimmed reads moved to the temp1 folder: 
```
mkdir temp1
mv NS.LH00147_0009.001.IDT_i7_18---IDT_i5_18.BJE1507-pool_trim_R*.fq.gz temp1
mv NS.LH00147_0009.001.IDT_i7_19---IDT_i5_19.Cas260421_female_trim_R*.fq.gz temp1/
mv NS.LH00147_0009.001.IDT_i7_20---IDT_i5_20.AMNH17295_female_trim_R*.fq.gz temp1/
mv NS.LH00147_0009.001.IDT_i7_21---IDT_i5_21.Z23339_male_trim_R*.fq.gz temp1/
mv NS.LH00147_0009.001.IDT_i7_30---IDT_i5_30.Z23341_female_trim_R*.fq.gz temp1/
mv NS.LH00147_0009.001.IDT_i7_31---IDT_i5_31.AMNH17292_female_trim_R*.fq.gz temp1/
mv NS.LH00147_0009.001.IDT_i7_32---IDT_i5_32.Cas262486_male_trim_R*.fq.gz temp1/
mv  NS.LH00147_0009.001.IDT_i7_33---IDT_i5_33.BJE1506-pool_trim_R*.fq.gz temp1/
mv  NS.LH00147_0009.001.IDT_i7_42---IDT_i5_42.Cas262409_male_trim_R*.fq.gz temp1/
mv NS.LH00147_0009.001.IDT_i7_43---IDT_i5_43.Z23350_male_trim_R*.fq.gz temp1/
```

**mapping script**
```
#!/bin/sh
#SBATCH --job-name=bwa_align
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=120:00:00
#SBATCH --mem=32gb
#SBATCH --output=bwa_align.%J.out
#SBATCH --error=bwa_align.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2020_align_paired_fq_to_ref.sh pathandname_of_ref path_to_paired_fq_filez
# sbatch 2020_align_paired_fq_to_ref.sh /home/ben/projects/rrg-ben/ben/2018_Austin_XB_genome/Austin_genome/Xbo.v1.fa.gz pathtofqfilez

module load bwa/0.7.17
module load samtools/1.10


for file in ${2}/*_trim_R2.fq.gz ; do         # Use ./* ... NEVER bare *    
    if [ -e "$file" ] ; then   # Check whether file exists.
        bwa mem ${1} ${file::-14}_trim_R1.fq.gz ${file::-14}_trim_R2.fq.gz -t 16 | samtools view -Shu - | samtools sort - -o ${file::-14}_sorted.bam
        samtools index ${file::-14}_sorted.bam
  fi
done
```

**mapping script executed by**
```
sbatch ../../ben_scripts/2020_align_paired_fq_to_ref.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa .
```
```.``` /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data/temp1

all codes are equivalent:

```
sbatch ../../ben_scripts/2020_align_paired_fq_to_ref.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa .

sbatch ../../ben_scripts/2020_align_paired_fq_to_ref.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data/temp2

sbatch /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/ben_scripts/2020_align_paired_fq_to_ref.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data/temp2
```

movement the other 10 trimmed files:

from directory
`/home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data`
```
ls *trim_R*
```

copy and paste first 20 files to Gibhub, remove enters, and add `mv` at the very beginning and `temp2` as a destination folder. All files will be moved to the `temp2` folder
```
mv NS.LH00147_0009.001.IDT_i7_44---IDT_i5_44.XEN170_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_44---IDT_i5_44.XEN170_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_45---IDT_i5_45.Z23340_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_45---IDT_i5_45.Z23340_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_54---IDT_i5_54.BJE1508-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_54---IDT_i5_54.BJE1508-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_55---IDT_i5_55.CAS260423-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_55---IDT_i5_55.CAS260423-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_56---IDT_i5_56.Cas262487_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_56---IDT_i5_56.Cas262487_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_66---IDT_i5_66.Z23342_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_66---IDT_i5_66.Z23342_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_67---IDT_i5_67.AMNH17293_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_67---IDT_i5_67.AMNH17293_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_68---IDT_i5_68.Z23337-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_68---IDT_i5_68.Z23337-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_6---IDT_i5_6.Cas260390_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_6---IDT_i5_6.Cas260390_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_78---IDT_i5_78.Cas260422_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_78---IDT_i5_78.Cas260422_female_trim_R2.fq.gz temp2/
```

similarly for the `temp3` folder etc. to the `tempX`
```
ls *trim_R*
```
```
mv NS.LH00147_0009.001.IDT_i7_79---IDT_i5_79.CAS260425-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_79---IDT_i5_79.CAS260425-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_7---IDT_i5_7.Z23349_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_7---IDT_i5_7.Z23349_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_80---IDT_i5_80.Cas262488_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_80---IDT_i5_80.Cas262488_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_8---IDT_i5_8.Cas260426_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_8---IDT_i5_8.Cas260426_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_90---IDT_i5_90.BJE1509-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_90---IDT_i5_90.BJE1509-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_91---IDT_i5_91.AMNH17294_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_91---IDT_i5_91.AMNH17294_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_92---IDT_i5_92.Z23338_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_92---IDT_i5_92.Z23338_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_9---IDT_i5_9.BJE1505-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_9---IDT_i5_9.BJE1505-pool_trim_R2.fq.gz NS.LH00147_0009.002.IDT_i7_18---IDT_i5_18.BJE1507-pool_trim_R1.fq.gz NS.LH00147_0009.002.IDT_i7_18---IDT_i5_18.BJE1507-pool_trim_R2.fq.gz NS.LH00147_0009.002.IDT_i7_19---IDT_i5_19.Cas260421_female_trim_R1.fq.gz NS.LH00147_0009.002.IDT_i7_19---IDT_i5_19.Cas260421_female_trim_R2.fq.gz temp3/
```

### 5) merge 3 bam files in 1

```
mkdir bams_combined
```

go to the `ben_scripts` folder and check if there are 3 files with the name AMNH17292:
`ls ../raw_data/*/*AMNH17292*` shows 12 files: 3 R1_trimmed, 3 R2_trimmed, 3 sorted.bam, and 3 sorted.bam.bai
`ls ../raw_data/*/*AMNH17292*bam*` shows 6 files: 3 sorted.bam, and 3 sorted.bam.bai
`ls ../raw_data/*/*AMNH17292*bam` shows 3 sorted.bam files
`/*/` involves all folders in the `raw_data` folder (all `temp` folders)

a script for merging:
```
#!/bin/sh
#SBATCH --job-name=merge
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=4
#SBATCH --time=32:00:00
#SBATCH --mem=16gb
#SBATCH --output=merge.%J.out
#SBATCH --error=merge.%J.err
#SBATCH --account=def-ben

# run by passing an argument like this (in the directory with the files)
# sbatch 2021_samtools_merge3.sh merged_path_and_outfile in1_path_and_file in2_path_and_file in3_path_and_file in4_path_and_file
module load samtools/1.10

samtools merge ${1} ${2} ${3} ${4}
samtools index ${1)
```

the script execution:
```
sbatch 2021_samtools_merge3.sh ../bams_combined/AMNH17292_female_sorted.bam ../raw_data/temp1/NS.LH00147_0009.001.IDT_i7_31---IDT_i5_31.AMNH17292_female_sorted.bam ../raw_data/temp4/NS.LH00147_0009.002.IDT_i7_31---IDT_i5_31.AMNH17292_female_sorted.bam ../raw_data/temp7/NS.LH00147_0009.003.IDT_i7_31---IDT_i5_31.AMNH17292_female_sorted.bam
```
by the same way to merge all sorted.bam files

### 6) read groups using picard

a script for read groups:

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
# sbatch 2022_picard_add_read_groups.sh bamfile_prefix

module load picard/2.23.3

    java -jar $EBROOTPICARD/picard.jar AddOrReplaceReadGroups I=${1}.bam O=${1}_rg.bam RGID=4 RGLB=${1} RGPL=ILLUMINA RGPU=${1} RGSM=${1}

module load StdEnv/2020 samtools/1.12
samtools index ${1}_rg.bam
```

go to the `ben_scripts` folder and every single script is separetely executed:
```
sbatch 2022_picard_add_read_groups_and_index.sh ../bams_combined/AMNH17292_female_sorted
sbatch 2022_picard_add_read_groups_and_index.sh ../bams_combined/AMNH17293_male_sorted
sbatch 2022_picard_add_read_groups_and_index.sh ../bams_combined/AMNH17294_male_sorted
```
etc till 
```
sbatch 2022_picard_add_read_groups_and_index.sh ../bams_combined/Z23350_male_sorted
```

### 7) call haplotype

a script for read groups:
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

module load nixpkgs/16.09 gatk/4.1.0.0

for file in ${2}*_rg.bam
do
    gatk --java-options -Xmx24G HaplotypeCaller  -I ${file} -R ${1} -L ${3} -O ${file}_${3}.g.vcf -ERC GVCF
done
```

from the `ben_scripts/` folder:
```
sbatch 2021_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ../bams_combined/temp1/ Chr1L

sbatch 2021_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ../bams_combined/temp1/ Chr2L

sbatch 2021_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ../bams_combined/temp1/ Chr3L

sbatch 2021_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ../bams_combined/temp1/ Scaffolds

sbatch 2021_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ../bams_combined/temp2/ Chr1L
```
for faster coding, use loop. This is an example of haplotype caller execution for Chr1L-8L in temp8 folder:
```
for x in {1..6}; do sbatch 2021_HaplotypeCaller.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ../bams_combined/to_complete_temp2II/ Chr$x\L; done
```

there are three files in each temp folder. Haplotype callet has to be done for each chromosome and scaffolds in each folder 

### 8) CombineGVCFs

```
#!/bin/sh
#SBATCH --job-name=CombineGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=12gb
#SBATCH --output=CombineGVCFs.%J.out
#SBATCH --error=CombineGVCFs.%J.err
#SBATCH --account=def-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:
# sbatch ../../ben_scripts/2021_CombineGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ Chr1L

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx10G CombineGVCFs -R ${1}"
for file in ${2}*${3}*g.vcf
do
    commandline+=" -V ${file}"
done

commandline+=" -L ${3} -O ${2}allsites_${3}.g.vcf.gz"

${commandline}
```

### OR, if this does not work (which it didn't for everyone except fischbergi, probably because of many scaffolds in the reference genome) use GenomicsDBImport to combine

```
#!/bin/sh
#SBATCH --job-name=GenomicsDBImport
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=60:00:00
#SBATCH --mem=24gb
#SBATCH --output=GenomicsDBImport.%J.out
#SBATCH --error=GenomicsDBImport.%J.err
#SBATCH --account=def-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:
# sbatch 2021_GenomicsDBImport.sh /home/ben/projects/rrg-ben/ben/2020_XL_v9.2_re
fgenome/XENLA_9.2_genome.fa /home/ben/projects/rrg-ben/ben/2020_GBS_muel_fish_al
lo_cliv_laev/raw_data/cutaddapted_by_species_across_three_plates/clivii/vcf/ chr
1L temp_path_and_dir db_path_and_dir_prefix

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx20G GenomicsDBImport -R ${1}"
for file in ${2}*g.vcf
do
    commandline+=" -V ${file}"
done

commandline+=" -L ${3} --tmp-dir=${4} --batch-size 50 --genomicsdb-workspace-pat
h ${5}_${3}"

${commandline}
```

I skipped this step for our largeni, clivii, and pygmaeus samples and performed the following one:

### 9) If gvcfs were combined with CombineGVCFs, use this:

```
#!/bin/sh
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=12gb
#SBATCH --output=GenotypeGVCFs.%J.out
#SBATCH --error=GenotypeGVCFs.%J.err
#SBATCH --account=def-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:
# sbatch ../../ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr1L.g.vcf.gz Chr1L

# ./ = /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/bams_combined/combined_GVCFs

# or loop for chromosomes 1-8L: for x in {1..8}; do sbatch ../../ben_scripts/2021_GenotypeGVCFs.sh /home/knedlo/projects/rrg-ben/knedlo/laevis_genome/2021_XL_v10_refgenome/XL_v10.1_concatenatedscaffolds.fa ./ allsites_Chr$x\L.g.vcf.gz Chr$x\L; done

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx10G GenotypeGVCFs -R ${1} -V ${2}${3} -L ${
4} -O ${2}${3}_${4}.vcf.gz"

${commandline}
```

### If gvcfs were combined using GenomicsDBImport use this:

```
#!/bin/sh
#SBATCH --job-name=GenotypeGVCFs
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=12:00:00
#SBATCH --mem=23gb
#SBATCH --output=GenotypeGVCFs.%J.out
#SBATCH --error=GenotypeGVCFs.%J.err
#SBATCH --account=def-ben


# This script will read in the *.g.vcf file names in a directory, and 
# make and execute the GATK command "GenotypeGVCFs" on these files. 

# execute like this:
# sbatch 2021_GenotypeGVCFs_DB.sh ref DB_path chr1L

module load nixpkgs/16.09 gatk/4.1.0.0

commandline="gatk --java-options -Xmx18G GenotypeGVCFs -R ${1} -V gendb://${2}_$
{3} -O ${2}_${3}_out.vcf"

${commandline}
```

### 10) VariantFiltration

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
# for example: sbatch ../../ben_scripts/2021_VariantFiltration.sh ./allsites_Chr2L.g.vcf.gz_Chr2L.vcf.gz
# ./ = /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/bams_combined/variant_Filtration


module load nixpkgs/16.09 gatk/4.1.0.0

gatk --java-options -Xmx8G VariantFiltration -V ${1}\
    -filter "QD < 2.0" --filter-name "QD2" \
    -filter "QUAL < 30.0" --filter-name "QUAL30" \
    -filter "SOR > 3.0" --filter-name "SOR3" \
    -filter "FS > 60.0" --filter-name "FS60" \
    -filter "MQ < 40.0" --filter-name "MQ40" \
    -filter "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
    -filter "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
    -O ${1}_filtered.vcf.gz
```

### 11) SelectVariants

```
#!/bin/sh
#SBATCH --job-name=SelectVariants
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=36:00:00
#SBATCH --mem=10gb
#SBATCH --output=SelectVariants.%J.out
#SBATCH --error=SelectVariants.%J.err
#SBATCH --account=def-ben


# This script will execute the GATK command "SelectVariants" on a file

# execute like this:
# sbatch 2021_SelectVariants.sh pathandfile
# example of loop: for x in {1..8}; do sbatch ../../ben_scripts/2021_SelectVariants.sh ./allsites_Chr$x\L.g.vcf.gz_Chr$x\L.vcf.gz_filtered.vcf.gz; done
# no loop: sbatch ../../ben_scripts/2021_SelectVariants.sh ./allsites_Chr9_10L.g.vcf.gz_Chr9_10L.vcf.gz_filtered.vcf.gz


module load nixpkgs/16.09 gatk/4.1.0.0

gatk --java-options -Xmx8G SelectVariants \
        --exclude-filtered \
        -V ${1} \
        -O ${1}_filtered_removed.vcf
```

### 12) to make tab delimited files out of the genotype files (gvcf files)

```
module load StdEnv/2020 vcftools/0.1.16
zcat file.vcf.gz | vcf-to-tab > out.tab
```

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
# sbatch 2024_vcf-to-tab.sh 

module load StdEnv/2020 vcftools/0.1.16

# loop:
for x in {1..8}; do sbatch ../../ben_scripts/2023_vcf-to-tab.sh zcat allsites_Chrx\L.g.vcf.gz_Chrx\L.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chrx\L.tab; done

for x in {1..8}; do sbatch ../../ben_scripts/2023_vcf-to-tab.sh zcat allsites_Chrx\S.g.vcf.gz_Chrx\S.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chrx\S.tab; done

# no loop: sbatch ../../ben_scripts/2021_SelectVariants.sh zcat allsites_Chr9_10L.g.vcf.gz_Chr9_10L.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chr9_10L.tab

zcat allsites_Chr9_10L.g.vcf.gz_Chr9_10L.vcf.gz_filtered.vcf.gz_filtered_removed.vcf.gz | vcf-to-tab > Chr9_10L.tab
```

### 13) Searching for sex-specific heterozygous positions

For sex-linked regions, we expect divergence in the heterogametic sex and not the homogametic sex. We can identify these sites and map them.  

Here is a perl script that will look for these sites in a tab delimited file:
```perl #!/usr/bin/env perl
#!/usr/bin/env perl
use strict;
use warnings;
use lib qw(~/perl_modules);
use List::MoreUtils qw/ uniq /;

# Prepare input file
# module load tabix
# bgzip -c file.vcf > file.vcf.gz
# tabix -p vcf file.vcf.gz
#  Now use vcftools to make a tab delimited file:
# module load StdEnv/2020 vcftools/0.1.16
# zcat file.vcf.gz | vcf-to-tab > out.tab

# on computecanada:
# module load perl/5.30.2
# module load gcc/9.3.0
# cpan
# install List::MoreUtils

#  This program reads in a tab delimited genotype file generated
#  by vcftools (vcf2tab) and searches for sites that are homozygous
#  in one sex for one SNP and at least partially heterozygous in the other sex

# to execute type Parse_tab.pl inputfile.tab 1111100110000111100011100110010100002200 interesting_sites.out proportion
# where 1111100110000111100011100110010100002200 refers to whether or not each individual in the ingroup 
# in the vcf file is (0) male, (1) female, and or (2) skipped

# proportion is the proportion of genotyped alleles in the heterogametic sex that are required to be
# different from the homogametic sex in order for the position to be reported.  This is a way to reduce reporting
# of low frequency polymorphisms (which are unlikely to be sex-linked but likely to have one sex all homozygous).
# the proportion parameter should be less than or equal to 0.5 

# if it is 0.5, this means all females are heterozygous and all males are homozygous (for positions with only 2 variants)

# we will also use this proportion to be a requirement for male-specific or female-specific SNPs, meaning at least
# this proportion of the individuals within each sex is required to have a genotype.

# het_sites.out sex_specific_sites.out diverged_sites.out are the output files that have the positions and chr of interesting sites

# example for clivii
# perl Parse_tab.pl clivii_unfiltered_removed_allchrs.vcf.tab 111111111111111111111111110000000000000000000 interesting_sites.out 0.35
# include only Eritrea:
# 222222221111111111111112222222222200000000222

# exclude Eritrea:
# 111111112222222222222221110000000022222222000

# Example for XB_WGS
# perl Parse_tab.pl XB_WGS_not_filtered_allchrs.vcf.gz.tab 100110011101010000102222 interesting_sites.out 0.5

my $inputfile = $ARGV[0];
my $input2 = $ARGV[1];
my $outputfile1 = $ARGV[2];
my $proportion = $ARGV[3];

print "hello ",$proportion,"\n";

unless (open DATAINPUT, $inputfile) {
	print "Can not find the input file.\n";
	exit;
}

unless (open(OUTFILE1, ">$outputfile1"))  {
	print "I can\'t write to $outputfile1\n";
	exit;
}
print "Creating output file: $outputfile1\n";
print OUTFILE1 "CHR\tPOS\tTYPE\tCATEGORY\tn_FEMs\tn_MALS\n";




my @sexes = split("",$ARGV[1]);

my @males=();
my @females=();
my @temp;
my @unique_male_nucleotides;
my @unique_female_nucleotides;
my $y;
my $x;
my $counter=0;
my $diverged=0;
my $number_of_male_individuals_genotyped=0;
my $number_of_female_individuals_genotyped=0;

for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 0){
		$number_of_male_individuals_genotyped +=1;
	}	
}	
for ($y = 0 ; $y <= $#sexes ; $y++ ) {
	if($sexes[$y] == 1){
		$number_of_female_individuals_genotyped +=1;
	}	
}
print "This includes ",$number_of_female_individuals_genotyped," female(s) and  ", $number_of_male_individuals_genotyped," males\n";

while ( my $line = <DATAINPUT>) {
	chomp($line);
	@temp=split /[\t\/]/,$line;
	if($temp[0] ne '#CHROM'){
		if($#temp ne (($#sexes+1)*2)+2){
			print "The number of individuals in the input line does not match the number of individuals genotyped ",
			$temp[0],"\t",$temp[1],"\t",$#temp," ",(($#sexes+1)*2)+2,"\n";
		}

		# parse the bases in all genotypes in each sex
		@males=();
		@females=();
		$counter=0;
		for ($y = 3 ; $y <= $#temp; $y=$y+2 ) {
			if(($temp[$y] ne ".")&&($temp[$y+1] ne ".")){
				if($sexes[$counter] == 0){
						push(@males, $temp[$y]);
						push(@males, $temp[$y+1]);
				}
				elsif($sexes[$counter] == 1){
					push(@females, $temp[$y]);
					push(@females, $temp[$y+1]);
				}	
			}
			$counter+=1;
		}	
		# OK I should have all the bases loaded for non-missing genotypes for each male and each female
		
		@unique_male_nucleotides = uniq @males;
		@unique_female_nucleotides = uniq @females;
		#print @females," ",@males,"\n";
		#print $#unique_male_nucleotides," ",$#unique_female_nucleotides,"\n";
		# looks fine
		if(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides != -1)){
			# we can compare homoz and het genotypes because both sexes have data
			if(($#unique_male_nucleotides == 0)&&($#unique_female_nucleotides > 0)){
				# all males are homoz but at least some females are hets or homoz for another SNP
				# check if the proportion of divergent positions in females is high enough
				$diverged=0;
				for ($x = 0 ; $x <= $#females ; $x++ ) {
					if($females[$x] ne $males[0]){
						$diverged+=1;
					}
				}
				if($diverged > $proportion*($#females+1)){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_SNP\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_SNP
					# Category of 1 means ZW or a female-specific SNP
				}
				# now check for the extreme case where all females are heterozygous and all males are homoz
				# this is rare because we expect some genotypes in females to be undercalled, even in sex-linked regions
				$diverged=0;
				for ($x = 0 ; $x <= $#females ; $x=$x+2 ) {
					if($females[$x] ne $females[$x+1]){
						$diverged+=1;
					}
				}
				if($diverged == ($#females+1)/2){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_heterozygosity\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_heterozygosity
					# 1 means all females are hets and all males are homoz
				}					
			}
			elsif(($#unique_male_nucleotides > 0)&&($#unique_female_nucleotides == 0)){
				# all females are homoz but at least some males are hets or homoz for another SNP
				# check if the proportion of divergent positions in males is high enough
				$diverged=0;
				for ($x = 0 ; $x <= $#males ; $x++ ) {
					if($males[$x] ne $females[0]){
						$diverged+=1;
					}
				}
				if($diverged > $proportion*($#males+1)){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_SNP\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_SNP
					# -1 means XY or a male-specific SNP
				}	
				# now check for the extreme case where all females are heterozygous and all males are homoz
				# this is rare because we expect some genotypes in females to be undercalled, even in sex-linked regions
				$diverged=0;
				for ($x = 0 ; $x <= $#males ; $x=$x+2 ) {
					if($males[$x] ne $males[$x+1]){
						$diverged+=1;
					}
				}
				if($diverged == ($#males+1)/2){	
					print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_heterozygosity\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
					# Sex_specific_heterozygosity
					# -1 means all males are hets and all females are homoz
				}
			}
			elsif(($#unique_male_nucleotides == 0)&&
			($#unique_female_nucleotides == 0)&&
			($unique_male_nucleotides[0] ne $unique_female_nucleotides[0])){
				# males are homoz, females are homoz, but fixed for a different diverged nucleotide
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tFixed_divergence\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
				# Fixed_divergence
				# 1 means diverged
			}
		}
		elsif(($#unique_male_nucleotides != -1)&&($#unique_female_nucleotides == -1)){
			# females have no data
			# could be male-specific
			if((($#males +1)/2) > $proportion*$number_of_male_individuals_genotyped){
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_nucleotides\t-1\t",($#females+1)/2,"\t",($#males+1)/2,"\n";
				# Sex_specific_nucleotides
				# -1 means male specific or male specific
			}	
		}
		elsif(($#unique_male_nucleotides == -1)&&($#unique_female_nucleotides != -1)){
			# males have no data
			# could be female-specific
			if((($#females +1)/2) > $proportion*$number_of_female_individuals_genotyped){
				print OUTFILE1 $temp[0],"\t",$temp[1],"\tSex_specific_nucleotides\t1\t",($#females+1)/2,"\t",($#males+1)/2,"\n"; 
				# Sex_specific_nucleotides
				# 1 means fem specific or female specific
			}	
		}
	}
} # end while	
close OUTFILE1;
```

# Plot Output with R

```R
setwd("/Users/Shared/Previously\ Relocated\ Items/Security/projects/2023_clivii_largeni_pygmaeus/parsetab")
library(ggplot2)
library(plyr)
library(viridis)
options(scipen=999)
# Parsetab collects stats on these site patterns:
# [1] "Fixed_divergence"  # each sex is homoz for a different nucleotide       
# [2] "Sex_specific_heterozygosity" # this is the most stringent and should be used - all individuals
        # from one sex are heteroz and all individuals of the other sex are homoz
        # 1 means female-specific heterozygosity; -1 means male-specific heterozygosity
# [3] "Sex_specific_SNP"  # a heterozygous genotype is found only in one sex          
# [4] "Sex_specific_nucleotides" # I think this menas that a nucleotide is found in only one sex, 
                                 # and that it can be homozygous or heterozygous


# cliv WGS - this inputfile has only "Sex_specific_heterozygosity"
dat <-read.table("allcliv_out_Sex_specific_heterozygosity.txt",header=F)
dat <-read.table("alllarg_out_Sex_specific_heterozygosity.txt",header=F)
dat <-read.table("allpygm_out_Sex_specific_heterozygosity.txt",header=F)


colnames(dat) <- c("CHR","POS","TYPE","CATEGORY","n_FEMs","n_MALS")

head(dat)



# Make a column that summarizes the number of fems and male
dat$FandM <- paste(dat$n_FEMs,dat$n_MALS,sep="_")

dat$POS <- as.numeric(dat$POS)
dat$CATEGORY <- as.numeric(dat$CATEGORY)
dat$n_FEMs <- as.numeric(dat$n_FEMs)
dat$n_MALS <- as.numeric(dat$n_MALS)

unique(dat$FandM)

summary <- count(dat, "FandM");summary

# remove rows where FandM is "1_1" - this means only 2 individuals are genotyped
dat <- dat[dat$FandM != "1_1", ]
summary <- count(dat, "FandM");summary

only_55_dat <- dat[dat$FandM == "5_5", ]

# make separate dataframes for female-specific and male-specific heterozygosity
fem_specific <- dat[dat$CATEGORY == 1, ]
mal_specific <- dat[dat$CATEGORY == -1, ]
summary <- count(fem_specific, "FandM");summary
summary <- count(mal_specific, "FandM");summary

fem_specific_noscaf <- fem_specific[fem_specific$CHR != "Scaffolds", ]

fplot<-ggplot(fem_specific, aes(x = POS/1000000, fill=FandM)) +
    geom_density(linewidth=0.15)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

fplot_noscaf<-ggplot(fem_specific_noscaf, aes(x = POS/1000000, fill=FandM)) +
    geom_density(linewidth=0.15)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

mplot<-ggplot(mal_specific, aes(x = POS/1000000, fill=FandM)) +
    geom_density(linewidth=0.15)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

# make the combined plot more stringent
fem_specific_noscaf <- fem_specific_noscaf[(fem_specific_noscaf$FandM != "1_2")&
                                               (fem_specific_noscaf$FandM != "2_1")&
                                               (fem_specific_noscaf$FandM != "1_3")&
                                               (fem_specific_noscaf$FandM != "3_1")&
                                               (fem_specific_noscaf$FandM != "2_2"), ]

fplot_allcombined<-ggplot(fem_specific_noscaf, aes(x = POS/1000000)) +
    geom_density(linewidth=0.15)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./pygm_female_specific_heterozygosity.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    fplot
dev.off()

pdf("./pygm_male_specific_heterozygosity.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    mplot
dev.off()

pdf("./pygm_female_noscaf_specific_heterozygosity.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    fplot_noscaf
dev.off()

pdf("./larg_female_specific_allcombined.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    fplot_allcombined
dev.off()



# make separate 55 dataframes for female-specific and male-specific heterozygosity
fem_55_specific <- dat[only_55_dat$CATEGORY == 1, ]
mal_55_specific <- dat[only_55_dat$CATEGORY == -1, ]
summary <- count(fem_55_specific, "FandM");summary
summary <- count(mal_55_specific, "FandM");summary


f55plot<-ggplot(fem_55_specific, aes(x = POS/1000000, fill=FandM)) +
    geom_density(linewidth=0.15)+
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlim(0,250)+
    xlab("Position") + ylab("Density") +
    scale_fill_viridis(discrete = TRUE) +
    theme_bw()

pdf("./larg_55_female_specific.pdf",w=8, h=30.0, version="1.4", bg="transparent")
    f55plot
dev.off()

# need to add endpoints to prevent densities from tanking
dat[nrow(dat) + 1,] <- c("Chr1L","232529967","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr1S","196169796","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr2L","184566229","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr2S","167897111","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr3L","145564449","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr3S","127416162","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr4L","156120765","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr4S","131359388","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr5L","174499024","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr5S","139053354","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr6L","157843502","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr6S","137668413","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr7L","136892544","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr7S","105895006","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr8L","123836259","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr8S","105436522","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr9_10L","135078614","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Chr9_10S","110702964","Sex_specific_heterozygosity","1","5","5")
dat[nrow(dat) + 1,] <- c("Scaffolds","42147407","Sex_specific_heterozygosity","1","5","5")
#View(all_SL)






View(all_SL)
# write.csv(all_SL,"all_SL_mappedtoXL.csv", row.names = FALSE)
other_types_of_SL <- dat[((dat$TYPE == "Sex_specific_SNP")|(dat$TYPE == "Fixed_divergence")|(dat$TYPE == "Sex_specific_nucleotides")) &
                  (as.numeric(dat$n_FEMs) >= 10) &  
                  (as.numeric(dat$n_MALS) >= 10),]
View(other_types_of_SL)
high_het_SL <- all_SL[(as.numeric(all_SL$POS) >=140372057)&(as.numeric(all_SL$POS) <=140433397),]
View(high_het_SL) # lots!
the_rest_of_dmrt1L <- all_SL[(as.numeric(all_SL$POS) >140449273)&(as.numeric(all_SL$POS) <=140462368),]
View(the_rest_of_dmrt1L) # nothing!

high_het_other_SL <- other_types_of_SL[(as.numeric(other_types_of_SL$POS) >=140372057)&(as.numeric(other_types_of_SL$POS) <=140433397),]
View(high_het_other_SL) # nothing!
the_other_rest_of_dmrt1L <- other_types_of_SL[(as.numeric(other_types_of_SL$POS) >140449273)&(as.numeric(other_types_of_SL$POS) <=140462368),]
View(the_other_rest_of_dmrt1L) # nothing!

#    geom_rect(aes(xmin=140372057,xmax=140433397,ymin=-2,ymax=Inf),color="blue",alpha=0) + # high heteroz
#    geom_rect(aes(xmin=140462145,xmax=140462368,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex2
#    geom_rect(aes(xmin=140449273,xmax=140449553,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex3
#    geom_rect(aes(xmin=140395623,xmax=140395753,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex4
#    ex5 missing
#    geom_rect(aes(xmin=140449273,xmax=140449319,ymin=-2,ymax=Inf),color="red",alpha=0) + # XB dmrt1L ex6    
    
all_SL_ZW <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
              (dat$n_FEMs == "5") &  
              (dat$n_MALS == "5") &
              (dat$CATEGORY == "1")   ,]

all_SL_XY <- dat[(dat$TYPE == "Sex_specific_heterozygosity") &
                     (dat$n_FEMs == "9") &  
                     (dat$n_MALS == "11") &
                     (dat$CATEGORY == "-1")   ,]

dim(all_SL_ZW)
dim(all_SL_XY)
View(all_SL)
View(all_SL_XY)





library(plyr)
summary <- count(all_SL_ordered, "CHR")
chr_sizes <- c(232529967,196169796,184566229,167897111,145564449,127416162,156120765,
               131359388,174499024,139053354,157843502,137668413,136892544,105895006,123836259,
               105436522,135078614,110702964,42147407)
summary <- cbind(summary,chr_sizes)
summary$proportions <- summary$freq/summary$chr_sizes
summary[summary$proportions == max(summary$proportions),]


all_SL_ordered <- all_SL[order(all_SL$CHR, all_SL$POS),]

#View(temp)
temp <-c("red","blue")

p<-ggplot(all_SL_ordered, aes(x = POS/1000000)) +
    geom_density()+
    #geom_vline(xintercept=49000000)+ 
    facet_wrap(~CHR, ncol = 1, scales = "free_x") +
    xlab("Position") + ylab("Density") +
    ylim(0,1) +
    theme_bw()


pdf("./Clivi_female_specific_heterozygosity.pdf",w=10, h=30.0, version="1.4", bg="transparent")
p
dev.off()


# filter data to retain only sites with a reasonable amount of
# genotypes in both sexes

subset_data <- dat[(dat$n_FEMs + dat$n_MALS) >= 10,]
dim(subset_data)


p<-ggplot(dat, aes(x=POS, y=CATEGORY*(n_FEMs+n_MALS))) + 
    # add points
    geom_point(size=2, alpha = 0.7 ) +
    # add loess line
    # geom_smooth() +
    # color the stuff the way I want
    facet_wrap(~CHR, ncol = 1) +
    # get rid of gray background
    theme_bw()
pdf("./out_0.5.pdf",w=18, h=18.0, version="1.4", bg="transparent")
    p
dev.off()


```
