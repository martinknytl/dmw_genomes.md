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
