### 1) dowload sequencing files to /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data**

each sample has R1 and R2 raw reads = 28 samples = 56 files
```
mkdir 2023_clivii_largeni_pygmaeus
mkdir raw_data
cp ~/projects/rrg-ben/for_martin/2023_clivii_largeni_pygmaeus/raw_data/* .
```
```.``` /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data

### 2) Run trimmomatic
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

### 3) map trimmed sequences to X. laevis genome

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

copy and paste first 20 files to Gibhub, remove enters, and add `mv` at the very beginning and `temp2` as a destination folder
```
mv NS.LH00147_0009.001.IDT_i7_44---IDT_i5_44.XEN170_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_44---IDT_i5_44.XEN170_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_45---IDT_i5_45.Z23340_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_45---IDT_i5_45.Z23340_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_54---IDT_i5_54.BJE1508-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_54---IDT_i5_54.BJE1508-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_55---IDT_i5_55.CAS260423-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_55---IDT_i5_55.CAS260423-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_56---IDT_i5_56.Cas262487_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_56---IDT_i5_56.Cas262487_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_66---IDT_i5_66.Z23342_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_66---IDT_i5_66.Z23342_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_67---IDT_i5_67.AMNH17293_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_67---IDT_i5_67.AMNH17293_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_68---IDT_i5_68.Z23337-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_68---IDT_i5_68.Z23337-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_6---IDT_i5_6.Cas260390_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_6---IDT_i5_6.Cas260390_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_78---IDT_i5_78.Cas260422_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_78---IDT_i5_78.Cas260422_female_trim_R2.fq.gz temp2/
```
these files move to the `temp2` folder

```
ls *trim_R*
```
mv NS.LH00147_0009.001.IDT_i7_79---IDT_i5_79.CAS260425-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_79---IDT_i5_79.CAS260425-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_7---IDT_i5_7.Z23349_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_7---IDT_i5_7.Z23349_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_80---IDT_i5_80.Cas262488_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_80---IDT_i5_80.Cas262488_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_8---IDT_i5_8.Cas260426_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_8---IDT_i5_8.Cas260426_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_90---IDT_i5_90.BJE1509-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_90---IDT_i5_90.BJE1509-pool_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_91---IDT_i5_91.AMNH17294_male_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_91---IDT_i5_91.AMNH17294_male_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_92---IDT_i5_92.Z23338_female_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_92---IDT_i5_92.Z23338_female_trim_R2.fq.gz NS.LH00147_0009.001.IDT_i7_9---IDT_i5_9.BJE1505-pool_trim_R1.fq.gz NS.LH00147_0009.001.IDT_i7_9---IDT_i5_9.BJE1505-pool_trim_R2.fq.gz NS.LH00147_0009.002.IDT_i7_18---IDT_i5_18.BJE1507-pool_trim_R1.fq.gz NS.LH00147_0009.002.IDT_i7_18---IDT_i5_18.BJE1507-pool_trim_R2.fq.gz NS.LH00147_0009.002.IDT_i7_19---IDT_i5_19.Cas260421_female_trim_R1.fq.gz NS.LH00147_0009.002.IDT_i7_19---IDT_i5_19.Cas260421_female_trim_R2.fq.gz temp3/
```



