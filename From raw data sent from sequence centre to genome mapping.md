### 1) dowload sequencing files to /home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/raw_data**

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
cp ~/projects/rrg-ben/for_martin/2023_clivii_largeni_pygmaeus/ben_scripts/TruSeq2_and_3-PE-2.fa .
