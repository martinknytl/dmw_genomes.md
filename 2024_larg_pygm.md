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
# sbatch ./2020_trimmomatic.sh ../raw_data/plate1


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

**mapping script**
