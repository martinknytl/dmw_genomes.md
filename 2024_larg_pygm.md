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

or even sequences of one individual can be moved to one folder and thus the script can run for each sequence separately (this was preformed in `2024_larg_pygm` 

**mapping script**
