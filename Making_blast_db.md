## connect to compute Canada
```
ssh username@graham.computecanada.ca
```
## change directory

first, I must go to the folder "projects", than to the "rrg-ben/knedlo"
```
cd /home/knedlo/projects
```
```
cd rrg-ben/knedlo
```
## make new directory
```
mkdir borealis_genome
```
## copy genome assembly from the folder "for_martin" to "borealis_genome
```
cp ../../for_martin/XB_genome_concat_scafs/* .
```
## download genome
```
wget https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
```
## unzip a genome assembly 
```
gunzip XENLA_10.1_genome.fa.gz
```
## loding modules on computecanada
```
module load 'blast+/2.13.0' StdEnv/2020 gcc/9.3.0
```

## check what dependences are needed for a module
```
module spider 'blast+/2.13.0'
```

## check jobs that are running
```
squeue -u name
```
## cancel a job
```
scancel jobID
```
## start a job
```
sbatch script.sh
```
## making XBO, XLA, and XTR blast databases
```
makeblastdb -in xxx.fa -dbtype nucl -out xxx.fa_blastable
```
my XLA example:
```
makeblastdb -in XENLA_10.1_genome.fa -dbtype nucl -out XENLA_10.1_genome.fa_blastable
```
## make and edit text files for each sequence separately
```
emacs -nw example.txt
```
my example:
```
emacs -nw XBO_SOX3L.fa
```

`control` + `x`

`control` + `c`

confirm by `y` like yes

## blast a query
X. borealis ndufs1L to XBO genome:
```
[knedlo@gra-login2 knedlo]$ blastn -query XB_cytogenetic_probes/XBO_ndufs1L.fa -db borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out XB_cytogenetic_probes/XBO_ndufs1L_to_XB
```
X. borealis ndufs1L to XLA genome:
```
[knedlo@gra-login2 knedlo]$ blastn -query XB_cytogenetic_probes/XBO_ndufs1L.fa -db XL_v10.1_genome/XENLA_10.1_genome.fa_blastable -outfmt 6 -out XB_cytogenetic_probes/XBO_ndufs1L_to_XL
```
## open mapped table
```
cd XB_cytogenetic_probes/
```
```
[knedlo@gra-login2 XB_cytogenetic_probes]$ more XBO_ndufs1L_to_XB
```
```
[knedlo@gra-login2 XB_cytogenetic_probes]$ more XBO_ndufs1L_to_XL 
```
## save all blast results to Google Drive (last step for all mapped sequences)
```
rsync -axvH --no-g --no-p knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/XB_cytogenetic_probes/* .
```
```*``` means the path ```/drives/g/My Drive/pracovni slozka/vyzkum/moje publikace/rozpracovane/Xenopus borealis/blast```
