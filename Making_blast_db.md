# connect to compute Canada
```
ssh username@graham.computecanada.ca
```
# change directory

first, I must to the folder "projects", than to the "rrg-ben/knedlo"
```
cd /home/knedlo/projects
```
```
cd rrg-ben/knedlo
```
# make new directory
```
mkdir borealis_genome
```
# copy genome assembly from the folder "for_martin" to "borealis_genome
```
cp ../../for_martin/XB_genome_concat_scafs/* .
```
# print working diretory
```
pwd
```

# list files

```
ls
```

# detailed list files
```
ls -lh
```
or 
```
ls -l
```
# make and edit text files
```
emacs -nw example.txt
```
my example:
```
emacs -nw XBO_SOX3L.fa
```

# loding moduals on computecanada
```
module load 'blast+/2.13.0' StdEnv/2020 gcc/9.3.0
```

# check what dependences are needed for a module
```
module spider 'blast+/2.13.0'
```

# check jobs that are running
```
squeue -u name
```
# cancel a job
```
scancel jobID
```
# start a job
```
sbatch script.sh
```
# making blast database
working with Ben to make blast database using the borealis genome
```
makeblastdb -in xxx.fa -dbtype nucl -out xxx.fa_blastable
```
my example:
```
makeblastdb -in XENLA_10.1_genome.fa -dbtype nucl -out XENLA_10.1_genome.fa_blastable
```
# download genome

```
wget https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
```
# unzip a genome assembly 
```
gunzip XENLA_10.1_genome.fa.gz
```
# blast a query
```
blastn -query XBO_SOX3L.fa -db ../borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out XBO_SOX3L_to_XB
```
# generate table format 6
```
more XBO_SOX3L_to_XL
```
