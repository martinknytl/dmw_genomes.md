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
# download genome
```
wget https://ftp.xenbase.org/pub/Genomics/JGI/Xenla10.1/XENLA_10.1_genome.fa.gz
```
# unzip a genome assembly 
```
gunzip XENLA_10.1_genome.fa.gz
```
# make and edit text files
```
emacs -nw example.txt
```
my example:
```
emacs -nw XBO_SOX3L.fa
```

# loding modules on computecanada
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
# making XBO, XLA, and XTR blast databases
```
makeblastdb -in xxx.fa -dbtype nucl -out xxx.fa_blastable
```
my XLA example:
```
makeblastdb -in XENLA_10.1_genome.fa -dbtype nucl -out XENLA_10.1_genome.fa_blastable
```
# blast a query
```
blastn -query XBO_SOX3L.fa -db ../borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out XBO_SOX3L_to_XB
```
# convert mapping to a file format
```
more XBO_SOX3L_to_XB
```
# save all blast results to Google Drive
```
rsync -axvH --no-g --no-p knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/XB_cytogenetic_probes/* .
```
```*``` means the path ```/drives/g/My Drive/pracovni slozka/vyzkum/moje publikace/rozpracovane/Xenopus borealis/blast```
