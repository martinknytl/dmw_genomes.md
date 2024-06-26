## identification of coding sequences for variable positions that have the best hit 4/4 and higher

print rows with the hits 4/4, 4/5, and 5/4

```
awk '$(NF-1) >= 4 && $NF >= 4 {print}' all_larg_Sex_specific_heterozygosity.txt > all_larg_Sex_specific_heterozygosity_4_4_and_higher.txt
```

copy gff table (annotation file with exons, modified for synteny mapping) from computecanada to my Googledisk folder, cut 1-8 columns and convert to txt

```
cut -f1,2,3,4,5,6,7,8 XENLA_10.1_Xenbase_longest_CDSonly.gff > XENLA_10.1_Xenbase_longest_CDSonly_cut.txt
```
maybe tab was added, but not sure

then these two tables were transferred to RStudio and using script `sex_linked_largeni.R` I added the column 'exon' in which numbers 0 and 1 occur (0 = SNP is located in intronic locus, 1 = SNP is located in exonic locus). Table saved as all_larg_Sex_specific_heterozygosity_4_4_and_higher_exons.txt

THIS DOES NOT WORK: from the txt file (all_larg_Sex_specific_heterozygosity_5_5_only.txt) print rows, which has value in the column 2 in between values of columns 4 and 5 in the gff file (XENLA_10.1_Xenbase_longest_CDSonly.gff)

THIS WAS NOT USED:
```
awk -F $'\t' ' { if ($4 > $5) {t = $4; $4 = $5; $5 = t; print; } } ' OFS=$'\t' XENLA_10.1_Xenbase_longest_CDSonly.gff  > XENLA_10.1_Xenbase_longest_CDSonly_swap.txt
awk -F $'\t' ' { if ($4 < $5) {print; } } ' OFS=$'\t' XENLA_10.1_Xenbase_longest_CDSonly.gff  > XENLA_10.1_Xenbase_longest_CDSonly_nonswap.txt
awk '{print}' XTlongCDS_to_XL_Ssubgenome_nonswap.txt XTlongCDS_to_XL_Ssubgenome_swap.txt > XTlongCDS_to_XL_Ssubgenome_final.txt
<xxx.fa>: The desired output filename.

awk 'FNR==NR {min=$4; max=$5; next} $2 >= min && $2 <= max {print}' XENLA_10.1_Xenbase_longest_CDSonly.gff all_larg_Sex_specific_heterozygosity_5_5_only.txt > exons.txt

awk 'NR==FNR {min[$1]=$4; max[$1]=$5; next} $2 >= min[$2] && $2 <= max[$2]' XENLA_10.1_Xenbase_longest_CDSonly.gff all_larg_Sex_specific_heterozygosity_5_5_only.txt > exons.txt

awk -F $'\t' ' { if ($4 > $5) {t = $4; $4 = $5; $5 = t; print; } } ' OFS=$'\t' XENLA_10.1_Xenbase_longest_CDSonly_cut.txt  > 1_Xenbase_longest_CDSonly_cut_swap.txt
```

## Xenopus laevis

I selected a positions in largeni genome (all_larg_Sex_specific_heterozygosity.txt; for exonic regions ans SNPs with hits 4 and higher: all_larg_Sex_specific_heterozygosity_4_4_and_higher_exons.txt) that is heterozygous in all females and homozygous in all males in exonic region. If there are at least 5 variable positions within a locus

I Added candidate coordinates to the laevis genome on the web site https://www.xenbase.org/xenbase/displayJBrowse.do?data=data/xl10_1 and "go"

Now I see the selected region and can to distinguish if it is exon, intron, UTR. I primarily focused on exons. Then using a "copy" and "paste" fasta sequence of features from Xenbase to Geneious. Put annotations as much as it possible.

The interesting region using coordinates from table (all_larg_Sex_specific_heterozygosity.txt) was extracted within the Bash interpreter from downloaded genome using the `blastdbcmd` 

one of candidates regions is 136504257-136504516 of slc26a7.L gene located on Chr6L of X. largeni. This region corresponds to the 5'UTR and start of coding region of slc26a7.L in X. laevis. Using the code:

```
blastdbcmd -entry all -range 136504257-136504516 -db ../../laevis_genome/XENLA_10.1_genome.fa_blastable -out candidate_sex_locus_slc26a7.L+UTR_exon1.fa
```

the code gave me extractions from each chromosome. So I opend the output file in vi using `vi candidate_sex_locus_slc26a7.L+UTR_exon1.fa` and removed all sequences but Chr6L. Sequence can be directly copied from vi to Geneious, annotate exons and UTRs, and design primers in the conserved regions for both species.

## Xenopus tropicalis

I used the sequence extracted from the laevis genome including nucleotides from both sides to have a sequence of app 1 kb. Then I copied this sequece to the Blast on Xenbase. Click on the best hit ('Overview of Results'). This redirected me to the trop database where I can see the structure of a gene and trop coordinates. Then I used the coordinates for extraction of the sequence using this code:

```
blastdbcmd -entry all -range 2881130-2881853 -db /home/knedlo/projects/rrg-ben/knedlo/tropicalis_genome/XENTR_10.0_genome.fa_blastable -out candidate_sex_locus_LOC121395205_trop_exon.fa
```

the code gave me extractions from each chromosome. So I opend the output file in vi using `vi candidate_sex_locus_LOC121395205_trop_exon.fa` and removed all sequences but Chr6.

Sequence can be directly copied to Geneious, annotate exons and UTRs, and design primers in the conserved regions for both species.


Another way is:
I used the same laevis coordinates for the trop genome but the corresponding region has no gene there. slc26a7 gene is 10 Mbp upstream. Therefore I used `blastn` to search homologous trop sequence:

```
blastn -query candidate_sex_locus_slc26a7.L+UTR_exon1.fa -db /home/knedlo/projects/rrg-ben/knedlo/tropicalis_genome/XENTR_10.0_genome.fa_blastable -outfmt 6 -out candidate_sex_locus_slc26a7.L+UTR_exon1.fa_to_XT
```

I opened the output result using `more candidate_sex_locus_slc26a7.L+UTR_exon1.fa_to_XT` and see this: `:136504257-136504516	Chr6	84.277	159	18	4	1	159	126806460	126806611	3.00e-34	148`. trop coordinates are 126806460-126806611.

I went to xenbase and added 126806460-126806611 for search and the result showed me the trop 5'UTR region

I extracted this trop sequence including 1000 pb on each side using the code:

```
blastdbcmd -entry all -range 126805460-126807611 -db /home/knedlo/projects/rrg-ben/knedlo/tropicalis_genome/XENTR_10.0_genome.fa_blastable -out candidate_sex_locus_slc26a7_extracted_from_XT
```

files with extravted sequences converted (remaned) to the '.fa' 

files downloaded to Google disk (or copied directly to Geneios)

```
scp knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/largeni_candidate_sex_loci/slc26a7/ .
```

`.` is `/Users/martinknytl/Google Drive/My Drive/pracovni_slozka/vyzkum/xenopus/clivii_largeni_pygmaeus/heterozygous_positions/largeni_candidate_sex_loci/candicate_loci_from_Xenbase`

## Xenopus borealis

borealis genome has not graphical search such as trop and laevis and thus an approach of sequence extraction is a bit different than in trop and laevis

blastn -query xlae_slco3a1.S.fa -db ../../borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out xlae_slco3a1.S_to_xb

Search of the slco3a1.S gene in borealis genome:

```
blastn -query xlae_slco3a1.S.fa -db ../../borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable -outfmt 6 -out xlae_slco3a1.S_to_xb
```

The best bit score is a region 8S:15341892-15343716 (1829 bp). I need to extract this region from blastable database in fasta. This region is candidate region for sex determination in X. largeni and it is gene slco3a1.S

```
blastdbcmd -entry all -range 15341892-15343716 -db ../../borealis_genome/Xbo.v1_chrs_and_concatscafs_blastable -out xxx.fa
```

<all>: The identifier of the entire sequence you want to extract.
<15341892> and <15343716>: The start and end positions of the region you wish to extract.
<Xbo.v1_chrs_and_concatscafs_blastable>: The name of your database.

The command extracted 15341892-15343716 positions from each chromosomes. So then open `vi xxx.fa` and remove all sequences that we do not need

file renamed as using vi `xbo_slco3a1.S.fa`
