I selected a position in largeni genome (all_larg_Sex_specific_heterozygosity.txt) that is heterozygous in all females and homozygous in all males

I Added candidate coordinates to the laevis genome on the web site https://www.xenbase.org/xenbase/displayJBrowse.do?data=data/xl10_1 and "go"

Now I see the selected region and can to dostinguis if it is exon, intron, UTR. I primarily focused on exons.

If the region is interesting, I extracted the refion within the Bash interpreter from downloaded genome using the `blastdbcmd` 

one of candidates regions is 136504257-136504516 of slc26a7.L gene located on Chr6L of X. largeni. This region corresponds to the 5'UTR and start of coding region of slc26a7.L in X. laevis. Using the code:

```
blastdbcmd -entry all -range 136504257-136504516 -db ../../laevis_genome/XENLA_10.1_genome.fa_blastable -out candidate_sex_locus_slc26a7.L+UTR_exon1.fa
```

the code gave me extractions from each chromosome. So I opend the output file in vi using `vi candidate_sex_locus_slc26a7.L+UTR_exon1.fa` and removed all sequences but Chr6L. 

I used the same laevis coordinates for the trop genome but the corresponding region has no gene there. slc26a7 gene is 10 Mbp upstream. Therefore I used `blastn` to search homologous trop sequence:

```
blastn -query candidate_sex_locus_slc26a7.L+UTR_exon1.fa -db /home/knedlo/projects/rrg-ben/knedlo/tropicalis_genome/XENTR_10.0_genome.fa_blastable -outfmt 6 -out candidate_sex_locus_slc26a7.L+UTR_exon1.fa_to_XT
```

I opened the output result using `more candidate_sex_locus_slc26a7.L+UTR_exon1.fa_to_XT` and see this: `:136504257-136504516	Chr6	84.277	159	18	4	1	159	126806460	126806611	3.00e-34	148`. trop coordinates are 126806460-126806611.

I went to xenbase and added 126806460-126806611 for search and the result whowed me the trop 5'UTR region

I extracted this trop sequence including 1000 pb on each side using the code:

```
blastdbcmd -entry all -range 126805460-126807611 -db /home/knedlo/projects/rrg-ben/knedlo/tropicalis_genome/XENTR_10.0_genome.fa_blastable -out candidate_sex_locus_slc26a7_extracted_from_XT
```

files with extravted sequences converted (remaned) to the '.fa' 

files downloaded to Google disk 

```
scp knedlo@graham.computecanada.ca:/home/knedlo/projects/rrg-ben/knedlo/2023_clivii_largeni_pygmaeus/largeni_candidate_sex_loci/slc26a7/ .
```

`.` is `/Users/martinknytl/Google Drive/My Drive/pracovni_slozka/vyzkum/xenopus/clivii_largeni_pygmaeus/heterozygous_positions/largeni_candidate_sex_loci/candicate_loci_from_Xenbase`

trop and laevis sequences aligned in Geneious, primer designed in the conserved regions for both species
