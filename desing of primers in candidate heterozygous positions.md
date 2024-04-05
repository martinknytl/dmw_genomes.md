I selected a position in largeni genome (all_larg_Sex_specific_heterozygosity.txt) that is heterozygous in all females and homozygous in all males

I Added candidate coordinates to the laevis genome on the web site https://www.xenbase.org/xenbase/displayJBrowse.do?data=data/xl10_1 and "go"

Now I see the selected region and can to dostinguis if it is exon, intron, UTR. I primarily focused on exons.

If the region is interesting, I extracted the refion within the Bash interpreter from downloaded genome using the `blastdbcmd` 

one of candidates regions is 136504257-136504516 of slc26a7.L gene located on Chr6L of X. largeni. This region corresponds to the 5'UTR and start of coding region of slc26a7.L in X. laevis. Using the code:

```
blastdbcmd -entry all -range 136504257-136504516 -db ../../laevis_genome/XENLA_10.1_genome.fa_blastable -out candidate_sex_locus_slc26a7.L+UTR_exon1.fa
```

the code gave me extractions from each chromosome. So I opend the output file in vi using `vi candidate_sex_locus_slc26a7.L+UTR_exon1.fa` and removed all sequences but Chr6L. 



one of candidates regions is 42296253-42296435 of slco3a1.S gene located on Chr3S of X. largeni. This region corresponds to the 5'UTR and start of coding region of slc26a7.L in X. laevis. Using the code:


blastdbcmd -entry all -range 42296253-42296435 -db ../../laevis_genome/XENLA_10.1_genome.fa_blastable -out candidate_sex_locus.fa
