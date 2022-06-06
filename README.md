# Training for construction of polygenic risk score using PRSice2
This tutorial is in development and is not yet finalized.
In different sources, the terms ‘polygenic score (PGS)’, ‘polygenic risk scores (PRS)’, and ‘genetic risk score (GRS)’ are used interchangeably. All refer to the same score where “[multi-locus profiles of genetic risk](https://pubmed.ncbi.nlm.nih.gov/23701538/), so-called genetic risk scores, can be used to translate discoveries from genome-wide association studies (GWAS) into tools for population health research” . Evident from the explanation, construction of a PRS is dependent on findings from GWAS.

This [weblink](https://www.genome.gov/Health/Genomics-and-Medicine/Polygenic-risk-scores) gives a very nice overview of PRS for readers who might need an update on their understanding of genetic variations and disease development and how complex diseases are different from single-gene (Mendelian) diseases. 

The outcome of interest in this training is Anorexia Nervosa (AN) and we use genetics data from 1000 Genomes phaase 3 release. The outcome (phenotype) data (AN) is simulated.

## Installation of packages
The following two packages are essential for this training: [PRSice2 package (v2.3.5)](http://www.prsice.info/), and [PLINK v.1.90b6.2](https://www.cog-genomics.org/plink/). Optionally, PRSice2 is also able to produce graphs using R. If interested, R version 4.0.0 is recommended.
```bash
# Make a directory and download the packages 
mkdir ~/prstrain
cd prstrain
wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip PRSice_linux.zip
wget https://s3.amazonaws.com/plink1-assets/plink_mac_20220402.zip
# test the binary files for executability
# use chmod for correctly setting the executability if required.
./plink --help
./PRSice_linux --help
```

## Download genetic dataset
We ues publicly available datasets from [1000 Genomes phase 3 release](https://www.internationalgenome.org/data-portal/data-collection/phase-3).<br/>
We use mostly unrelated individuals, and use SNPs in common with HapMap3 and UK Biobank. This dataset is provided by Florian Privé and is available [online](https://figshare.com/articles/dataset/1000_genomes_phase_3_files_with_SNPs_in_common_with_HapMap3/9208979).
```bash
# Make a directory for 1000 genome data and download the genetics files 
mkdir ~/prstrain/1000G
cd ~/prstrain/1000G
wget https://figshare.com/ndownloader/articles/9208979/versions/4
unzip 1000G_phase3_common_norel.zip
# You should have binary PLINK files
```
* Spend some minutes exploring the .fam and .bim file. Note .bed is not human-readable.

## Find and download GWAS summary results for Anorexia Nervosa (AN) and do initial QC
Get yourself familarize with [GWAS Catalog](https://www.ebi.ac.uk/gwas/), and search in PubMed for major GWAS in European populations for AN. Try to find the latest, and largest (highest sample size) GWAS for AN. Find if the GWAS summary results are publicly available and where you can download it.
<details>
<summary>Find the answer here</summary>
The largest GWAS for AN in European populatiosn as of June 2022 was published by [Watson et al](https://pubmed.ncbi.nlm.nih.gov/31308545/). Their summary results could be downloaded from the [PGS website](https://www.med.unc.edu/pgc/download-results/).
</details>

```bash
# Make a directory for GWAS summary results, also known as "base" data, and download the result
mkdir ~/prstrain/base
cd ~/prstrain/base
wget https://figshare.com/ndownloader/articles/14671980/versions/1
unzip 14671980.zip
# The .vcf file has a header. In Readme file of the donwloaded dataset, it has given an R code to remove the header.
# Remove the header of .vcf file and save it
##===R code to read in the TSV version of the VCF
##library(data.table)
##AN_basegwas.txt <-fread(file="pgcAN2.2019-07.vcf.tsv.gz", skip="CHROM\tPOS",stringsAsFactors=FALSE, data.table=FALSE)
# Compress back the file 
tar -czvf AN_basegwas.txt.gz  AN_basegwas.txt
# Do some QC of base data as descibed in PRSice2 basic tutorial (https://choishingwan.github.io/PRS-Tutorial/base/)
# Check if imputation quality is above 0.8. The base data does not have minor allele frequency (MAD) information to check, but the Readme file states MAF > 0.01 
gunzip -c AN_basegwas.txt.gz |\
awk 'NR==1 || ($10 > 0.8) {print}' |\
gzip > AN_basegwas.gz
# No SNP dropped due to low imputation quality (8219103 SNPs).
# Check for duplicate SNPs
gunzip -c AN_basegwas.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip - > AN_basegwas.nodup.gz
# After the duplicate SNPs are dropped, we now have 8125798 SNPs.
# Removing ambigeous SNPs
gunzip -c AN_basegwas.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > AN_basegwas.QC.gz
# After the ambigeous SNPs are dropped, we now have 7002697 SNPs.
    
```







