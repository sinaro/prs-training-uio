# Training for construction of polygenic risk score using PRSice2
This tutorial is in development and is not yet finalized. <br/>
In different sources, the terms ‘polygenic score (PGS)’, ‘polygenic risk scores (PRS)’, and ‘genetic risk score (GRS)’ are used interchangeably. All refer to the same score where “[multi-locus profiles of genetic risk](https://pubmed.ncbi.nlm.nih.gov/23701538/), so-called genetic risk scores, can be used to translate discoveries from genome-wide association studies (GWAS) into tools for population health research” . Evident from the explanation, construction of a PRS is dependent on findings from GWAS.

This [weblink](https://www.genome.gov/Health/Genomics-and-Medicine/Polygenic-risk-scores) gives a very nice overview of PRS for readers who might need an update on their understanding of genetic variations and disease development and how complex diseases are different from single-gene (Mendelian) diseases. 

The outcome of interest in this training is Anorexia Nervosa (AN) and we use genetics data from 1000 Genomes phaase 3 release. The outcome (phenotype) data (AN) is simulated.

## Downloading packages
The following two packages are essential for this training: [PRSice2 package (v2.3.5)](http://www.prsice.info/), and [PLINK v.1.90b6.2](https://www.cog-genomics.org/plink/). Optionally, PRSice2 is also able to produce graphs using R. If interested, R version 4.0.0 is recommended.
Downloading packages
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

## Find and download GWAS summary results for Anorexia Nervosa (AN) and perform QC of base data
Get yourself familarize with [GWAS Catalog](https://www.ebi.ac.uk/gwas/), and search in PubMed for major GWAS in European populations for AN. Try to find the latest, and largest (highest sample size) GWAS for AN. Find if the GWAS summary results are publicly available and where you can download it.
<details>
<summary>Find the answer here</summary>
The largest GWAS for AN in European populatiosn as of June 2022 was published by [Watson et al](https://pubmed.ncbi.nlm.nih.gov/31308545/). Their summary results could be downloaded from the [PGS website](https://www.med.unc.edu/pgc/download-results/).
</details>

* QC of base GWAS summary data
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

## Download genetics dataset
We ues publicly available datasets from [1000 Genomes phase 3 release](https://www.internationalgenome.org/data-portal/data-collection/phase-3).<br/>
We use mostly unrelated individuals, and use SNPs in common with HapMap3 and UK Biobank. This dataset is provided by Florian Privé and is available [online](https://figshare.com/articles/dataset/1000_genomes_phase_3_files_with_SNPs_in_common_with_HapMap3/9208979).
* Downloading genetics dataset
```bash
# Make a directory for 1000 genome data and download the genetics files 
mkdir ~/prstrain/1000G
cd ~/prstrain/1000G
!wget --content-disposition https://figshare.com/ndownloader/files/17838962
unzip 1000G_phase3_common_norel.zip
# You should have binary PLINK files
```
* Spend some minutes exploring the .fam and .bim file. Note .bed is not human-readable.

## QC of genetics data and SNP clumping using PLINK
The genetics dataset should not contain any duplicate SNP. Othwerwise, construction of PRS might run into trouble.
We are performing clumping of SNPs using PLINK. Although PRSice2 package is also able to perform clumping, we saw it was not optimized for our dataset.
* Removing duplicated SNP from genetics data
```bash
# See which SNP(s) is a duplicate in the .bim file
cd  ~/prstrain/1000G
cut -f 2 1000G_phase3_common_norel.bim | sort | uniq -d > 1.dups
# You see one SNP was duplicated. You can remove this SNP pair from the genetics data using PLINK.
./plink --bfile 1000G_phase3_common_norel --exclude 1.dups --make-bed --out 1000G_phase3_common_norel.nodup
```
* Clumping of SNPs
```bash
# You are recommended to write a bash script. For users running on HPC clusters, it is recommened to use a job scheduler (eg. Slurm).


```bash
# Example bash script to run on a personal device:
#!/bin/bash

./plink \
        --bfile 1000G_phase3_common_norel.nodup \
        --clump ~/prstrain/base/AN_basegwas.QC.gz \
        --clump-p1 1.0 \
        --clump-kb 250 \
        --clump-r2 0.1 \
        --clump-snp-field ID \
        --clump-field PVAL \
        --out 1000G_phase3_common_norel.nodup.clumped
        
# Example Slurm script to run on an HPC cluster:
#!/bin/bash
#SBATCH --jobname=plink_prstrain_20220925
#SBATCH --account=[to be assigned]
#SBATCH --partition=[accel]
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=12G
#SBATCH --cpus-per-task=6

# Set up job environment
module purge
module load plink/1.90b6.2

./plink \
        --bfile 1000G_phase3_common_norel.nodup \
        --clump ~/prstrain/base/AN_basegwas.QC.gz \
        --clump-p1 1.0 \
        --clump-kb 250 \
        --clump-r2 0.1 \
        --clump-snp-field ID \
        --clump-field PVAL \
        --out 1000G_phase3_common_norel.nodup.clumped

# It reports 177330 clumps formed frin 1259600 top variants, and results written to .clumped file.

```


