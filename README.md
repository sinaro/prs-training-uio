# Training for construction of polygenic risk score using PRSice2 (September 26th, 2022)
[Last update: 23.09.2022. This tutorial is in development and is not yet finalized]. <br/>
In different sources, the terms ‘polygenic score (PGS)’, ‘polygenic risk scores (PRS)’, and ‘genetic risk score (GRS)’ are used interchangeably. All refer to the same score where “[multi-locus profiles of genetic risk](https://pubmed.ncbi.nlm.nih.gov/23701538/), so-called genetic risk scores, can be used to translate discoveries from genome-wide association studies (GWAS) into tools for population health research”. It is evident from the explanation, construction of a PRS is dependent on findings from GWAS.

This [weblink](https://www.genome.gov/Health/Genomics-and-Medicine/Polygenic-risk-scores) gives a very nice overview of PRS for readers who might need an update on their understanding of genetic variations and disease development and how complex diseases are different from single-gene (Mendelian) diseases. 

To get familiar with PRS, read this [tutorial paper](https://pubmed.ncbi.nlm.nih.gov/32709988/).

The outcome of interest in this training is Anorexia Nervosa (AN). "[AN](https://www.mayoclinic.org/diseases-conditions/anorexia-nervosa/symptoms-causes/syc-20353591) is an eating disorder characterized by an abnormally low body weight, an intense fear of gaining weight and a distorted perception of weight. People with anorexia place a high value on controlling their weight and shape, using extreme efforts that tend to significantly interfere with their lives".
We use genetics data from 1000 Genomes phase 3 release. The outcome (phenotype) data (AN) is simulated.

## Source materials for practical
* Source materials for this tutorial can be found in the [Github repository](https://github.com/sinaro/prs-training-uio).
* all code lines are compatible with bash, unless otherwise stated for R.

## Downloading programs
* The following two packages are essential for this training: [PRSice2 package (v2.3.5)](http://www.prsice.info/), and [PLINK v.1.90b6.2](https://www.cog-genomics.org/plink/). Optionally, PRSice2 is also able to produce graphs using R. If interested, R version 4.0.0 is recommended. R, or another statistical package might be needed to prepare datasets for training. <br/>
* This training was tested on Ubuntu 20.04.4 LTS, Red Hat Enterprise Linux Server release 7.9 (Maipo), and Mac OS (Catalina, version 10.15.7).
* The tutorial assumes that you have root/admin privileges on your device. On an HPC system, this also works but it might differ in some steps.

* Downloading and preparing packages:
```bash
# Make a directory and download the packages 
mkdir ~/prstrain
cd prstrain
## Download PRSice2 and unzip it
wget https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip
unzip PRSice_linux.zip
# Note "wget" might not be available by default on Mac. You might want to test "curl" instead. Check their availability by "which wget" and "which curl".
# If "curl" is available. You can try downloading the link by:
# curl -o ./prs "https://github.com/choishingwan/PRSice/releases/download/2.3.5/PRSice_linux.zip"
# If this solution did not work as well, you need to directly download it and then transfer it to your work directory.
# You might need to go to security on Mac and choose "allow it anyway".
## Download PLINK and unzip it
wget https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip
unzip plink_linux_x86_64_20220402.zip
# test the binary files for executability
./plink --help
./PRSice_linux --help
# Optional: you can put the binaries in the PATH so it will be accessible system-wide (rather than in a specific directory).
```

## Find and download GWAS summary results for Anorexia Nervosa (AN) and perform QC of base data
* Get yourself familiarize with [GWAS Catalog](https://www.ebi.ac.uk/gwas/), and search in [PubMed](https://pubmed.ncbi.nlm.nih.gov/) for major GWAS in European populations for AN. <br/> 
* Is Pubmed able to interpret the phenotype "AN", and "GWAS" correctly? <br/>
* Try to find the latest, and largest (highest sample size) GWAS for AN. Find if the GWAS summary results are publicly available and where you can download it.
<details>
<summary>Click and find the answer here</summary>
The largest GWAS for AN in European population as of June 2022 was published by [Watson et al](https://pubmed.ncbi.nlm.nih.gov/31308545/). Their summary results could be downloaded from the [PGS website](https://www.med.unc.edu/pgc/download-results/).
</details>


* QC of base GWAS summary data
```bash
# Make a directory for GWAS summary results, also known as "base" data, and download the summary results
mkdir ~/prstrain/base
cd ~/prstrain/base
wget --content-disposition https://figshare.com/ndownloader/files/28169271
# Compare MD5 checksum between the downloaded file and the one mentioned. They should report the same string.
md5sum pgcAN2.2019-07.vcf.tsv.gz
# The .vcf file has a header.
# Remove the header of .vcf file and save it
##===R code to read in the TSV version of the VCF
install.packages("data.table")
library(data.table)
AN_basegwas <- fread(file="pgcAN2.2019-07.vcf.tsv.gz", skip="CHROM\tPOS",stringsAsFactors=FALSE, data.table=FALSE)
# 8219102 SNPs in the file.
### Do some QC of base data as described in PRSice2 basic tutorial (https://choishingwan.github.io/PRS-Tutorial/base/)
## general QC information mentioned in Readme that is of interest.
# genotyping rate > 0.99 (here a call rate >= 98%)
# sample missingness < 0.02 (ok)
# Hardy-Weinberg equilibrium (HWE) (p > 1 * 10-6)
# the human genome build is "GRCh37". This is important as both base and target data should be from the same build (or should become comparable with LiftOver).
## common QC steps for base file.
# Check Imputation quality and MAF. 
# Check if imputation quality is above 0.8 in R. The base data does not have minor allele frequency (MAF) information to check, but the Readme file states MAF > 0.01 
library(data.table)
install.packages("dplyr")
library(dplyr)
AN_basegwas_imp <- subset(AN_basegwas, IMPINFO > 0.8)
fwrite(AN_basegwas_imp, "AN_basegwas_imp.gz", sep="\t")
# 7818559 SNPs left after imputation quality check.
# Check for duplicate SNPs
gunzip -c AN_basegwas_imp.gz |\
awk '{seen[$3]++; if(seen[$3]==1){ print}}' |\
gzip > AN_basegwas.nodup.gz
# After that the duplicate SNPs are dropped, we now have 7732135 SNPs.
## Removing ambigeous SNPs
gunzip -c AN_basegwas.nodup.gz |\
awk '!( ($4=="A" && $5=="T") || \
        ($4=="T" && $5=="A") || \
        ($4=="G" && $5=="C") || \
        ($4=="C" && $5=="G")) {print}' |\
    gzip > AN_basegwas.QC.gz
# After the ambigueous SNPs are dropped, we now have 6663432 SNPs.
## Mismatching SNPs: this is taken care of by the program. The program also reports this in a file.
## Check the presence of effect and non-effect allele. Here we have REF and ALT allele. If this was not directly mentioned, this should be asked from the GWA authors. Otherwise, the association will be in the opposite direction.
## SNPs on sex chromosome. There are models that can use such information but we are working with autosomal SNPs.
less AN_basegwas.QC.gz | awk '{print $1}' | sort -n | uniq -c
## Here you see you only have SNPs on autosomal chromosomes (ch 1-22).
```
* LDlink interactive tool. Visit https://ldlink.nci.nih.gov/, choose LDpair tool. As an example, choose the first SNP in the base file "rs62513865". Evaluate LD pair with two other variants in the base file "rs16898460" and "rs2507789".


## Download 1000 Genomes genetics dataset and prepare it for analysis
* We ues publicly available datasets from [1000 Genomes phase 3 release](https://www.internationalgenome.org/data-portal/data-collection/phase-3).<br/>
* We use mostly unrelated individuals, and use SNPs in common with HapMap3 and UK Biobank. This dataset is provided by Florian Privé and is available [online](https://figshare.com/articles/dataset/1000_genomes_phase_3_files_with_SNPs_in_common_with_HapMap3/9208979).
* QC is done both at SNP level and at sample/population level.
* Downloading genetics dataset
```bash
# Make a directory for 1000 genome data, download the binary plink files, and unzip it.
mkdir ~/prstrain/1000G
cd ~/prstrain/1000G
wget --content-disposition https://figshare.com/ndownloader/files/17838962
unzip 1000G_phase3_common_norel.zip 
# on Mac OS, you can try: <curl -O http://figshare.com/ndownloader/files/17838962 --location-trusted>
# You should now have binary PLINK files
```
* Spend some minutes exploring the .fam and .bim file. Note .bed is not human-readable. Note the number of individuals in the .fam file is 2490.
Have a look at .fam2 file in the repository. Take a look at the populations. Not all of them are European. We generally would like to construct PRS in specific populations (in our tutorial, only European ancestry). 

* Visualize the population ancestry structure by their first two PCs

```bash
# In this tutorial, we will use the QC'ed SNPs and individuals to get eigenvec of population.
# To get familiar with PCA (a dimensionality reduction method), please follow a video by (StatQuest)[https://www.youtube.com/watch?v=FgakZw6K1QQ)].
# PC1 and PC2 explain the the first and second highest variance, respectively.
# We LD-prune the SNPs first with plink
./plink --bfile 1000G_phase3_common_norel --indep-pairwise 50 10 0.1 --out pcaprun  #1455110 of 1664850 variants removed.
# Getting 10 first PCs using the pruned SNP file.
./plink --bfile 1000G_phase3_common_norel --extract pcaprun.prune.in --pca 10 --out pcaprun
# Loading the eigenvec file in R and using ggplot to visualize the population substructure.
# Note we use .fam2 file to color code points by population ancestry.
install.packages("dplyr")
install.packages("ggplot2")
library(dplyr)
library(ggplot2)
pcaprun <- read.table("~/prstrain/1000G/pcaprun.eigenvec", quote="\"", comment.char="")
fam2 <- read.delim("~/prstrain/1000G/1000G_phase3_common_norel.fam2")
colnames(pcaprun)[2] <- "sample.ID"
pcaprunfam <- left_join(pcaprun, fam2, by="sample.ID")
gpcaprunfam <- ggplot(pcaprunfam, aes(x = V3, y = V4, color = Super.Population)) +
        geom_point() +
        xlab("PCA1") +
        ylab("PCA2") +
        coord_fixed()
plot(gpcaprunfam)
# you can clearly see that populations have a distinct cluster for PC1 and PC2, and that members of the European ancestry tend to cluster nicely (no outlier).
```

* Get a list of IDs for participants of European ancestry in the dataset.

```bash
# Here we provide it using R.
family.fam <- read.delim("~/prstrain/1000G/1000G_phase3_common_norel.fam2")
# We use R dplyr package to manipulate the dataset easier
library(dplyr)
# Restrict it to European samples only
eurfamily.fam <- filter(family.fam, Super.Population == "EUR") #503 samples of the total 2490 are European
# Drop columns we no longer need, and some steps to make it look like a .fam file format
eurfamily.fam <- eurfamily.fam %>% 
                  select(-c("Population", "Population.Description", "Super.Population")) %>% 
                  rename(IID=sample.ID) %>% 
                  mutate(FID = IID) %>% 
                  select(FID, IID, sex)
# save the covariate file (this can be useful when we want to compare the the null model explained variance (sex, age, PCs, etc) to that of PRS. Keep header of the file).
utils::write.table(eurfamily.fam, "eurfamily.cov", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# save the list of IDs of participants of European ancestry.
eurfamily.fam <- select(eurfamily.fam, IID)
utils::write.table(eurfamily.fam, "eurfamily.txt", sep="\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
# alternatively you can download it from the Github of this tutorial.
```


## QC of genetics data and SNP clumping using PLINK
* The target sample is recommended to be at least 100, this is the case here even after selection for participants of European ancestry. <br/>
* The human genome build should be similar between the base and target data. The [target data has used GRCh37 reference genome] (https://www.internationalgenome.org/data-portal/data-collection/phase-3). This is also the case with that of the base data (as we previously discussed in Readme of the base file) <br/>
* The genetics dataset should not contain any duplicate SNP. Otherwise, construction of PRS might run into trouble.
We are performing clumping of SNPs using PLINK. Although PRSice2 package is also able to perform clumping, we saw it was not optimized for our dataset.
* Removing duplicated SNP from genetics data
```bash
# See which SNP(s) is a duplicate in the .bim file
cd  ~/prstrain/1000G
cut -f 2 1000G_phase3_common_norel.bim | sort | uniq -d > 1.dups
# You see one SNP was duplicated. You can remove this SNP pair from the genetics data using PLINK. Use "out" to create new PLINK binary files. 
# Easier to copy plink binary to this directory to work with.
# Note that the plink binary files need to be treated all together with plink. You cannot, and should not, for instance, treat .bim file as a text file to manipulate!
./plink --bfile 1000G_phase3_common_norel --exclude 1.dups --make-bed --out 1000G_phase3_common_norel.nodup
```
* Clumping of SNPs <br/>
You are recommended to write a bash script. For users running on HPC clusters, it is recommended to use a job scheduler (eg. Slurm).

```bash
# Example bash script to run on a personal device (create a file with Nano or Vim text editor called "clump.sh" and execute it with ./clump.sh):
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
# It reports 168483 clumps formed from 1232246 top variants, and results written to .clumped file
```

```bash
# Example Slurm script to run on an HPC cluster:
#!/bin/bash
#SBATCH --jobname=plink_prstrain_20220926
#SBATCH --account=[to be assigned]
#SBATCH --partition=[accel]
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1

# Set up job environment
module purge
module load plink/1.90b6.2

./plink \
        # The rest of arguments (bfile etc) is the same as above.

# It reports 168483 clumps formed from 1232246 top variants, and results written to .clumped file.
```

## QC relevant for base and target data
* Human genome build comparability was already checked. <br/>
* Sample overlap between base and target data: this is not an issue here as we worked with 1000 genome data which was not part of the base data. However, this should be checked. Othwerwise, there will be inflation of the association of PRS with the target data. Practically, overlapping samples should be removed and base GWAS should be re-calculated. <br/>
* Relatedness: ideally no relatedness of second degree or closer within base, within target, and between base and target. R file for the target has used a King threshold of 0.0884. Therefore, individuals with second degree or closer familiar relationship have been removed. The base data has also mentioned PiHat > 0.2 which has to do with removal of related individuals. <br/>
* Similar ancestry: base and target data should be from the same ancestry which is the case here (European). <br/>

## Downloading phenotype (AN) data
* Find the phenotype dataset on Github and download it as a "raw" file in the main prstrain directory:
```bash
# Download phenotype file
cd ~/prstrain
wget https://raw.githubusercontent.com/sinaro/prs-training-uio/main/ANpheno_1000G_EUR
# Take a look at the phenotype file and check if the coding of phenotype file corresponds to that of documentation of plink .fam file.
```

## Running PRSice2 to construct PRS

```bash
# Example bash script to run on a personal device (create a file called "prsiceR_AN.sh" in the main prstrain directory with the content indicated below and execute it in bash):
# To check which commands are available in PRSice2 and their description, check: http://www.prsice.info/command_detail/
# You can use --fastscore to only calculate threshold stated in --bar-levels. If --fast-score is not used, it creates PRS for a range stated from --lower, to --upper, with the specified --interval.

# Run the prsice analysis with the Rscript
#!/bin/bash
Rscript PRSice.R \
        --prsice PRSice_linux \
        --base ~/prstrain/base/AN_basegwas.QC.gz \
        --snp ID \
        --a1 REF \
        --a2 ALT \
        --beta \
        --stat BETA \
        --pvalue PVAL \
        --target ~/prstrain/1000G/1000G_phase3_common_norel.nodup \
        --extract ~/prstrain/1000G/1000G_phase3_common_norel.nodup.clumped.clumped \
        --keep ~/prstrain/1000G/eurfamily.txt \
        --bar-levels 1e-08,1e-05,0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
        --lower 5e-08 \
        --interval 5e-05 \
        --upper 0.5 \
        --num-auto 22 \
        --no-clump \
        --pheno ~/prstrain/ANpheno_1000G_EUR \
        --pheno-col caco \
        --binary-target T \
        --seed 4059203949 \
        --out AN_1000GEUR_phen \
        --thread 1 \
        --device pdf \
        --quantile 10 \
        --print-snp \
        --ignore-fid \
        --all-score
```

* Check and discuss the output of Prsice2, and the generated files.

## Running PRSice2 on HPC (eg. TSD at UiO)
* a brief overview of offline and online HPC system, transferring of data/environment to an offline system, login node and submit node, modules, and Slurm.
* you need to be on the compute node (not the login node) [on TSD, "ssh -Y p****-submit"]
* you can search for available modules [module spider plink]
* if available, you can load the module and use it [module load plink/1.90b6.2]
* however, we do not have PRsice2 available as a module on TSD. The installation can be "ordered" to the helpdesk, or the binary file can be transferred and directly used.
* This requires tsd-api-client to be installed and registered on an online environment (given the offline nature of TSD), to upload the file to the tsd project (please get in touch if you need assistance with this).


```bash
# Example Slurm script to run on an HPC cluster:
#!/bin/bash
#SBATCH --jobname=prscie_prstrain_20220926
#SBATCH --account=[to be assigned]
#SBATCH --partition=[accel]
#SBATCH --time=00:15:00
#SBATCH --mem-per-cpu=1G
#SBATCH --cpus-per-task=1

# Set up job environment
module purge
module load R-bundle-Bioconductor/3.11.foss-2020a-R-4.0.0

# Run the prsice analysis with the Rscript
Rscript PRSice.R \
        # The rest of arguments (base etc) is the same as above.
```

## Questions and comments
Please do no hesitate to get in touch in case you have any questions or comments (sina.rostami@farmasi.uio.no). A 30 to 45 minute Zoom meeting is scheduled for Friday, October 7th at 13.00 for those who would like to discuss anything related to this session (if interested, send an email).

## Acknowledgements
I would like to thank Nasimeh Naseri who contributed to this document by testing compatibility with Mac OS. 

