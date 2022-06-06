# Training for construction of polygenic risk score using PRSice2
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

## Find and download GWAS summary results for Anorexia Nervosa (AN)
Get yourself familarize with [GWAS Catalog](https://www.ebi.ac.uk/gwas/), and search in PubMed for major GWAS in European populations for AN. Try to find the latest, and largest (highest sample size) GWAS for AN. Find if the GWAS summary results are publicly available and where you can download it.
<details>
<summary>Find the answer here</summary>

The largest GWAS for AN in European populatiosn as of June 2022 was published by [Watson et al](https://pubmed.ncbi.nlm.nih.gov/31308545/). Their summary results could be downloaded from the [PGS website](https://www.med.unc.edu/pgc/download-results/).

</details>
The largest GWAS for AN in European populatiosn as of June 2022 was published by [Watson et al](https://pubmed.ncbi.nlm.nih.gov/31308545/). Their summary results could be downloaded from the [PGS website](https://www.med.unc.edu/pgc/download-results/).


```bash
```







