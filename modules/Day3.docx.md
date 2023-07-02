# Advanced Polygenic Risk Score Analyses

## Day 3 - Polygenic Risk Score Analyses Workshop 2023

## Table of Contents

  1. [Key Learning Outcomes](#key-learning-outcomes)
  2. [Resources you will be using](#resources-you-will-be-using)
  3. [Datasets](#data-sets) 
  4. [Exercise 1: Estimating R2 in case/control studies](#Exercise-1:-Estimating-R2-in-case/control-studies)

## Key Learning Outcomes
After completing this practical, you should be able to:
  1. know how to adjust for ascertainment bias in case-control analysis
  2. Know how over-fitting aﬀects PRS results and how to handle it 
  3. understand distribution of PRS
  4. understand diﬀerent file formats involved in gene-set analysis
  5. understand diﬀerence between self-contained and competitive gene-set analyses
  6. Calculate pathway based PRS

## Resources you will be using 
To perform PRS analyses, summary statistics from Genome-Wide Association Studies (GWAS) are required. In this workshop, the following summary statistics are used:

|**Phenotype**|**Provider**|**Description**|**Download Link**|
|:---:|:---:|:---:|:---:|
|Height|[GIANT Consortium](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files)|GWAS of height on 253,288 individuals| [Link](https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz)|
|Coronary artery disease (CAD)|[CARDIoGRAMplusC4D Consortium](http://www.cardiogramplusc4d.org/)|GWAS on 60,801 CAD cases and 123,504 controls| [Link](http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.additive.Oct2015.pub.zip)|

## Additional Resources

|**Data Set**|**Description**|**Download Link**|
|:---:|:---:|:---:|
|Ensembl Human Genome GTF file|A file containing the coordinates for genes in the human genome. Used by PRSice to map the SNPs onto genic regions| [Link](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/) |
|MSigDB Gene Sets | File containing the gene-set information. *Free registration required.*| [Download here after registration](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.1/h.all.v6.1.symbols.gmt)|

In this practical, we will explore how to perform gene-set based PRS analyses. To perform this analysis, gene-set information and coordinates for the genic regions are required. These information can be obtained from the following database:


## Data Sets
You will find all practical materials in the **PRS_Workshop/Day_3** directory. Relevant materials that you should see there at the start of the practical are as follows:

 📂: Base_Data
  - GIANT_Height.txt,
  - cad.add.txt,
  - cad.add.readme.

 📂: Target_Data
  - TAR.fam
  - TAR.bim
  - TAR.bed
  - TAR.height
  - TAR.cad 
  - TAR.covariate
    
  📁: Reference files
   - Homo_sapiens.GRCh38.86.gt
   - Sets.gmt
     
 🛠️: Software
  - plink_mac
  - plink_linux
  - plink.exe
  - PRSice.R 
  - PRSice_mac
  - PRSice_linux
  - PRSice_win64.exe
    
---
> 
> ‼️ All target phenotype data in this worshop are **simulated**. They have no specific biological meaning and are for demonstration purposes only. 
> 
---
<a href="#top">[Back to Top](#table-of-contents)</a>

## Exercise 1: Estimating R2 in case/control studies
Bias in R2 estimation caused by ascertained case/control samples can be adjusted using the equation proposed by Lee et al (2011), which requires the sample prevalence (case/control ratio) and population prevalence as parameters. This function is implemented in PRSice and the adjustment can be performed by providing the population prevalence to the command --prevalence.

Residuals of logistic regression is not well defined, and in PRS analyses, Nagelkerke R2 is usually used to represent the model R2 (this is the default of PRSice). However, this R2 does not account for the diﬀerence between sample prevalence (i.e. case-control ratio) and population prevalence, which can lead to bias in the reported R2 (Fig.1.1a).

>
   **Figure 1.1 Performance of diﬀerent R2 when the study contains equal portion of cases and controls**
   ** (a) Nagelkerke R2  **
![Figure 1.1a](images/day3/images-004.png)
---

Bias in R2 estimation caused by ascertained case/control samples can be adjusted using the equation proposed by Lee et al. 2012 (Fig.2.1b), which requires the sample
prevalence (case/control ratio) and population prevalence as parameters. This function is implemented in PRSice and the adjustment can be performed by providing the population prevalence to the command --prevalence.



Now, account for the ascertainment of the case/control sample by including the population prevalence (let’s assume e.g. 5% here) in the PRSice command to obtain the adjusted (Lee) R2 :

