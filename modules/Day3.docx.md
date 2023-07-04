# Advanced Polygenic Risk Score Analyses

## Day 3 - Polygenic Risk Score Analyses Workshop 2023

## Table of Contents

  1. [Key Learning Outcomes](#key-learning-outcomes)
  2. [Resources you will be using](#resources-you-will-be-using)
  3. [Datasets](#data-sets) 
  4. [Exercise 1 Estimating R<sup>2</sup> in case and control studies](#exercise-1-estimating-r2-in-case-and-control-studies)
  5. [Exercise 2 Overfitting caused by model optimisation](#exercise-2-Overfitting-caused-by-model-optimisation)
  6. [Exercise 3 Distribution of PRS](#exercise-3-distribution-of-prs)

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

In this practical, we will explore how to perform gene-set based PRS analyses. To perform this analysis, gene-set information and coordinates for the genic regions are required. These information can be obtained from the following database:

|**Data Set**|**Description**|**Download Link**|
|:---:|:---:|:---:|
|Ensembl Human Genome GTF file|A file containing the coordinates for genes in the human genome. Used by PRSice to map the SNPs onto genic regions| [Link](https://ftp.ensembl.org/pub/release-109/gtf/homo_sapiens/) |
|MSigDB Gene Sets | File containing the gene-set information. *Free registration required.*| [Download here after registration](http://software.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/6.1/h.all.v6.1.symbols.gmt)|


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


## Exercise 1 Estimating R<sup>2</sup> in case and control studies
Bias in R<sup>2</sup> estimation caused by ascertained case/control samples can be adjusted using the equation proposed by Lee et al (2011), which requires the sample prevalence (case/control ratio) and population prevalence as parameters. This function is implemented in PRSice and the adjustment can be performed by providing the population prevalence to the command --prevalence.

Residuals of logistic regression is not well defined, and in PRS analyses, Nagelkerke R<sup>2</sup> is usually used to represent the model R<sup>2</sup> (this is the default of PRSice). However, this R<sup>2</sup> does not account for the diﬀerence between sample prevalence (i.e. case-control ratio) and population prevalence, which can lead to bias in the reported R<sup>2</sup> (Figure 1.1a). 
>
  **Figure 1.1 Performance of diﬀerent R<sup>2</sup> when the study contains equal portion of cases and controls**
>
  **(a) Nagelkerke R<sup>2</sup>**
![Figure 1.1a](/images/day3/images004.png)
---

Bias in R<sup>2</sup> estimation caused by ascertained case/control samples can be adjusted using the equation proposed by Lee et al. 2012 (Figure 1.1b), which requires the sample
prevalence (case/control ratio) and population prevalence as parameters. This function is implemented in PRSice and the adjustment can be performed by providing the population prevalence to the command --prevalence.
>
  **Figure 1.1 Performance of diﬀerent R<sup>2</sup> when the study contains equal portion of cases and controls**
>
  **(b) Lee adjusted R<sup>2</sup>**
![Figure 1.1b](/images/day3/images006.png)
---


Now, account for the ascertainment of the case/control sample by including the population prevalence (let’s assume e.g. 5% here) in the PRSice command to obtain the adjusted (Lee) R<sup>2</sup> :

```
Rscript ./Software/PRSice.R \
    --prsice Software/PRSice_linux \
    --base  Base_Data/cad.add.txt \
    --target Target_Data/TAR \
    --snp markername \
    --A1 effect_allele \
    --A2 noneffect_allele \
    --chr chr \
    --bp bp_hg19 \
    --stat beta \
    --beta \
    --pvalue p_dgc \
    --pheno Target_Data/TAR.cad \
    --prevalence 0.493 \
    --binary-target T \
--out Results/CAD.highres.LEER2
```
Check the *.summary file in the Results folder where you will find the usual (Nagelkerke) R<sup>2</sup> and the adjusted (Lee) R<sup>2</sup>.

>
  **Plot of CAD Lee R<sup>2</sup>**
>
> <img src="https://github.com/WCSCourses/prs_2023/blob/main/images/day3/CADLEER2HIGHRESPLOT.png"/>

---

---
>
> 
> 📌 To speed up the practical, we have generated a smaller gene-set file. If you want the full gene-set file, you can download it from the link above.
> 
> 
> 📌 All target phenotype data in this workshop are simulated. While they reflect the corresponding trait data, they have no specific biological meaning and are for demonstration purposes only.
---

---
> 
>❓Has accounting for the population prevalence aﬀected the R<sup>2</sup>?
> 
---
> 
>❓Would you expect a diﬀerence between the Nagelkerke R<sup>2</sup> and the Lee adjusted R<sup>2</sup> if the case/control ratio in the target sample reflects the disease prevalence in the population?
> 
---
<a href="#top">[Back to Top](#table-of-contents)</a>

## Exercise 2 Overfitting caused by model optimisation

In PRS analyses, the shrinkage or tuning parameter is usually optimized across a wide range of parametric space (e.g. P -value threshold, proportion of causal SNPs). When both optimisation and association testing are performed on the target data, over-fitted results will be obtained. The accuracy and predictive power of over-fitted results are likely to diminish when replicated in an independent data set.

A simple solution is to perform permutation to obtain an empirical P -value for the association model, which is implemented in PRSice. Briefly, permutation is performed as follows:
1) Compute the P -value in your original data, denoted as obs.p, at the "best" threshold.
2) Then shuﬄe the phenotype and obtain the P -value of the "best" threshold for this null phenotype, denoted as null.p
3) Repeat 2) N times
4) Calculate the empirical P-value as:
>
 $` Pemp = (\sum(obs.p > null.pi + 1) / (N + 1) `$
You will have to specify the number of permutation (N ) to perform by providing --perm N as a parameter to PRSice.

```
Rscript ./Software/PRSice.R \
    --prsice Software/PRSice_mac \
    --base  Base_Data/GIANT_Height.txt \
    --target Target_Data/TAR \
    --snp MarkerName \
    --A1 Allele1 \
    --A2 Allele2 \
    --stat b \
    --beta \
    --pvalue p \
    --pheno Target_Data/TAR.height \
    --binary-target F \
    --cov Target_Data/TAR.covariate \
    --cov-col Sex \
    --perm 1000 \
    --out Results/Height.perm
```
---
>
> 📝 **10000 permutations typically provide empirical P -values with high accuracy to the second decimal place (eg. 0.05), but smaller empirical P -values should be considered approximate.**
>
---
>
> ❓ What is the smallest possible empirical P-value when 10000 permutation are performed?
>
> ❓ Is the height PRS significantly associated with height after accounting for the over-fitting implicit in identifying the best-fit PRS? How about CAD?
>
---

The best way to avoid having results that are over-fit is to perform validation on an independent validation data set. We can perform validation of the previous height + covariate analysis with PRSice, using the independent VAL target sample as validation data and the "best" P-value threshold predicted in the VAL samples:

```
Rscript ./Software/PRSice.R \
    --prsice Software/PRSice_mac \
    --base  Base_Data/GIANT_Height.txt \
    --target Target_Data/VAL \
    --snp MarkerName \
    --A1 Allele1 \
    --A2 Allele2 \
    --stat b \
    --beta \
    --pvalue p \
    --pheno Target_Data/VAL.height \
    --binary-target F \
    --no-full \
    --bar-levels 0.0680001 \
    --fastscore \
    --cov Target_Data/VAL.covariate \
    --cov-col Sex \
    --out Results/Height.val
```
---
>
> ❓ Why do we use --bar-levels 0.0680001 --no-full and --fastscore in this script?
>
> ❓ How does the PRS R2 and P -value for the validation data set compare to the analysis on the TAR target data? Is this what you would expect? Why?
>
---

<a href="#top">[Back to Top](#table-of-contents)</a>

## Exercise 3 Distribution of PRS

Many PRS study publications include quantile plots that show an exponential increase in phenotypic value or / Odd Ratios (OR) among the top quantiles (e.g. an S-shaped quantile plot, e.g. Figure 1.3). 
>
  **Figure 1.3 An example of a S-shaped quantile plot**
>
![Figure 1.3](/images/day3/images021.png)
---

This might lead us to believe that individuals with PRS values in the top quantiles have a distinctly diﬀerent genetic aetiology compared to the rest of the sample, or that there is epistasis/interactions causing there substantially higher risk. However, when we plot a normally distributed variable (e.g. a PRS) as quantiles on the X-axis then we expect to observe this exponential pattern even when the X variable only has a linear eﬀect on the Y variable. This is because the top (and bottom) quantiles are further away from each other on the absolute scale of the variable and so the diﬀerences in their eﬀects are larger than between quantiles in the middle of the distribution.

To understand this more, we will perform a simple simulation using R:
```
R
# First, we define some simulation parameters
n.sample <- 10000
PRS.r2 <- 0.01
# Then, we simulate PRS that follow a random normal distribution
prs <- rnorm(n.sample)
# We can then simulate the phenotype using the following script
pheno <- prs + rnorm(n.sample,mean=0, sd=sqrt(var(prs)*(1-PRS.r2)/(PRS.r2)))
# We can examine the relationship between the phenotype and prs 
# using linear regression
summary(lm(pheno~prs))
# Which shows that we have the expected PRS R2
# Group the phenotype and PRS into a data.frame
info <- data.frame(SampleID=1:n.sample, PRS=prs, Phenotype=pheno)
# Then we can generate the quantile plot. 
# To save time, we will load in the quantile plot script from Software
source("./Software/Quantile.R")
# Then we can plot the quantile plot using quantile_plot function
quantile_plot(info, "Results/Height", 100)
```
