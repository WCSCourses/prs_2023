# Advanced Polygenic Risk Score Analyses

## Day 3 - Polygenic Risk Score Analyses Workshop 2023

## Table of Contents

  1. [Key Learning Outcomes](#key-learning-outcomes)
  2. [Resources you will be using](#resources-you-will-be-using)
  3. [Datasets](#data-sets) 
  4. [Exercise 1 Estimating R<sup>2</sup> in case and control studies](#exercise-1-estimating-r2-in-case-and-control-studies)
  5. [Exercise 2 Overfitting caused by model optimisation](#exercise-2-Overfitting-caused-by-model-optimisation)
     1. [Out of Sample Validation](#out-of-sample-validation)
  7. [Exercise 3 Distribution of PRS](#exercise-3-distribution-of-prs)
  8. [Gene Set Analysis](#gene-set-analysis)
     1. [Molecular Signatures Database MSigDB](#molecular-signatures-Database-msigdb)
     2. [General Transfer Format file](#general-transfer-format-file)
     3. [Browser Extensible Data BED](#browser-extensible-data-bed)
  9. [Gene Set Enrichment Analysis](#gene-set-enrichment-analysis)
  10. [Exercise 4 Gene Set Based PRS Analysis](#exercise-4-gene-set-based-prs-analysis)
      

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

You will need to download the required files for this tutorial.       

[Download Backup WeTransfer](https://we.tl/t-YO8YbDjDK7)        

The data will be downloaded into your "Downloads" folder. You will need to move it to right directory, using the following commands.



```
cd data
mv ~/Downloads/Day_3.zip .
unzip Day_3.zip
cd Day_3
```

You will find all practical materials in the **data/Day_3** directory. Relevant materials that you should see there at the start of the practical are as follows:

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
  - PRSice.R 
  - PRSice_linux
  - nagelkerke.R
  - Quantile.R
    
---
> 
> ‼️ All target phenotype data in this worshop are **simulated**. They have no specific biological meaning and are for demonstration purposes only. 
> 
---
<a href="#top">[Back to Top](#table-of-contents)</a>


## Exercise 1 Estimating R<sup>2</sup> in case and control studies
Bias in R<sup>2</sup> estimation caused by ascertained case/control samples can be adjusted using the equation proposed by **Lee et al (2011)**, which requires the sample prevalence (case/control ratio) and population prevalence as parameters. This function is implemented in PRSice and the adjustment can be performed by providing the population prevalence to the command **--prevalence**.

Residuals of logistic regression is not well defined, and in PRS analyses, Nagelkerke R<sup>2</sup> is usually used to represent the model R<sup>2</sup> (this is the default of PRSice). However, this R<sup>2</sup> does not account for the diﬀerence between sample prevalence (i.e. case-control ratio) and population prevalence, which can lead to bias in the reported R<sup>2</sup> (Figure 1.1a). 
>
  **Figure 1.1: Performance of diﬀerent R<sup>2</sup> when the study contains equal portion of cases and controls**
>
  **(a) Nagelkerke R<sup>2</sup>** 
![Figure 1.1a](https://drive.google.com/uc?id=1r3LV442RQT3CWGSAnZsU6QT1DqHtRKkd)
---

Bias in R<sup>2</sup> estimation caused by ascertained case/control samples can be adjusted using the equation proposed by **Lee et al. 2011 (Figure 1.1b)**, which requires the sample
prevalence (case/control ratio) and population prevalence as parameters. This function is implemented in PRSice and the adjustment can be performed by providing the population prevalence to the command **--prevalence**.
>
  **Figure 1.1: Performance of diﬀerent R<sup>2</sup> when the study contains equal portion of cases and controls**
>
  **(b) Lee adjusted R<sup>2</sup>** 
![Figure 1.1b](https://drive.google.com/uc?id=19ACkFb2gr7EU4dfcsTDPtprwtYPDQufo)
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
--prevalence 0.05 \ 
--binary-target T \
--out Results/CAD.highres.LEER2
```
The results are written to the "Results" directory. Examine the results folder and each file that was generated. For more information about each file type, see  [here](https://choishingwan.github.io/PRSice/step_by_step/).

---
>
> ⭐ Check the *.summary file in the Results folder where you will find the usual (Nagelkerke) R<sup>2</sup> and the adjusted (Lee) R<sup>2</sup>.
>

>
  **Figure 1.2: Barplot of CAD Lee R<sup>2</sup>** 
![Figure 1.2](https://drive.google.com/uc?id=1ktv_UydzcZXUdup9oTB1LKkFyhCT9WdY)
---

---
>
> 
> 📌 To speed up the practical, we have generated a smaller gene-set file. If you want the full gene-set file, you can download it from the link above.
> 
> 📌 All target phenotype data in this workshop are simulated. While they reflect the corresponding trait data, they have no specific biological meaning and are for demonstration purposes only.
---

---
>
> 
>❓Has accounting for the population prevalence aﬀected the R<sup>2</sup>?
><details>
> <summary>Solution</summary>     
>
> Yes, the adjusted R<sup>2</sup> = 0.0521524 and default R<sup>2</sup> = 0.0442664      
>
></details>
>
> 
---
>
> 
>❓Would you expect a diﬀerence between the Nagelkerke R<sup>2</sup> and the Lee adjusted R<sup>2</sup> if the case/control ratio in the target sample reflects the disease prevalence in the population?
>
><details>
> <summary>Solution</summary>     
>
> No, the R<sup>2</sup> will be the same because the prevalence of the target is a true representation of the population prevalence.      
>
></details>
---
<a href="#top">[Back to Top](#table-of-contents)</a>

## Exercise 2 Overfitting caused by model optimisation

In PRS analyses, the shrinkage or tuning parameter is usually optimized across a wide range of parametric space (e.g. P -value threshold, proportion of causal SNPs). When both optimisation and association testing are performed on the target data, over-fitted results will be obtained. The accuracy and predictive power of over-fitted results are likely to diminish when replicated in an independent data set.

A simple solution is to perform permutation to obtain an empirical P -value for the association model, which is implemented in PRSice. Briefly, permutation is performed as follows:
1) Compute the P -value in your original data, denoted as obs.p, at the "best" threshold.
2) Then shuﬄe the phenotype and obtain the P -value of the "best" threshold for this null phenotype, denoted as null.p
3) Repeat 2) N times
4) Calculate the empirical P-value as:

 $` Pemp = (\sum(obs.p > null.pi + 1) / (N + 1) `$
---

You will have to specify the number of permutation (N ) to perform by providing --perm N as a parameter to PRSice.

```
Rscript ./Software/PRSice.R \
    --prsice Software/PRSice_linux \
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
  **Figure 1.3: Barplot of Height using 1000 permutations**
>
![Figure 1.3](/images/day3/Height.perm_BARPLOT_2023-06-30.png)
>
---

---
>
> 📝 **10000 permutations typically provide empirical P-values with high accuracy to the second decimal place (eg. 0.05), but smaller empirical P-values should be considered approximate.**
>
---
>
> ❓ What is the smallest possible empirical P-value when 10000 permutation are performed? 
>
><details>
> <summary>Solution</summary>     
>
> 1.5 X 10<sup>-34</sup>.      
>
></details>
>
> ❓ Is the height PRS significantly associated with height after accounting for the over-fitting implicit in identifying the best-fit PRS? How about CAD?
>
><details>
> <summary>Solution</summary>     
>
> Yes, the height PRS is significantly associated with height. After accounting for the over-fitting implicit in identifying the best-fit PRS, the emprical p-value is 0.000999001.
>
></details>
---
### Out of Sample Validation

The best way to avoid having results that are over-fit is to perform validation on an independent validation data set. We can perform validation of the previous height + covariate analysis with PRSice, using the independent VAL target sample as validation data and the "best" P-value threshold predicted in the VAL samples:

```
Rscript ./Software/PRSice.R \
    --prsice Software/PRSice_linux \
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
  **Figure 1.4: Barplot of Height validation dataset** 
![Figure 1.4](/images/day3/Height.val_BARPLOT_2023-06-30.png)
---
---
>
> ❓ Why do we use --bar-levels 0.0680001 --no-full and --fastscore in this script?
>
> ❓ How does the PRS R2 and P -value for the validation data set compare to the analysis on the TAR target data? Is this what you would expect? Why?
>
---

<a href="#top">[Back to Top](#table-of-contents)</a>

## Exercise 3 Distribution of PRS

Many PRS study publications include quantile plots that show an exponential increase in phenotypic value or / Odd Ratios (OR) among the top quantiles (e.g. an S-shaped quantile plot, e.g. Figure 1.6). 
>
  **Figure 1.5: An example of density plot for PRS**
>
![Figure 1.5](/images/day3/images020.png)
---
>
  **Figure 1.6: An example of a S-shaped quantile plot**
>
![Figure 1.6](/images/day3/images021.png)
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
>
  **Figure 1.7: The resulting quantile plot**
>
![Figure 1.7](/images/day3/Height_QUANTILES_PLOT_2023-07-05.png)
---

---
>
> ❓ What is the shape of the resulting quantile plot?
>
> ❓ Try plotting the densities of the height or CAD PRS in R * - do they look normally distributed? Why? (*Hint: You can generate a density plot for the
PRS in R using plot(density(x)) where x is a vector of the PRS values in the sample).
>
---

<a href="#top">[Back to Top](#table-of-contents)</a>

## Gene Set Analysis
Currently, most PRS analyses have been performed on a genome-wide scale, disregarding the underlying biological pathways. Udler et al. 2018 suggest that grouping Single Nucleotide Polymorphisms (SNPs) into biological functional groups can lead to PRS that are more relevant to clinical risk. In this practical, we will go through some common file formats for gene-set analysis and will then calculate some gene-set (or pathway) based PRS.

### Molecular Signatures Database MSigDB
The Molecular Signatures Database (MSigDB) oﬀers an excellent source of gene-sets, including the hallmark genes, gene-sets of diﬀerent biological processes, gene-sets of diﬀerent oncogenic signatures etc. All gene-sets from MSigDB follows the Gene Matrix Transposed file format (GMT), which consists of one line per gene-set, each containing at least 3 column of data:

| | | | | |
|:---:|:---:|:---:|:---:|:---:|
|Set A| Description | Gene 1 | Gene 2 | ...
|Set A| Description | Gene 1 | Gene 2 | ...

---
> ** Have a look at the Reference/Sets.gmt file. **
>
> ❓ How many gene-sets are there in the Reference/Sets.gmt file? 
>
> ❓ How many genes does the largest gene-set contain?
>
---
>
> 💬 While you can read the GMT file using Excel. You should be aware that Excel has a tendency to convert gene names into dates (e.g. SEPT9 to Sep-9)
>
---

As GMT format does not contain the chromosomal location for each individual gene, an additional file is required to provide the chromosoaml location such that SNPs can be map to genes.

### General Transfer Format file
The General Transfer Format (GTF) file contains the chromosomal coordinates for each gene. It is a **tab** separated file and all but the final field in each feature line must contain a value. "Empty" columns should be denoted with a ‘.’. You can read the full format specification here. One column that might be of particular interest is column 3: **feature**, which indicates what feature that line of GTF represents. This allows us to select or ignore features that are of interest.

You can find the description of each feature [here](http://www.sequenceontology.org/browser/obob.cgi).


### Browser Extensible Data BED
Browser Extensible Data (BED) file (diﬀerent to the binary ped file from PLINK) is a file format to define genetic regions. It contains 3 required fields per line (chromosome, start coordinate and end coordinate) together with 9 additional optional field. A special property of BED is that it is a 0-based format, i.e. chromosome starts at 0, as opposed to the usual 1-based format such as the PLINK format. For example, a SNP on chr1:10000 will be represented as:

| | | |
|:---:|:---:|:---:|
|**1**|**9999**|**10000**|

---
>
> ❓ How should we represent the coordinate of rs2980300 (chr1:785989) in BED format?
>
---

<a href="#top">[Back to Top](#table-of-contents)</a>


## Gene Set Enrichment Analysis

Now we have gone through all the files involved in gene-set analysis, we should consider one of the most important aspects of gene-set (or pathway) enrichment
analyses, which is the diﬀerent types of testing that we can perform when doing them:

### Self-Contained vs Competitive Testing
The null-hypothesis of self-contained and competitive test statistics is diﬀerent:
  – **Self-Contained** - None of the genes within the gene-set are associated with the phenotype
  – **Competitive** - Genes within the gene-set are no more associated with the phenotype than genes outside the gene-set
Therefore, a bigger gene-set will have a higher likelihood of having a significant P -value from self-contained test, which is not desirable.


## Exercise 4 Gene Set Based PRS Analysis

Having learnt about the basics of gene-set analyses, we are now ready to perform gene-set association analyses using PRSet.

To perform the PRSet analysis and obtain the set based PRS and competitive P-value, simply provide the GTF file and the GMT file to PRSice and specify the number of permutation for competitive P-value calculation using the --set-perm option.

```
Rscript ./Software/PRSice.R \
    --prsice Software/PRSice_linux  \
    --base Base_Data/GIANT_Height.txt \
    --target Target_Data/TAR \
    --A1 Allele1 \
    --A2 Allele2 \
    --snp MarkerName \
    --pvalue p \
    --stat b \
    --beta \
    --binary-target F \
    --pheno Target_Data/TAR.height \
    --cov Target_Data/TAR.covariate \
    --out Results/Height.set \
    --gtf Reference/Homo_sapiens.GRCh38.86.gtf \
    --wind-5 5kb \
    --wind-3 1kb \
    --msigdb Reference/Sets.gmt \
    --multi-plot 10 \
    --set-perm 1000
```

>
  **Figure 1.8: An example of the multi-set plot. Sets are sorted based on their self-contained R2 . Base is the genome wide PRS**
>
![Figure 1.8](/images/day3/Height.set_MULTISET_BARPLOT_2023-06-30.png)
---
>
> 📌 If the --wind-5 and --wind-3 flag is not specified, PRSet will use the exact coordinates of each gene as the boundary. By specifying eg. --wind-5 5kb and --wind-3 1kb then the boundary of each gene will be extended 5 kb towards the 5’ end and 1 kb towards the 3’ end so that regulatory elements of the gene can be included.
>
> 📍 By default, when calculating set based PRS, PRSet will not perform P -value thresholding. This is because the aim in gene-set analyses is to assess the overall signal in each gene-set, and compare which is most enriched for signal, rather than optimise predictive power as typically desirable for genome-wide PRS. Providing any of the following commands will activate P -value thresholding for set based PRS calculation: --lower, --upper, --inter, --bar-levels, --fastscore
>
---

>
> ❓ Can you plot the relationship between the gene-set R2 and the number of SNPs in each gene-set? What general trend can be seen?
>
> ❓ Considering the plot, what gene-sets do you think are most interesting and why?
>
> ❓ Why is it useful to have polygenic scores measured across gene-sets (or pathways) for individuals? Isn’t it suﬃcient to just obtain a ranking of gene-sets according to GWAS-signal enrichment?
>
---
