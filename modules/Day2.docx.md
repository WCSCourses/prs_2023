# Introduction to Polygenic Risk Scores

 # Day 2 Timetable


|**TIME**   |   **TITLE** |  **PRESENTER**     |
|      :---: |    :---:      |      :---:  |
|8:30 - 9:00| Recap and day overview  |                                  |
|9:00 - 10:30 | ***Lecture***: Introduction to PRS I |Carene Ndong Sima |
|             | ***Practical***: Performing QC + Computing PRS|  Carene + Marion            | 
|10:30 - 11:00 | Coffee break and Q&A |
|11:00 - 12:00| ***Practical***: PRS I - continued |Carene + Marion | 
|12:00 - 13:00| Seminar | Abdel |
|13:00 - 14:00| Lunch |
|14:00 - 15:30| ***Lecture***: Introduction to PRS II| Carene Ndong Sima|
|             | ***Practical***: PRS prediction + cross-trait analyses, visualizing results| Carene + Marion | 
|15:30 - 16:00| Break|
|16:00 - 17:00| ***Practical***: Intro to PRS II - continued| Carene + Marion |
|17:00 - 17:30| Seminar: PRS ethics consideration| |
|17:30 - 18:00| Feedback & reflecxion|

## Table of Contents

   1. [Introduction to Polygenic Score Analyses](#Introduction-to-polygenic-score-analyses)
        1. [Key Learning Outcomes](#key-learning-outcomes)
        2. [Resources you will be using](#resources-you-will-be-using)
        3. [Data Structure](#data-structure) 
        4. [Introduction](#introduction) 
        5. [Understanding GWAS Summary Statistics](#understanding-gwas-summary-statistics) 
        6. [Matching the Base and Target Data sets](#matching-the-base-and-target-data-sets) 
        7. [Linkage Disequilibrium in PRS Analyses](#linkage-disequilibrium-in-prs-analyses) 
           1. [Performing Clumping](#performing-clumping) 
       8. [P-Value Thresholding](#p-value-thresholding) 
           1. [Height PRS using GW-significant SNPs
            only](#height-prs-using-gw-significant-snps-only) 
           2. [Height PRS across multiple P-value
            thresholds](#height-prs-across-multiple-p-value-thresholds)   
          3. [High Resolution Scoring](#high-resolution-scoring)
       9. [Stratifying Samples by PRS](#stratifying-samples-by-prs) 
      10. [Case Control Studies](#case-control-studies) 
      11. [Cross-Trait Analysis](#cross-trait-analysis)

---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


# Introduction to Polygenic Score Analyses

## Key Learning Outcomes
After completing this practical, you should be able to:

1.  Perform basic Polygenic Risk Score (PRS) analyses using PRSice
    (Euesden, Lewis & O'Reilly 2015; Choi & O'Reilly 2019)
2.  Interpret results generated from PRS analyses
3.  Customise visualisation of results

## Resources you will be using
To perform PRS analyses, summary statistics from Genome-Wide Association Studies (GWAS) are required. In this workshop, the following summary statistics are used:

|**Phenotype**|**Provider**|**Description**|**Download Link**|
|:---:|:---:|:---:|:---:|
|Height|GIANT Consortium|GWAS of height on 253,288 individuals ([**wood_defining_2014**](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files))| [Download](https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz) [Link](https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz)|
|Coronary artery disease (CAD)|CARDIoGRAM plus C4D Consortium|GWAS on 60,801 CAD cases and 123,504 controls ([**consortium_comprehensive_2015**](http://www.cardiogramplusc4d.org/))|[Download](http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.additive.Oct2015.pub.zip) [Link](http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.additive.Oct2015.pub.zip)|

## Data Structure
You will find all practical materials in the **PRS_Workshop/Day_2** directory. Relevant materials that you should see there at the start of the practical are as follows:

 :file_folder: Base_Data
  - GIANT_Height.txt,
  - cad.add.txt,
  - cad.add.readme.

 :file_folder: Target_Data
  - TAR.fam
  - TAR.bim
  - TAR.bed
  - TAR.height
  - TAR.cad 
  - TAR.covariate

 :hammer_and_wrench: Software
  - plink_mac
  - plink_linux
  - plink.exe
  - PRSice.R 
  - PRSice_mac
  - PRSice_linux
  - PRSice_win64.exe

 >
 > ‼️ All target phenotype data in this worshop are **simulated**. They have no specific biological meaning and are for demonstration purposes only.
 >
---
<a href="#top">Back to top</a>


## Introduction
A PRS is a (usually weak) estimate of an individual's genetic propensity to a phenotype, calculated as a sum of their genome-wide genotypes weighted by corresponding genotype effect sizes obtained from GWAS summary statistics. In the next section we will consider what the effect size means and how it is used in computing PRS.

## Understanding GWAS Summary Statistics

 When GWAS are performed on a quantitative trait, the effect size is
 typically given as a beta coefficient (𝛽) from a linear regression
 with Single Nucleotide Polymorphism (SNP) genotypes as predictor of
 phenotype. The 𝛽 coefficient estimates the increase in the phenotype
for each copy of the *effect allele*. For example, if the effect
 allele of a SNP is **G** and the non-effect allele is **A** , then the
 genotypes **AA**, **AG** and **GG** will be coded as 0, 1 and 2
 respectively. In this scenario, the 𝛽 coefficient reflects how much
 the phenotype changes for each **G** allele present (NB. The 𝛽 can be
 positive or negative - so the 'effect allele' is simply the allele
 that was coded in the regression, not necessarily the allele with a
 positive effect).

 When a GWAS is performed on a binary trait (e.g. case-control study),
 the effect size is usually reported as an Odd Ratios (OR). Using the
 same example, if the OR from the GWAS is 2 with respect to the **G**
 allele, then the OR of **AG** relative to **AA** is 2, and the OR of
 **GG** relative to **AA** is 4. So an individual with the **GG**
 genotype are estimated\* to be 4 times more likely to be a case than
 someone with the **AA** genotype (\*an Odds Ratio is itself an
 estimate of a Risk Ratio, which cannot be calculated from a
 case/control study)

 >
 > 📜 The relationship between the 𝛽 coefficient from the logistic regression and the OR is: 𝑂𝑅 = 𝑒<sup>𝛽</sup> and 𝑙𝑜𝑔<sub>𝑒</sub>(𝑂𝑅) = 𝛽.
 >  While GWAS usually convert from the 𝛽 to the OR when reporting results, most PRS software convert OR back to 𝛽's(𝑙𝑜𝑔<sub>𝑒</sub>(𝑂𝑅)) to allow simple addition.
 > 
---
 >
 > 📜 Column names are not standardised across reported GWAS results, thus it is important to check which column is the effect (coded) allele and which is the non-effect allele. For example, in the height GWAS conducted by the GIANT consortium, the effect allele is in the column Allele1, while Allele2 represents the non-effect allele.
 >
 > ![](/images/Day2.docs_folder/images-006.png)
 > Let us open the Height GWAS file (**GIANT_Height.txt**) and inspect the SNPs at the top of the file. If we only consider SNPs *rs4747841* and *rs878177*, what will the ‘PRS’ of an individual with genotypes **AA** and **TC**, respectively, be? And what about for an individual with **AG** and **CC**, respectively? (Careful these are not easy to get correct! This shows how careful PRS algorithms/code need to be).
 > 
 > ❓❓What do these PRS values mean in terms of the height of those individuals?
 >
---
<a href="#top">Back to top</a>


## Matching the Base and Target Data sets

 The first step in PRS calculation is to ensure consistency between the GWAS summary statistic file (*base data*) and the target genotype file (*target data*). Since the base and target data are generated independently, they often relate to different SNPs and so the first job is to identify the overlapping SNPs across the two data sets and remove non-overlapping SNPs (this is usually done for you by PRS software). If the overlap is low then it would be a good idea to perform imputation on your target data to increase the number of SNPs that overlap between the data sets.
 
 The next, more tricky issue, is that the genotype encoding between the data sets may differ. For example, while the effect allele of a SNP is **T** in the base data, the effect allele in the target might be **G** instead. When this occurs, *allele flipping* should be performed,where the genotype encoding in the target data is reversed so that **TT**, **TG** and **GG** are coded as 2, 1 and 0. Again, this is usually performed automatically by PRS software.

![](media/image5.jpeg){width="0.3229166666666667in"
height="0.3229166666666667in"}

---
<a href="#top">Back to top</a>


## Linkage Disequilibrium in PRS Analyses

 GWAS are typically performed one-SNP-at-a-time, which, combined with the strong correlation structure across the genome (Linkage Disequilibrium (LD)), makes identifying the independent genetic effects (or their best proxies if these are not genotyped/imputed) challenging. There are two main options for approximating the PRS that would have been generated from full conditional GWAS: 1. SNPs are *clumped* so that the retained SNPs are largely independent of each other, allowing their effects to be summed, assuming additive effects, 2. all SNPs are included and the LD between them is accounted for.

 While option 2 is statistically appealing, option 1 has been most adopted in PRS studies so far, most likely due to its simplicity and the similarity of results of methods using the different options to date ([**mak_polygenic_2017**](https://doi.org/10.1002/gepi.22050)). In this workshop we will consider option 1, implemented in PRSice, but if you are interested in how LD can be incorporated as a parameter in PRS calculation then see the LDpred ([**vilhjalmsson_modeling_2015**](https://doi.org/10.1016/j.ajhg.2015.09.001))and lassosum ([**mak_polygenic_2017**](https://doi.org/10.1002/gepi.22050)) papers.

---
<a href="#top">Back to top</a>


### Performing Clumping

 *Clumping* is the procedure where a SNP data set is 'thinned' by removing SNPs across the genome that are correlated (in high LD) with a nearby SNP that has a smaller association 𝑃 -value. 

 SNPs are first sorted (i.e. ranked) by their 𝑃 -values. Then, starting from the most significant SNP (denoted as the *index SNP*), any SNPs in high LD (eg. 𝑟<sup>2</sup>>0.1, with 𝑟<sup>2</sup> typically calculated from *phased haplotype* data) with the index SNP are removed. To reduce computational burden, only SNPs that are within e.g. 250 kb of the *index SNP* are 𝑐𝑙𝑢𝑚𝑝𝑒𝑑. This process is continued until no *index SNPs* remain.

Use the command below to perform clumping of the Height GWAS data using PLINK([**chang_second_2015**](https://doi.org/10.1186/s13742-015-0047-8)). First, you will have to navigate to the right folder where the data are stored using the terminal. Open the terminal and type the command below at the terminal prompt:
```
 cd \~/Desktop/PRS\\\_Workshop/

#Next type the following command (NB. See warning below):

  ./Software/plink_linux \
  --bfile Target_Data/TAR \
  --clump Base_Data/GIANT_Height.txt \
  --clump-p1 1 \
  --clump-snp-field MarkerName \
  --clump-field p \
  --clump-kb 250 \
  --clump-r2 0.1 \
  --out Results/Height
```
![](media/image5.jpeg){width="0.3229155730533683in"
height="0.3229166666666667in"}

The command above performs clumping on the height GWAS using LD calculated based on the **TAR** genotype file. SNPs that have 𝑟<sup>2</sup>>0.1 within a 250 kb window of the index SNP are removed. This will generate the **Height.clumped** file, which contains the SNPs retained after clumping.

![](media/image7.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

---
<a href="#top">Back to top</a>

## P-Value Thresholding

 Deciding which SNPs to include in the calculation of PRS is one of the major challenges in the field. A simple and popular approach is to include SNPs according to their GWAS association 𝑃-value. For example, we may choose to include only the genome-wide significant SNPs from the GWAS because those are the SNPs with significant evidence for association. In the next subsection you will compute PRS from GW-significant SNPs only, and then in the subsequent subsection you will generate multiple PRSs using different 𝑃-value thresholds.

### Height PRS using GW-significant SNPs only

Use the commands below to run PRSice with GIANT Height GWAS as base data and the height phenotype as target data. PRSice will calculate Height PRS in the target data and then perform a regression of the Height PRS against the target individual's true height values. From the **PRS_Workshop/Day_2** directory, run the following command in the terminal:
```
Rscript ./Software/PRSice.R \
--prsice Software/PRSice_linux \
--base Base_Data/GIANT_Height.txt \
--target Target_Data/TAR \
--snp MarkerName \
--A1 Allele1 \
--A2 Allele2
--stat b \
--beta \
--pvalue p \
--pheno Target_Data/TAR.height \
--binary-target F \
--bar-levels 5e-8 \
--no-full \
--fastscore \
--out Results/Height.gws

```

This command takes the Height GWAS summary statistic file (\--base), informs PRSice of the column name for the column containing the SNP ID(\--snp), the effect allele (--A1), the non-effect allele (\--A2),the effect size (\--stat) and the 𝑃-value (\--pvalue). We also inform PRSice that the effect size is a 𝛽 coefficient (\--beta) instead of an OR. The \--binary-target F command informs PRSice that the target phenotype is a quantitative trait and thus linear regression should be performed. In addition, we ask PRSice not to perform high-resolution scoring over multiple thresholds (\--fastscore), and to compute the PRS using only those SNPs with 𝑃-value \< 5*×*10<sup>−8</sup>.

![](media/image6.jpeg){width="0.3281244531933508in"
height="0.3281244531933508in"}

PRSice performs strand flipping and clumping automatically and generates the **Height.gws.summary** file, together with other output that we will look into later in the practical. The summary file contains the following columns:

1.  **Phenotype** - Name of Phenotype.
2.  **Set** - Name of Gene Set. Default is *Base*
3.  **Threshold** - Best P-value Threshold
4.  **PRS.R2** - Variance explained by the PRS
5.  **Full.R2** - Variance explained by the full model (including the covariates)
6.  **Null.R2** - Variance explained by the covariates (none provided here)
7.  **Prevalence** - The population disease prevalence as indicated by the user (not provided here due to testing continuous trait)
8.  **Coefficient** - The 𝛽 coefficient corresponding to the effect estimate of the best-fit PRS on the target trait in the regression. A one unit increase in the PRS increases the outcome by 𝛽
9.  **Standard.Error** - The standard error of the best-fit PRS 𝛽 coefficient (see above)
10. **P** - The 𝑃 -value relating to testing the null hypothesis that the best-fit PRS 𝛽 coefficient is zero.
11. **Num_SNP** - Number of SNPs included in the best-fit PRS
12. **Empirical-P** - Only provided if permutation is performed. This is the empirical 𝑃 -value corresponding to the association test of the best-fit PRS - this controls for the over-fitting that occurs when multiple thresholds are tested.

For now, we can ignore most columns and focus on the **PRS.R2** and the **P** column, which provide information on the model fit.

![](media/image7.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

---
<a href="#top">Back to top</a>

### Height PRS across multiple P-value thresholds

A disadvantage of using only genome-wide significant SNPs is that there are likely to be many true signals among SNPs that did not reach genome-wide significance. However, since we do not know what 𝑃-value threshold provides the \"best\" prediction for our particular data, then we can calculate the PRS under several 𝑃-value thresholds and test their prediction accuracy to identify the \"best\" threshold (NB. See [**dudbridge_power_2013**](https://doi.org/10.1371/journal.pgen.1003348) for theory on factors affecting the best-fit PRS)

> []{#_bookmark13 .anchor}Figure 1.1: BARPLOT generated by PRSice

![](media/image8.png){width="2.7440616797900264in"
height="2.761874453193351in"}

This process is implemented in PRSice and can be performed automatically as follows:
```
Rscript ./Software/PRSice.R \
--dir . \
--prsice Software/PRSice_linux \
--base Base_Data/GIANT_Height.txt \
--target Target_Data/TAR \
--snp MarkerName \
--A1 Allele1 \
--A2 Allele2 \
--stat b \
--beta \
--pvalue p \
--pheno Target_Data/TAR.height \
--binary-target F \
--fastscore \
--out Results/Height.fast 
```

By removing the \--bar-levels and \--no-full command, we ask PRSice to perform PRS calculation with a number of predefined thresholds (0.001,0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1). The **.prsice** file is very similar to the **.summary** file, the only difference is that **.prsice** file reports the results of model fits for **all thresholds** instead of the most predictive threshold. This allow us to observe the change in model fitting across different thresholds, visualized by the BARPLOT (fig. [1.1)](#_bookmark13) generated by PRSice

> ![](media/image7.jpeg){width="0.3541655730533683in"
> height="0.3593744531933508in"}

---
<a href="#top">Back to top</a>


### High Resolution Scoring

If we limit ourselves to a small number of 𝑃-value thresholds, we might \"miss\" the most predictive threshold. In order to identify this \"best\" threshold, we will need \"high-resolution scoring\", that is, to test the predictive power of PRS generated under a larger number of p-value thresholds. We can achieve that by simply removing the \--fastscore command from the PRSice script:
```
 Rscript ./Software/PRSice.R
--dir .
--prsice Software/PRSice_linux
--base Base_Data/GIANT_Height.txt
--target Target_Data/TAR
--snp MarkerName
--A1 Allele1
--A2 Allele2
--stat b
--beta
--pvalue p
--pheno Target_Data/TAR.height
--binary-target F
--out Results/Height.highres
```

 When PRSice performs high-resolution scoring, it will generate a plot (fig. [1.2](#_bookmark15)) presenting the model fit of PRS calculated at all P-value thresholds.

> []{#_bookmark15 .anchor}
> Figure 1.2: High Resolution Plot generated by PRSice

![](media/image9.jpeg){width="2.463332239720035in"
height="2.455in"}![](media/image7.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

---
<a href="#top">Back to top</a>


**Accounting for Covariates**

When performing PRS, one might want to account for covariates. Based on user inputs, PRSice can automatically incorporate covariates into its model. For example, the following commands will include sex into the regression model as a covariate:
```
Rscript ./Software/PRSice.R
--prsice Software/PRSice_linux
--base Base_Data/GIANT_Height.txt
--target Target_Data/TAR
--snp MarkerName
--A1 Allele1
--A2 Allele2
--stat b
--beta
--pvalue p
--pheno Target_Data/TAR.height
--binary-target F
--cov Target_Data/TAR.covariate
--cov-col Sex
--out Results/Height.sex
```
When covariates are included in the analysis, PRSice will use the model fit of **only PRS** for all its output. This 𝑃𝑅𝑆.𝑅<sup>2</sup> is calculated by substracting the 𝑅<sup>2</sup> of the null model (e.g. 𝐻𝑒𝑖𝑔ℎ𝑡 *∼* 𝑆𝑒𝑥) from the 𝑅<sup>2</sup> of the full model (e.g.𝐻𝑒𝑖𝑔ℎ𝑡 *∼* 𝑆𝑒𝑥 + 𝑃𝑅𝑆)

> ![](media/image10.jpeg){width="0.34781167979002625in"
> height="0.34781167979002625in"}![](media/image7.jpeg){width="0.3541666666666667in"
> height="0.3593744531933508in"}


Usually, with categorical variables, dummy variables have to be generated to represent the different categories. Alternatively, residualized phenotype, generated by regresing the covariates against the phenotype, can be used for downstream analyses. A useful feature of PRSice is to automatically generate the dummy variable for users. This can be achieved with the following command:
```
Rscript.exe ./Software/PRSice.R
--prsice Software/PRSice_linux
--base Base_Data/GIANT_Height.txt
--target Target_Data/TAR
--snp MarkerName
--A1 Allele1
--A2 Allele2
--stat b
--beta
--pvalue p
--pheno Target_Data/TAR.height
--binary-target F
--cov Target_Data/TAR.covariate
--cov-col Sex
--cov-factor Sex
--out Results/Height.sex
```
> []{#_bookmark16 .anchor}
> Figure 1.3: Example of a quantile plot generated by PRSice

![](media/image11.png){width="2.4641655730533683in"
height="2.4674989063867017in"}

---
<a href="#top">Back to top</a>


## Stratifying Samples by PRS

 An interesting application of PRS is to test whether samples with higher PRS have higher phenotypic values. This can be nicely visualized using the quantile plot (fig. [1.3).](#_bookmark16)

To generate quantile plots in PRSice, simply add *\--quantile 10* option.

![](media/image6.jpeg){width="0.3281244531933508in"
height="0.3281244531933508in"}
```
Rscript ./Software/PRSice.R
--prsice Software/PRSice_linux
--base Base_Data/GIANT_Height.txt
--target Target_Data/TAR
--snp MarkerName
--A1 Allele1
--A2 Allele2
--stat b
--beta
--pvalue p
--pheno Target_Data/TAR.height
--binary-target F
--cov Target_Data/TAR.covariate
--cov-col Sex
--plot
--quantile 10
--out Results/Height.sex
```

![](media/image10.jpeg){width="0.34781167979002625in"
height="0.34781167979002625in"}
Figure 1.4: Example of a strata plot generated by PRSice

![](media/image12.jpeg){width="2.4641655730533683in"
height="2.4674989063867017in"}

 A disadvantage of the quantile plot is that it only separate samples into quantiles of equal size. However, it is sometimes interesting to investigate whether a specific strata (e.g. top 5% of samples),contain a higher PRS than the reference strata. For example, ([**mavaddat_prediction_2015**](https://doi:10.1093/jnci/djv036)) found that samples in the highest 1% of PRS distribution have a 2.81 increased OR of breast cancer when comparing to samples at the middle quantiles (40th to 60th percentile). We can mimic their table by using *\--quant-break*, which represents the upper bound of each strata, and *\--quant-ref*, which represents the upper bound of the reference quantile:
```
Rscript ./Software/PRSice.R
--prsice Software/PRSice_linux
--base Base_Data/GIANT_Height.txt
--target Target_Data/TAR
--snp MarkerName
--A1 Allele1
--A2 Allele2
--stat b
--beta
--pvalue p
--pheno Target_Data/TAR.height
--binary-target F
--cov Target_Data/TAR.covariate
--cov-col Sex
--plot
--quantile 100
--quant-break 1,5,10,20,40,60,80,90,95,99,100
--quant-ref 60
--out Results/Height.sex
```
![](media/image6.jpeg){width="0.3281244531933508in"
height="0.3281244531933508in"}

---
<a href="#top">Back to top</a>


## Case Control Studies

 In the previous exercises, we have performed PRS analyses on height, which is a quantitative trait. For binary phenotypes (e.g case-control) there are a number of differences in the analysis:

1.  Logistic regression has to be performed instead of linear regression
2.  ORs are usually provided and need to be converted to 𝛽's when
    constructing PRS

 Here, we will use CAD as an example. You will find the summary statistic under *Base_Data* (**cad.add.txt**) and the phenotype file (**TAR.cad**) under *Target_Data*. You will also need to specify *\--binary-target T* in the PRSice command to indicate that the phenotype is binary.

![](media/image10.jpeg){width="0.34781167979002625in"
height="0.34781167979002625in"}
```
Rscript ./Software/PRSice.R
--prsice Software/PRSice_linux
--base Base_Data/cad.add.txt
--target Target_Data/TAR
--snp markername
--A1 effect_allele
--A2 noneffect_allele
--chr chr
--bp bp_hg19
--stat beta
--beta
--pvalue p_dgc
--pheno Target_Data/CAD.pheno
--binary-target T
--out Results/CAD.highres
```
![](media/image10.jpeg){width="0.34781167979002625in"
height="0.34781167979002625in"}

> ![](media/image13.jpeg){width="3.5316666666666667in"
> height="2.066457786526684in"}
>
> []{#_bookmark19 .anchor}
> Figure 1.5: Plot taken from Ruderfer et al. 2014

![](media/image7.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

## Cross-Trait Analysis

 A popular application of PRS is in performing cross-trait analyses. This allows some interesting analyses such as those performed by ([**ruderfer_polygenic_2014**](https://doi:10.1038/mp.2013.138)) (fig. [1.5),](#_bookmark19) which used the bipolar PRS to predict into different clinical dimensions of schizophrenia.

 In this practical, we will perform cross-trait analyses between CAD and Height, using height as the base and CAD as the target.

![](media/image10.jpeg){width="0.34781167979002625in"
height="0.34781167979002625in"}
```
Rscript ./Software/PRSice.R
--prsice Software/PRSice_linux
--base Base_Data/GIANT_Height.txt
--target Target_Data/TAR
--snp MarkerName
--A1 Allele1
--A2 Allele2
--stat b
--beta
--pvalue p
--pheno Target_Data/CAD.pheno
--binary-target T
--out Results/Cross.highres
```
![](media/image7.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

---
<a href="#top">Back to top</a>

