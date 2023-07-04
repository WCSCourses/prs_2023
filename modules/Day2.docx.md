# Introduction to Polygenic Risk Scores

## Table of Contents
      
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


## Key Learning Outcomes
After completing this practical, you should be able to:

1.  Perform basic Polygenic Risk Score (PRS) analyses using PRSice:
    ([Euesden, Lewis & O'Reilly 2015](https://doi:10.1093/bioinformatics/btu848);[Choi & O'Reilly 2019](https://doi:10.1093/gigascience/giz082))
2.  Interpret results generated from PRS analyses
3.  Customise visualisation of results

## Resources you will be using
To perform PRS analyses, summary statistics from Genome-Wide Association Studies (GWAS) are required. In this workshop, the following summary statistics are used:

|**Phenotype**|**Provider**|**Description**|**Download Link**|
|:---:|:---:|:---:|:---:|
|Height|[GIANT Consortium](https://portals.broadinstitute.org/collaboration/giant/index.php/GIANTconsortium_data_files)|GWAS of height on 253,288 individuals| [Link](https://portals.broadinstitute.org/collaboration/giant/images/0/01/GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz)|
|Coronary artery disease (CAD)|[CARDIoGRAMplusC4D Consortium](http://www.cardiogramplusc4d.org/)|GWAS on 60,801 CAD cases and 123,504 controls| [Link](http://www.cardiogramplusc4d.org/media/cardiogramplusc4d-consortium/data-downloads/cad.additive.Oct2015.pub.zip)|

## Data Structure
You will find all practical materials in the **/home/manager/data/Day2_Target_Data/Day2_Target_Data** directory. Relevant materials that you should see there at the start of the practical are as follows:

> First, download the Base datasets using this [link](https://www.dropbox.com/s/w7g75bo069tupuq/Day2_Base_Data.zip?dl=0).
> The data will be downloaded into your "Downloads" folder. You will need to move it to right directory, using the following command.
> ```
> cd ~/data/Day2_Target_Data/Day2_Target_Data
> ```
> ```
> mv ~/Downloads/Day2_Base_Data.zip ~/data/Day2_Target_Data/Day2_Target_Data
> ```
> Now, your Day2_Base_dat.zip has been moved to your Day2_Target_Data directory. However, the file is zipped so you will need to unzip it:
> ```
> unzip Day2_Base_Data.zip
> ```
>
> For clarity purposes, rename your Day2_Base_data to Base_data and make a directory for your Target data called "Target_data" in which you will move all your target datasets.
> ```
> mv Day2_Base_Data Base_Data
> ```
> ```
> mkdir Target_Data
> ```
> ```
> mv TAR.* Target_Data
> ```
>
> You also want to create a "Results" directory where all your results will be saved
> ```
> mkdir Results
> ```
>So now in your current working directory (for Day2), you should have 3 directories (in blue) called Base_data, Target_data, and Results
>

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
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


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

---
>
> 📜 The relationship between the 𝛽 coefficient from the logistic regression and the OR is: 
> **OR = *e*<sup>𝛽</sup>** and 
> **log<sub>*e*</sub>(OR)** **=** 𝛽.
> 
>  While GWAS usually convert from the 𝛽 to the OR when reporting results, most PRS software convert OR back to 𝛽's(𝑙𝑜𝑔<sub>𝑒</sub>(𝑂𝑅)) to allow simple addition.
>
> 📜 Column names are not standardised across reported GWAS results, thus it is important to check which column is the effect (coded) allele and which is the non-effect allele. For example, in the height GWAS conducted by the GIANT consortium, the effect allele is in the column Allele1, while Allele2 represents the non-effect allele.
> 
---
>
> 🔎 Let us open the Height GWAS file (**GIANT_Height.txt**) and inspect the SNPs at the top of the file. If we only consider SNPs *rs4747841* and *rs878177*, what will the ‘PRS’ of an individual with genotypes **AA** and **TC**, respectively, be? And what about for an individual with **AG** and **CC**, respectively? (Careful these are not easy to get correct! This shows how careful PRS algorithms/code need to be).
><details>
>   <summary>Answer</summary>     
>
> In the GIANT_Height.txt, SNPs effect sizes are reported as 𝛽 coefficient and measures the effect of the 'effect allele'. So, first check identify which allele is the effect allele for the SNPs of interest
>```
>grep -E 'rs4747841|rs878177' ./Base_data/GIANT_Height.txt
>```
> rs4747841 A G 0.551 -0.0011 0.0029 0.7 253213
>
> rs878177 T C 0.3 0.14 0.0031 8.2e-06 251271      
>      
> Second, multiply the weight of each risk allele by their dosage     
> Individual_1: PRS=?,  
> Individual_2: PRS=?
>      
> </details>

---
>
>❓What do these PRS values mean in terms of the height of those individuals?
>
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


## Matching the Base and Target Data sets

 The first step in PRS calculation is to ensure consistency between the GWAS summary statistic file (*base data*) and the target genotype file (*target data*). Since the base and target data are generated independently, they often relate to different SNPs and so the first job is to identify the overlapping SNPs across the two data sets and remove non-overlapping SNPs (this is usually done for you by PRS software). If the overlap is low then it would be a good idea to perform imputation on your target data to increase the number of SNPs that overlap between the data sets.
 
 The next, more tricky issue, is that the genotype encoding between the data sets may differ. For example, while the effect allele of a SNP is **T** in the base data, the effect allele in the target might be **G** instead. When this occurs, *allele flipping* should be performed,where the genotype encoding in the target data is reversed so that **TT**, **TG** and **GG** are coded as 2, 1 and 0. Again, this is usually performed automatically by PRS software.

---
>
‼️For SNPs that have complementary alleles, e.g. **A/T**, **G/C**, we cannot be certain that the alleles referred to in the target data correspond to those of the base data or whether they are the 'other way around' due to being on the other DNA strand (unless the same genotyping chip was used for all data). These SNPs are known as ***ambiguous SNPs***, and while allele frequency information can be used to match the alleles, we remove ambiguous SNPs in PRSice to avoid the possibility of introducting unknown bias.
>  
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


## Linkage Disequilibrium in PRS Analyses

 GWAS are typically performed one-SNP-at-a-time, which, combined with the strong correlation structure across the genome (Linkage Disequilibrium (LD)), makes identifying the independent genetic effects (or their best proxies if these are not genotyped/imputed) challenging. There are two main options for approximating the PRS that would have been generated from full conditional GWAS: 1. SNPs are *clumped* so that the retained SNPs are largely independent of each other, allowing their effects to be summed, assuming additive effects, 2. all SNPs are included and the LD between them is accounted for.

 While option 2 is statistically appealing, option 1 has been most adopted in PRS studies so far, most likely due to its simplicity and the similarity of results of methods using the different options to date ([**Mak, T et al., 2017**](https://doi.org/10.1002/gepi.22050)). In this workshop we will consider option 1, implemented in PRSice, but if you are interested in how LD can be incorporated as a parameter in PRS calculation then see the LDpred ([**Vilhjalmsson, B et al., 2015**](https://doi.org/10.1016/j.ajhg.2015.09.001))and lassosum ([**Mak, T et al., 2017**](https://doi.org/10.1002/gepi.22050)) papers.

---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


### Performing Clumping

 *Clumping* is the procedure where a SNP data set is 'thinned' by removing SNPs across the genome that are correlated (in high LD) with a nearby SNP that has a smaller association 𝑃 -value. 

 SNPs are first sorted (i.e. ranked) by their 𝑃 -values. Then, starting from the most significant SNP (denoted as the *index SNP*), any SNPs in high LD (eg. 𝑟<sup>2</sup>>0.1, with 𝑟<sup>2</sup> typically calculated from *phased haplotype* data) with the index SNP are removed. To reduce computational burden, only SNPs that are within e.g. 250 kb of the *index SNP* are 𝑐𝑙𝑢𝑚𝑝𝑒𝑑. This process is continued until no *index SNPs* remain.

Use the command below to perform clumping of the Height GWAS data using PLINK([**Chang, C et al., 2015**](https://doi.org/10.1186/s13742-015-0047-8)). First, you will have to navigate to the right folder where the data are stored using the terminal. Open the terminal and type the command below at the terminal prompt:
```
  plink \
  --bfile Target_Data/TAR \
  --clump Base_Data/GIANT_Height.txt \
  --clump-p1 1 \
  --clump-snp-field MarkerName \
  --clump-field p \
  --clump-kb 250 \
  --clump-r2 0.1 \
  --out Results/Height
```
---
>
> ‼️You can copy & paste code from this document directly to the terminal, but this can cause problems (e.g. when opened by Preview in Mac) and distort the code. Try using Adobe Reader or first copy & pasting to a text editor (e.g. notepad) or use the script file provided that contains all the commands.
> 
---

The command above performs clumping on the height GWAS using LD calculated based on the **TAR** genotype file. SNPs that have 𝑟<sup>2</sup>>0.1 within a 250 kb window of the index SNP are removed. This will generate the **Height.clumped** file, which contains the SNPs retained after clumping.

---
>
>❓How many SNPs were in the GIANT_Height.txt file before clumping?
><details>
> <summary>Answer</summary>     
>
> 2,183,049 variants     
>
></details>
>
>❓How many SNPs remain after clumbing?
><details>
> <summary>Answer</summary>     
>
> 109,496 variants     
>
></details>
>
>❓If we change the r<sup>2</sup> threshold to 0.2, how many SNPs remain? Why are there, now, more SNPs remaining?
><details>
> <summary>Answer</summary>     
>
> 162,275 variants     
>
></details>
> 
>❓Why is clumping performed for calculation of PRS? (in the standard approach).
> 
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>

## P-Value Thresholding

 Deciding which SNPs to include in the calculation of PRS is one of the major challenges in the field. A simple and popular approach is to include SNPs according to their GWAS association 𝑃-value. For example, we may choose to include only the genome-wide significant SNPs from the GWAS because those are the SNPs with significant evidence for association. In the next subsection you will compute PRS from GW-significant SNPs only, and then in the subsequent subsection you will generate multiple PRSs using different 𝑃-value thresholds.

### Height PRS using GW-significant SNPs only

Use the commands below to run PRSice with GIANT Height GWAS as base data and the height phenotype as target data. PRSice will calculate Height PRS in the target data and then perform a regression of the Height PRS against the target individual's true height values. From the **/home/manager/data/Day2_Target_Data/Day2_Target_Data** directory, run the following command in the terminal:
```
Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/GIANT_Height.txt --target Target_Data/TAR --snp MarkerName --A1 Allele1 --A2 Allele2 --stat b --beta --pvalue p --pheno Target_Data/TAR.height --binary-target F --bar-levels 5e-8 --no-full --fastscore --out Results/Height.gws
```

This command takes the Height GWAS summary statistic file (\--base), informs PRSice of the column name for the column containing the SNP ID(\--snp), the effect allele (--A1), the non-effect allele (\--A2),the effect size (\--stat) and the 𝑃-value (\--pvalue). We also inform PRSice that the effect size is a 𝛽 coefficient (\--beta) instead of an OR. The \--binary-target F command informs PRSice that the target phenotype is a quantitative trait and thus linear regression should be performed. In addition, we ask PRSice not to perform high-resolution scoring over multiple thresholds (\--fastscore), and to compute the PRS using only those SNPs with 𝑃-value \< 5*×*10<sup>−8</sup>.

---
> 
> 📜The default of PRSice is to perform clumping with r<sup>2</sup> threshold of 0.1 and a window size of 250kb.
> 
> To see a full list of command line options available in PRSice, type: 
> 
>  ```
> ~/PRSice/PRSice_linux -h
>  ```
>  
>  Take some time to have a look through some of these user options. By looking at the user options, work out which user option or options were used to ensure that the command above only calculated 1 PRS at genome-wide significance of \< 5*×*10<sup>−8</sup>.
>  
---

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

![Figue](/images/Day2.docxfolder/Height.gws_BARPLOT_2023-07-03.png)

---
>
> ❓What is the R<sup>2</sup> for the PRS constructed using only genome-wide significant SNPs?
><details>
> <summary>Answer</summary>     
>
> 0.071     
>
></details>
>
> ❓What is the P-value for the association between the PRS and the outcome? Is this significant? (explain your answer)
><details>
> <summary>Answer</summary>     
>
> 1.36e-09     
>
></details>
>
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>

### Height PRS across multiple P-value thresholds

A disadvantage of using only genome-wide significant SNPs is that there are likely to be many true signals among the SNPs that did not reach genome-wide significance. However, since we do not know what 𝑃-value threshold provides the \"best\" prediction for our particular data, then we can calculate the PRS under several 𝑃-value thresholds and test their prediction accuracy to identify the \"best\" threshold (NB. See [**Dudbridge, F 2013**](https://doi.org/10.1371/journal.pgen.1003348) for theory on factors affecting the best-fit PRS)
>
>

This process is implemented in PRSice and can be performed automatically as follows:
```
Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/GIANT_Height.txt --target Target_Data/TAR --snp MarkerName --A1 Allele1 --A2 Allele2 --stat b --beta --pvalue p --pheno Target_Data/TAR.height --binary-target F --fastscore --out Results/Height.fast 
```

By removing the \--bar-levels and \--no-full command, we ask PRSice to perform PRS calculation with a number of predefined thresholds (0.001,0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1). The **.prsice** file is very similar to the **.summary** file, the only difference is that **.prsice** file reports the results of model fits for **all thresholds** instead of the most predictive threshold. This allow us to observe the change in model fitting across different thresholds, visualized by the BARPLOT (<a href="#top">[Figure 1.1](#_Figure-1.1)</a>) generated by PRSice

>
>   **Figure 1.1: BARPLOT generated by PRSice**   
![Figure 1.1](/images/Day2.docxfolder/images-018.png)

---
>
>❓Which is the most predictive threshold?
><details>
>  <summary>Answer</summary>     
>
> 0.05     
>
></details>
>
>❓What is the R<sup>2</sup> of the most predictive threshold and how does it compare to PRS generated using genome-wide significant SNPs?
><details> 
> <summary>Answer</summary>     
>
> 0.26523     
>
></details>
>
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


### High Resolution Scoring

If we limit ourselves to a small number of 𝑃-value thresholds, we might \"miss\" the most predictive threshold. In order to identify this \"best\" threshold, we will need \"high-resolution scoring\", that is, to test the predictive power of PRS generated under a larger number of p-value thresholds. We can achieve that by simply removing the \--fastscore command from the PRSice script:
```
 Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/GIANT_Height.txt --target Target_Data/TAR --snp MarkerName --A1 Allele1 --A2 Allele2 --stat b --beta --pvalue p --pheno Target_Data/TAR.height --binary-target F --out Results/Height.highres
```

 When PRSice performs high-resolution scoring, it will generate a plot ([Figure 1.2.1](#_bookmark15)) presenting the model fit of PRS calculated at all P-value thresholds.
>
   **Figure 1.2.1: High resolution Plot generated by PRSice**
![Figure 1.2.1](/images/Day2.docxfolder/images-022.png)

   **Figure 1.2.2: BARPLOT generated by PRSice**
   
![Figure 1.2.2](/images/Day2.docxfolder/Height.highresBARPLOT2023-06-23.png) 
 
---
>
>❓Which is the most predictive threshold?
> <details>
>  <summary>Answer</summary>     
>
> 0.058   
>
></details>
>
>❓How much better is the threshold identified using high-resolution scoring, in terms of model R<sup>2</sup>?
><details>
> <summary>Answer</summary>     
>
> The R<sup>2</sup> using high-resolution scoring is 0.268237. In this case, there is not significant difference between the two results.        
> 
>
></details>
>
>❓How does running fastscore vs high-resolution change the most predictive threshold identified?
>
>❓The default of PRSice is to iterate the P-value threshold from \< 5*×*10<sup>−8</sup> to 0.5 with a step size of 0.00005, and to inlcude the P-value threshold of 1. Can you identify the commands controlling these parameters?
> <details>
>  <summary>Hint</summary>     
>
> ./Software/PRSice_linux -h
>
></details>
>
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


**Accounting for Covariates**

When performing PRS, one might want to account for covariates. Based on user inputs, PRSice can automatically incorporate covariates into its model. For example, the following commands will include sex into the regression model as a covariate:
```
Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/GIANT_Height.txt --target Target_Data/TAR --snp MarkerName --A1 Allele1 --A2 Allele2 --stat b --beta --pvalue p --pheno Target_Data/TAR.height --binary-target F --cov Target_Data/TAR.covariate --cov-col Sex --out Results/Height.sex
```
When covariates are included in the analysis, PRSice will use the model fit of **only PRS** for all its output. This 𝑃𝑅𝑆.𝑅<sup>2</sup> is calculated by substracting the 𝑅<sup>2</sup> of the null model (e.g. 𝐻𝑒𝑖𝑔ℎ𝑡 *∼* 𝑆𝑒𝑥) from the 𝑅<sup>2</sup> of the full model (e.g.𝐻𝑒𝑖𝑔ℎ𝑡 *∼* 𝑆𝑒𝑥 + 𝑃𝑅𝑆)

---
> 
> ℹ️ Usually, Principal Components (PCs) are also included as covariates in the analysis to account for population structure and it can be tedious to type all 20 or 40 PCs (e.g. *PC1*, *PC2*,...,*PC20*). In PRSice, you can add @ in front of the *--cov-col* string to activate the automatic substitution of numbers. If @ is found in front of the *--cov-col* string, any numbers within **\[** and **\]** will be parsed. E.g. @PC\[1-3] will be read as *PC1*, *PC2*, *PC3*. Discontinuous input are also supported: @cov\[1.3-5] will be parsed as *cov1*, *cov3*, *cov4*, *cov5*. You can also mix it up, e.g. @PC\[1-3]. Sex will be interpreted as *PC1*, *PC2*, *PC3*, *Sex* by PRSice.
>
---
>
>❓How does the inclusion of sex as covariate change the results?
> <details>
>  <summary>Answer</summary>     
>
> The inclusion of sex as a covariate increased the phenotypic variance explained by both PRS and sex (FULL.R2 = 0.47)  
>
></details>
>
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


## Stratifying Samples by PRS

 An interesting application of PRS is to test whether samples with higher PRS have higher phenotypic values. This can be nicely visualized using the quantile plot ([Figure 1.3](#_bookmark16))

To generate quantile plots in PRSice, simply add *\--quantile 10* option.

---
>
>📜 We can skip the PRS calculation using the --plot option, which will use previoulsy calculated PRS to generate the plots
>
---

```
Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/GIANT_Height.txt --target Target_Data/TAR --snp MarkerName --A1 Allele1 --A2 Allele2 --stat b --beta --pvalue p --pheno Target_Data/TAR.height --binary-target F --cov Target_Data/TAR.covariate --cov-col Sex --plot --quantile 10 --out Results/Height.sex
```

 **Figure 1.3: Example of a quantile plot generated by PRSice**
![Figure 1.3](/images/Day2.docxfolder/images-030.png)

---
>
> ℹ️ The --plot option tells PRSice to generate the plots without re-running the whole PRSice analysis. This is handy when you want to change some of the parameters for plotting e.g. the number of quantiles. **Try running the previous command with 20 quantiles, and again with 50 quantiles (checking the quantile plot each time**.
> 
---

 A disadvantage of the quantile plot is that it only separate samples into quantiles of equal size. However, it is sometimes interesting to investigate whether a specific strata (e.g. top 5% of samples),contain a higher PRS than the reference strata. For example, ([**Mavaddat, N et al., 2015**](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4754625/?report=reader)) found that samples in the highest 1% of PRS distribution have a 2.81 increased OR of breast cancer when comparing to samples at the middle quantiles (40th to 60th percentile). We can mimic their table by using *\--quant-break*, which represents the upper bound of each strata, and *\--quant-ref*, which represents the upper bound of the reference quantile:
```
Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/GIANT_Height.txt --target Target_Data/TAR --snp MarkerName --A1 Allele1 --A2 Allele2 --stat b --beta --pvalue p --pheno Target_Data/TAR.height --binary-target F --cov Target_Data/TAR.covariate --cov-col Sex --plot --quantile 100 --quant-break 1,5,10,20,40,60,80,90,95,99,100 --quant-ref 60 --out Results/Height.sex
```

 **Figure 1.4: Example of a strata plot generated by PRSice**
![Figure 1.4](/images/Day2.docxfolder/images-036.png)


---
>
>📜 See the quantile results in Table from in the *_QUANTILES_* file, and the plots in the *_QUANTILES_PLOT_* file. Due to the small sample size of the target data the results here are underwhelming, but with high power we may observe strong deviation in the extreme quantiles.
>
---
<a href="#top">[Back to Table of Contents](#table-of-contents)</a>


## Case Control Studies

 In the previous exercises, we have performed PRS analyses on height, which is a quantitative trait. For binary phenotypes (e.g case-control) there are a number of differences in the analysis:

1.  Logistic regression has to be performed instead of linear regression
2.  ORs are usually provided and need to be converted to 𝛽's when
    constructing PRS

 Here, we will use CAD as an example. You will find the summary statistic under *Base_Data* (**cad.add.txt**) and the phenotype file (**TAR.cad**) under *Target_Data*. You will also need to specify *\--binary-target T* in the PRSice command to indicate that the phenotype is binary.

---
>
>ℹ️ GWAS summary statistics for binary traits tend to report the OR instead of the 𝛽 coefficient, in which case the --or should be used. However, CARDIoGRAM plus C<sub>4</sub>D consortium provided the 𝛽 coefficient, therefore we will include --beta in our code and specify --binary-target T to indicate that the phenotype is binary.
>
---

```
Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/cad.add.txt --target Target_Data/TAR --snp markername --A1 effect_allele --A2 noneffect_allele --chr chr --bp bp_hg19 --stat beta --beta --pvalue p_dgc --pheno Target_Data/CAD.pheno --binary-target T --out Results/CAD.highres
```
---
>
>ℹ️ --chr and --bp inform PRSice the columns containing the chromosomal coordinates. This enables PRSice to check whether the SNPs in the Base and Target data have the same chromosomal coordinate.
>
---

   **Figure 1.5: BARPLOT generated from PRSice** 
![Figure 1.5](/images/Day2.docxfolder/CAD.highresBARPLOT2023-06-23.png)

---
>
>❓What is the R<sup>2</sup> and P-value of the best-fit PRS?
><details>
> <summary>Answer</summary>     
>
> 0.39 and 0.001, respectively
>
></details>
>
>❓Does this suggest that there is a significant association between the CAD PRS and CAD status in the target sample?
><details>
> <summary>Answer</summary>     
>
> The p-value does not suggest a significant association between the CAD PRS and CAD status.
>
></details>
>
---
## Cross-Trait Analysis

 A popular application of PRS is in performing cross-trait analyses. This allows some interesting analyses such as those performed by ([**Ruderfer, D et al., 2014**](https://www.nature.com/articles/mp2013138)) which used the bipolar PRS to predict into different clinical dimensions of schizophrenia.

 In this practical, we will perform cross-trait analyses between CAD and Height, using height as the base and CAD as the target.

---
>
>ℹ️ We will only focus on the simplest form of cross-trait analysis in this practical. To perform the multi-phenotype cross-trait analysis similar to that of ([**Ruderfer, D et al., 2014**](https://www.nature.com/articles/mp2013138)), you can use the --pheno-col to include multiple target phenotype into the analysis.
> 
---
```
Rscript ~/PRSice/PRSice.R --prsice ~/PRSice/PRSice_linux --base Base_Data/GIANT_Height.txt --target Target_Data/TAR --snp MarkerName --A1 Allele1 --A2 Allele2 --stat b --beta --pvalue p --pheno Target_Data/CAD.pheno --binary-target T --out Results/Cross.highres
```
---

**Figure 1.6: BARPLOT generated from PRSice** 
![Figure 1.6](/images/Day2.docxfolder/Cross.highBARPLOT2023-06-23.png)


>❓What is the R<sup>2</sup> for the most predictive threshold when using height as the base phenotype and CAD as the target phenotype?
><details>
> <summary>Answer</summary>     
>
> 0.032
>
></details>
>
>❓Now try using CAD as the base to predict height as the target trait? What is the PRS R<sup>2</sup> for that?
><details>
> <summary>Answer</summary>     
>
> 
>
></details>
>
---
<a href="#top">Back to top</a>

