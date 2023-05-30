> **Polygenic Risk Score Analyses Workshop 2023**

![](media/image1.jpeg){width="2.94in" height="2.8874989063867016in"}

> **Day 1: Introduction to GWAS & PRS**

> **Day 1 Timetable**

> 9:00 - 9:15 Welcome Address Dr Daneshwar and Dr Baichoo
>
> 9:15 - 9:30 Opening Speech from Organisers Dr Segun Fatumo and
>
> Dr Nicki Tiﬃn
>
> 9:30 - 10:30 [Lecture]{.underline}: Background to PRS: GWAS &
>
> relevant Statistics
>
> Dr Paul O'Reilly
>
> 10:30 - 11:00 Coﬀee Break and Q&A -
>
> 11:00 - 12:00 [Practical]{.underline}: Introduction to Bash and R Dr
> Paul O'Reilly & Tu-
>
> tors
>
> 12:00 - 13:30 Lunch -
>
> 13:30 - 15:00 [Practical]{.underline}: Introduction to PLINK I -
>
> Basics
>
> Dr Conrad Iyegbe & Tutors
>
> 15:00 - 15:30 Coﬀee Break and Q&A -
>
> 15:30 - 16:30 [Practical]{.underline}: Introduction to PLINK II -
>
> QC & GWAS
>
> Dr Conrad Iyegbe & Tutors
>
> **Contents**

[Day 1 Timetable](#day-1-timetable) 1

[D](#day-1-timetable)[ay 1 Timetable](#time-title-presenter) . . . . . .
. . . . . . . . . . . . . . . . . . . . . . . . . 1

1.  [Introduction to PLINK I: Basics](#introduction-to-plink-i-basics) 3

    1.  [Key Learning Outcomes](#key-learning-outcomes) . . . . . . . .
        . . . . . . . . . . . . . . . . 3

    2.  [Introduction](#introduction) . . . . . . . . . . . . . . . . .
        . . . . . . . . . . . . . . 3

    3.  [Command line basics](#command-line-basics) . . . . . . . . . .
        . . . . . . . . . . . . . . . . 3

    4.  [Exploring Data Sets](#exploring-data-sets) . . . . . . . . . .
        . . . . . . . . . . . . . . . . 4

    5.  [Recoding alleles as counts](#recoding-alleles-as-counts) . . .
        . . . . . . . . . . . . . . . . . . . . 5

    6.  [PLINK website](#plink-website) . . . . . . . . . . . . . . . .
        . . . . . . . . . . . . . 6

    7.  [Write SNP list and extract
        SNPs](#write-snp-list-and-extract-snps) . . . . . . . . . . . .
        . . . . . . . 6

2.  [Introduction to PLINK II: Performing QC &
    > GWAS](#introduction-to-plink-ii-performing-qc-gwas) 8

    1.  [Key Learning Outcomes](#key-learning-outcomes-1) . . . . . . .
        . . . . . . . . . . . . . . . . . 8

    2.  [Generate summaries to perform
        QC](#generate-summaries-to-perform-qc) . . . . . . . . . . . . .
        . . . . . 8

        1.  [Individual missingness](#individual-missingness) . . . . .
            . . . . . . . . . . . . . . . . 8

        2.  [SNP Missingness](#snp-missingness) . . . . . . . . . . . .
            . . . . . . . . . . . . 9

        3.  [Hardy-Weinberg Equilibrium](#hardy-weinberg-equilibrium) .
            . . . . . . . . . . . . . . . . . 9

        4.  [Allele frequencies](#allele-frequencies) . . . . . . . . .
            . . . . . . . . . . . . . . . 9

    3.  [Apply QC filters](#apply-qc-filters) 10

        1.  [Apply individual missingness
            thresholds](#apply-individual-missingness-thresholds) 10

        2.  [Apply SNP missingness and MAF
            thresholds](#apply-snp-missingness-and-maf-thresholds) 11

        3.  [Apply Hardy-Weinberg
            thresholds](#apply-hardy-weinberg-thresholds) 11

    4.  [Perform GWAS](#perform-gwas) 12

        1.  [Case/Control GWAS - no
            covariates](#casecontrol-gwas---no-covariates) 12

        2.  [Case/Control GWAS - with
            covariates](#casecontrol-gwas---with-covariates) 12

[License](#license) 14

# Practical 1: Introduction to PLINK I basics

## Key Learning Outcomes

> After completing this practical, you should be able to:

1.  Explore and generate genetic data sets needed for GWAS

2.  Recode and reorder allelic data

3.  Use the PLINK website

4.  Select and exclude lists of samples and SNPs
----
⚠️ All data used in this workshop are **simulated**. They have no specific biological meaning and are for demonstration purposes only.

----
## Introduction

PLINK is the most popular software program for performing genome-wide association analyses it is extremely extensive, allowing a huge
number of analyses to be performed. It also includes many options for reformatting your data and provides useful data summaries. Software
packages are usually best learnt by having a go at running some of their basic applications and progressing from there (rather than
reading the entire user manual first!) - so we begin by running some basic PLINK commands and then work steadily towards performing more
sophisticated analyses through these PLINK tutorials.

## Command line basics

> In all of the instructions below, *italics* indicate commands - they can be directly copy and pasted. 
> Anything in between the symbols *\<\>* needs to be changed in some way. For example, \<file_name\>
> indicates that you should replace that entire statement (including the *\<\>* symbols) with the appropriate file name. 
> **Bold** indicates non- command-line instructions (e.g. **right-click**)

### Let's begin

1.  Open up a terminal

2.  Navigate to the Practical/ folder and then the Day 1 working
    directory (cd \<directory_name\>)

3.  List all files in this directory by typing *ls*

4.  Test PLINK with no input by typing *./Software/plink*

5.  Note that you can see lots of PLINK options by using the built-in
    help function: > *./Software/plink --help*
----    
 📝 Calling PLINK with no output will test if PLINK is installed and available in the directory, because you should see some 
 output showing the PLINK license and some commands. If you do not see this, then please ask for help now!
 
 ----


## Exploring Data Sets

1.  Open an Explorer window ('Finder' on a Mac) and navigate to your PLINK working directory.

----
📝 An explorer window should show the same files as the *ls* command

----

2.  Open the file called 'D1D.map' with a Text Editor e.g. by typing
    **right-click \> Open**.

3.  Open the file 'D1D.ped'. Note this is a large file - if it will not
    open or is very slow, skip this step.

4.  Go to the PLINK website and investigate the format of the MAP/PED
    files (<http://zzz.bwh.harvard.edu/plink/download.shtml)>
    (Look in the blue column on the left side)
    
 ----  
 󠁑
❓ What are the 4 columns in the map file?

   What are the first 6 columns in a ped file?
   
   What information is in the remaining columns of the ped file?
   
----

5.  Create 'binary' format PLINK files using the recode command:

> 1 *./Software/plink*

> 2 *--file Data/D1D*

> 3 *--make-bed*

> 4 *--out Data/D1D*

6.  List files (*ls*) and check which new files have appeared

7.  Open and examine files ending .bim and .fam. Do not open the .bed
    file.

8.  Open and skim the '.log' file.

----
❓ How is the fam file similar to the ped file? How is it different? 

   How is the bim file similar to the map file? How is it different?
   
  (Use the PLINK website if necessary)
  
----
## Recoding alleles as counts

Genotype data in allele count format is very useful, for example to use in regression modelling in statistical software such as R.
Generate the D1D data in allele count format:

> 1 *./Software/plink
> 
> 2 *--bfile Data/D1D
> 
> 3 *--recodeA
> 
> 4 *--out Data/D1D_AC*

----
📝 There are several options for recoding SNPs in different ways - more information on the PLINK website (see next section).
    Again note that a log file was created - skim the log file or screen output
    
----

## PLINK website

Go to <http://zzz.bwh.harvard.edu/plink/download.shtml> and skim through the front page to get an idea of PLINK's functionality. Note
the list of clickable links on the left side of the website. 

Under 'Data Management' (click the heading on the left) and read the list of the diﬀerent ways you may want to recode and reorder data sets. Don't attempt to read much further as this is a very large and detailed section - a useful future resource but too much for today. 

Under 'Data Management', click 'Write SNP list' and read the instructions there to write SNP lists.

## Write SNP list and extract SNPs

You will now use the information that you found on the PLINK website to create a command to extract a list of SNPs. Below is a list of requirements - try to do this before you go to the end of this section, where the full command is given and explained.

1.  Set the D1D binary file as input

2.  Set MAF threshold to 0.05

3.  Set SNP missingness threshold to 0.05

4.  Add the appropriate command to write out a snp list containing only
    those SNPs with MAF above 0.05 and missingness below 0.05

5.  Use 'D1D_snps' as the output file name

6.  After the command has run, check the output for your SNP list and
    look at it with the default viewer.

You will now use the SNP list that you have created to extract those SNPs and create a new set of data files in a single command.

1.  Use the D1D binary file set as input

2.  Find the command for extracting a set of SNPs listed in a file (hint: Data Management section) and combine it with a command that
    you learned above to create binary files

3.  Use the output file name 'D1D_MAF_MISS'

**Solutions:**

> 1 *./Software/plink
>
> 2  *--bfile Data/D1D
>
> 3 *--maf 0.05
>
> 4 *--geno 0.05
>
> 5 *--write-snplist
>
> 6 *--out Data/D1D_snps*


> 1 *./Software/plink
>
> 2 *--bfile Data/D1D
>
> 3 *--extract Data/D1D_snps.snplist
>
> 4 *--make-bed
>
> 5 *--out Data/D1D_MAF_MISS*

# Introduction to PLINK II: Performing QC & GWAS

## Key Learning Outcomes

> After completing this practical, you should be able to:

1.  Generate summaries of the data needed for QC

2.  Apply QC thresholds

3.  Perform GWAS

## Generate summaries to perform QC

> There are many kinds of summaries of the data that can generated in
> PLINK in order to perform particular quality control (QC) steps, which
> help to make our data more reliable. Some of these involve summaries
> in relation to the individuals (e.g. individual missingness,
> sex-check) and some relate to summaries of SNP data (e.g. MAF,
> Hardy-Weinburg Equilibrium). Over the next few sub-sections you will
> go through some examples of generating summary statistics that can be
> used to perform QC.

### Individual missingness

1.  Use the D1D binary files to generate files containing missingness
    information (*\--*missing). Use the output file name 'D1D_miss'

2.  Open the 2 files that were generated (lmiss & imiss).

![](media/image4.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### SNP Missingness

1.  Use the D1D binary files to generate files containing missingness
    information (*\--*missing). Use the output file name 'D1D_miss'

2.  Look inside the file containing SNP missingness information:
    D1D_miss.lmiss.

![](media/image4.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

### Hardy-Weinberg Equilibrium

1.  Generate HWE statistics using the *\--*hardy option. Use output file
    name D1D_hardy.

2.  Open and examine results.

![](media/image4.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Allele frequencies

1.  Generate allele frequencies using the command *\--*freq. Use
    D1D_freq as the output name.

2.  Examine the output.

![](media/image4.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

## Apply QC filters

#### There are diﬀerent strategies for performing QC on your data: {#there-are-diﬀerent-strategies-for-performing-qc-on-your-data .unnumbered}

(a) create lists of SNPs and individuals and use *\--*remove,
    *\--*extract, *\--*exclude,

> *\--*include to create new file sets (good for documentation,
> collaboration)

(b) apply thresholds one at a time and generate new bed/bim/fam files
    (good for applying sequential filters)

(c) use options (e.g. *\--*maf ) in other commands (e.g. *\--*assoc) to
    remove SNPs or samples at required QC thresholds during analysis.

![](media/image3.jpeg){width="0.3281244531933508in"
height="0.3281244531933508in"}

### Apply individual missingness thresholds

1.  Generate new binary file sets (*\--*make*-*bed) from the 'D1D'
    binary file set, removing individuals with missingness greater than
    3% using a single command (hint: In the 'Inclusion thresholds'
    section, see the 'Missing/person' sub-section). Use the output file
    name 'D1D_imiss3pc'

2.  Examine the output files (no need to open, and remember the bed file
    can not be read) and the log file

> ![](media/image4.jpeg){width="0.3541655730533683in"
> height="0.3593744531933508in"}

### Apply SNP missingness and MAF thresholds

1.  Create new binary file sets from the 'D1D_imiss3pc' binary file set
    (NOT the original D1D files) by setting MAF threshold to 0.05 and
    SNP missingness threshold to 0.02 (See 'Inclusion thresholds' to
    obtain the correct threshold flags). Use the output file name
    'D1D_imiss3pc_lmiss2pc_maf5pc

2.  Examine the output files and the log file

![](media/image4.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Apply Hardy-Weinberg thresholds

1.  Generate a new binary file set called 'D1D_QC' from the
    D1D_imiss3pc_lmiss2pc_maf5pc file, applying a HWE threshold of
    0.0001.

2.  This is our final, QC'ed file set.

3.  Examine log and output files.

![](media/image4.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

## Perform GWAS

### Case/Control GWAS - no covariates

> Run the following code, which performs a genetic association study
> using logistic regression on some case/control data:
>
> 1 ./Software/plink
>
> 2 \--plink
>
> 3 \--bfile D1D_QC
>
> 4 \--logistic
>
> 5 \--adjust
>
> 6 \--pheno D1D.pheno1
>
> 7 \--out Results/D1D_CC

![](media/image4.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Case/Control GWAS - with covariates

> Here we repeat the previous analysis but this time including some
> covariates. The file D1D.pcs1234 contains the first 4 principal
> components from a PCA on the ge- netic data.
>
> 1\. Run the analysis specifying the covariates file:
>
> 1 ./Software/plink
>
> 2 \--plink
>
> 3 \--bfile D1D_QC
>
> 4 \--logistic
>
> 5 \--adjust
>
> 6 \--pheno D1D.pheno1
>
> 7 \--covar D1D.pcs.1234
>
> 8 \--out Results/D1D_CC_PCadj

![](media/image4.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

# License {#license .unnumbered}

> This work is licensed under Creative Commons
> Attribution-NonCommercial-ShareAlike
>
> 4.0 International Public License and the below text is a summary of
> the main terms of the full Legal Code (the full licence) available at
> [https://creativecommons.](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode)
> [org/licenses/by-nc-sa/4.0/legalcode](https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode).

#### You are free to: {#you-are-free-to .unnumbered}

-   **Share** --- copy and redistribute the material in any medium or
    > format

-   **Adapt** --- remix, transform, and build upon the material

> The licensor cannot revoke these freedoms as long as you follow the
> license terms.

#### Under the following terms: {#under-the-following-terms .unnumbered}

-   **Attribution** --- You must give appropriate credit, providea link
    > to the license, and indicate if changes were made. You may do so
    > in any reasonable manner, but not in any way that suggests the
    > licensor endorses you or your use.

-   **NonCommercial** --- You may not use the material for commercial
    > purposes.

-   **ShareAlike** --- If you remix, transform, or build upon the
    > material, you must distribute your contributions under the same
    > license as the original.

> No additional restrictions --- You may not apply legal terms or
> technological mea- sures that legally restrict others from doing
> anything the license permits.

#### Notices: {#notices .unnumbered}

> You do not have to comply with the license for elements of the
> material in the public domain or where your use is permitted by an
> applicable exception or limitation.
>
> No warranties are given. The license may not give you all of the
> permissions neces- sary for your intended use. For example, other
> rights such as publicity, privacy, or moral rights may limit how you
> use the material.
