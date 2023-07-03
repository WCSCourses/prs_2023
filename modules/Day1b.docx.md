**Polygenic Risk Score Analyses Workshop 2023**
# Practical 1
## Introduction to PLINK I: basics

## Key Learning Outcomes

After completing this practical, you should be able to:

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

In all of the instructions below:
- Anything in between the symbols *\<\>* needs to be changed in some way. For example, \<file_name\>
  indicates that you should replace that entire statement (including the *\<\>* symbols) with the appropriate file name. 
- **Bold** indicates non- command-line instructions (e.g. **right-click**)

### Let's begin

1.  Open up a terminal

2.  Navigate to the Day 1b working directory
    
            cd ~/data/Data_Day1b/
            
3.  List all files in this directory by typing 
   
            ls

4.  Test PLINK with no input by typing 

         plink

5.  Note that you can see lots of PLINK options by using the built-in
    help function: 
    
         plink --help
         
----    
 📝 Calling PLINK with no output will test if PLINK is installed and available in the directory, because you should see some 
 output showing the PLINK license and some commands. If you do not see this, then please ask for help now!
 
 ----

## Exploring Data Sets

1.  Open an Explorer window ('Finder' on a Mac) and navigate to your PLINK working directory.

----
📝 An explorer window should show the same files as the *ls* command

----

2.  Open the file called 'D1D.map' with a Text Editor e.g. by typing **right-click \> Open**.

3.  Open the file 'D1D.ped'. Note this is a large file - if it will not open or is very slow, skip this step.

4.  Go to the PLINK website http://zzz.bwh.harvard.edu/plink/download.shtml and investigate the format of the MAP/PED files 
    (Look in the blue column on the left side)

What do you observe?
- What are the 4 columns in the map file?
- What are the first 6 columns in a ped file?
- What information is in the remaining columns of the ped file?
   
5.  Create 'binary' format PLINK files using the recode command:

         plink --file D1D --make-bed --out D1D

6.  List files (*ls*) and check which new files have appeared

7.  Open and examine files ending .bim and .fam. Do not open the .bed
    file.

8.  Open and skim the '.log' file.

What do you observe?

- How is the fam file similar to the ped file?
- How is it different?
- How is the bim file similar to the map file?
- How is it different?
   (Use the PLINK website if necessary)
  
## Recoding alleles as counts

Genotype data in allele count format is very useful, for example to use in regression modelling in statistical software such as R.
Generate the D1D data in allele count format:

       plink --bfile D1D --recodeA --out D1D_AC

📝 There are several options for recoding SNPs in different ways - more information on the PLINK website (see next section).
    Again note that a log file was created - skim the log file or screen output

Look inside the .raw file. 

- What do you think the 0/1/2 represent?
- Do there appear to be more 0s or 2s?
- Why might this be?
    
## PLINK website

Go to http://zzz.bwh.harvard.edu/plink/download.shtml and skim through the front page to get an idea of PLINK's functionality. Note
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

----
📝 Log files are uselful to check that the number of SNPs and samples is as expected. Always check your your log files to ensure that they are sensible.
SNP lists can also be used to EXCLUDE SNPs - select 'exclude' above instead of 'extract'. 
Sample ID lists can also be used to 'keep' or 'remove' individuals in the same 'filter' window. Note that both sample IDs (FID IID,separated by a space are required in the sample file list.

----

**Solution 1:**

TO BE REVEALED LATER!!

**Solution 2:**

TO BE REVEALED LATER!!
       
# Practical 2 

## Introduction to PLINK II: Performing QC & GWAS

## Key Learning Outcomes

After completing this practical, you should be able to:

1.  Generate summaries of the data needed for QC

2.  Apply QC thresholds

3.  Perform GWAS

## Generate summaries to perform QC

There are many kinds of summaries of the data that can generated in PLINK in order to perform particular quality control (QC) steps, which
help to make our data more reliable. Some of these involve summaries in relation to the individuals (e.g. individual missingness, sex-check) and some 
relate to summaries of SNP data (e.g. MAF, Hardy-Weinburg Equilibrium). Over the next few sub-sections you will go through some examples of generating 
summary statistics that can be used to perform QC.

### Individual missingness

1.  Use the D1D binary files to generate files containing missingness information (--missing). Use the output file name 'D1D_miss'

2.  Open the 2 files that were generated (lmiss & imiss).

- What do the two output files contain?
- In the imiss file, what is the meaning of the data in the column headed "F_MISS"?


### SNP Missingness

1.  Use the D1D binary files to generate files containing missingness information (--missing). Use the output file name 'D1D_miss'

2.  Look inside the file containing SNP missingness information: D1D_miss.lmiss.

- What is the meaning of the value under F_MISS?
- What does the command --test-missing do and why might it be useful?

### Hardy-Weinberg Equilibrium

1.  Generate HWE statistics using the --hardy option. Use output file name D1D_hardy.

2.  Open and examine results.

- Why are there multiple rows for each SNP and what does each mean?
- Which of the rows do you think should be used to exclude SNPs from the subsequent analysis (if any) for failing the HWE test? Why?

### Allele frequencies

1.  Generate allele frequencies using the command *\--*freq. Use
    D1D_freq as the output name.

2.  Examine the output.

- What is the heading of the column that tells you which nucleotide is the minor allele?
  
 📝 **This information is important to remember as many PLINK files use this notation. The minor allele is always labeled the same way** 


## Apply QC filters

#### There are diﬀerent strategies for performing QC on your data:

(a) Create lists of SNPs and individuals and use --remove, --extract, --exclude, --include to create new file sets (good for documentation, collaboration)

(b) Apply thresholds one at a time and generate new bed/bim/fam file (good for applying sequential filters)

(c) Use options (e.g. --maf ) in other commands (e.g. --assoc) to remove SNPs or samples at required QC thresholds during analysis.

----
📝 We have already seen how to select or exclude individuals or SNPs by first creating lists (a), so in this section we will set thresholds to generate new files 
   sets in a single command. However, it is useful to have lists of all SNPs and individuals excluded pre-analysis, according to the reason for exclusion, so 
   generating and retaining such files using the techniques that we used before for good practice.

----

### Apply individual missingness thresholds

1.  Generate new binary file sets (--make-bed) from the 'D1D' binary file set, removing individuals with missingness greater than 3% using a single command 
    (hint: In the 'Inclusion thresholds' section, see the 'Missing/person' sub-section). Use the output file name 'D1D_imiss3pc'

2.  Examine the output files (no need to open, and remember the bed file cannot be read) and the log file

- How many individuals were in the original file?
- How many individuals were removed?
- How many males and females were left after screening?

### Apply SNP missingness and MAF thresholds

1.  Create new binary file sets from the 'D1D_imiss3pc' binary file set (NOT the original D1D files) by setting MAF threshold to 0.05 and
    SNP missingness threshold to 0.02 (See 'Inclusion thresholds' to obtain the correct threshold flags). Use the output file name'D1D_imiss3pc_lmiss2pc_maf5pc

2.  Examine the output files and the log file

- How many SNPs were in the original files?
- How many SNPs were removed for low minor allele frequency?
- How many SNPs were removed for missingness?

### Apply Hardy-Weinberg thresholds

1.  Generate a new binary file set called 'D1D_QC' from the D1D_imiss3pc_lmiss2pc_maf5pc file, applying a HWE threshold of 0.0001.

2.  This is our final, QC'ed file set.

3.  Examine log and output files.

-How many SNPs were removed for HWE *p-values* below the threshold?

📝 **It is useful to know how to do this, but be careful about setting this threshold - strong association signals can cause departure from HWE and you may remove great results! Use a lenient threshold and apply to controls only to avoid this problem. HWE can also be checked post-hoc for each SNP.**

## Perform GWAS

### Case/Control GWAS - no covariates

Run the following code, which performs a genetic association study using logistic regression on some case/control data:

        plink --bfile D1D_QC --logistic --adjust --pheno D1D.pheno1 --out D1D_CC

- What are the raw and Bonferroni-adjusted p-values for the top hit?
- What does this mean - is there a significant association?
- Are there any other significant associations?

### Case/Control GWAS - with covariates

Here we repeat the previous analysis but this time including some covariates. The file D1D.pcs1234 contains the first 4 principal components from a PCA on the genetic data.

Run the analysis specifying the covariates file:

  plink --bfile D1D_QC --logistic --adjust --pheno D1D.pheno1 --covar D1D.pcs.1234 --out D1D_CC_PCadj

- What are the raw and Bonferroni-adjusted p-values for the top hit?
- What does this mean - is there a significant association?
- Suggest a reason for the different results when adjusting for the 4 PCs?

# License

This work is licensed under Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International Public License and the below text is a summary of
the main terms of the full Legal Code (the full licence) available at https://creativecommons.org/licenses/by-nc-sa/4.0/legalcode

#### You are free to:

-   **Share** --- copy and redistribute the material in any medium or format

-   **Adapt** --- remix, transform, and build upon the material

The licensor cannot revoke these freedoms as long as you follow the license terms.

#### Under the following terms: {#under-the-following-terms .unnumbered}

-   **Attribution** --- You must give appropriate credit, providea link to the license, and indicate if changes were made. You may do so in any reasonable manner, but not in any way that suggests the licensor endorses you or your use.

-   **NonCommercial** --- You may not use the material for commercial purposes.

-   **ShareAlike** --- If you remix, transform, or build upon the material, you must distribute your contributions under the same license as the original.

No additional restrictions --- You may not apply legal terms or technological measures that legally restrict others from doing anything the license permits.

#### Notices:

You do not have to comply with the license for elements of the material in the public domain or where your use is permitted by an applicable exception or limitation. No warranties are given. The license may not give you all of the permissions necessary for your intended use. For example, other rights such as publicity, privacy, or moral rights may limit how you use the material.
