 # Polygenic Risk Score Analyses Workshop 2023


 # Day 1a: GWAS & relevant Statistics


## Introduction to Bash

 Most software in Bioinformatics and Statistical Genetics need to be
 run in a Unix environment (e.g. Linux or Mac OS) and most
 high-performance computer clusters run Unix systems. Therefore,
 although there are alternatives available on Windows (command line,
 Linux subsystems or Virtual Machines), it will be highly beneficial to
 become familiar with performing research in a Unix-only environment.


## Moving around the File System

 To begin our practical, please open up a \"terminal\" on your computer
 (on a Mac this is stored in Applications/Utilities/).

 We can change our directory using the following command:

       cd \<Path\

 where *\<Path\* is the path to the target directory.


 Some common usage of cd includes

 

       cd ~/ # will bring you to your home directory

       cd ../ # will bring you to the parent directory (up one level)

       cd XXX # will bring you to the XXX directory, so long as it is in the current directory

 As an example, we can move to the **PRS_Workshop** directory by
 typing:

       cd ~/Desktop/PRS_Workshop/

 [PRACTICAL 1. INTRODUCTION TO BASH]{.underline}

## Looking at the Current Directory

 Once we have moved into the **PRS_Workshop** folder, we can list out
 the folder content by typing:

        ls

 For ls, there are a number of additional Unix command options that you
 can append to it to get additional information, for example:

        ls -l \# shows files as list

        ls -lh \# shows files as a list with human readable format

        ls -lt \# shows the files as a list sorted by time-last-edited

        ls -lS \# shows the files as a list sorted by size



## Counting Number of Lines in File

 We can also count the number of lines in a file with the following
 command (where

 *\<file\* is the file of interest):

        wc -l \<file\

 1.4. SEARCH FILE CONTENT

 Often we would like to store the output of a command, which we can do
 by *redirecting* the output of the command to a file. For example, we
 can redirect the count of the **GIANT_Height.txt** to **giant_count**
 using the following command:

     wc -l GIANT_Height.txt \ giant_count.txt

**Search File Content**

 Another common task is to search for specific words or characters in a
 file (e.g. does this file contain our gene of interest?). This can be
 performed using the \"grep\" command as follows:

       grep \<string\ \<file\

 For example, to check if the [Single Nucleotide Polymorphism
 (SNP)](#_bookmark0) *rs10786427* is present in **GIANT_Height.txt**,
 we can do:

       grep rs10786427 GIANT_Height.txt

 In addition, grep allows us to check if patterns contained in one file
 can be found in another file. For example, if we want to extract a
 subset of samples from the phenotype file (e.g. extract the list of
 samples in **Data/Select.sample**), we can do:

        grep -f Select.sample TAR.height

 An extremely useful feature of the terminal is chaining multiple
 commands into one command, which we call ***piping*** .

 For example, we can use piping to count the number of samples in
 **Select.sample**

 that were found in **TAR.height** in a single command, as follows:

   grep -f Select.sample TAR.height | wc -l


## Filtering and Reshuﬄing Files

 A very powerful feature of the terminal is the **awk** programming
 language, which allows us to extract subsets of a data file, filter
 data according to some criteria or perform arithmetic operations on
 the data. awk manipulates a data file by per- forming operations on
 its **columns** - this is extremely useful for scientific data sets
 because typically the columns features or variables of interest.

 For example, we can use awk to produce a new results file that only
 contains SNP rsIDs (column 1), allele frequencies (column 4) and *P*
 -values (column 7) as follows:

      awk '{ print $1,$4,$7}' GIANT_Height.txt < GIANT_Height_3cols.txt

 We can also use a \"conditional statement\" in awk to extract all
 *significant [SNPs](#_bookmark0)*

 from the results file, using the following command:

     awk '{if($7 < 5e-8) { print } }' GIANT_Height.txt < Significant_SNPs.txt

# Or the short form:

     awk '$7 < 5e-8{ print}' GIANT_Height.txt < Significant_SNPs.txt

 "if($7<5e-8)" and "$7 < 5e-8" tell awk to extract any rows
 with column 7 (the column containing *P* -value) with a value of
 smaller than 5e-8 and {print} means that we would like to print the
 entire row when this criterion is met.


# Introduction to R

**R** is a useful programming language that allows us to perform a
variety of statis- tical tests and data manipulation. It can also be
used to generate fantastic data visualisations. Here we will go
through some of the basics of **R** so that you can better understand
the practicals throughout the workshop.


## Basics

 If you are not using R Studio then you can type **R** in your terminal
 to run **R** in the terminal.

 **Working Directory**

 When we start **R**, we will be working in a specific folder called
 the **working direc- tory**. We can check the current/working
 directory we are in by typing:

       getwd()

 And we can change our working directory to the **Practical** folder by

     setwd("~/Desktop/PRS_Workshop/Day_1a")

### Libraries

 Most functionality of **R** is organised in \"packages\" or
 \"libraries\". To access these functions, we will have to install and
 \"load\" these packages. Most commonly used packages are installed
 together with the standard installation process. You can install a new
 library using the install.packages function.

 For example, to install *ggplot2*, run the command:

        install.packages("ggplot2")


 After installation, you can load the library by typing

        library(ggplot2)

 Alternatively, we can import functions (e.g. that we have written)
 from an R script file on our computer. For example, you can load the
 Nagelkerke *R*^2^ function by typing

         source(\"Software/nagelkerke.R\")

 And you are now able to use the NagelkerkeR2 function (we will use
 this function at the end of this worksheet).

### Variables in R

 You can assign a value or values to any variable you want using \<-.
 e.g

 # Assign a number to a

      a <- 1

# Assign a vector containing a,b,c to b

      v1 <- c("a", "b","c")

### Functions

 You can perform lots of operations in **R** using diﬀerent built-in R
 functions. Some examples are below:

 ### Assign number of samples

        nsample <- 10000

 ### Generate nsample random normal variable with mean = 0 and sd = 1

        normal <- rnorm(nsample, mean=0,sd=1)

        normal.2 <- rnorm(nsample, mean=0,sd=1)

 ### We can examine the first few entries of the result using head

      head(normal)

 ### And we can obtain the mean and sd using

 9 mean(normal)

 10 sd(normal)

 11 \# We can also calculate the correlation between two variables
 using cor

 12 cor(normal, normal.2)

## Plotting

 While **R** contains many powerful plotting functions in its base
 packages, customisa- tion can be diﬃcult (e.g. changing the colour
 scales, arranging the axes). **ggplot2**

 [2.3. REGRESSION MODELS]{.underline}

 is a powerful visualization package that provides extensive
 flexibility and customi- sation of plots. As an example, we can do the
 following

 1 \# Load the package

 2 library(ggplot2)

 3 \# Specify sample size

 4 nsample\<-1000

 5 \# Generate random grouping using sample with replacement

 6 groups \<- sample(c(\"a\",\"b\"), nsample, replace=T)

 7 \# Now generate the data

 8 dat \<- data.frame(x=rnorm(nsample), y=rnorm(nsample), groups)

 9 \# Generate a scatter plot with diﬀerent coloring based on group

 10 ggplot(dat, aes(x=x,y=y,color=groups))+geom_point()

2.  []{#_bookmark15 .anchor}**Regression Models**

 In statistical modelling, regression analyses are a set of statistical
 techniques for estimating the relationships among variables or
 features. We can perform regression analysis in **R**.

 Use the following code to perform linear regression on simulated
 variables \"x\" and \"y\":

 1 \# Simulate data

 2 nsample \<- 10000

 3 x \<- rnorm(nsample)

 4 y \<- rnorm(nsample)

 5 \# Run linear regression

 6 lm(y\~x)

 7 \# We can store the result into a variable

 8 reg \<- lm(y\~x)

 9 \# And get a detailed output using summary

 10 summary(lm(y\~x))

 11 \# We can also extract the coeﬃcient of regression using

 12 reg\$coeﬃcient

 13 \# And we can obtain the residuals by

 14 residual \<- resid(reg)

 15 \# Examine the first few entries of residuals

 16 head(residual)

 17 \# We can also include covariates into the model

 18 covar \<- rnorm(nsample)

 19 lm(y\~x+covar)

 20 \# And can even perform interaction analysis

 21 lm(y\~x+covar+x∗covar)

 Alternatively, we can use the glm function to perform the regression:

 1 glm(y\~x)

 For binary traits (case controls studies), logistic regression can be
 performed using

 1 \# Simulate samples

 [PRACTICAL 2. INTRODUCTION TO R]{.underline}

 2 nsample\<- 10000

 3 x \<- rnorm(nsample)

 4 \# Simulate binary traits (must be coded with 0 and 1)

 5 y \<- sample(c(0,1), size=nsample, replace=T)

 6 \# Perform logistic regression

 7 glm(y\~x, family=binomial)

 8 \# Obtain the detail output

 9 summary(glm(y\~x, family=binomial))

 10 \# We will need the NagelkerkeR2 function

 11 \# to calculate the pseudo R2 for logistic model

 12 source(\"Software/nagelkerke.R\")

 13 reg \<- glm(y\~x, family=binomial)

 14 \# Calculate the Nagelkerke R2 using the NagelkerkeR2 function

 15 NagelkerkeR2(reg)
