> **Polygenic Risk Score Analyses Workshop 2022**
>
> **Day 2: Introduction to Polygenic Score Analyses**
>
> **University of Mauritius**
>
> **Contents**

1.  [Descriptive analyses of the 1000 Genomes dataset](#_bookmark0) . .
    . . . . . . . . ii

2.  [Cross-ancestry GWAS](#_bookmark1) . . . . . . . . . . . . . . . . .
    . . . . . . . . x

    1.  [Using PRS-CSx](#_bookmark2) . . . . . . . . . . . . . . . . . .
        . . . . . . . x

3.  1.  
    2.  

[Appendices](#appendices) xiv[Downloading and preparing 1000Genomes
data](#downloading-and-preparing-1000genomes-data) xiv[Investigation of
relatedness within the
sample](#investigation-of-relatedness-within-the-sample) xv

i

1.  []{#_bookmark0 .anchor}**Descriptive analyses of the 1000 Genomes
    dataset**

### Introduction to the multi-ancestry 1000Genomes dataset {#introduction-to-the-multi-ancestry-1000genomes-dataset .unnumbered}

> These data are the product of a whole-genome sequencing initiative
> completed in 2013. Individuals from 26 different source populations
> from around the world are included. For simplicity these divergent
> populations have been collapsed into 5 continental
> *super-populations*. The scripts used to download and process the
> 1000Genomes data for the purposes of this course can be found in the
> appendix. The cleaned Plink binary files should be located inside your
> home directory

### Sample sizes {#sample-sizes .unnumbered}

> First of all we would like to know the number of individuals within
> each super- population. Type the following command to query the number
> of European ancestry individuals in the downloaded dataset
>
> 1 grep -F \'EUR\' all_phase3.king.psam \| wc -l
>
> Repeat the same command for East Asian, African, South Asian and
> Amerindian
>
> *superpopulations*, by inserting the relevant ancestry codes (EAS,
> AFR, SAS, AMR).

![](media/image1.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Number of genetic variants {#number-of-genetic-variants .unnumbered}

> We do not need to use the full genome-wide data for this tutorial,
> only a small frac- tion of the 80 million total available variants.
> This provides a reliable approximation for the genomic analyses in
> this tutorial and importantly, reduces the computation time required
> to complete the tutorial. The following command derives the number of
> genetic variants on chromosomes 1 to chromsome 22 by counting the
> number of lines in the relevant (.bim) file, which contains a single
> variant per line.
>
> 1 wc chr1-22.bim -l
>
> **How many genetic variants (SNPs) are there in this dataset?** To
> quantify the number of single nucleotide polymorphisms (SNPs) we can
> ask plink to write a list of SNPs
>
> 1 ./plink \--bfile chr1-22 \--snps-only \--write-snplist
>
> See the output file plink.snplist, which contains a list of all the
> SNPs in the dataset
>
> ![](media/image1.jpeg){width="0.3541655730533683in"
> height="0.3593744531933508in"}

### Quantification of variable SNPs {#quantification-of-variable-snps .unnumbered}

> The rate at which a genetic variant occurs in a population is also
> known as its allelic frequency. Allele frequencies are shaped by
> evolutionary forces over a long period of time and hence can vary.
> This has implications for PRS research as the allelic frequency
> distribution of a disease or trait may vary between populations. It is
> possible to generate allele frequency statistics for each SNP in a
> given population, using the population information in the file
> pop_info.pheno.
>
> 1 ./plink \--bfile chr1-22 \--snps-only \--freq \--within
> pop\\\_info.pheno
>
> Population-stratified allele frequency results can be found in the
> output file **plink.frq.strat**. For each population, print the
> numbers of total SNPs to screen, as follows:
>
> 1 grep -F \'AFR\' plink.frq.strat \| wc -l
>
> Compare the totals against number of SNPs which have minor allele
> frequencies greater than 0 (and hence are useful for statistical
> analysis). Do this for all 5 populations (EAS, EUR, SAS, EUR and AFR),
> using the code given below:
>
> 1 grep -F \'AFR\' plink.frq.strat \> freq_report.afr
>
> 2 grep -F \'AMR\' plink.frq.strat \> freq_report.amr
>
> 3 grep -F \'EUR\' plink.frq.strat \> freq_report.eur
>
> 4 grep -F \'EAS\' plink.frq.strat \> freq_report.eas
>
> 5 grep -F \'SAS\' plink.frq.strat \> freq_report.sas
>
> 6 grep -F \'AFR\' plink.frq.strat \| awk \'\$6 \>0\' freq_report.afr
> \| wc -l
>
> 7 grep -F \'EUR\' plink.frq.strat \| awk \'\$6 \>0\' freq_report.eur
> \| wc -l
>
> 8 grep -F \'EAS\' plink.frq.strat \| awk \'\$6 \>0\' freq_report.eas
> \| wc -l
>
> 9 grep -F \'AMR\' plink.frq.strat \| awk \'\$6 \>0\' freq_report.amr
> \| wc -l
>
> 10 grep -F \'SAS\' plink.frq.strat \| awk \'\$6 \>0\' freq_report.sas
> \| wc -l
>
> Having compared the number of SNPs that show variation in each
> population:

![](media/image1.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Investigating missingness {#investigating-missingness .unnumbered}

> Genotype missingness, caused by genotyping failure can potentially
> lead to biased allele frequency estimation. Therefore missingness
> needs to be excluded as a possible source of bias when calculating
> allele frequency differences.
>
> 1 ./plink \--bfile chr1-22 \--missing \--within pop_info.pheno
>
> The output file **plink.lmiss** provides a variant-based missing data
> report). Use the following code to query the number of genotyping
> failures based on the missingness information in the **NMISS** column:
>
> 1 awk \'\$4 \> 0\' plink.lmiss \| wc -l

![](media/image1.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

### Cross-population allele frequency comparisons {#cross-population-allele-frequency-comparisons .unnumbered}

> Here we compare profiles of allele frequency across the five ancestral
> populations. To do this we will use the previously-generated output on
> minor allele frequencies per ancestry group (the file
> \"plink.frq.strat\").

### In R: {#in-r .unnumbered}

> 1 library(dplyr)
>
> 2 library(ggplot2)
>
> 3 freq \<-read.table(\"plink.frq.strat\", header =T)
>
> 4 plotDat \<- freq %\>%
>
> 5 mutate(AlleleFrequency = cut(MAF, seq(0, 1, 0.25))) %\>%
>
> 6 group_by(AlleleFrequency, CLST) %\>%
>
> 7 summarise(FractionOfSNPs = n()/nrow(freq) ∗ 100)
>
> 8
>
> 9 ggplot(na.omit(plotDat),
>
> 10 aes(AlleleFrequency, FractionOfSNPs, group = CLST, col = CLST)) +
>
> 11 geom_line() +
>
> 12 scale_y_continuous(limits = c(0, 12)) +
>
> 13 ggtitle(\"Distribution of allele frequency across genome\")

![](media/image1.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Calculation of Fst {#calculation-of-fst .unnumbered}

> Fst is a formal way to estimate divergence between populations, using
> information derived from a set of genome-wide and mutually independent
> SNPs. *Fst* is estimated efficiently in **Plink-2**, which calculates
> *Fst* for all pairwise ancestry combinations, (the same application in
> Plink-1.9 has more limited capabilities).
>
> Use the following command to calculate Fst:
>
> 1 ./plink2 \--bfile chr1-22 \--fst POP \--pheno pop_info.pheno
>
> Check the output file **plink2.fst.summary**
>
> ![](media/image1.jpeg){width="0.3541655730533683in"
> height="0.3593744531933508in"}

### Linkage disequilibrium versus genomic distance, across populations {#linkage-disequilibrium-versus-genomic-distance-across-populations .unnumbered}

> We will now perform pairwise LD comparisons between genome-wide snps
> in order to show cross-populations relationships between genomic
> distance and LD strength.
>
> -------- We derive information on pairwise R2 between all SNPs:
>
> 1 ./plink \--bfile chr1-22 \\
>
> 2 \--keep-cluster-names AFR \\
>
> 3 \--within pop_info.pheno \\
>
> 4 \--r2 \\
>
> 5 \--ld-window-r2 0 \\
>
> 6 \--ld-window 999999 \\
>
> 7 \--ld-window-kb 2500 \\
>
> 8 \--threads 30 \\
>
> 9 \--out chr1-22.AFR
>
> Repeat this step for all five 1000Genomes populations. Output files
> containing LD info for all pairwise SNPs, have a *'.ld'* suffix Next,
> create a summary file contain- ing the base-pair distance between each
> pair and the corresponding *r*2 value. The following example shows
> this for **AFR** and **EUR** populations only, as just these
> populations will be used in the plot.
>
> 1 cat chr1-22.AFR.ld \| sed 1,1d \| awk -F \" \" \'function abs(v)
> {return v \< 0 ? -v : \\
>
> v}BEGIN{OFS=\"\\t\"}{print abs(\$5-\$2),\$7}\' \| sort -k1,1n \> \\
>
> chr1-22.AFR.ld.summary
>
> 2 cat chr1-22.EUR.ld \| sed 1,1d \| awk -F \" \" \'function abs(v)
> {return v \< 0 ? -v : \\
>
> v}BEGIN{OFS=\"\\t\"}{print abs(\$5-\$2),\$7}\' \| sort -k1,1n \> \\
>
> chr1-22.EUR.ld.summary
>
> --------

### LD decay versus chromosomal distance {#ld-decay-versus-chromosomal-distance .unnumbered}

> In order to visualise LD behaviour as a function of chromosomal
> distance we need to add additional functionality to be able to carry
> out the necessary data transformation (*dplyr* ) and manipulation of
> character strings (*stringr* )
>
> 1 R
>
> 2 install.packages(\"dplyr\")
>
> 3 install.packages(\"stringr\")
>
> 4 install.packages(\"ggplot2\")
>
> 5 library(dplyr)
>
> 6 library(stringr)
>
> 7 library(ggplot2)
>
> The next piece of code does the following **4 steps**:

1.  Loads the previously-generated information on pairwise LD.

2.  Categorises distances into intervals of fixed length (100 Kb).

3.  Computes mean and median *r*2 within blocks.

4.  Obtains mid-points for each distance interval.

> 1
> dfr\<-read.delim(\"chr1-22.AFR.ld.summary\",sep=\"\",header=F,check.names=F,\\newline
> \\
>
> stringsAsFactors=F)\\newline
>
> 2 colnames(dfr)\<-c(\"dist\",\"rsq\")
>
> 3
> dfr\$distc\<-cut(dfr\$dist,breaks=seq(from=min(dfr\$dist)-1,to=max(dfr\$dist)+1,by=100000))
>
> 4 dfr1\<-dfr %\>% group_by(distc) %\>% \\
>
> summarise(mean=mean(rsq),median=median(rsq))
>
> 5 dfr1 \<- dfr1 %\>% \\
>
> mutate(start=as.integer(str_extract(str_replace_all(distc,\"\[\\\\(\\\\)\\\\\[\\\\\]\]\",\"\"),\"\^\[0-9-e+.\]+\")),
>
> 6 \\
>
> end=as.integer(str_extract(str_replace_all(distc,\"\[\\\\(\\\\)\\\\\[\\\\\]\]\",\"\"),\"\[0-9-e+.\]+\$\")),
>
> 7 mid=start+((end-start)/2))
>
> Steps 1-4 should be repeated for the file chr1-22.\_EUR.ld.summary.
>
> The output object *dfr1* on lines 4 and 5 should be renamed *dfr2* to
> prevent the object *df1* being over-written. Finally, we plot LD decay
> for **AFR** and **EUR** reference populations in a single graph:
>
> 1 ggplot()+
>
> 2
> geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour=\"grey20\")+
>
> 3
> geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour=\"grey40\")+
>
> 4 labs(x=\"Distance (Megabases)\",y=expression(LD\~(r\^{2})))+
>
> 5 \\
>
> scale_x_continuous(breaks=c(0,2∗10\^6,4∗10\^6,6∗10\^6,8∗10\^6),labels=c(\"0\",\"2\",\"4\",\"6\",\"8\"))+
>
> 6 theme_bw()

![](media/image1.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Distribution of LD-block length {#distribution-of-ld-block-length .unnumbered}

> The next set of scripts will allow us to visualise the distribution of
> LD block length across the different 1000Genomes populations.
>
> 1 ./plink \--bfile chr1-22 \--keep-cluster-names AFR \--blocks
> no-pheno-req \\
>
> no-small-max-span \--blocks-max-kb 250 \--within pop_info.pheno
> \--threads 30 \\
>
> \--out AFR
>
> The "--block\" flag estimates haplotype blocks using the same block
> definition imple- mented by the software Haploview. The default
> setting for the flag *\--blocks-max-kb* only considers pairs of
> variants that are within **200 kilobases** of each other. The output
> files from the above command is a ***.blocks*** file. Use the same
> code to gen- erate output for *EUR, EAS, SAS* and *AMR* populations
> (as it is not possible to generate population-specific information
> using the \--within flag).
>
> Then, in R:
>
> 1 R
>
> 2 dfr.afr \<- \\
>
> read.delim(\"AFR.blocks.det\",sep=\"\",header=T,check.names=F,stringsAsFactors=F)
>
> 3 colnames(dfr.afr) \<- tolower(colnames(dfr.afr))
>
> Load each of the 5 datasets and set column names to lower case. Now
> plot the data
>
> 1 plot (density(dfr.afr\$kb), main=\"LD block length distribution\",
> \\
>
> ylab=\"Density\",xlab=\"LD block length (Kb)\" )
>
> 2 lines (density(dfr.eur\$kb), col=\"blue\")
>
> 3 lines (density(dfr.eas\$kb), col=\"red\")
>
> 4 lines (density(dfr.amr\$kb), col=\"purple\")
>
> 5 lines (density(dfr.sas\$kb), col=\"green\")
>
> 6 legend(\"topright\",c(\"AFR\",\"EAS\",\"EUR\",\"SAS\",\"AMR\"), \\
>
> fill=c(\"black\",\"red\",\"blue\",\"green\",\"purple\"))

![](media/image1.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}

### Principle Component Analysis {#principle-component-analysis .unnumbered}

> In the case where you have genetic data for *N* individuals whose
> genetic ancestry status is unknown possible strategies that can narrow
> the bounds of this uncertainty include Principle Components Analysis
> **(PCA)**; a visual inquiry that can be used to determine the spatial
> distance between genetic samples in question and reference
> populations, of which 1000 genomes is one. Each sample is projected
> onto principle components that reflect different linear combinations
> of conditionally independent SNPs ). The steps involved are as
> follows.
>
> 1 ./plink \--bfile chr1-22 \--indep-pairwise 250 25 0.1 \--maf 0.1
> \--threads 30 \--out \\
>
> chr1-22.ldpruned.all.1kgv2
>
> 2 ./plink \--bfile chr1-22 \--extract
> chr1-22.ldpruned.all.1kgv2.prune.in \--pca \--threads 30

### In R: {#in-r-1 .unnumbered}

> 1 require(\'RColorBrewer\')
>
> 2 options(scipen=100, digits=3)
>
> Read in the eigenvectors, produced in PLINK

### In R: {#in-r-2 .unnumbered}

> 1 eigenvec \<- read.table(\'plink.eigenvec\', header = F, skip=0, sep
> = \' \')
>
> 2 rownames(eigenvec) \<- eigenvec\[,2\]
>
> 3 eigenvec \<- eigenvec\[,3:ncol(eigenvec)\]
>
> 4 colnames(eigenvec) \<- paste(\'Principal Component \', c(1:20), sep
> = \'\')
>
> First we read in the PED data
>
> 1 PED \<- read.table(\"all_phase3.king.psam\", header = TRUE, skip =
> 0, sep = \'\\t\')
>
> 2 PED \<- PED\[which(PED\$IID %in% rownames(eigenvec)), \]
>
> 3 PED \<- PED\[match(rownames(eigenvec), PED\$IID),\]
>
> Then we set up the colour scheme we are going to use
>
> 1 require(\'RColorBrewer\')
>
> 2 PED\$Population \<- factor(PED\$Population, levels=c(
>
> 3 \"ACB\",\"ASW\",\"ESN\",\"GWD\",\"LWK\",\"MSL\",\"YRI\",
>
> 4 \"CLM\",\"MXL\",\"PEL\",\"PUR\",
>
> 5 \"CDX\",\"CHB\",\"CHS\",\"JPT\",\"KHV\",
>
> 6 \"CEU\",\"FIN\",\"GBR\",\"IBS\",\"TSI\",
>
> 7 \"BEB\",\"GIH\",\"ITU\",\"PJL\",\"STU\"))
>
> 8
>
> 9 col \<- colorRampPalette(c(
>
> 10
> \"yellow\",\"yellow\",\"yellow\",\"yellow\",\"yellow\",\"yellow\",\"yellow\",
>
> 11 \"forestgreen\",\"forestgreen\",\"forestgreen\",\"forestgreen\",
>
> 12 \"grey\",\"grey\",\"grey\",\"grey\",\"grey\",
>
> 13
> \"royalblue\",\"royalblue\",\"royalblue\",\"royalblue\",\"royalblue\",
>
> 14 \\
>
> \"black\",\"black\",\"black\",\"black\",\"black\"))(length(unique(PED\$Population)))\[factor(PED\$Population)\]
>
> Finally we generate our PCA plots
>
> 1 project.pca \<- eigenvec
>
> 2 par(mar = c(5,5,5,5), cex = 2.0,
>
> 3 cex.main = 7, cex.axis = 2.75, cex.lab = 2.75, mfrow = c(1,2))
>
> 4
>
> 5 plot(project.pca\[,1\], project.pca\[,2\],
>
> 6 type = \'n\',
>
> 7 main = \'A\',
>
> 8 adj = 0.5,
>
> 9 xlab = \'First component\',
>
> 10 ylab = \'Second component\',
>
> 11 font = 2,
>
> 12 font.lab = 2)
>
> 13 points(project.pca\[,1\], project.pca\[,2\], col = col, pch = 20,
> cex = 2.25)
>
> 14 legend(\'bottomright\',
>
> 15 bty = \'n\',
>
> 16 cex = 3.0,
>
> 17 title = \'\',
>
> 18 c(\'AFR\', \'AMR\', \'EAS\',
>
> 19 \'EUR\', \'SAS\'),
>
> 20 fill = c(\'yellow\', \'forestgreen\', \'grey\', \'royalblue\',
> \'black\'))
>
> 21
>
> 22 plot(project.pca\[,1\], project.pca\[,3\],
>
> 23 type=\"n\",
>
> 24 main=\"B\",
>
> 25 adj=0.5,
>
> 26 xlab=\"First component\",
>
> 27 ylab=\"Third component\",
>
> 28 font=2,
>
> 29 font.lab=2)
>
> 30 points(project.pca\[,1\], project.pca\[,3\], col=col, pch=20,
> cex=2.25)

![](media/image1.jpeg){width="0.3541655730533683in"
height="0.3593744531933508in"}

2.  []{#_bookmark1 .anchor}**Cross-ancestry GWAS**

    1.  []{#_bookmark2 .anchor}**Using PRS-CSx**

> PRS-CSx is a new software that allows integration of GWAS summary
> statistics from multiple populations to improve cross-population
> polygenic prediction. It does this by modelling how GWAS effects at
> each SNP are linked across different ancestry groups. It can be
> applied in a variety of situations incorporating different genetic
> architectures and levels of cross-population genetic overlap. In this
> practical we will be familiarising ourselves with running the software
> and interpreting the output.
>
> Simulated dataset are well-suited for these purposes and also allow us
> to vary dif- ferent parameters of the underlying genetic architecture,
> such as the heritability of the phenotype used for analysis. Hence we
> will use simulated data to explore the PRScsx software, using PRSice
> as a benchmark for the conventional clumping and thresholding method.
> As PRSice is designed for modelling PRS risk in the con- text of a
> single homogenous population group, we expect there to be a shortfall
> in performance, compared to PRS-CSx.

### Getting started {#getting-started .unnumbered}

> In this practical session we will be using PRS-CSx and PRSice
> interchangebly. To start please navigate to the PRSice directory and
> check for the following objects, consisting of:
>
> 1 Simulated summary statistics based on a cohort of \~250,000 European
> ancestry \\
>
> individuals and \~4,800 African ancestry individuals for BMI-derived
> traits.
>
> 2 AFR-BMI-ukb-simulated.sumstats.final.prscsx
>
> 3 EUR-BMI-ukb-simulated.sumstats.final.prscsx
>
> A set of Plink files for each target ancestry group
>
> 1 EUR_1kg-updated2.bim, EUR_1kg-updated2.bed, EUR_1kg-updated2.fam
>
> 2 AFR_1kg-updated2. bim, EUR_1kg-updated2.bed, EUR_1kg-updated2 .fam
>
> A list of pre-clumped SNPs for AFR and EUR
>
> 1 AFR.valid.snps \\newline
>
> 2 EUR.valid.snps
>
> A file containing the phenotypes of interest. Multiple variants of the
> main GWAS trait are included. The additional phenotypes will be used
> to demonstrate how heritability influences the predictability of a
> trait.
>
> 1 bmi_afr_1kg.sim_pheno.final
>
> 2 bmi_eur_1kg.sim_pheno.final

### Creating benchmarks for comparison using PRSice {#creating-benchmarks-for-comparison-using-prsice .unnumbered}

> First we will set benchmarks for the PRS-CSx method by using PRSice to
> apply clumping and threshold-derived PRS scores to assess within and
> cross-ancestry per- formance. The following script performs standard
> clumping and thresholding on single ancestry base and target samples.
> **Within-AFR BMI**
>
> 1 Rscript ./PRSice.R \\
>
> 2 \--prsice ./PRSice_linux \\
>
> 3 \--extract PRSice_input/AFR.valid.snps \\
>
> 4 \--base ./PRSice_input/AFR-BMI-ukb-simulated.sumstats.final \\
>
> 5 \--A1 A1 \\
>
> 6 \--pvalue P \\
>
> 7 \--stat BETA \\
>
> 8 \--snp ID \\
>
> 9 \--no-clump \\
>
> 10 \--target ./PRSice_input/AFR_1kg-updated2 \\
>
> 11 \--binary-target F,F,F,F,F \\
>
> 12 \--pheno ./PRSice_input/bmi_afr_1kg.sim_pheno.final \\
>
> 13 \--pheno-col pheno10,pheno20,pheno33,pheno50,pheno100 \\
>
> 14 \--thread 30 \\
>
> 15 \--score-sum \\
>
> 16 \--out ./PRSice_output_practical_day4/afr
>
> \[language=bash\]

### Within-EUR BMI {#within-eur-bmi .unnumbered}

> 1 Rscript ./PRSice.R \\
>
> 2 \--prsice ./PRSice_linux \\
>
> 3 \--extract PRSice_input/EUR.valid.snps \\
>
> 4 \--base ./PRSice_input/EUR-BMI-ukb-simulated.sumstats.final \\
>
> 5 \--A1 A1 \\
>
> 6 \--pvalue P \\
>
> 7 \--stat BETA \\
>
> 8 \--snp ID \\
>
> 9 \--no-clump \\
>
> 10 \--target ./PRSice_input/EUR_1kg-updated2 \\
>
> 11 \--binary-target F,F,F,F,F \\
>
> 12 \--pheno ./PRSice_input/bmi_eur_1kg.sim_pheno.final \\
>
> 13 \--pheno-col pheno10,pheno20,pheno33,pheno50,pheno100 \\
>
> 14 \--thread 30 \\
>
> 15 \--score-sum \\
>
> 16 \--out ./PRSice_output_practical_day4/eur
>
> Next we will set benchmarks for cross-ancestry prediction of the
> BMI-derived
>
> trait textbfAFR into EUR
>
> 1 Rscript ./PRSice.R \\
>
> 2 \--prsice ./PRSice_linux \\
>
> 3 \--base PRSice_input/AFR-BMI-ukb-simulated.sumstats.final \\
>
> 4 \--no-clump \\
>
> 5 \--extract PRSice_input/AFR.valid.snps \\
>
> 6 \--A1 A1 \\
>
> 7 \--pvalue P \\
>
> 8 \--stat BETA \\
>
> 9 \--snp ID \\
>
> 10 \--target PRSice_input/EUR_1kg-updated2 \\
>
> 11 \--binary-target F,F,F,F,F \\
>
> 12 \--pheno PRSice_input/bmi_eur_1kg.sim_pheno.final \\
>
> 13 \--pheno-col pheno10,pheno20,pheno33,pheno50,pheno100 \\
>
> 14 \--thread 30 \\
>
> 15 \--score-sum \\
>
> 16 \--out ./PRSice_output_practical_day4/afr-into-eur
>
> textbfEUR into AFR
>
> 1 Rscript ./PRSice.R \\
>
> 2 \--prsice ./PRSice_linux \\
>
> 3 \--base PRSice_input/EUR-BMI-ukb-simulated.sumstats.final \\
>
> 4 \--no-clump \\
>
> 5 \--extract PRSice_input/EUR.valid.snps \\
>
> 6 \--A1 A1 \\
>
> 7 \--pvalue P \\
>
> 8 \--stat BETA \\
>
> 9 \--snp ID \\
>
> 10 \--target PRSice_input/AFR_1kg-updated2 \\
>
> 11 \--binary-target F,F,F,F,F \\
>
> 12 \--pheno PRSice_input/bmi_afr_1kg.sim_pheno.final \\
>
> 13 \--pheno-col pheno10,pheno20,pheno33,pheno50,pheno100 \\
>
> 14 \--thread 30 \\
>
> 15 \--score-sum \\
>
> 16 \--out ./PRSice_output_practical_day4/eur-into-afr
>
> Navigate to the output directory and then save a shortened summary of
> each set of results, by substituting-in the relevant summary filename
> and output file to the following code:
>
> 1 awk \'{print \$1,\$4,\$8,\$10}\' summary \> output_file
>
> 2 \\begin{end}
>
> 3 \\begin{questions}
>
> 4 What pattern do you notice in the way that the performance metric is
> varying across \\
>
> the different analyses?
>
> 5 \\end{questions}
>
> 6
>
> 7 \\noindent\\textbf{Running PRS-CSx} \\newline
>
> 8 We now proceed to analyses in which you\'ll be using PRS-CSx (you
> will need to switch \\
>
> over the directory of the same name). When you navigate to the PRS-CSx
> \\
>
> directory you should see the following set of 5 PRScsx run files:
>
> 9 \\begin{lstlisting}\[language=bash\]
>
> 10 PRScsx.py
>
> 11 prscx_run.sh
>
> 12 gigrnd.py
>
> 13 mcmc_gtb.py
>
> 14 parse_genet.py
>
> In addition to the above the following summary files contain the same
> informa-
>
> tion as those used in Plink but in a much more condensed form (as per
> PRS-CSx requirements):
>
> 1 EUR-BMI-ukb-simulated.sumstats.final.prscsx
>
> 2 AFR-BMI-ukb-simulated.sumstats.final.prscsx
>
> The file textitpop_info.pheno contains population information for the
> 1000Genomes individuals belonging tos this dataset. An additional
> folder containing a list of HAPMAP SNPs, (snpinfo_mult_1kg_hm3) at MAF
> threshold of 0.01 can be found in the same location. The same folder
> also contains pre-calculated LD profiles for HAPMAP SNPs in the 1000
> genomes sample. For each input GWAS used, PRS- CSx will write the
> posterior SNP effect size estimates for each chromosome to the
> user-specified directory. The output file contains chromosome, rsID,
> base position, A1, A2 and posterior effect size estimate for each SNP.
>
> Type the following code to start modelling coupling between effect
> sizes at indi- vidual SNPs across ancestries. PRS-CSx implements a
> Bayesian MCMC solution based on maximum likelihood to estimate values
> of the global shrinkage parameter, textitphi. For datasets in which
> the target population contains a mix of ancestries, or individuals of
> mixed ancestral heritage, the value of the shrinkage parameter can be
> learned using autonomous methods.
>
> For AFR target population:
>
> 1 python ./PRScsx.py \--ref_dir=1KG_reference_files \\
>
> 2 \--bim_prefix=AFR_1kg_wg_updated.ncbi.csx \\
>
> 3 \\
>
> \--sst_file=EUR-BMI-ukb-simulated.sumstats.final.prscsx,AFR-BMI-ukb-simulated.sumstats.final.prscsx
> \\
>
> \\
>
> 4 \--n_gwas=257000,4855 \\
>
> 5 \--phi=1e-4 \\
>
> 6 \--pop=EUR,AFR \\
>
> 7 \--out_dir=output2 \\
>
> 8 \--out_name=AFR_1kg.ncbi.csx
>
> 9 For EUR target population:
>
> 10 python ./PRScsx.py \--ref_dir=1KG_reference_files \\
>
> 11 \--bim_prefix=EUR_1kg_wg_updated.ncbi.csx \\
>
> 12 \\
>
> \--sst_file=EUR-BMI-ukb-simulated.sumstats.final.prscsx,AFR-BMI-ukb-simulated.sumstats.final.prscsx
> \\

+--------+-------+----------------------------------------------------+
|        | \\    |                                                    |
+========+=======+====================================================+
| > 13   |       | > \--n_gwas=257000,4855 \\                         |
+--------+-------+----------------------------------------------------+
| > 14   |       | > \--phi=1e-4 \\                                   |
+--------+-------+----------------------------------------------------+
| > 15   |       | > \--pop=EUR,AFR \\                                |
+--------+-------+----------------------------------------------------+
| > 16   |       | > \--out_dir= out-afr-as-target \\                 |
+--------+-------+----------------------------------------------------+
| > 17   |       | > \--out_name=EUR_1kg.ncbi.csx                     |
+--------+-------+----------------------------------------------------+

> **Notes** The snpset we are using have been filtered down to the same
> subset used for ancestry-specific pre-calculation of LD, on which the
> method relies. The PRS-CSx authors recommend setting phi values of
> 1e-2 (for highly polygenic traits, or GWAS with limited sample sizes )
> or 1e-4 (for less polygenic traits). The ideal would be a small-scale
> grid search across different values of phi (for example, phi=1e-6,
> 1e-4 and 1e-2, 1) to find the optimal phi value in subsequent
> validation datasets, as this can improve predictive performance
> overall. The goal here would be to implement feature selection via
> cross-validation to determine the optimal linear combination of
> ancestry specific score weights. A next step is to then validate the
> utility of the chosen model to improve the precision of cross-ancestry
> trait prediction in a specific research scenario.

# Appendices

## Downloading and preparing 1000Genomes data

> The instructions in this first appendix can be implemented as a script
> by typing vi download.sh and copying and pasting the code that
> follows. When finished type esc, full-colon, then w, q and return, (in
> sequence) to save. Once this is done enter chmod 777 download.sh to
> make the file executable. Then type ./download.sh to run the commands
> automatically.
>
> Download and install Plink v1.9 and Plink v2.0 to your home directory
>
> 1 wget
> https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip
>
> 2 wget
> https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20220519.zip
>
> 3 unzip plink_linux_x86_64_20220402.zip
>
> 4 unzip plink2_linux_avx2_20220519.zip
>
> Download and decompress 1000 Genomes phase 3 data
>
> 1 refdir=\'{Path/to/home/directory}\' #Enter path to your home
> directory path here
>
> 2
>
> 3 cd \$refdir
>
> 4
> pgen=h[ttps://www.dropb](http://www.dropbox.com/s/e5n8yr4n7y91fyp/all_hg38.pgen.zst?dl=1)o[x.com/s](http://www.dropbox.com/s/e5n8yr4n7y91fyp/all_hg38.pgen.zst?dl=1)/e5n8y[r4n7y91fyp/all_hg38.pgen.zst?dl=1](http://www.dropbox.com/s/e5n8yr4n7y91fyp/all_hg38.pgen.zst?dl=1)
>
> 5
> pvar=h[ttps://www.dropb](http://www.dropbox.com/s/cy46f1c8yutd1h4/all_hg38.pvar.zst?dl=1)o[x.com/s/cy46f1c8yutd1h4/all_hg38.p](http://www.dropbox.com/s/cy46f1c8yutd1h4/all_hg38.pvar.zst?dl=1)v[ar.zst?dl=1](http://www.dropbox.com/s/cy46f1c8yutd1h4/all_hg38.pvar.zst?dl=1)
>
> 6
> sample=h[ttps://www](http://www.dropbox.com/s/3j9zg103fi8cjfs/hg38_corrected.psam?dl=1).dropbo[x.com/s/3j9zg103fi8cjfs/hg38_corrected.psam?dl=1](http://www.dropbox.com/s/3j9zg103fi8cjfs/hg38_corrected.psam?dl=1)
>
> 7
>
> 8 wget \$pgen
>
> 9 mv \'all_hg38.pgen.zst?dl=1\' all_phase3.pgen.zst
>
> 10 ./plink2 \--zst-decompress all_phase3.pgen.zst \> all_phase3.pgen
>
> 11
>
> 12 wget \$pvar
>
> 13 mv \'all_hg38.pvar.zst?dl=1\' all_phase3.pvar.zst
>
> 14
>
> 15 wget \$sample
>
> 16 mv \'hg38_corrected.psam?dl=1\' all_phase3.psam
>
> Download list of 1st and 2nd degree relatiives from above plink page
> (based on KING output)
>
> 1 wget
> h[ttps://www.dropb](http://www.dropbox.com/s/129gx0gl2v7ndg6/deg2_hg38.king.cutoff.out.id?dl=1)o[x.com/s/129gx0gl2v7ndg6/](http://www.dropbox.com/s/129gx0gl2v7ndg6/deg2_hg38.king.cutoff.out.id?dl=1)deg2_h[g38.king.cutoff.out.id?dl=1](http://www.dropbox.com/s/129gx0gl2v7ndg6/deg2_hg38.king.cutoff.out.id?dl=1)
>
> 2 mv deg2_hg38.king.cutoff.out.id?dl=1 deg2_hg38.king.cutoff.out
>
> 3
>
> 4 ./plink2 \\
>
> 5 \--pfile \$refdir/all_phase3 vzs \\
>
> 6 \--allow-extra-chr \\
>
> 7 \--remove deg2_hg38.king.cutoff.out \\
>
> 8 \--make-pgen \\
>
> 9 \--out \$refdir/all_phase3.king
>
> Convert 1000 Genomes phase 3 data to plink 1 binary format
>
> 1 ./plink2 \\
>
> 2 \--pfile \$refdir/all_phase3.king \\
>
> 3 \--allow-extra-chr \\
>
> 4 \--max-alleles 2 \\
>
> 5 \--make-bed \\
>
> 6 \--out \$refdir/all_phase3
>
> 7 mv \$refdir/all_phase3.log \$refdir/log
>
> The following selects the correct versions of chromosomes 1 -- 22.
>
> 1 for i in {1..22}; do
>
> 2 plink \--bfile all_phase3 \--chr \${i} \--make-bed
> \--allow-extra-chr \--threads 50 \--out \\
>
> chr\${i}

+----+----------+------------------------------------------------------+
| 3  | > done   |                                                      |
+====+==========+======================================================+
| 4  |          |                                                      |
+----+----------+------------------------------------------------------+
| 5  | > for i  | > {1..22}; do                                        |
|    | > in     |                                                      |
+----+----------+------------------------------------------------------+
| 6  |          | > echo chr\${i} \>\> mergelist.txt;                  |
+----+----------+------------------------------------------------------+
| 7  | > done   |                                                      |
+----+----------+------------------------------------------------------+

> 8
>
> 9 plink \--merge-list mergelist.txt \--make-bed \--threads 50 \--out
> chr1-22 #This step takes \\
>
> a while
>
> Cleanup
>
> 1 for i in {1..22}; do
>
> 2 rm chr\${i}.bed chr\${i}.fam chr\${i}.bim chr\${i}.log
>
> 3 done

## Investigation of relatedness within the sample

> PiHat is a metric used to estimate relatedness due to IBD. LD pruning
> is the first prerequisite first step
>
> 1 ./plink \--bfile chr1-22 \--keep-cluster-names EUR \--within
> pop_info.pheno \\
>
> \--indep-pairwise 250 25 0.1 \--maf 0.1 \--threads 30 \--out
> chr1-22_pruned_EUR
>
> *\--indep-pairwise* has 3 parameters: (1) a window size in variant
> count or kilobase (if the *kb* modifier is present) units, (2)a
> variant count to shift the window at the end of each step (3) a
> pairwise *r*2 threshold. Briefly, one of each pair of variants within
> an active window, with *r*2 values in excess of specific cut-offs, are
> noted and variants are pruned until no pairs remain. Run the
> Identity-By-Descent (IBD) analysis for one ancestry group, before
> repeating for the remaining ancestries:
>
> 1 ./plink \--bfile chr1-22 \--keep-cluster-names EUR \--within
> pop_info.pheno \--extract \\
>
> chr1-22_pruned_EUR.prune.in \--genome \--threads 30 \--out
> chr1-22_pruned_EUR

### In R: {#in-r-3 .unnumbered}

> 1 R
>
> 2 library(ggplot2)
>
> 3 library(reshape2)
>
> 4 pihat_eas\<-read.table(\"chr1-22_pruned_EAS.genome\", header =T)
>
> 5 pihat_afr\<-read.table(\"chr1-22_pruned_AFR.genome\", header =T)
>
> Load all 5 population file as shown below, taking care to ensure that
> the name of each imported object is informative.
>
> 1 pihat_eas\<-read.table(\"chr1-22_pruned_EAS.genome\", header =T)
>
> 2 pihat_afr\<-read.table(\"chr1-22_pruned_AFR.genome\", header =T)
>
> Extract data for each population
>
> 1 eas=pihat_eas\$PI_HAT
>
> 2 afr=pihat_afr\$PI_HAT
>
> 3 n\<-max(length(eas),length(afr))
>
> 4 length(eur)\<-n
>
> 5 length(eas)\<-n
>
> Repeat the last step for the 3 other populations. After this the pihat
> info for all five populations can be merged together and formatted.
>
> 1 merge\<-data.frame(cbind(eur,afr,eas,sas,amr))
>
> 2 data\<-melt(merge)
>
> Plotting the data:
>
> 1 ggplot(data,aes(x=value, fill=variable)) +
>
> 2 geom_density(alpha=1.00)+xlim(-0.01,0.04) +
>
> 3 ggtitle(\"Pi Hat vs density\") +
>
> 4 xlab(\"Pi Hat\") +
>
> 5 ylab(\"Density\") +
>
> 6 labs(caption = \"(Pi Hat = 0.04 is less than a 4th degree
> relative)\")

![](media/image1.jpeg){width="0.3541666666666667in"
height="0.3593744531933508in"}
