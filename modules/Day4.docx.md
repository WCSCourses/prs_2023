## Introduction to the multi-ancestry 1000Genomes dataset 
These data are the product of a whole-genome sequencing initiative completed in 2013. Individuals from 26 different source populations from around the world are included. For simplicity these divergent populations have been collapsed into 5 continental super-populations. The scripts used to download and process the 1000Genomes data for the purposes of this course can be found in the appendix. The cleaned Plink binary files should be located inside your home directory

#### Sample sizes  
First of all we would like to know the number of individuals within each super-population. Type the following command to query the number of European ancestry individuals in the downloaded dataset.
```sh
grep -F 'EUR' all_phase3.king.psam | wc -l
```
Repeat the same command for East Asian, African, South Asian and Amerindian
superpopulations, by inserting the relevant ancestry codes (EAS, AFR, SAS, AMR).
#### Number of genetic variants
We do not need to use the full genome-wide data for this tutorial, only a small frac- tion of the 80 million total available variants. This provides a reliable approximation for the genomic analyses in this tutorial and importantly, reduces the computation time required to complete the tutorial. The following command derives the number of genetic variants on chromosomes 1 to chromsome 22 by counting the number of lines in the relevant (.bim) file, which contains a single variant per line.
```sh
wc chr1-22.bim -l
```
**Question: _How many genetic variants (SNPs) are there in this dataset?_**
To quantify the number of single nucleotide polymorphisms (SNPs) we can ask plink to write a list of SNPs:
```sh
./plink --bfile chr1-22 --snps-only --write-snplist
```
See the output file plink.snplist, which contains a list of all the SNPs in the dataset
#### Quantification of variable SNPs
The rate at which a genetic variant occurs in a population is also known as its allelic frequency. Allele frequencies are shaped by evolutionary forces over a long period of time and hence can vary. This has implications for PRS research as the allelic frequency distribution of a disease or trait may vary between populations. It is possible to generate allele frequency statistics for each SNP in a given population, using the population information in the file pop_info.pheno.
```sh
./plink --bfile chr1-22 --snps-only --freq --within pop_info.pheno
```
Population-stratified allele frequency results can be found in the output file plink.frq.strat. For each population, print the numbers of total SNPs to screen, as follows:
```sh
grep -F 'AFR' plink.frq.strat | wc -l
```
Compare the totals against number of SNPs which have minor allele frequencies greater than 0 (and hence are useful for statistical analysis). Do this for all 5 populations (EAS, EUR, SAS, EUR and AFR), using the code given below:
```sh
grep -F 'AFR' plink.frq.strat > freq_report.afr
grep -F 'AMR' plink.frq.strat > freq_report.amr
grep -F 'EUR' plink.frq.strat > freq_report.eur
grep -F 'EAS' plink.frq.strat > freq_report.eas
grep -F 'SAS' plink.frq.strat > freq_report.sas
grep -F 'AFR' plink.frq.strat | awk '$6 >0' freq_report.afr | wc -l
grep -F 'EUR' plink.frq.strat | awk '$6 >0' freq_report.eur | wc -l
grep -F 'EAS' plink.frq.strat | awk '$6 >0' freq_report.eas | wc -l
grep -F 'AMR' plink.frq.strat | awk '$6 >0' freq_report.amr | wc -l
grep -F 'SAS' plink.frq.strat | awk '$6 >0' freq_report.sas | wc -l
```
**_Question: Having compared the number of variable SNPs sites across populations:_**
* Which populations have the largest number (density) of SNP sites that can actually be considered polymorphic?
* What do you think is significant about the observed population ordering?
#### Interrogation of missingness
Genotype missingness, caused by genotyping failure can potentially lead to biased allele frequency estimation. Therefore missingness needs to be excluded as a possible source of bias when calculating allele frequency differences.
```sh
./plink --bfile chr1-22 --missing --within pop_info.pheno
```
The output file plink.lmiss provides a variant-based missing data report). Use the following code to query the number of genotyping failures based on the missingness information in the NMISS column:
```sh
awk '$4 > 0' plink.lmiss | wc -l
```
#### Cross-population allele frequency comparisons 
Here we compare profiles of allele frequency across the five ancestral populations. To do this we will use the previously-generated output on minor allele frequencies per ancestry group (the file "**_plink.frq.strat_**"). 
_In R_:
```sh
library(dplyr)
library(ggplot2)
freq <-read.table("plink.frq.strat", header =T)
plotDat <- freq %>%
mutate(AlleleFrequency = cut(MAF, seq(0, 1, 0.25))) %>%
group_by(AlleleFrequency, CLST) %>%
summarise(FractionOfSNPs = n()/nrow(freq) ∗ 100)
ggplot(na.omit(plotDat),
aes(AlleleFrequency, FractionOfSNPs, group = CLST, col = CLST)) +
geom_line() +
scale_y_continuous(limits = c(0, 12)) +
ggtitle("Distribution of allele frequency across genome")
```
**Question: _How are the allele frequencies in AFR distinguishable from the other global reference groups?_**
#### Calculation of Fst 
Fst is a formal way to estimate divergence between populations, using information derived from a set of genome-wide and mutually independent SNPs. Fst is estimated efficiently in Plink-2, which calculates Fst for all pairwise ancestry combinations, (the same application in Plink-1.9 has more limited capabilities).

Use the following command to calculate Fst:
```sh
./plink2 --bfile chr1-22 --fst POP --pheno pop_info.pheno
```
Check the output file **_plink2.fst.summary_**

**Questions**
* For which populations pairs is Fst greatest?
* For which populations is Fst smallest?
* What factors have contributed to the genetic divergence between human populations we see today?
  
#### LD structural differences across populations
Linkage disequilibrium (or LD) is used to define regions of the genome at which the correlation in allelic variability, spanning two or more SNPs, exceeds a threshold value. It is an important consideration in GWAS given that it provides an im- portant means of capturing genetic variation linked to disease remotely, (i.e. even when the causal variants in question are not directly typed by the genotyping array used). As well as being one of the critical factors that determines whether individual causal variants may be identified in a given population at GWAS, structural differences in LD between populations are also relevant in determining whether collective risk variants (encapsulated by polygenic risk scores) are portable across ancestrally-diverse populations. We will investigate cross-population LD profiles using the 1000Genomes dataset.

#### Linkage disequilibrium versus genomic distance
We will now perform pairwise LD comparisons between genome-wide snps in order to show cross-populations relationships between genomic distance and LD strength.We derive information on pairwise _R^2_ between all SNPs:
```sh
./plink --bfile chr1-22 \
--keep-cluster-names AFR \
--within pop_info.pheno \
--r2 \
--ld-window-r2 0 \
--ld-window 999999 \
--ld-window-kb 2500 \
--threads 30 \
--out chr1-22.AFR
```
Repeat this step for all five 1000Genomes populations. Output files containing LD info for all pairwise SNPs, have a '_.ld_' suffix Next, create a summary file contain- ing the base-pair distance between each pair and the corresponding r2 value. The following example shows this for AFR and EUR populations only, as just these populations will be used in the plot.
```sh
cat chr1-22.AFR.ld |sed 1,1d|awk -F " " 'function abs(v){return v < 0 ? -v:v}BEGIN{OFS="\t"}{print abs($5-$2),$7}' |sort -k1,1n > chr1-22.AFR.ld.summary
cat chr1-22.EUR.ld |sed 1,1d|awk -F " " 'function abs(v){return v < 0 ? -v:v}BEGIN{OFS="\t"}{print abs($5-$2),$7}'|sort -k1,1n > chr1-22.EUR.ld.summary
```

#### LD decay versus chromosomal distance
In order to visualise LD behaviour as a function of chromosomal distance we need to add additional functionality to be able to carry out the necessary data transformation (_dplyr_) and manipulation of character strings (_stringr_).
```sh
R
install.packages("dplyr")
install.packages("stringr")
install.packages("ggplot2")
library(dplyr)
library(stringr)
library(ggplot2)
```
The next piece of code does the following 4 steps:
* Loads the previously-generated information on pairwise LD.
* Categorises distances into intervals of fixed length (100 Kb).
* Computes mean and median r2 within blocks.
* Obtains mid-points for each distance interval.
```sh
dfr<-read.delim("chr1-22.AFR.ld.summary",sep="",header=F,check.names=F, stringsAsFactors=F)
colnames(dfr)<-c("dist","rsq")
dfr$distc<-cut(dfr$dist,breaks=seq(from=min(dfr$dist)-1,to=max(dfr$dist)+1,by=100000))
dfr1<-dfr %>% group_by(distc) %>% summarise(mean=mean(rsq),median=median(rsq))
dfr1 <- dfr1 %>% mutate(start=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"^[0-9-e+.]+")),
                        end=as.integer(str_extract(str_replace_all(distc,"[\\(\\)\\[\\]]",""),"[0-9-e+.]+$")),
                        mid=start+((end-start)/2))
```
Steps 1-4 should be repeated for the file _chr1-22._EUR.ld.summary_. The output object "_dfr1_" on lines 4 and 5 should be renamed "_dfr2_" to prevent the object _df1_ being over-written. 

Finally, we plot LD decay for AFR and EUR reference populations in a single graph:
```sh
ggplot()+
  geom_point(data=dfr1,aes(x=start,y=mean),size=0.4,colour="grey20")+
  geom_line(data=dfr1,aes(x=start,y=mean),size=0.3,alpha=0.5,colour="grey40")+
  labs(x="Distance (Megabases)",y=expression(LD~(r^{2})))+
  scale_x_continuous(breaks=c(0,2*10^6,4*10^6,6*10^6,8*10^6),labels=c("0","2","4","6","8"))+
  theme_bw()
```
**Questions**
(i) _What differences do you observe in terms of LD decay between AFR and EUR genomes?_
(ii) _How does this impact the portability of PRS performance between these two populations?_

#### Distribution of LD-block length
The following scripts allow us to visualise the distribution of LD block length across 1000Genomes populations.
```sh
./plink --bfile chr1-22 --keep-cluster-names AFR --blocks no-pheno-req no-small-max-span --blocks-max-kb 250 --within pop_info.pheno  --threads 30 --out AFR
./plink --bfile chr1-22 --keep-cluster-names EUR --blocks no-pheno-req no-small-max-span --blocks-max-kb 250 --within pop_info.pheno  --threads 30 --out EUR
```
The “_–block_" flag estimates haplotype blocks using the same block definition implemented by the _Haploview_ software. The default setting for the flag _--blocks-max-kb_ only considers pairs of variants within 200 kilobases of each other. The output files from the above command is a _.blocks_ file. The same code may be used to generate output for _EUR, EAS, SAS_ and _AMR_ populations (as it is not possible to generate population-specific information using the _--within_ flag).
Then, in R:
```sh
R
dfr.afr <-read.delim("AFR.blocks.det",sep="",header=T,check.names=F,stringsAsFactors=F) 
colnames(dfr.afr) <- tolower(colnames(dfr.afr))
```
Use the above code to load each of the 5 datasets and set column names to lower case. To plot the data do the following:
```sh
plot (density(dfr.afr$kb), main="LD block length distribution", ylab="Density",xlab="LD block length (Kb)" )
lines (density(dfr.eur$kb), col="blue")
lines (density(dfr.eas$kb), col="red")
lines (density(dfr.amr$kb), col="purple")
lines (density(dfr.sas$kb), col="green")
legend("topright",c("AFR","EAS","EUR","SAS","AMR"), fill=c("black","red","blue","green","purple"))
```





