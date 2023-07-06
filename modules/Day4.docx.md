## Introduction to Cross-Ancestry PRS analysis
Before starting the practical the following commands will need to be run from within your virtual machine:

(1) conda create -n "PRScsx" python=3.7

(2) conda activate PRScsx

(3) pip install scipy

(4) pip install h5py

 The goal of this practical is to provide you with basic understanding and experience of running the PRS-CSx software. After completing this practical, you should:
* Be able to perform cross-population descriptives.
* Be familiar with running cross-ancestry PRS analyses using PRS-CSx.
* Understand how to evaluate linear models using Akaike’s Information Criterion

#### 1. The 1000 Genomes dataset
The data we will be working with are coming from the 1000 Genomes Project reference panel. The data relates to individuals from 26 different source populations around the world. For simplicity, the populations have been collapsed into 5 broader continental super-populations: East Asian, European, South Asian, Amerindian, African ((EAS, EUR, SAS, EUR and AFR)). The scripts used to download and process the 1000Genomes data for the purposes of this course will be provided in the course appendix at the end of this week. 

#### 2. Cross-population allele frequency
Genetic variation is conveyed using allelic frequencies. Allele frequency is shaped by evolutionary forces and drift.  Here we compare profiles of allele frequency across the five ancestral populations. Global differences in allelic frequency has important implications for the portability of PRS across populations. Using plink it is possible to generate allele frequency statistics for each SNP, across populations, using the annotations provided in the file pop_info.pheno. In _/home/manager/data/Data_Day4_:
```sh
./software/plink_linux --bfile ./data/chr1-22 --freq --within ./data/pop_info.pheno
```
Population-stratified allele frequencies are reported in the output file plink.frq.strat.
For each population, print the numbers of total SNPs to screen, as follows:
```sh
grep -F 'AFR' plink.frq.strat | wc -l
```
From there we can print the number of SNPs with minor allele frequencies greater than 0 (and are hence potentially available for genomic analyes).
```sh
grep -F 'EAS' plink.frq.strat | awk '$6 >0' | wc -l
```
Recycle the code above to find the number of available SNPs in each of the 4 other global populations (EUR, AFR, SAS, AMR).

#### **Questions**
##### (i) Which population contains the most SNPs?
##### (ii) What  is the significance of the observed population order?  
&nbsp;
#### 3. Distribution of allele frequencies
```sh
R
library(dplyr)
library(ggplot2)
freq <-read.table("plink.frq.strat", header =T)
plotDat <- freq %>%
  mutate(AlleleFrequency = cut(MAF, seq(0, 1, 0.25))) %>%
  group_by(AlleleFrequency, CLST) %>%
  summarise(FractionOfSNPs = n()/nrow(freq) * 100)

ggplot(na.omit(plotDat),
       aes(AlleleFrequency, FractionOfSNPs, group = CLST, col = CLST)) +
  geom_line() +
  scale_y_continuous(limits = c(0, 12)) +
  ggtitle("Distribution of allele frequency across genome")
```
#### **Questions**
##### (i) Which population has the most SNPs?
##### (ii) What  is the significance of the observed population ordering?
##### (iii) What is the reason behind these two features?
&nbsp;
#### 4. Calculation of Fst
Fst is a formal metric which is used to convey the level of genetic divergence between populations (on a scale between 0 and 1), using information derived from a set of genome-wide and mutually independent SNPs. Fst between parwise populations is estimated efficiently in Plink-2. So first we need to download Plink-2 using the command below:
```sh
sudo apt install plink2
```
A higher Fst corresponds to greater divergence between populations. Use the following command to calculate Fst:
```sh
plink2 --bfile ./data/chr1-22 --fst POP --pheno ./data/pop_info.pheno
```
Check the output file plink2.fst.summary
#### **Questions**
##### (i) Which population pairs have the highest Fst ?
##### (ii) For which populations is Fst smallest??
&nbsp;
## Introduction to PRS-CSx
#### 5. Background to PRS-CSX
PRS-CSx is a Python based command line tool that integrates GWAS summary statistics and external LD reference panels from multiple populations to improve cross-population polygenic prediction. We will be using simulated trait data pertaininng to systolic blood pressure (SBP) to explore PRS performance using 2 target populations that consist of 650 individuals of African ancestry and 500 individuals of European ancestry. Please note when running PRSice that the object of the flag "--prsice" will change according to whether plink is being called within the linux-like environment of the virtual machine (PRSice_linux) or a mac (PRSice_mac). Both executables can be found in the _/home/manager/data/Data_Day4_ directory. 

#### 6. Extraction of SNPs
PRS-CSx uses only HAPMAP3 SNPs therefore we produce a set of plink files containing this SNP subset.
```sh
./software/plink_linux --bfile ./data/1kg.eur.dbSNP153.hg38 --extract ./data/csxsnp --make-bed --out ./data/EUR_1kg.hm3.only.csx

./software/plink_linux --bfile ./data/1kg.afr.dbSNP153.hg38 --extract ./data/csxsnp --make-bed --out ./data/AFR_1kg.hm3.only.csx
```

#### 7. Running PRS-CSx
To model the coupling of effect sizes at individual SNPs across ancestries PRS-CSx uses an MCMC (Bayesian) sampling algorithm to determine values of the global shrinkage parameter ("phi") by Maximum likelihood. For samples of mixed or admixed genetic ancestry (which ours are not) the optimal value of the shrinkage parameter is estimated autonomously from the data. Here we use the value of phi (1e-4), which is suggested by the software authors, given that our trait is simulated to have a relatively small number (N=110) causal variants, distributed genome-wide.
To save time, we will be running the analyses across chromosome 15, rather than the entire genome. The commands needed to run PRS-CSx are contained in two script files, located in /home/manager/data/Data_Day4/scripts. The file run_prscsx_afr-target.sh is used to estimate optimal SNP weights for predicting into the African target population. 

NB - move the snp file to the reference folder

```
cp /home/manager/PRScsx/snpinfo_mult_1kg_hm3 /home/manager/data/Data_Day4/reference/csx    
```

NB - the ld block files must be moved to /home/manager/data/Data_Day4/reference/csx to /home/manager/PRScsx/

```
mv /home/manager/PRScsx/*.tar.gz /home/manager/PRScsx/
```
then extracted with tar

```
tar -xvfz ld...tar.gz
```
be careful of space, your vm has 100GB cap. remove the .gz files once extracted


**Script contents**:
```sh
python /home/manager/data/Data_Day4/software/PRScsx.py \
--ref_dir=/home/manager/data/Data_Day4/reference/csx \
--bim_prefix=/home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
--sst_file=/home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx,/home/manager/data/Data_Day4/data/AFR-SBP-simulated.sumstats.prscsx \
--n_gwas=25732,4855 \
--pop=EUR,AFR \
--chrom=15 \
--phi=1e-4 \
--out_dir=/home/manager/data/Data_Day4/out/csx \
--out_name=afr.target.csx
```
The file run_prscsx_eur-target.sh is used to estimate optimal SNP weights for predicting into the European target population:
**Script contents**:
```sh
python /home/manager/data/Data_Day4/software/PRScsx.py \
--ref_dir=/home/manager/data/Data_Day4/reference/csx \
--bim_prefix=/home/manager/data/Data_Day4/data/EUR_1kg.hm3.only.csx \
--sst_file=/home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx,/home/manager/data/Data_Day4/data/AFR-SBP-simulated.sumstats.prscsx \
--n_gwas=25732,4855 \
--pop=EUR,AFR \
--chrom=15 \
--phi=1e-4 \
--out_dir=/home/manager/data/Data_Day4/out/csx \
--out_name=eur.target.csx
```
Prior to running each script you will need to update the first line with a reference (i.e the pathname) to your home directory. From _/home/manager/data/Data_Day4_ the scripts can then be run as follows:
```sh
./scripts/run_prscsx_afr-target.sh
./scripts/run_prscsx_eur-target.sh
```
#### **Questions**
##### (i) How many results files do you see in the output directory?
##### (ii) What does each file correspond to?
##### (ii) In which column are the adjusted SNP weights contained?      
&nbsp;
#### 8. Processing
The next step would be to collate adjusted SNP weight information from multiple chromosomes and relabel the consolidated files for clarity. The following code works regardless of whether the preceding PRS-CSx analysis was performed for single or multiple chromosomes. In _/home/manager/data/Data_Day4/out/csx_
```sh
for file in afr.target.csx_AFR_*; do
cat $file >> posteriors.afr.by.afr
done

for file in afr.target.csx_EUR_*; do
cat $file >> posteriors.afr.by.eur
done

for file in eur.target.csx_AFR_*; do
cat $file >> posteriors.eur.by.afr
done

for file in eur.target.csx_EUR_*; do
cat $file >> posteriors.eur.by.eur
done
```

#### 9. Data processing in R
The next stage of data processing, done in R will create a new set of summary statistics files, containing the adjusting SNP weights from PRS-CSx. Staying in _/home/manager/data/Data_Day4/out/csx_ :
```sh
R
# Load posteriors effect sizes from PRS-CSx output

afr.afr<-read.table("posteriors.afr.by.afr", col.names=c("chr","SNP","bp","a1","a2","post.afr.afr"))
afr.by.eur<-read.table("posteriors.afr.by.eur", col.names=c("chr","SNP","bp","a1","a2","post.afr.by.eur"))
eur.eur<-read.table("posteriors.eur.by.eur", col.names=c("chr","SNP","bp","a1","a2","post.eur.eur"))
eur.by.afr<-read.table("posteriors.eur.by.afr", col.names=c("chr","SNP","bp","a1","a2","post.eur.by.afr"))

# Load original summary statistics
afr.SBP<-read.table("/home/manager/data/Data_Day4/data/AFR-SBP-simulated.sumstats.prscsx", header=T)
eur.SBP<-read.table("/home/manager/data/Data_Day4/data/EUR-SBP-simulated.sumstats.prscsx", header=T)

# Combine the posterior derived weights and summary statistics
eur.SBP.merge1<-merge(x=eur.SBP, y=eur.eur[c(2,6)], by="SNP", all.y=T)
eur.SBP.merge<-merge(x=eur.SBP.merge1, y=afr.by.eur[c(2,6)], by="SNP", all.y=T)
afr.SBP.merge1<-merge(x=afr.SBP, y=afr.afr[c(2,6)], by="SNP", all.y=T)
afr.SBP.merge<-merge(x=afr.SBP.merge1, y=eur.by.afr[c(2,6)], by="SNP", all.y=T)

# Save files
write.table(afr.SBP.merge, "/home/manager/data/Data_Day4/afr.SBP.posterior.sumstats", quote=F, row.names=F)
write.table(eur.SBP.merge, "/home/manager/data/Data_Day4/eur.SBP.posterior.sumstats", quote=F, row.names=F)

q()
```

#### 10. Use PRSice to apply the adjusted SNP effects to target phenotypes
We first need to create a list of SNPs to be used as input, based on the previous PRS-CSx analysis. This is done in location _/home/manager/data/Data_Day4_.
```sh
awk 'NR>1 {print $1}' afr.SBP.posterior.sumstats > snps.afr.posterior
awk 'NR>1 {print $1}' eur.SBP.posterior.sumstats > snps.eur.posterior
```
The next series of commands implement in-ancestry and cross-ancestry prediction of the systolic blood pressure phenotype using the PRS-CSx optimised SNP weights. **Do not forget** to exchange _PRSice_linux_ for _PRSice_mac_ if running the commands in a Mac environment. The phenotypic files sbp_afr_1kg.sim_pheno and sbp_eur_1kg.sim_pheno contain data on systolic blood pressure with simulated heritabilities that vary from 10%, 20%, 33%, 50% and 100%. To compensate for the fact that the adjusted weights generated are based on chromosome 15, rather than the entire genome, we will be using the **_100%_** trait version (pheno100) for this analysis.

##### 11. Predicting from African training to African target data
```sh
Rscript /home/manager/data/Data_Day4/software/PRSice.R \
--prsice /home/manager/data/Data_Day4/software/PRSice_linux \
--base /home/manager/data/Data_Day4/afr.sbp.posterior.sumstats \
--extract /home/manager/data/Data_Day4/snps.afr.posterior \
--A1 A1 \
--pvalue P \
--no-clump \
--stat post.afr.afr \
--beta \
--snp SNP \
--score sum \
--target /home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
--binary-target F \
--pheno /home/manager/data/Data_Day4/data/sbp_afr_1kg.sim_pheno \
--pheno-col pheno100 \
--thread 8 \
--out /home/manager/data/Data_Day4/out/prscsx_prsice/SBP.afr.afr
```
##### 12. Predicting from European training to African target data
```sh
Rscript /home/manager/data/Data_Day4/software/PRSice.R \
--prsice /home/manager/data/Data_Day4/software/PRSice_linux \
--base /home/manager/data/Data_Day4/eur.SBP.posterior.sumstats \
--extract /home/manager/data/Data_Day4/snps.eur.posterior \
--A1 A1 \
--pvalue P \
--no-clump \
--stat post.afr.by.eur \
--beta \
--snp SNP \
--score sum \
--target /home/manager/data/Data_Day4/data/AFR_1kg.hm3.only.csx \
--binary-target F \
--pheno /home/manager/data/Data_Day4/data/sbp_afr_1kg.sim_pheno \
--pheno-col pheno100 \
--thread 8 \
--out /home/manager/data/Data_Day4/out/prscsx_prsice/SBP.afr.by.eur
```
##### 13. Predicting from European training to European target data
```sh
Rscript /home/manager/data/Data_Day4/software/PRSice.R \
--prsice /home/manager/data/Data_Day4/software/PRSice_linux \
--base /home/manager/data/Data_Day4/eur.SBP.posterior.sumstats \
--extract /home/manager/data/Data_Day4/snps.eur.posterior \
--A1 A1 \
--pvalue P \
--no-clump \
--stat post.eur.eur \
--beta \
--snp SNP \
--score sum \
--target /home/manager/data/Data_Day4/data/EUR_1kg.hm3.only.csx \
--binary-target F \
--pheno /home/manager/data/Data_Day4/data/sbp_eur_1kg.sim_pheno \
--pheno-col pheno100 \
--thread 8 \
--out /home/manager/data/Data_Day4/out/prscsx_prsice/SBP.eur.eur
```
##### 14. Predicting from African training to European target data
```sh
Rscript /home/manager/data/Data_Day4/software/PRSice.R \
--prsice /home/manager/data/Data_Day4/software/PRSice_linux \
--base /home/manager/data/Data_Day4/afr.SBP.posterior.sumstats \
--extract /home/manager/data/Data_Day4/snps.afr.posterior \
--A1 A1 \
--pvalue P \
--no-clump \
--stat post.eur.by.afr \
--beta \
--snp SNP \
--score sum \
--target /home/manager/data/Data_Day4/data/EUR_1kg.hm3.only.csx \
--binary-target F \
--pheno /home/manager/data/Data_Day4/data/sbp_eur_1kg.sim_pheno \
--pheno-col pheno100 \
--thread 8 \
--out /home/manager/data/Data_Day4/out/prscsx_prsice/SBP.eur.by.afr
```
#### 15. Create summary files
From location _/home/manager/data/Data_Day4/out/prscsx_prsice_
```sh
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.afr.afr.summary > quicksum.afr-afr
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.afr.by.eur.summary > quicksum.afr-by-eur
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.eur.by.afr.summary > quicksum.eur-by-afr
awk '{print $1,"\t",$4,"\t",$8,"\t",$10}' SBP.eur.eur.summary > quicksum.eur-eur
```
#### **Questions**
##### **(i) What information is being summarized in the output files??** 
##### **(ii) For each prediction model what is the R2 and  value???** 
&nbsp;

#### 16. Complete the remaining PRS analysis in R
##### Step 1: Load and prepare data
```sh 
R
# Load libraries
sudo apt install cmake
install.packages("AICcmodavg")
install.packages("fmsb")
library(AICcmodavg)
library(fmsb)

afr.afr<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.afr.afr.best", header=T)
afr.eur<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.afr.by.eur.best", header=T)
eur.eur<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.eur.eur.best", header=T)
eur.afr<-read.table("/home/manager/data/Data_Day4/out/prscsx_prsice/SBP.eur.by.afr.best", header=T)

colnames(afr.afr)[4]<-"afr.afr"   # african prediction using african PRS
colnames(afr.eur)[4]<-"afr.eur"   # african prediction using european PRS
colnames(eur.eur)[4]<-"eur.eur" # european prediction using european PRS
colnames(eur.afr)[4]<-"eur.afr" # european prediction using african PRS

# Merge PRSs according to target ancestry and source population
combined.afr<-merge(x=afr.afr, y=afr.eur[c(2,4)], by="IID", all.x=T)
combined.eur<-merge(x=eur.eur, y=eur.afr[c(2,4)], by="IID", all.x=T)

# Rescale scores to have mean = ‘0’ and standard deviation = ‘1’
combined.afr[4:5] <- as.data.frame(scale(combined.afr[4:5]))
combined.eur[4:5] <- as.data.frame(scale(combined.eur[4:5]))

# Load phenotype data
sbp.eur<-read.table("/home/manager/data/Data_Day4/data/SBP_eur_1kg.sim_pheno", header=T)
sbp.afr<-read.table("/home/manager/data/Data_Day4/data/SBP_afr_1kg.sim_pheno", header=T)

# Merge phenotypes and PRS
combined.eur.pheno<-merge(x=combined.eur[-c(2,3)], y=sbp.eur[,c(2,3)], by="IID", all.x=T)
combined.afr.pheno<-merge(x=combined.afr[-c(2,3)], y=sbp.afr[,c(2,3)], by="IID", all.x=T)
```
##### Step 2: Model selection using Akaike’s Information Criterion
In the final step we want to determine whether the use of the multi-ancestry PRS formulation generated by PRS-CSx performs better at predicting systolic blood pressure, compared to the single-ancestry PRS formulation. To do this we use (AIC Akaike’s Information Criterion). AIC is used to assess the performance of a competing set of regression models. We use it to compare the performance of models that contain different numbers of predictors, given that R2 is inflated by the inclusion of redundant terms in a model. 

##### 16b. Model evaluation
```sh 
# EUROPEAN
model1.eur <- glm(pheno100 ~ eur.eur, data = combined.eur.pheno, family=gaussian)
model2.eur <- glm(pheno100 ~ eur.afr, data = combined.eur.pheno, family=gaussian)
model3.eur <- glm(pheno100 ~ eur.eur + eur.afr, data = combined.eur.pheno, family=gaussian)

models.eur <- list(model1.eur, model2.eur, model3.eur) # define model set
mod.names.eur <- c('eur.eur', 'eur.afr', 'eur.combined') # add model names
aictab(cand.set = models.eur, modnames = mod.names.eur) # calculate model AICs model

# AFRICAN
model1.afr <- glm(pheno100 ~ afr.afr, data = combined.afr.pheno, family=gaussian)
model2.afr <- glm(pheno100 ~ afr.eur, data = combined.afr.pheno, family=gaussian)
model3.afr <- glm(pheno100 ~ afr.afr + afr.eur, data = combined.afr.pheno, family=gaussian)

models.afr <- list(model1.afr, model2.afr, model3.afr) # define model set
mod.names.afr <- c('afr.afr', 'afr.eur', 'afr.combined') # add model names
aictab(cand.set = models.afr, modnames = mod.names.afr) # calculate model AICs

# NagelkerkeR2 calculations in R
# African target
NagelkerkeR2(model1.afr)
NagelkerkeR2(model2.afr)
NagelkerkeR2(model3.afr)

# European target
NagelkerkeR2(model1.eur)
NagelkerkeR2(model2.eur)
NagelkerkeR2(model3.eur)

For each ancestry which predictive model performs best and worst?
Does combining the two sets of adjusted scores perform consistently better than modelling each one separately?
```
#### **Questions**
##### **(i) For each ancestry which model performs best and which performs worst?**   
##### **(ii) Does linearly combining adjusted scores (European plus African) always result in better performance compared to single-ancestry prediction?**   


