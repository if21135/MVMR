### MVMR analyses exploring the direct effects of nicotine and cigarette smoke on 
#sleep behaviours using data excluding UK Biobank and 23andme

#install.packages
library(TwoSampleMR)
library(ieugwasr)
library(googleAuthR)
library(tidyverse)
library(stringr)
library(dplyr)
library(forestplot)
library(plyr)
library(gtable)
library(psych)
install_github("WSpiller/RadialMR")
library(RadialMR)
library(reshape)
library(gplots)
require(ggplot2)
library(ggplot2)
library(gridExtra)
library(grid)
library(extrafont)
library(plotly)
library(data.table)
library(curl)
install.packages("jsonlite", type = "source")
install_github("WSpiller/MVMR", build_opts = c("--no-resave-data", "--no-manual"), build_vignettes = TRUE)
library(MVMR)
library(MendelianRandomization)
library(simex)

#Set directory and memory
setwd("~/Desktop")

#1. Read in exposures
nmr_cpd_noukb <- read.csv("exp1_exp2_noukb.csv", header = TRUE)

#2. Check for overlapping SNPs
n_occur<- data.frame(table(nmr_cpd_noukb$SNP)) #59
n_occur[n_occur$Freq >1, ]
nmr_cpd_noukb[nmr_cpd_noukb$SNP %in% n_occur$Var1[n_occur$Freq >1], ]

#3. Read in SNP lists for EXP1, EXP2 and both
EXP1_SNPlist <- read_csv("MP1/Project portal/data/NMR and CPD exposure data for MVMR/SNP lists/SNPlist_NMR.txt", col_names = FALSE)
EXP2_SNPlist <- read_csv("MP1/Project portal/data/NMR and CPD exposure data for MVMR/SNP lists/SNPlist_CPD.txt", col_names = FALSE)
EXP1_EXP2_SNPlist <- read_csv("MP1/Project portal/data/NMR and CPD exposure data for MVMR/SNP lists/NMR_CPD_SNPlist.txt", col_names = FALSE)

colnames(EXP1_SNPlist)[colnames(EXP1_SNPlist) == "X1"] <- "SNP"
colnames(EXP2_SNPlist)[colnames(EXP2_SNPlist) == "X1"] <- "SNP"
colnames(EXP1_EXP2_SNPlist)[colnames(EXP1_EXP2_SNPlist) == "X1"] <- "SNP"
EXP1_EXP2_SNPlist <- EXP1_EXP2_SNPlist[-1, ]

#4. If one exposure is stronger than another, you may want to only drop SNPs from the stronger instrument during the clumping process. 
#To do so, save the p-values for the weaker exposure in another variable, then set the p-value column to the lowest possible value.

#4a. save old p-values 
nmr_cpd_noukb$oldpvalues <-nmr_cpd_noukb$pval.exposure

#4b. Change all p-values for EXP1 to 1e-200 for clumping so that none are dropped
nmr_cpd_noukb<- nmr_cpd_noukb %>% 
  mutate(pval.exposure = if_else(nmr_cpd_noukb$SNP %in% EXP1_SNPlist, 1e-201, pval.exposure))

#5. Clump
nmr_cpd_noukb$id.exposure[nmr_cpd_noukb$id.exposure == "2"] <- "1"
nmr_cpd_noukb<- clump_data(nmr_cpd_noukb, clump_kb=500, clump_r2=0.1) 
str(nmr_cpd_noukb)

#6. Keep only snps that are present across both exposures - they would have frequency 1 if only available in one dataset
n_occur <- data.frame(table(nmr_cpd_noukb$SNP))
n_occur[n_occur$Freq == 2,]
nmr_cpd_noukb<- nmr_cpd_noukb[nmr_cpd_noukb$SNP %in% n_occur$Var1[n_occur$Freq == 2],]
str(nmr_cpd_noukb)

#7. Add ID's back
nmr_cpd_noukb$id.exposure[nmr_cpd_noukb$samplesize.exposure<6000] <- "1"
nmr_cpd_noukb$id.exposure[nmr_cpd_noukb$samplesize.exposure>6000] <- "2"

#8. Revert all p-values for EXP1 from 1e-200 back to their original value if changed above.
nmr_cpd_noukb$pval.exposure<-nmr_cpd_noukb$oldpvalues
nmr_cpd_noukb<-select(nmr_cpd_noukb,-c(oldpvalues))

#9. Split to harmonise based on exposure id
EXP1 = split(nmr_cpd_noukb, nmr_cpd_noukb$id.exposure)[['1']]
EXP2 = split(nmr_cpd_noukb, nmr_cpd_noukb$id.exposure)[['2']]

#10. Harmonise EXP1 on EXP2
names(EXP1) = gsub( "exposure", "outcome", names(EXP1))
EXP1_EXP2 = harmonise_data(EXP2, EXP1)

#11. Keep only snps where mr_keep= TRUE#
EXP1_EXP2 = EXP1_EXP2[EXP1_EXP2$mr_keep== TRUE, ]
str(EXP1_EXP2) 

#12a. Split the tables#EXP2#
EXP2_H<- subset(EXP1_EXP2, id.exposure== id.exposure[2], select= c(SNP, exposure, id.exposure, effect_allele.exposure, other_allele.exposure, beta.exposure, se.exposure, pval.exposure, eaf.exposure))

#12b. Split the tables#EXP1#
EXP1_H<- subset(EXP1_EXP2, id.outcome== id.outcome[1], select= c(SNP, outcome, id.outcome, effect_allele.outcome, other_allele.outcome, beta.outcome, se.outcome, pval.outcome, eaf.outcome))

#13. Turn EXP1 to exposure to merge the datasets#
names(EXP1_H) <- gsub("outcome", "exposure", names(EXP1_H))
Exposures_H<- merge(EXP1_H, EXP2_H, all= TRUE)
Exposures_H["Phenotype"]<- NA
Exposures_H$Phenotype[Exposures_H$id.exposure == 1] <- "EXP1"
Exposures_H$Phenotype[Exposures_H$id.exposure == 2] <- "EXP2"
str(Exposures_H)

setwd("/Users/if21135/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Chronotype")

#14. Select SNP-outcome summary data
#a. Current smokers
memory.limit(size = 80000)
chronotype_current <- read_outcome_data(
  "GWASofChronotypeInCurrentSmokers_imputed.txt",
  snps = Exposures_H$SNP,
  sep="\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

#Remember to reset wd when returning to this point
setwd("~/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Chronotype")

#b. Ever smokers
memory.limit(size = 80000)
chronotype_ever <- read_outcome_data(
  "GWASofChronotypeInEverSmokers_imputed.txt",
  snps = Exposures_H$SNP,
  sep="\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)


#c. Former smokers
setwd("~/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Chronotype")
memory.limit(size = 80000)
chronotype_former <- read_outcome_data(
  "GWASofChronotypeInFormerSmokers_imputed.txt",
  snps = Exposures_H$SNP,
  sep="\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)

#d. Never smokers
setwd("~/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Chronotype")
memory.limit(size = 80000)
chronotype_never <- read_outcome_data(
  "GWASofChronotypeInNeverSmokers_imputed.txt",
  snps = Exposures_H$SNP,
  sep="\t",
  snp_col = "SNP",
  beta_col = "BETA",
  se_col = "SE",
  effect_allele_col = "ALLELE1",
  other_allele_col = "ALLELE0",
  eaf = "A1FREQ",
  pval_col = "P_BOLT_LMM_INF",
  samplesize_col = 9859,
  min_pval = 1e-200,
  log_pval = FALSE
)


#15. #Harmonise with outcome (choose from a, b, c or d depending on what your group is)
#a.
mvdat <- mv_harmonise_data(Exposures_H, chronotype_current, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, chronotype_current)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#b.
mvdat <- mv_harmonise_data(Exposures_H, chronotype_ever, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, chronotype_ever)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#c.
mvdat <- mv_harmonise_data(Exposures_H, chronotype_former, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, chronotype_former)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#d.
mvdat <- mv_harmonise_data(Exposures_H, chronotype_never, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, chronotype_never)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#16. Find proxies to add to missing outcome SNPs - if needed see code from Jaz (again choose from a, b, c or d depending on what your group is)
#a.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, chronotype_current$SNP))

#b.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, chronotype_ever$SNP))

#c.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, chronotype_former$SNP))

#d.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, chronotype_never$SNP))


#17. Run MVMR analyses#

#17a. IVW#
bX1<- c(mvdat_1$beta.exposure[mvdat_1$id.exposure== 1])
bX2<- c(mvdat_1$beta.exposure[mvdat_1$id.exposure== 2])
bY<- c(mvdat_1$beta.outcome[mvdat_1$id.exposure== 1])
bYse<- c(mvdat_1$se.outcome[mvdat_1$id.exposure== 1])

set.seed(1234)
mod.MVMR<-lm(bY~bX1+bX2-1, weights=bYse^-2)
se_theta1MI.random = summary(lm(bY~bX1+bX2-1, weights=bYse^-2))$coef[1,2]/
  min(summary(lm(bY~bX1+bX2-1, weights=bYse^-2))$sigma,1)

mod<- summary(mod.MVMR)

mod_or <- coef(summary(mod.MVMR))
colnames(mod_or) <- c("b", "se", "t", "p")
mod_or<-as.data.frame(mod_or)
mod_or<-generate_odds_ratios(mod_or)

#17b. Orientation EXP1#
clist<-c("bX2","bY")
for (var in clist){
  eval(parse(text=paste0(var,"<-ifelse(bX1>0,",var,",",var,"*-1)")))
}
bX1<-abs(bX1)

#17c. MVMR Egger#
mod.MVMRME <- summary(lm(bY~bX1+bX2, weights=bYse^-2))
se_theta1ME.random = summary(lm(bY~bX1+bX2, weights=bYse^-2))$coef[2,2]/
  min(summary(lm(bY~bX1+bX2, weights=bYse^-2))$sigma,1)
mod_ME<- summary(mod.MVMRME)

mod_ME_or <- data.frame(mod.MVMRME[["coefficients"]])
colnames(mod_ME_or) <- c("b", "se", "t", "p")
mod_ME_or<-as.data.frame(mod_ME_or)
mod_ME_or<-generate_odds_ratios(mod_ME_or)

#17d. Orientation EXP2#
clist<-c("bX1","bY")
for (var in clist){
  eval(parse(text=paste0(var,"<-ifelse(bX2>0,",var,",",var,"*-1)")))
}
bX2<-abs(bX2)

#17d. MVMR Egger#
mod.MVMRME_2 <- summary(lm(bY~bX1+bX2, weights=bYse^-2))
se_theta1ME.random = summary(lm(bY~bX1+bX2, weights=bYse^-2))$coef[2,2]/
  min(summary(lm(bY~bX1+bX2, weights=bYse^-2))$sigma,1)
mod_ME_2<- summary(mod.MVMRME_2)

mod_ME_2_or <- data.frame(mod.MVMRME_2[["coefficients"]])
colnames(mod_ME_2_or) <- c("b", "se", "t", "p")
mod_ME_2_or<-as.data.frame(mod_ME_2_or)
mod_ME_2_or<-generate_odds_ratios(mod_ME_2_or)

#Saving MVMR Egger model with ORs
#a.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Chronotype/Current")
write.csv(mod_ME_2_or,"No_UKB_Chronotype_current_MVMR_Egger_noukb.csv", row.names = TRUE)

#b.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Chronotype/Ever")
write.csv(mod_ME_2_or,"No_UKB_Chronotype_ever_MVMR_Egger.csv", row.names = TRUE)

#c.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Chronotype/Former")
write.csv(mod_ME_2_or,"No_UKB_Chronotype_former_MVMR_Egger.csv", row.names = TRUE)

#d.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Chronotype/Never")
write.csv(mod_ME_2_or,"No_UKB_Chronotype_never_MVMR_Egger.csv", row.names = TRUE)


#17f. Format to analyse with MVMR package# 
bX1<- c(mvdat_1$beta.exposure[mvdat_1$id.exposure== 1])
bX2<- c(mvdat_1$beta.exposure[mvdat_1$id.exposure== 2])
bY<- c(mvdat_1$beta.outcome[mvdat_1$id.exposure== 1])
bYse<- c(mvdat_1$se.outcome[mvdat_1$id.exposure== 1])
bXse1<- c(mvdat_1$se.exposure[mvdat_1$id.exposure== 1])
bXse2<- c(mvdat_1$se.exposure[mvdat_1$id.exposure== 2])
df<- data.frame(bX1, bXse1, bX2, bXse2, bY, bYse)

library(MVMR)
df_mvmr <- format_mvmr(df[, c(1, 3)], df[,5], df[,c(2, 4)], df[, 6])

#17g. cross check result with mvmr package
res<- ivw_mvmr(df_mvmr)

#Tried creating CIs but kept saying the replacement has 0 rows so can't do it. Decided to use CIs from mod_or instead as presumably if cross-checked results are fine then this is okay?
#df_or<-as.data.frame(df)
#df_or <- generate_odds_ratios(df_or)

#Attempting to add CI columns from mod_or to res...
res_w_CIs <- cbind(res, mod_or$lo_ci, mod_or$up_ci)
colnames(res_w_CIs)[5] = "lo_ci"
colnames(res_w_CIs)[6] = "up_ci"

#Saving csv files for current, ever, former and never smokers (again choose from a, b, c or d depending on what your group is)
#a.
setwd("~/Desktop/MP1/Project portal/results/MVMR analysesMVMR no UKB/Chronotype/Current")
write.csv(res_w_CIs,"No_UKB_Chronotype_current_MVMR_noukb.csv", row.names = TRUE)

#b.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Chronotype/Ever")
write.csv(res_w_CIs,"No_UKB_Chronotype_ever_MVMR.csv", row.names = TRUE)

#c.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Chronotype/Former")
write.csv(res_w_CIs,"No_UKB_Chronotype_former_MVMR.csv", row.names = TRUE)

#d.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Chronotype/Never")
write.csv(res_w_CIs,"No_UKB_Chronotype_never_MVMR.csv", row.names = TRUE)

#Once group analysis is complete, return to point 14 (selecting outcome data) and rerun down to here

#17. Calculate F-statistic. An F-stat >10 is strong
# find the correlation between EXP1 and EXP2 e.g., in Buchwald et al.the correlation between nmr and smoking heaviness  = -0.019
cov <- matrix(c(1,-0.019,-0.019,1), nrow=2, ncol=2)
Xcovmat <-phenocov_mvmr(cov, df_mvmr[,c(6,7)] )
Fstat<- strength_mvmr(df_mvmr, gencov = Xcovmat)

#As before, choose from the four options depending on smoking status you are looking at:
write.csv(Fstat,"No_UKB_Chronotype_current_MVMR_FStat_noukb.csv", row.names = TRUE)
write.csv(Fstat,"No_UKB_Chronotype_ever_MVMR_FStat.csv", row.names = TRUE)
write.csv(Fstat,"No_UKB_Chronotype_former_MVMR_FStat.csv", row.names = TRUE)
write.csv(Fstat,"No_UKB_Chronotype_never_MVMR_FStat.csv", row.names = TRUE)

#17. Test for horizontal pleiotropy - Q should be greater than the number of SNPs included
ptr<-pleiotropy_mvmr(df_mvmr, gencov=Xcovmat)

#Again, choose from the four options depending on smoking status you are looking at:
write.csv(ptr,"No_UKB_Chronotype_current_MVMR_Qstat_noukb.csv", row.names = TRUE)
write.csv(ptr,"No_UKB_Chronotype_ever_MVMR_Qstat.csv", row.names = TRUE)
write.csv(ptr,"No_UKB_Chronotype_former_MVMR_Qstat.csv", row.names = TRUE)
write.csv(ptr,"No_UKB_Chronotype_never_MVMR_Qstat.csv", row.names = TRUE)