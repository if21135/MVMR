# MVMR No UKB: Sleep duration
# Running straight after prev. scripts so exposure data is already ready

setwd("/Users/if21135/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Sleep duration")

#14. Select SNP-outcome summary data
#a. Current smokers
memory.limit(size = 80000)
sleep_duration_current <- read_outcome_data(
  "GWASofSleepDurationInCurrentSmokers_imputed.txt",
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
setwd("~/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Sleep duration")

#b. Ever smokers
memory.limit(size = 80000)
sleep_duration_ever <- read_outcome_data(
  "GWASofSleepDurationInEverSmokers_imputed.txt",
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
setwd("~/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Sleep duration")
memory.limit(size = 80000)
sleep_duration_former <- read_outcome_data(
  "GWASofSleepDurationInFormerSmokers_imputed.txt",
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
setwd("~/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Sleep duration")
memory.limit(size = 80000)
sleep_duration_never <- read_outcome_data(
  "GWASofSleepDurationInNeverSmokers_imputed.txt",
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
mvdat <- mv_harmonise_data(Exposures_H, sleep_duration_current, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, sleep_duration_current)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#b.
mvdat <- mv_harmonise_data(Exposures_H, sleep_duration_ever, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, sleep_duration_ever)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#c.
mvdat <- mv_harmonise_data(Exposures_H, sleep_duration_former, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, sleep_duration_former)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#d.
mvdat <- mv_harmonise_data(Exposures_H, sleep_duration_never, harmonise_strictness = 2)
mvdat_1<- harmonise_data(Exposures_H, sleep_duration_never)
mvdat_1<- mvdat_1[mvdat_1$mr_keep== TRUE, ]
str(mvdat_1)

#16. Find proxies to add to missing outcome SNPs - if needed see code from Jaz (again choose from a, b, c or d depending on what your group is)
#a.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, sleep_duration_current$SNP))

#b.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, sleep_duration_ever$SNP))

#c.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, sleep_duration_former$SNP))

#d.
proxy_needed3 <- data.frame(setdiff(EXP1_H$SNP, sleep_duration_never$SNP))


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
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Current")
write.csv(mod_ME_2_or,"No_UKB_Sleep_duration_current_MVMR_Egger.csv", row.names = TRUE)

#b.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Ever")
write.csv(mod_ME_2_or,"No_UKB_Sleep_duration_ever_MVMR_Egger.csv", row.names = TRUE)

#c.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Former")
write.csv(mod_ME_2_or,"No_UKB_Sleep_duration_former_MVMR_Egger.csv", row.names = TRUE)

#d.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Never")
write.csv(mod_ME_2_or,"No_UKB_Sleep_duration_never_MVMR_Egger.csv", row.names = TRUE)


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
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Current")
write.csv(res_w_CIs,"No_UKB_Sleep_duration_current_MVMR.csv", row.names = TRUE)

#b.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Ever")
write.csv(res_w_CIs,"No_UKB_Sleep_duration_ever_MVMR.csv", row.names = TRUE)

#c.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Former")
write.csv(res_w_CIs,"No_UKB_Sleep_duration_former_MVMR.csv", row.names = TRUE)

#d.
setwd("~/Desktop/MP1/Project portal/results/MVMR analyses/MVMR no UKB/Sleep duration/Never")
write.csv(res_w_CIs,"No_UKB_Sleep_duration_never_MVMR.csv", row.names = TRUE)

#Once group analysis is complete, return to point 14 (selecting outcome data) and rerun down to here

#17. Calculate F-statistic. An F-stat >10 is strong
# find the correlation between EXP1 and EXP2 e.g., in Buchwald et al.the correlation between nmr and smoking heaviness  = -0.019
cov <- matrix(c(1,-0.019,-0.019,1), nrow=2, ncol=2)
Xcovmat <-phenocov_mvmr(cov, df_mvmr[,c(6,7)] )
Fstat<- strength_mvmr(df_mvmr, gencov = Xcovmat)

#As before, choose from the four options depending on smoking status you are looking at:
write.csv(Fstat,"No_UKB_Sleep_duration_current_MVMR_FStat.csv", row.names = TRUE)
write.csv(Fstat,"No_UKB_Sleep_duration_ever_MVMR_FStat.csv", row.names = TRUE)
write.csv(Fstat,"No_UKB_Sleep_duration_former_MVMR_FStat.csv", row.names = TRUE)
write.csv(Fstat,"No_UKB_Sleep_duration_never_MVMR_FStat.csv", row.names = TRUE)

#17. Test for horizontal pleiotropy - Q should be greater than the number of SNPs included
ptr<-pleiotropy_mvmr(df_mvmr, gencov=Xcovmat)

#Again, choose from the four options depending on smoking status you are looking at:
write.csv(ptr,"No_UKB_Sleep_duration_current_MVMR_Qstat.csv", row.names = TRUE)
write.csv(ptr,"No_UKB_Sleep_duration_ever_MVMR_Qstat.csv", row.names = TRUE)
write.csv(ptr,"No_UKB_Sleep_duration_former_MVMR_Qstat.csv", row.names = TRUE)
write.csv(ptr,"No_UKB_Sleep_duration_never_MVMR_Qstat.csv", row.names = TRUE)
