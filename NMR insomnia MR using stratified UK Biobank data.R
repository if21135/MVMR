#########################################################################################
# Assessing the effect of NMR on sleep behaviours using stratified data from UK Biobank #
#########################################################################################

#install.packages("devtools")
#library(devtools)
#devtools::install_github("MRCIEU/TwoSampleMR") #to update the package
#devtools::install_github("MRCIEU/MRInstruments")
#install.packages("plyr")
#install.packages("ggplot2")
#install.packages("png")
#install.packages("ggrepel")
#install.packages("ggthemes")
#install.packages("calibrate")
library(calibrate)
library(ggrepel)
library(ggthemes)
library(TwoSampleMR)
library(MRInstruments)
library(plyr) 
library(ggplot2)
library(png)

# Set working directory
setwd("/Users/if21135/Desktop/MP1/Project portal/data/Unzipped stratified UK Biobank data for MVMR/Insomnia")

##############
# Insomnia #
##############

######################################
#1. Select SNP-exposure summary data #
######################################
nmr <- read_exposure_data(filename="nmr_dat_mr.csv", sep = ",", clump = FALSE, phenotype_col="phenotype", snp_col = "SNP", beta_col = "beta.exposure", se_col = "se.exposure", eaf_col="eaf.exposure", effect_allele_col = "effect_allele.exposure", other_allele_col = "other_allele.exposure", pval_col = "pval.exposure")
#SNP causing noise already removed

#####################################
#2. Select SNP-outcome summary data #
#####################################

#As there are four groups of smoking status, have copied and pasted code for each group. Don't run whole script at once
#but instead run the appropriate line for whichever group is being analysed

insomnia_current <- read_outcome_data(filename="GWASofInsomniaInCurrentSmokers_imputed.txt", sep="\t", snp_col = "SNP", chr_col = "CHR", pos_col = "BP", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf = "A1FREQ", beta_col = "BETA", se_col = "SE", pval_col = "P_BOLT_LMM_INF")
insomnia_ever <- read_outcome_data(filename="GWASofInsomniaInEverSmokers_imputed.txt", sep="\t", snp_col = "SNP", chr_col = "CHR", pos_col = "BP", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf = "A1FREQ", beta_col = "BETA", se_col = "SE", pval_col = "P_BOLT_LMM_INF")
insomnia_former <- read_outcome_data(filename="GWASofInsomniaInFormerSmokers_imputed.txt", sep="\t", snp_col = "SNP", chr_col = "CHR", pos_col = "BP", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf = "A1FREQ", beta_col = "BETA", se_col = "SE", pval_col = "P_BOLT_LMM_INF")
insomnia_never <- read_outcome_data(filename="GWASofInsomniaInNeverSmokers_imputed.txt", sep="\t", snp_col = "SNP", chr_col = "CHR", pos_col = "BP", effect_allele_col = "ALLELE1", other_allele_col = "ALLELE0", eaf = "A1FREQ", beta_col = "BETA", se_col = "SE", pval_col = "P_BOLT_LMM_INF")

########################
#3. Harmonise datasets #
########################
dat <- harmonise_data(nmr, insomnia_current, action = 2)
dat <- harmonise_data(nmr, insomnia_ever, action = 2)
dat <- harmonise_data(nmr, insomnia_former, action = 2)
dat <- harmonise_data(nmr, insomnia_never, action = 2)

#Harmonise negative effect alleles (negative exposure effects need to be made positive for some analyses e.g., MR Egger)
for(i in 1:length(dat$beta.exposure)){
  if(dat$beta.exposure[i]<0){dat$beta.outcome[i] <- -1*dat$beta.outcome[i] }
  if(dat$beta.exposure[i]<0){dat$beta.exposure[i] <- -1*dat$beta.exposure[i] }
}

######################################################
#4. Estimate the causal effect of NMR on insomnia #
######################################################
mr_results <- mr(dat, method_list=c("mr_ivw","mr_egger_regression","mr_weighted_median", "mr_weighted_mode"))
#Add CIs (ignore or, or_lci95 and or_uci95)
mr_results<-generate_odds_ratios(mr_results)

#Export results
setwd("/Users/if21135/Desktop/MP1/Project portal/results/MR analyses/UVMR stratified/NMR/Insomnia/Current")
setwd("/Users/if21135/Desktop/MP1/Project portal/results/MR analyses/UVMR stratified/NMR/Insomnia/Ever")
setwd("/Users/if21135/Desktop/MP1/Project portal/results/MR analyses/UVMR stratified/NMR/Insomnia/Former")
setwd("/Users/if21135/Desktop/MP1/Project portal/results/MR analyses/UVMR stratified/NMR/Insomnia/Never")
write.csv(mr_results,"./NMR_insomnia_current_stratified_UKB_results.csv")
write.csv(mr_results,"./NMR_insomnia_ever_stratified_UKB_results.csv")
write.csv(mr_results,"./NMR_insomnia_former_stratified_UKB_results.csv")
write.csv(mr_results,"./NMR_insomnia_never_stratified_UKB_results.csv")

# B.
# Run some sensitivity analyses

#Test instrument strength
F = dat$beta.exposure^2/dat$se.exposure^2
mF = mean(F)
write.csv(mF, "./NMR_insomnia_current_stratified_UKB_F_stat.csv")
write.csv(mF, "./NMR_insomnia_ever_stratified_UKB_F_stat.csv")
write.csv(mF, "./NMR_insomnia_former_stratified_UKB_F_stat.csv")
write.csv(mF, "./NMR_insomnia_never_stratified_UKB_F_stat.csv")

nmr_insomnia_uvmr_strat_heterogeneity<-mr_heterogeneity(dat)
write.csv(nmr_insomnia_uvmr_strat_heterogeneity,"./NMR_insomnia_current_stratified_UKB_heterogeneity.csv")
write.csv(nmr_insomnia_uvmr_strat_heterogeneity,"./NMR_insomnia_ever_stratified_UKB_heterogeneity.csv")
write.csv(nmr_insomnia_uvmr_strat_heterogeneity,"./NMR_insomnia_former_stratified_UKB_heterogeneity.csv")
write.csv(nmr_insomnia_uvmr_strat_heterogeneity,"./NMR_insomnia_never_stratified_UKB_heterogeneity.csv")

# C.
nmr_insomnia_uvmr_strat_pleiotropy<-mr_pleiotropy_test(dat)
write.csv(nmr_insomnia_uvmr_strat_pleiotropy,"./NMR_insomnia_current_stratified_UKB_pleiotropy.csv")
write.csv(nmr_insomnia_uvmr_strat_pleiotropy,"./NMR_insomnia_ever_stratified_UKB_pleiotropy.csv")
write.csv(nmr_insomnia_uvmr_strat_pleiotropy,"./NMR_insomnia_former_stratified_UKB_pleiotropy.csv")
write.csv(nmr_insomnia_uvmr_strat_pleiotropy,"./NMR_insomnia_never_stratified_UKB_pleiotropy.csv")

res_single <- mr_singlesnp(dat)
res_single
write.csv(res_single,"./NMR_insomnia_current_stratified_UKB_res_single.csv")
write.csv(res_single,"./NMR_insomnia_ever_stratified_UKB_res_single.csv")
write.csv(res_single,"./NMR_insomnia_former_stratified_UKB_res_single.csv")
write.csv(res_single,"./NMR_insomnia_never_stratified_UKB_res_single.csv")

#######################################################
#5. Visualise the causal effect of NMR on insomnia #
######################################################

# Generate a scatter plot comparing the different methods
png("./NMR_insomnia_current_stratified_UKB_scatter.png")
png("./NMR_insomnia_ever_stratified_UKB_scatter.png")
png("./NMR_insomnia_former_stratified_UKB_scatter.png")
png("./NMR_insomnia_never_stratified_UKB_scatter.png")
mr_scatter_plot(mr_results, dat)
dev.off()

# Generate a forest plot of each of the SNP effects, which are then meta-analysed using the IVW and MR-Egger methods
png("./NMR_insomnia_current_stratified_UKB_forest.png")
png("./NMR_insomnia_ever_stratified_UKB_forest.png")
png("./NMR_insomnia_former_stratified_UKB_forest.png")
png("./NMR_insomnia_never_stratified_UKB_forest.png")
mr_forest_plot(res_single)
dev.off()

# Generate a funnel plot to check asymmetry
png("./NMR_insomnia_current_stratified_UKB_funnel.png")
png("./NMR_insomnia_ever_stratified_UKB_funnel.png")
png("./NMR_insomnia_former_stratified_UKB_funnel.png")
png("./NMR_insomnia_never_stratified_UKB_funnel.png")
mr_funnel_plot(res_single)
dev.off()

# Run a leave-one-out analysis and generate a plot to test whether any one SNP is driving any pleiotropy or asymmetry in the estimates
res_loo <- mr_leaveoneout(dat)
png("./NMR_insomnia_current_stratified_UKB_loo.png")
png("./NMR_insomnia_ever_stratified_UKB_loo.png")
png("./NMR_insomnia_former_stratified_UKB_loo.png")
png("./NMR_insomnia_never_stratified_UKB_loo.png")
mr_leaveoneout_plot(res_loo)
dev.off()

#Tidy up
rm(nmr_insomnia_uvmr_strat_heterogeneity,nmr_insomnia_uvmr_strat_pleiotropy,dat,mr_results,insomnia_current,res_loo,res_single)
rm(nmr_insomnia_uvmr_strat_heterogeneity,nmr_insomnia_uvmr_strat_pleiotropy,dat,mr_results,insomnia_ever,res_loo,res_single)
rm(nmr_insomnia_uvmr_strat_heterogeneity,nmr_insomnia_uvmr_strat_pleiotropy,dat,mr_results,insomnia_former,res_loo,res_single)
rm(nmr_insomnia_uvmr_strat_heterogeneity,nmr_insomnia_uvmr_strat_pleiotropy,dat,mr_results,insomnia_never,res_loo,res_single)
