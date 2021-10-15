#knee:ebi-a-GCST007090
#hip:ebi-a-GCST007091 
exposure_dat<- extract_instruments("ebi-a-GCST007091")

#OA at any site
setwd("D:/R/a-code/PA/R")
exposure_dat <- extract_instruments("large")
file <- system.file("extdata","OAs.csv", package="TwoSampleMR")
exposure_dat <- read_exposure_data(
  filename = "OAs.csv",
  sep = ",",
  snp_col = "rsID",
  beta_col = "BETA",
  se_col = "SE",
  eaf_col = "WEAF",
  effect_allele_col = "EA",
  other_allele_col = "NEA",
  pval_col = "P",
) 
exposure_dat <- clump_data(exposure_dat)
 
#stroke 
stroke_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006906")
AIS_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006908")
LAS_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST005840")
CES_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST006910")
SVS_outcome_dat <- extract_outcome_data(snps=exposure_dat$SNP, outcomes="ebi-a-GCST005841")
#stroke--- ICH
setwd("D:/R/a-code/PA/R")
file <- system.file("extdata","ICH.csv", package="TwoSampleMR")
outcome_dat <- read_outcome_data(
  filename = "ICH.csv",
  snps=exposure_dat$SNP,
  sep = ",",
  snp_col = "SNP",
  beta_col = "Effect",
  se_col = "FreqSE",
  eaf_col = "Freq1",
  effect_allele_col = "Allele1",
  other_allele_col = "Allele2",
  pval_col = "P-value",
) 


#harmonise dat
harmonise_dat <- harmonise_data(exposure_dat, stroke_outcome_dat, action = 2)
dat <-harmonise_dat

#run MR
res <- mr(dat)
res <- generate_odds_ratios(res)
#run heterogeneity pleiotropy
mr_pleiotropy
mr_heterogeneity(dat)
mr_pleiotropy_test(dat)
#run mr_presso
mr_presso <- mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = dat, NbDistribution = 10000,  SignifThreshold = 0.05)
raps<- mr(dat, method_list = c("mr_raps"))
raps <- generate_odds_ratios(raps)
mr_presso_main<-mr_presso[["Main MR results"]]
names(mr_presso_main) <- c("Exposure","MR Analysis","b","se","T-stat","pval")
mr_presso_main<- generate_odds_ratios(mr_presso_main)

#run scatter plot 
p1 <- mr_scatter_plot(res, dat)
p1[[1]]
length(p1)
ggsave(p1[[1]], file="d://.tiff", width=6, height=6)
res_single <- mr_singlesnp(dat, all_method=c("mr_ivw"))

#run forest plot
p2 <- mr_forest_plot(res_single)
p2[[1]]
ggsave(p2[[1]], file="d://.tiff", width=9, height=9)
#run  leaveoneout
res_loo <- mr_leaveoneout(dat)
p3 <- mr_leaveoneout_plot(res_loo)
p3[[1]]
ggsave(p3[[1]], file="d://.tiff", width=6, height=6)
res_single <- mr_singlesnp(dat)
#run funnel plot
p4 <- mr_funnel_plot(res_single)
p4[[1]]
ggsave(p4[[1]], file="d://.tiff", width=6, height=6)
head(res) 
head(raps)
head(mr_presso_main)

write.csv(res,"D:/R/stroke and OA/OA exposure/allOA/ICH/res.csv")
write.csv(dat,"D:/R/stroke and OA/OA exposure/allOA/small/dat.csv")
write.csv(mr_presso_main,"D:/R/stroke and OA/OA exposure/allOA/small/presso.csv")
write.csv(mr_heterogeneity(dat),"D:/R/stroke and OA/OA exposure/allOA/small/heterogeneity.csv")
write.csv(mr_pleiotropy_test(dat),"D:/R/stroke and OA/OA exposure/allOA/small/pleiotropy.csv")
write.csv(raps,"D:/R/stroke and OA/OA exposure/allOA/small/raps.csv")
