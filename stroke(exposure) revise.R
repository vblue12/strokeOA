library(usethis)
library(devtools)
library(TwoSampleMR)
library(ieugwasr)
library(MRPRESSO)
library(ggplot2)
rm(list=ls())
 
#stroke 
stroke_exposure_dat<- extract_instruments("ebi-a-GCST006906")
IS_exposure_dat<- extract_instruments("ebi-a-GCST006908")
LAS_exposure_dat<- extract_instruments("ebi-a-GCST005840")
CES_exposure_dat<- extract_instruments("ebi-a-GCST006910")
SVS_exposure_dat<- extract_instruments("ebi-a-GCST005841")

exposure_dat<- extract_instruments("ebi-a-GCST007090")

# discover confounder
install.packages("phenoscanner")
H<-SVS_exposure_dat[,'SNP']
library(phenoscanner)
phe2<- phenoscanner(snpquery=H)
head(res$results)
res$snps

#OA knee:ebi-a-GCST007090
#OA hip:ebi-a-GCST007091
#OA hospital:ebi-a-GCST005814
outcome_dat <- extract_outcome_data(snps=SVS_exposure_dat$SNP, outcomes="ebi-a-GCST007090")
outcome_dat <- extract_outcome_data(snps=SVS_exposure_dat$SNP, outcomes="ebi-a-GCST007091")
outcome_dat <- extract_outcome_data(snps=SVS_exposure_dat$SNP, outcomes="ebi-a-GCST005814")

#harmonise dat
harmonise_dat <- harmonise_data(SVS_exposure_dat, outcome_dat, action = 2)
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

#print output
write.csv(res,"D:/R/stroke and OA/Stroke ??exposure??/small/hospital/res.csv")
write.csv(dat,"D:/R/stroke and OA/Stroke ??exposure??/small/hospital/dat.csv")
write.csv(mr_presso_main,"D:/R/stroke and OA/Stroke ??exposure??/small/hospital/presso.csv")
write.csv(mr_heterogeneity(dat),"D:/R/stroke and OA/Stroke ??exposure??/small/hospital/heterogeneity.csv")
write.csv(mr_pleiotropy_test(dat),"D:/R/stroke and OA/Stroke ??exposure??/small/hospital/pleiotropy.csv")
write.csv(raps,"D:/R/stroke and OA/Stroke ??exposure??/small/hospital/raps.csv")
