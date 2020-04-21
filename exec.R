load("f2GWAS.RData")
load("long.RData")
require(GenomicSEM)
#run the multivariate GWAS using parallel processing
pfactor <- commonfactorGWAS(covstruc = addlong, SNPs = p_sumstats_disease, estimation = "DWLS", cores = NULL, toler = FALSE, SNPSE = FALSE, parallel = TRUE, Output = NULL)
Effective_N<-(mean(((pfactor$Z_Estimate/pfactor$est)^2)/(2*pfactor$MAF*(1-pfactor$MAF))))
munge(pfactor,"w_hm3.noMHC.snplist",trait.names = "latent2",Effective_N)

ld <- "eur_w_ld_chr/"
wld <- "eur_w_ld_chr/"
traits<-c("systolic.sumstats.gz", "bmi.sumstats.gz", "cholesterol.sumstats.gz", "triglycerides.sumstats.gz",
          "coronary.sumstats.gz", "diabete.sumstats.gz", "kidney.sumstats.gz","latent2.sumstats.gz","ea.sumstats.gz")
sample.prev <- c(NA,NA,NA,NA,0.256,.068,.086,NA,NA)
population.prev <- c(NA,NA,NA,NA,0.064,.073,0.0114,NA,NA)
trait.names<-c("SBP","BMI","LDL","TRI","CAD","T2D","CKD","latent","EA")
ldsc_result<-ldsc(traits, sample.prev, population.prev, ld, wld, trait.names)
CommonFactor_DWLS<- commonfactor(ldsc_result, estimation="DWLS")
CommonFactor_DWLS