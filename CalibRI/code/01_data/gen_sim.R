# Simulated datasets from the RIbench parameter grid (Ammer et al., Clin Chem 2022)
spec <- read.csv(system.file("extdata", "SpecificationTestSets.csv", package = "RIbench")); spec <- spec[spec$N==1000,]
invBC <- function(x,l) if(abs(l)<1e-20) exp(x) else (l*x+1)^(1/l)
ANALYTES <- c("TSH","GGT","FT4","Hb"); NS <- c(300,500,1000,2000); FR <- c(0.1,0.3); REPS <- 200
set.seed(20261016)
meta <- list(); vals <- list(); id <- 0
for(a in ANALYTES) for(f in FR){
  rows <- spec[spec$Analyte==a & abs(spec$fractionPathol-f)<1e-9,]
  for(N in NS) for(r in 1:REPS){
    tmp <- rows[(r-1) %% nrow(rows) + 1,]
    phys <- invBC(rnorm(N*3,tmp$nonp_mu,tmp$nonp_sigma),tmp$nonp_lambda)+tmp$nonp_shift
    phys <- phys[is.finite(phys) & phys>=0]
    nph <- round(N*(1-f-tmp$bg_fraction)); phys <- sample(phys,nph)
    nl <- round(N*f*tmp$left_ratio/tmp$ratio_sum); nr <- round(N*f*tmp$right_ratio/tmp$ratio_sum)
    pl <- rnorm(N*3,tmp$left_mu,tmp$left_sigma); pl <- sample(pl[pl>=0],nl,replace=TRUE)
    pr <- rnorm(N*3,tmp$right_mu,tmp$right_sigma); pr <- sample(pr[pr>=0],nr,replace=TRUE)
    bg <- runif(max(1,round(N*tmp$bg_fraction)),tmp$bg_min,tmp$bg_max)
    x <- round(sample(c(phys,pl,pr,bg)),tmp$decimals)
    tLL <- invBC(tmp$nonp_mu+qnorm(.025)*tmp$nonp_sigma,tmp$nonp_lambda)+tmp$nonp_shift
    tUL <- invBC(tmp$nonp_mu+qnorm(.975)*tmp$nonp_sigma,tmp$nonp_lambda)+tmp$nonp_shift
    id <- id+1
    meta[[id]] <- data.frame(ds_id=id,source="sim",analyte=a,n=length(x),n_nominal=N,frac=f,rep=r,
                             spec_index=tmp$Index,trueLL=tLL,trueUL=tUL,decimals=tmp$decimals)
    vals[[id]] <- data.frame(ds_id=id,value=x)
  }
}
m <- do.call(rbind,meta); v <- do.call(rbind,vals)
dir.create("work/ribuild", showWarnings=FALSE, recursive=TRUE)
write.csv(m,"work/ribuild/sim_meta.csv",row.names=FALSE); write.csv(v,"work/ribuild/sim_values.csv",row.names=FALSE)
cat("sim datasets:",nrow(m),"values:",nrow(v),"\n"); print(aggregate(cbind(trueLL,trueUL)~analyte,m,unique))
