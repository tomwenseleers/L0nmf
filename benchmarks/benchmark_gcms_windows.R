# Run from the L0nmf repository root.
library(L0glm)
library(L0nmf)
results <- list()
for (name in c("GCMS_data1","GCMS_data2")) {
  data(list=name, package="L0nmf")
  Y <- get(name)
  time <- as.numeric(rownames(Y))
  X <- convmatrix(convkernel("gaussian",65L,
    list(sigma=if(name=="GCMS_data1")8 else 6)),nrow(Y))
  B <- if(name=="GCMS_data1") matrix(1,nrow(Y),1) else
    splines::bs(time,df=4L,intercept=TRUE)
  seconds <- system.time(fit <- splogconNMF(X,Y,background=B,
    family=quasipoisson("identity"),precluster=if(name=="GCMS_data1")16L else 32L,
    correlation="none"))[["elapsed"]]
  results[[name]] <- data.frame(dataset=name,n=nrow(Y),channels=ncol(Y),
    compounds=fit$rank,GIC=GIC(fit),seconds=seconds)
}
do.call(rbind,results)
