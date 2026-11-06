library(L0glm)
data(gcmsData2)
plotMat(gcmsData2,
      xlab = "Elution time (min)",
      ylab = "Mass (m/z)",
      main = "Section of beer headspace GC/MS profile")
matplot(x = as.numeric(rownames(gcmsData2)),
        y = gcmsData2^0.2, type = "l", lty = 1, col = rgb(0, 0, 0, 0.1),
        xlab = "Elution time (min)",
        ylab = "Intensity ^ 0.2",
        main = "Section of beer headspace GC/MS profile")
plot(x=as.numeric(names(attr(gcmsData2, "alkanes_TIC"))),
     y=attr(gcmsData2, "alkanes_TIC"),
     xlab="Retention time (min)", 
     ylab="Total ion intensity", 
     main="C7 to C22 alkane ladder standard", 
     type="l")
