library(L0glm)
data(gcmsData1)
plotMat(gcmsData1,
      xlab = "Elution time (min)",
      ylab = "Mass (m/z)",
      main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
matplot(x = as.numeric(rownames(gcmsData1)),
        y = gcmsData1^0.2, type = "l", lty = 1, col = rgb(0, 0, 0, 0.1),
        xlab = "Elution time (min)",
        ylab = "Intensity ^ 0.2",
        main = "Section of Cataglyphis iberica\nant queen cuticular GC/MS profile")
plot(x=as.numeric(names(attr(gcmsData1, "alkanes_TIC"))),
     y=attr(gcmsData1, "alkanes_TIC"),
     xlab="Retention time (min)", 
     ylab="Total ion intensity", 
     main="C9 to C40 alkane ladder standard", 
     type="l")
