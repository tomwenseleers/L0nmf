## available convolution kernels
cols <- 1:8
k1 <- convKernel("gaussian",    51, list(sigma = 4))
k2 <- convKernel("bigaussian",  51, list(sigma1 = 2, sigma2 = 6))
k3 <- convKernel("emg",         51, list(sigma = 2, tau = 6))
k4 <- convKernel("lorentzian",  51, list(gamma = 3))
k5 <- convKernel("laplace",     51, list(tau = 4))
k6 <- convKernel("exponential", 51, list(tau = 5))
k7 <- convKernel("gamma",       51, list(k = 9, tau = 0.5))
k8 <- convKernel("gammaexp",    51, list(k = 9, tau_g = 0.7, tau_e = 0.7))
matplot(cbind(k1,k2,k3,k4,k5,k6,k7,k8), type = "l", lty = 1, col = cols,
        xlab = "kernel index", ylab = "density",
        main = "Convolution kernels\n(point-spread functions)")
legend("topright",
       legend = c("Gaussian", "Bigaussian", "EMG", "Lorentzian",
                  "Laplace", "Exponential", "Gamma", "Gammaexp"),
       lty = 1, cex = 0.8, y.intersp = 0.8,
       col = cols, bty = "n")

## peak broadening as observed at higher concentrations in GC/MS analysis
cols=1:3
matplot(cbind(
  convKernel("bigaussian", 51,
             list(sigma1 = 2, sigma2 = 4), perc = 0.001),
  convKernel("bigaussian", 51,
             list(sigma1 = 2*2, sigma2 = 4*2), perc = 0.001),
  convKernel("bigaussian", 51,
             list(sigma1 = 2*3, sigma2 = 4*3), perc = 0.001)),
  xlab = "kernel index", ylab = "density",
  main = paste("Bigaussian point spread functions",
               "(sigma1 = 2, 4, 6; sigma2 = 2 * sigma1)",
               "perc = 0.001",
               sep = "\n"),
  col = cols, lty = 1,
  type="l")
