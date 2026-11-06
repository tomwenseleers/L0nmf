library(L0glm)
data(simulatedSpectra)
sim <- simulateGCMS(n = 10000,
                     k = 500,
                     S = simulatedSpectra,
                     peak_width_range = c(3, 5),
                     mean_beta  = 1e4,
                     sd_logbeta = 1,
                     width_height_coef = 3,
                     rt_height_coef    = 0.02,
                     sample_nr = 1,
                     rt_drift_coef = 0.001,
                     bg_scale = 100,
                     family = "poisson")
plotMat(sim$D,
        xlab = "Time (min)", 
        ylab = "Mass (m/z)", 
        main = "Simulated GC/MS chromatogram")
