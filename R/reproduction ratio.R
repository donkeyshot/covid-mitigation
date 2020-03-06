### Reproduction ratio for original model. In the publication, 
# there is no distinction anymore between hospitalisation and 
# self-quarantine. To keep that, make sure that reduction_self
# and reduction_hospital are always equal
covidR0 <- function(beta1 = 0.25, beta2 = 0.16, beta3 = 0.016,
                    sigma = 1, gamma1 = 0.2, gamma2 = 0.14, gamma3 = 0.14,
                    alpha = 0.5,
                    psevere = 0.5, phospital = 0.2,
                    reduction_severe = 0.8, reduction_self = 0.01, reduction_hosp = 0.01) {
  (beta1 / gamma1 + #from I1
     beta2 * (1 - psevere) / gamma2 + #from I2m
     beta3 * (1 - psevere) / gamma3 + #from I3m
     beta2 * reduction_severe * psevere / (gamma2 + alpha) + #from I2s
     reduction_severe * beta3 * psevere * (gamma2 / (gamma2 + alpha)) / gamma3 + #from I3s
     beta2 * reduction_hosp * psevere * (alpha / (gamma2 + alpha)) * phospital / gamma2 + 
     beta3 * reduction_hosp * psevere * (alpha / (gamma2 + alpha)) * phospital / gamma3 +
     beta2 * reduction_self * psevere * (alpha / (gamma2 + alpha)) * (1 - phospital) / gamma2 + 
     beta3 * reduction_self * psevere * (alpha / (gamma2 + alpha)) * (1 - phospital) / gamma3) 
}

### Making the plot: first calculation the R0 values for a range of
# 1/alpha and psevere
zvalues <- outer(seq(0.5,10,0.5), seq(0, 1, 0.05), 
                 function(X, Y) covidR0(alpha = 1/X, psevere = Y,
                                        beta1 = 0.25, beta2 = 0.16, beta3 = 0.016, gamma1 = 0.2))

### Then making the plot
nrz <- nrow(zvalues)
ncz <- ncol(zvalues)
surface.colors <- colorRampPalette( c("cyan", "red") )
# Generate the desired number of colors from this palette
nbcol <- 100
color <- surface.colors(nbcol)
zfacet <- zvalues[-1, -1] + zvalues[-1, -ncz] + zvalues[-nrz, -1] + zvalues[-nrz, -ncz]
facetcol <- cut(zfacet, nbcol)
par(cex = 1.2)
persp(x = seq(.5,10,0.5),
      y=seq(0,1,.05),z = zvalues, theta = 210, phi = 30,
      xlab = "Isolation delay (1/alpha)",
      ylab = "Proportion isolated", 
      zlab = "Reproductive number",
      xlim = c(0,10),
      zlim = c(1, max(zvalues)),
      ticktype = "detailed", col = color[facetcol],
      r = 10, d = 10)
tiff("Results/fig2c.tiff", units = "cm", width = 8, height = 8, res = 300)
par(cex = 1.2, mar = .5+c(0,0,0,0))
persp(x = seq(.5,10,0.5),
      y=seq(0,1,.05),z = zvalues, theta = 210, phi = 30,
      xlab = NA,
      ylab = NA, 
      zlab = NA,
      xlim = c(0,10),
      zlim = c(1, max(zvalues)),
      ticktype = "detailed", col = color[facetcol],
      r = 10, d = 10)
dev.off()





