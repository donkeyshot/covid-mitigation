### Reproduction ratio 
covidR0 <- function(beta1 = 0.25, beta2 = 0.16, beta3 = 0.016,
                    sigma = 1, gamma1 = 0.2, gamma2 = 0.14, gamma3 = 0.14,
                    alpha = 0.5, p = 0.5, b = 0.8, r = 0.01) {
  beta1 / gamma1 + #from I1
    beta2 * (1 - p) / gamma2 + #from I2m
    beta3 * (1 - p) / gamma3 + #from I3m
    beta2 * b * p / (gamma2 + alpha) + #from I2s
    b * beta3 * p * (gamma2 / (gamma2 + alpha)) / gamma3 + #from I3s
    beta2 * r * p * (alpha / (gamma2 + alpha)) / gamma2 + 
    beta3 * r * p * (alpha / (gamma2 + alpha)) / gamma3 
}

### Making the plot: first calculation the R0 values for a range of
# 1/alpha and psevere
zvalues <- outer(seq(0.5,10,0.5), seq(0, 1, 0.05), 
                 function(X, Y) covidR0(alpha = 1/X, p = Y,
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
tiff("Results/fig_appendix3.tiff", units = "cm", width = 8, height = 8, res = 300)
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





