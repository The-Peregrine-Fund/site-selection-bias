

simSiteBias <- function (nsites = 100, nsurveys = 1, nyears = 20,  
                          mean.beta = 0, sd.lam = c(0, 0),  
                          show.plot = TRUE, sd.N = 100) #change sd.N to adjust severity of yearly fluctuations
{
  year <- (1:nyears)-1
  alpha <- rgamma(nsites,shape=.1,rate=0.05) #change rate to adjust distance of outliers
  beta <- rnorm(nsites, mean.beta, sd.lam[2])
  lam <- N <- p <- ep <- array(NA, dim = c(nsites, nyears))
  for (t in (0:nyears)) { 
    lam[, t] <- exp(log(alpha) + beta * year[t])
  }
  
  
  for (t in (1:nyears)) { 
    N[, t] <- rnorm(nsites, lam[, t], sd.N)
  }
  N<-ifelse(N<0,0,N)
  
  totalN <- apply(N, 2, sum)
  if (show.plot) {
    op <- par(mfrow = c(2, 1), mar = c(5, 5, 5, 2), cex.lab = 1.5,
              cex.axis = 1.5, cex.main = 1.5)
    on.exit(par(op))
    ylim <- range(lam, N)
    matplot((1:nyears), t(lam), type = "l", lty = 1,
            lwd = 2, ylim = ylim, frame = FALSE, xlab = "Year",
            ylab = "lambda", main = "Expected abundance")
    matplot((1:nyears), jitter(t(N)), type = "l", lty = 1,
            lwd = 2, ylim = ylim, frame = FALSE, xlab = "Year",
            ylab = "N", main = "Realized abundance")
  }
  return(list(nsites = nsites, nsurveys = nsurveys, nyears = nyears,
              mean.beta = mean.beta, sd.lam = sd.lam,
              alpha = alpha, beta = beta,
              lam = lam, N = N, totalN = totalN, p = p, C = C,
              ep=ep))
  return(list(lam=lam))
}


simSiteBias()

