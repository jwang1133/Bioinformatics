###########################################################################################################
library(tidyverse)
library(broom)

##Read into the data
Data <- read_tsv("/XXX/file_name")

##Modeling the chrmosome index and standardized chromosome size with cubic distribution 
x <- Data2$Chr_index
y <- Data2$Standard_chr_len
x2 <- x * x
x3 <- x2 * x
##K_obs store the parameter estimates of cubic function
K_obs <- coefficients(nls(y ~ b0+ b1*x + b2*x2 + b3*x3, start = list(b0 = 0,b1 = 1, b2 = 1, b3 = 1)))

##Modeling the chrmosome index and standardized chromosome size with gamma distribution
##Function to estimate alpha, nls() is to use the nonlinear least square to estimate the shape parameter that best explain the data
IRWLS <- function (x, y, chro_num) {
  ReWeights <- rep (1, length(y)); 
  alpha_0 <- 6;
  repeat {
    NLS <- nls(y ~ qgamma(x, alpha, 1) / alpha, start = list(alpha = 6), weights = ReWeights);
    alpha <- summary(NLS)$coef[1];
    ReWeights = 1/(x*(1-x)*(gamma(alpha))^2/(chro_num*alpha^2*(qgamma(x, alpha, 1))^(2*alpha-2)*exp(-2*qgamma(x,alpha,1))));
    dif = abs(alpha - alpha_0);
    alpha_0 <- alpha;
    if(dif <= 1e-6) break;
  }
  alpha <- summary(NLS)$coef[1]
}

##Apply the function to estimate the alpha parameter of the gamma distribution from the data
alpha <- IRWLS(Data$Chr_index, Data$Standard_chr_len, Data$Total_chr_nums)

##plot a new figure for the distribution of chromosome index and standardized chromosome size
tiff("/XXX/Fig2.tiff", width = 8, height = 8.5, units = 'in', res = 400)

y_lims <- c(min(Data$Standard_chr_len), max(Data$Standard_chr_len))
x_lims <- c(0, 1)

plot( 0, 0,  xlab = "", ylab = "", ylim = y_lims, xlim = x_lims,  col = "white",  xaxt = "n", yaxt="n"); ###generate a big frame

xl <- seq(from = x_lims[1], to= x_lims[2], by= round((x_lims[2] - x_lims[1])/5, digits = 1))
axis(1, mgp=c(3, .1, 0), at=xl, labels=xl, cex.axis=1.2, tck=.01)
mtext("Chromosome index",side=1, line=1.4, cex=1.5,las=0)

yl <- seq(from = 0, to= 3.5, by= round((y_lims[2] - y_lims[1])/7, digits = 1))
axis(2, mgp=c(3, .1, 0), at=yl, labels=yl, cex.axis=1.2,tck=.01,las=0 )
mtext("Standardized chromosome size", side=2, line=1.4, cex=1.5,las=0)
points(Data2$Chr_index, Data2$Standard_chr_len, cex = 0.55, col = "dimgrey", pch = 19)

##Add the fitted line of cubic distribution

y_cubic <- K_obs[1] + K_obs[2] * x + K_obs[3] * x2 + K_obs[4] * x3
lines(x[order(x)], y_cubic[order(x)], col = "blue", lty = 2.5, lwd = 3 )

##Add the fitted line of gamma distribution
x <- Data$Chr_index
lines(x[order(x)], (qgamma(x,alpha,1)/alpha)[order(x)],  lwd = 2.5, col = "red")

##Add legend
legend(-0.001, 2.8, legend = c("Inverse of Gamma cdf (5.41, 1/5.41)", "Cubic function (0.26, 2.73, -4.47, 3.45)" ), bty = "n",
       lwd = 2.5, cex = 1.3, col = c("red", "blue"), lty = c(1, 2), horiz=F)
box(lwd = 2)

dev.off()


