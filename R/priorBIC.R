# return BIC (prior BIC)
# Author: tunenori
###############################################################################

priorBIC <- function(x, centers, q, ignore.covar){
	lnL0 <- lnL(x, centers, ignore.covar)
	bic <- -2 * lnL0$lnL + q * log(nrow(x)) # BIC
	# bic <- -2 * lnL0$lnL + q  # AIC
	list(lnL = lnL0$lnL, detVx = lnL0$detVx, bic = bic)
}
