# TODO: Add comment
# function for calculation of 
# t(xmu[i,]) %*% invVx %*% xmu[i,]
# Author: tunenori
###############################################################################
txInvVxX <- function(x, invVx){
	t(x) %*% invVx %*% x
}
