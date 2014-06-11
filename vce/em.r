y <- c(2.6,0.1,1.0,3.0,1.0)
n <- length(y)
s <- c(0,0,0,1,3,1,4,3)
d <- c(0,0,0,0,2,2,5,6)
Ainv <- quass(s,d)
# initial variance components
initE <- 0.4
initU <- 0.2
#alpha = initE/initU
# X matrix
X <- matrix(c(1,0,0,1,1,0,1,1,0,0), nrow=5, ncol=2)
# Z matrix
z1 <- matrix(0, ncol=3, nrow=5)
z2 <- matrix(c(1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1,0,0,0,0,0,1), ncol=5)
Z <- cbind(z1,z2)


## Given the components of MME, estimate variance components by EM Algorithm. 

## Arguments
## Ainv: an Inverse of additive relationship matrix
## y: a vector of response
## X: a incidence matrix for systematic effects
## Z: a incidence matrix for individual effects
## varE: initial value for the residual variance
## varA: initial value for the additive genetic variance

## Note 
## Unknown parents should be coded as zero.

## Suzuki, M. 2007. Applied Animal Breeding & Genetics. Course Notes. Obihiro University of Agriculture and Veterinary Medicine.

## Author: Gota Morota <morota at wisc dot edu>
## Create: 6-Apr-2010
## Last-Modified: 8-Apr-2010
## License: GPLv2 or later






# MME
# total sum of response per sex
Xpy <- t(X)%*%y
# reponse for each animal
Zpy <- t(Z)%*%y
XpX <- t(X)%*% X
XpZ <- t(X)%*%Z
ZpX <- t(Z)%*%X
ZpZ <- t(Z)%*%Z

#LHS <- rbind( cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*alpha))
RHS <- c(Xpy, Zpy)
oldE <- initE
oldU <- initU
rankX <- qr(X,LAPACK=TRUE)$rank
#rankZ <- qr(Z,LAPACK=TRUE)$rank
rankZ <- dim(Z)[2]
lhsRow <- length(RHS)
z <- length(RHS) - dim(ZpZ)[1] + 1

diff1 <- 1
diff2 <- 1
#B <- solve(LHS)%*%RHS
i <- 0
while (diff1 > 10E-6 & diff2 > 10E-6){
	i <- i+1
	print(i)
	alpha <- as.vector((oldE/oldU))
	LHS <- rbind( cbind(XpX, XpZ), cbind(ZpX, ZpZ+Ainv*alpha ) )
	B <- solve(LHS)%*%RHS
	e <-  y - cbind(X,Z)%*%B
	sigE <- (t(e)%*%y)/(n-rankX)
	c22 <- solve(LHS)[z:lhsRow, z:lhsRow]
	u <- B[z:length(B)]
	# sum(Ainv*c22) is same as sum(diag(Ainv%*%c22))
	sigU <- (t(u)%*%Ainv%*%u + sum(Ainv*c22)*oldE )/(rankZ)
	diff1 <- abs(sigE-oldE)
	diff2 <- abs(sigU-oldU)
	cat("oldE", oldE, '\n')
	cat("sigE", sigE, '\n')
	cat("oldU", oldU, '\n')
	cat("sigU", sigU, '\n')
	oldE <- sigE
	oldU <- sigU
}

sigE
sigU


