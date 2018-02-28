#
# Module that implements the automatic relevance determination technique for genetic network elucidation on a gene expression matrix.
#
# Author: Jose Miguel Hernandez Lobato and Daniel Hernandez Lobato
# Date: Feb 2009
#
# Reference:
#
# 	Tipping, M. E. and A. C. Faul (2003). Fast marginal likelihood maximisation for sparse Bayesian models.
#	In C. M. Bishop and B. J. Frey (Eds.), Proceedings of the Ninth International Workshop on Artificial Intelligence
#	and Statistics, Key West, FL, Jan 3-6.
#
#	Bishop, Christopher M. Pattern Recognition and Machine Learning. Springer, 2006.
#

source("updateQR.R")
source("cholUpdate.R")

gibbsSampling <- function(X, Y, beta0 = 1, p0 = 0.5, v0 = 1, burnin=1000,  nsamples=11000, smart.init=T, return.samples=F) {

	gc()

	# We compute the first point in the Markov chain and initialize some global vairables
        
	start <- if (smart.init) initializeChain(X, Y, p0, beta0, v0) else numeric(ncol(X))

	GlobalX <<- X
	GlobalY <<- Y

	# We perform some initializations

	samplesBeta <- rep(0, nsamples)
	samplesV <- rep(0, nsamples)
	samplesP <- rep(0, nsamples)
	samplesZ <- matrix(0, nsamples, ncol(X))
	samplesW <- matrix(0, nsamples, ncol(X))
	samplesZ[ 1 , ] <- start
	samplesBeta[ 1 ] <- beta0
	samplesV[ 1 ] <- v0
	samplesP[ 1 ] <- p0
	samplesW[ 1, ] <- sampleW(samplesZ[ 1, ], beta0, 1 / v0)

	# We do the sampling

	for (i in 2 : nsamples) {
            system.time({
		# We sample the weights

		samplesW[ i, ] <- sampleW(samplesZ[ i - 1, ], samplesBeta[ i - 1 ], samplesV[ i - 1 ]^-1)

		# We sample beta

		samplesBeta[ i ] <- sampleBeta(X, Y, samplesW[ i, ], beta0)

		# We sample v
	
		samplesV[ i ] <- sampleV(samplesW[ i, which(samplesZ[ i - 1, ] != 0)], v0)

		# We sample Z

		samplesZ[ i, ] <- sampleZ(samplesZ[ i - 1, ], samplesP[ i -1 ], samplesBeta[ i ], samplesV[ i ]^-1)

		# We sample p

		samplesP[ i ] <- sampleP(samplesZ[ i, ], p0)
            })
            cat("It:",i," #on:",sum(samplesZ[ i, ])," p0:",samplesP[i]," noisePrec:",samplesBeta[i]," slabVar:",samplesV[i],"\n")
	}

	samplesW <- samplesW[ (burnin+1) : nrow(samplesW), ]
	samplesZ <- samplesZ[ (burnin+1) : nrow(samplesZ), ]
	samplesBeta <- samplesBeta[ (burnin+1) : length(samplesBeta) ]
	samplesV <- samplesV[ (burnin+1) : length(samplesV) ]
	samplesP <- samplesP[ (burnin+1) : length(samplesP) ]
	
	results=list(m = apply(samplesW, 2, mean), v = apply(samplesW, 2, sd)^2, p = apply(samplesZ, 2, mean), beta = mean(samplesBeta),
		vSlab = mean(samplesV), p0 = mean(samplesP))
        if (return.samples){
            results$samplesW=samplesW
            results$samplesZ=samplesZ
            results$samplesBeta=samplesBeta
            results$samplesV=samplesV
            results$samplesP=samplesP
        }
        results
    }

sampleP <- function(z, p0) {

	a0 <- p0 * 3
	b0 <- (1 - p0) * 3

	a <- a0 + length(which(z == 1))
	b <- b0 + length(which(z == 0))

	rbeta(1, a, b)
}

##
# Function that samples the variance parameter v.
#
# @param 
# @param w	The current value for the regression coefficients.
#
# @return	A Gibbs sample of v.
#

sampleV <- function(w, v0) {

	a0 <- 0.5 * 3 
	b0 <- 0.5 * 3 * v0

	if (length(w != 0)) {
		a <- 0.5 * length(w) + a0
		b <- 0.5 * sum(w^2) + b0
	} else {
		a <- a0
		b <- b0
	}

	1 / rgamma(1, a, b)
}

##
# Function that samples the noise parameter beta.
#
# @param 
# @param X	n x d design matrix.
# @param y	n-dimensional vector of targets.
# @param w	The current value for the regression coefficients.
#
# @return	A Gibbs sample of beta.
#

sampleBeta <- function(X, y, w, beta0) {

	mPost <- X %*% w

	# Prior hyper-parameters

	a0 <- 0.5 * 3
	b0 <- 0.5 * 3 * beta0^-1
	
	# We compute the parameters of the gamma distribution

	b <- 0.5 * sum((y - mPost)^2) + b0
	a <- 0.5 * length(y) + a0

	# We sample beta

	rgamma(1, a, b)
}

initializeChain <- function(X, Y, p0, noisePrecision, v1) {

	size <- p0 * ncol(X)
	error <- rep(0, ncol(X))

	computeResidual <- function(selected) {
		if (length(selected) == 0) {
			Y
		} else {
			Xprima <- X[ , selected ]
			weights <- solve(diag(rep(v1^-1, length(selected)), length(selected), length(selected)) +
				noisePrecision * t(Xprima) %*% Xprima) %*% (noisePrecision * t(Xprima)) %*% Y
			Y - Xprima %*% weights
		}
	}

	selected <- c()
	residual <- Y
	for (i in 1 : size) {
		currentResidual <- computeResidual(selected)

		correlations <- rep(0, size)
		for (j in 1 : ncol(X)) {
			correlations[ j ] <- abs(cor(currentResidual, X[ , j ]))
		}
	
		newFeature <- which.max(correlations)
		selected <- c(selected, newFeature)

		print(i)
	}

	initialization <- rep(0, ncol(X))

	initialization[ selected ] <- 1

	initialization
}

sampleZ <- function(start, p0, noisePrecision, invV1) {

	GlobalColumns <- which(start == 1)
	GlobalPhi <- matrix(GlobalX[ , GlobalColumns ], nrow = nrow(GlobalX), ncol = length(GlobalColumns))
	invSigma <- invV1 * diag(length(GlobalColumns)) + noisePrecision * t.default(GlobalPhi) %*% GlobalPhi
	if (length(GlobalColumns) != 0) {
		GlobalCholInvSigma <- chol(invSigma)
		GlobalSigma <- chol2inv(GlobalCholInvSigma, LINPACK = T)
	} else {
		GlobalSigma <- NULL
	}
	GlobalSigmaChanged <- FALSE

	sampleZ <- start

	d <- length(sampleZ)

	order <- sample(1 : d)
	for (i in order) {

		phi <- GlobalX[ , i ]

		# We compute Qi and Si

		if (length(GlobalColumns) != 0) {
			aux <- ((t.default(phi) %*% GlobalPhi) %*% GlobalSigma) %*% t.default(GlobalPhi)
			Qi <- as.double((noisePrecision * t.default(phi) - noisePrecision^2 * aux) %*% GlobalY)
			Si <- as.double(noisePrecision * t.default(phi) %*% phi - noisePrecision^2 * aux %*% phi)
		} else {
			Qi <- noisePrecision * sum(GlobalY * phi)
			Si <- noisePrecision * sum(phi^2)
		}
	
		# We compute the probability of activation
	
		if (sampleZ[ i ] == 0) {
			qi <- Qi
			si <- Si
		} else {
			qi <- invV1 * Qi / (invV1 - Si)
			si <- invV1 * Si / (invV1 - Si)
		}
#		if (si < 0)
#                    cat("si: ",si,"\n")
		# We compute the probability of activation

		prob <- p0 / ((1 - p0) / (sqrt(invV1 / (invV1 + si)) * exp(0.5 * qi^2 / (invV1 + si))) + p0)

		# We sample

		if (runif(1) <= prob) {
			if (sampleZ[ i ] == 0) {
				sampleZ[ i ] <- 1
				GlobalColumns <- c(GlobalColumns, i)

				if (length(GlobalColumns) == 1)
					GlobalCholInvSigma <- chol(invV1 * diag(length(GlobalColumns)) + noisePrecision * t.default(phi) %*% phi)
				else
					GlobalCholInvSigma <-
						addColumn(GlobalCholInvSigma, sqrt(noisePrecision) * GlobalPhi, sqrt(noisePrecision) * phi, invV1)

				GlobalSigmaChanged <- TRUE
			}
		} else {
			if (sampleZ[ i ] == 1) {
				sampleZ[ i ] <- 0
				column <- which(GlobalColumns == i)
				GlobalColumns <- GlobalColumns[ -column ]

				if (length(GlobalColumns) == 0)
					GlobalCholInvSigma <- matrix(0,0,0)
				else
					GlobalCholInvSigma <- delColumn(GlobalCholInvSigma, column)

				GlobalSigmaChanged <- TRUE
			}
		}

		if (GlobalSigmaChanged) {
			GlobalPhi <- GlobalX[ , GlobalColumns ]
			if (length(GlobalColumns) == 0)
				GlobalSigma <- NULL
			else
				GlobalSigma <- chol2inv(GlobalCholInvSigma, LINPACK = T)
			GlobalSigmaChanged <- FALSE
		}
	}

	sampleZ
}

sampleW <- function(sampleZ, noisePrecision, invV1) {

	newSample <- rep(0, length(sampleZ))
	index <- which(sampleZ == 1)
	if (length(index) != 0) {

		GlobalColumns <- which(sampleZ == 1)
		GlobalPhi <- matrix(GlobalX[ , GlobalColumns ], nrow = nrow(GlobalX), ncol = length(GlobalColumns))
		invSigma <- invV1 * diag(length(GlobalColumns)) + noisePrecision * t.default(GlobalPhi) %*% GlobalPhi
		GlobalCholInvSigma <- chol(invSigma)
		GlobalSigma <- chol2inv(GlobalCholInvSigma, LINPACK = T)
		GlobalSigmaChanged <- FALSE

		indexColumns <- sort(GlobalColumns, index.return = T)$ix
		mean <- noisePrecision * (GlobalSigma %*% (t(GlobalPhi) %*% GlobalY))
		mean <- mean[ indexColumns ]
		newSample[ index ] <- mean + t(rnorm(length(mean))) %*% chol(GlobalSigma[ indexColumns, indexColumns ])
	}

	newSample
}
