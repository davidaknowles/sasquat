source("gibbsSampling.R")
N=100
P=1000
X = matrix(rnorm(N*P),N,P)
w = numeric(P)
w[1:10]=rnorm(10)
Y = X %*% w + rnorm(N)
p0=5.0/P
# p0 is the prior prob of being on
# noisePrec is the noise precision (inverse variance)
# slabVar is the variance of the gaussian coeffs
results=gibbsSampling(X,Y,p0=5.0/P,burnin=100,nsamples=200,smart.init=T,return.samples=T)
print(results$p[1:20])

require(BoomSpikeSlab)
bss=lm.spike(Y~X,niter=200)

bssm=colMeans(bss$beta[100:200,])

# error on first 100 coefs
myerr=mean(abs(w[1:100]-results$m[1:100]))
boomerr=mean(abs(w[1:100]-bssm[1:100]))

# logit test
y=runif(N) < 1/(1+exp(-X%*%w))
bssl=logit.spike(y~X,niter = 200)

bssm=colMeans(bssl$beta[100:200,])
boomerr=mean(abs(w[1:100]-bssm[1:100]))

require(glmnet)
cv=cv.glmnet(X[,1:20],y,family="binomial")
opt=glmnet(X[,1:20],y,family="binomial",lambda=cv$lambda.min)
coef(opt)[1:20]
