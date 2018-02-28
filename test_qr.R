N=5
P=3
x=matrix(rnorm(N*P),N,P)
lambdag=runif(P)
lambday=1.0
s=lambday * t(x) %*% x + diag(lambdag)
sp=s[1:(P-1),1:(P-1)]
up=chol(sp)
print(sp- (t(up)%*%up) )
source("updateQR.R")
u=addColumn(up,sqrt(lambday) * x[,1:(P-1)], sqrt(lambday) * x[,P], lambdag[P])
print(sum(abs(s-t(u)%*%u)))
u=chol(s)
print(sum(abs(s-t(u)%*%u)))

to.del=2
up=delColumn(u,2)
sp=s[c(1,3),c(1,3)]
print(sp- (t(up)%*%up) )
up=chol(sp)
print(sp- (t(up)%*%up) )
