## -----------------------------------------------------------------------------
X<-c(22.44,14.48,20.73,19.25,20.37,26.43,12.14,23.31,16.23,0.56,0.84,18.05,12.45,11.33)
Y<-c(2.40,2.98,2.06,1.09,1.96,1.55,2.16,1.60,0.80,1.94,3.00,0.28,0.84,1.80)
data<-data.frame(X,Y)
plot(X,Y,pch=23,col="red",bg="green")
LM<-lm(Y~X,data)
abline(LM)

## -----------------------------------------------------------------------------
generate.pareto<-function(n,a,b){
  u=runif(n)
  P=b*(1-u)^(-1/a)
  return (P)
}

set.seed(1000)
a=2
b=2
x=generate.pareto(1e4,a,b)

hist(x,prob=TRUE,breaks=50,main="Histogram of Pareto(2,2)")
y=sort(x)
fy=a*b^a*y^(-(a+1))
lines(y,fy,col="green")

## -----------------------------------------------------------------------------
generate.fe<-function(n){
  U1=runif(n,-1,1)
  U2=runif(n,-1,1)
  U3=runif(n,-1,1)
  U=ifelse((abs(U3)>=abs(U2) & abs(U3)>=abs(U1)),U2,U3)
  return(U)
}

set.seed(1000)
U=generate.fe(1e4)
hist(U,breaks=20,prob=TRUE,main=expression(paste("Histogram of ", f[e])))
lines(density(U),col="green")


## -----------------------------------------------------------------------------
generate.alt.pareto<-function(n,beta,r){
  U=runif(n)
  P=beta*(1/((1-U)^(1/r))-1)
  return(P)
}

set.seed(1000)
beta=2;r=4
x=generate.alt.pareto(1e4,beta,r)
hist(x,breaks=50,prob=TRUE,main="Histogram of Pareto(2,4)")
y=sort(x)
fy=r*beta^r*(beta+y)^(-(r+1))
lines(y,fy,col="green")

## ----warning=FALSE,message=FALSE----------------------------------------------
library(knitr)
UP = pi/3
U = runif(10000, 0, UP)
G = sin(U)
CONTROL = U - U^3/6
ECONTROL = pi/6 - pi^3/648
L = lm(G ~ CONTROL)
cstar = L$coefficients[2]
G_adjust = G - cstar*(CONTROL - ECONTROL)
result_adj = mean(G_adjust)*UP
result_0 = mean(G)*UP
var_0 = var(G)*UP^2
var_adjust = var(G_adjust)*UP^2
percentage_empirical = 1 - var_adjust/var_0
percentage_theorical = summary(L)$r.squared
result = as.matrix(c(result_0,result_adj))
result = cbind(result,c(sqrt(var_0),sqrt(var_adjust)))
result = cbind(result, c(0.5,0.5))
rownames(result) <- c('unadjust','adjust')
colnames(result) = c('estimate','standard error', 'True')
kable(result)
100*c(percentage_empirical,percentage_theorical) ->tb
names(tb) <-c('empirical', 'theorical')
tb = as.matrix(tb)
colnames(tb) <- c('percentage of variance reduced(control variate)')
kable(tb)

## ----warning=FALSE, message=FALSE---------------------------------------------
X = runif(10000)
Y = 1-X
G1 = exp(X)
G2 = exp(Y)
G = (G1+G2)/2
X2 = runif(10000)
G3 = exp(X2)
result_original = c(mean(c(G1,G3)), sd(c(G1,G3)))
result_anti = c(mean(G), sd(G))    
result =matrix(c(result_anti,result_original),2,2,byrow = T)
colnames(result) <- c('estimate', 'standard error')
rownames(result) <- c('anti', 'normal')
percentage_reduced = 1-var(G)/var(c(G1,G3))
kable(result)
print(percentage_reduced)

## -----------------------------------------------------------------------------
m<-10000
theta.hat<-se<-numeric(2)
g<-function(x){
  x^2*exp(-x^2/2)/sqrt(2*pi)*(x>1)
}

#f1
x<-rnorm(m)
fg<-g(x)/dnorm(x)
theta.hat[1]<-mean(fg)
se[1]<-sd(fg)

#f2
u<-runif(m)
x<-1/(1-u)
fg<-g(x)/x^2
theta.hat[2]<-mean(fg)
se[2]<-sd(fg)

rbind(theta.hat,se)

## -----------------------------------------------------------------------------
m<-1000
theta.hat<-se<-numeric(5)
g<-function(x){
  exp(-x-log(1+x^2))
}

for (i in 1:5) {
  u<-runif(m)
  x<--log(exp(-(i-1)/5)-(exp(-(i-1)/5)-exp(-i/5))*u)
  fg<-g(x)*(x>((i-1)/5))*(x<(i/5))/(exp(-x)/(exp(-(i-1)/5)-exp(-i/5)))
  theta.hat[i]<-mean(fg)
  se[i]<-sd(fg)
}
theta<-sum(theta.hat)
se<-sum(se)/5
rbind(theta,se)

## -----------------------------------------------------------------------------
calcCI <- function(n, alpha) {
  u <- rnorm(n, mean = 0, sd = 2)
  y<-exp(u)
  return(mean(y)-sd(y) * qt(alpha, df = n-1)/sqrt(n))
}
UCL <- replicate(1000, expr = calcCI(n = 20, alpha = .05)) 

## -----------------------------------------------------------------------------
n <- 20
alpha <- .05
UCL <- replicate(1000, expr = {
  x <- rnorm(n, mean = 0, sd = 2)
  y<-exp(x)
  mean(y)-sd(y) * qt(alpha, df = n-1)/sqrt(n)
})
#count the number of intervals that contain mu=0
sum(log(UCL)> 1)
#or compute the mean to get the confidence level
mean(log(UCL) > 1) 

## -----------------------------------------------------------------------------
set.seed(1020)
generate <- function(n, alpha){
  X <- rchisq(n,2)
  B <- abs(sqrt(var(X)/n)*qt(1-alpha/2, df=n-1))
  return(c(mean(X)-B,mean(X)+B))
}
RESULT <- replicate(1000, expr = generate(20,0.05))
REALmean <- 2
temp1 <- RESULT[1,]<2
temp2 <- RESULT[2,]>2
Proportion <- sum(temp1&temp2)/1000 * 100
# the proportion that confidence interval contains the real mean
Proportion

## ----warning=FALSE, message=FALSE---------------------------------------------
library(ggplot2)
data_0 = as.data.frame(cbind(t(RESULT),1:1000))
data_0$fac = as.factor(temp1&temp2)
ggplot(data_0) + geom_linerange(aes(x=V3, ymin=V1,ymax = V2,color=fac),size = 1) + scale_color_manual(values=c("#999999", "#E69F00"))

## -----------------------------------------------------------------------------
alpha <- 20
n <- c(10,20,30,50,100,500) #sample sizes
p.reject <- p.heavy <- cv <- xbar <-ybar <- m31 <- m32 <- 
  m21 <- m22 <- u <- v <- numeric(length(n)) 
#to store sim. results 
m <- 10000 #num. repl. each sim.
sktests <- heavy <- numeric(m) #test decisions 
set.seed(12345)
for (i in 1:length(n)) {
  for (j in 1:m) {
    cv[i] <- qnorm(.975, 0,sqrt(6*(n[i]-2) / ((n[i]+1)*(n[i]+3))))  
    #crit. values for each n
    x <- rbeta(n[i],alpha,alpha)
    y <- rt(n[i],2)
    xbar[i] <- mean(x)
    ybar[i] <- mean(y)
    m31[i] <- mean((x - xbar[i])^3)
    m32[i] <- mean((y - ybar[i])^3)
    m21[i] <- mean((x - xbar[i])^2)
    m22[i] <- mean((y - ybar[i])^2)
    u[i] <- m31[i] / ((m21[i])^1.5)
    v[i] <- m32[i] / ((m22[i])^1.5)
    sktests[j] <- as.integer(abs(u[i])>= cv[i] ) 
    heavy[j] <- as.integer(abs(v[i])>= cv[i] ) 
  }
  p.reject[i] <- mean(sktests) #proportion rejected 
  p.heavy[i] <- mean(heavy)
}
comparison <- data.frame(n,p.reject, p.heavy)
knitr::kable(comparison)

## -----------------------------------------------------------------------------
alpha <- 20
n <- c(1000,1500,2000) #sample sizes
p.reject <- p.heavy <- cv <- xbar <-ybar <- m31 <- m32 <- 
  m21 <- m22 <- u <- v <- numeric(length(n)) 
#to store sim. results 
m <- 10000 #num. repl. each sim.
sktests <- heavy <- numeric(m) #test decisions 
set.seed(1234)
for (i in 1:length(n)) {
  for (j in 1:m) {
    cv[i] <- qnorm(.975, 0, sqrt(6/n[i]))   
    #crit. values for each n
    x <- rbeta(n[i],alpha,alpha)
    y <- rt(n[i],2)
    xbar[i] <- mean(x)
    ybar[i] <- mean(y)
    m31[i] <- mean((x - xbar[i])^3)
    m32[i] <- mean((y - ybar[i])^3)
    m21[i] <- mean((x - xbar[i])^2)
    m22[i] <- mean((y - ybar[i])^2)
    u[i] <- m31[i] / ((m21[i])^1.5)
    v[i] <- m32[i] / ((m22[i])^1.5)
    sktests[j] <- as.integer(abs(u[i])>= cv[i] ) 
    heavy[j] <- as.integer(abs(v[i])>= cv[i] ) 
  }
  p.reject[i] <- mean(sktests) #proportion rejected 
  p.heavy[i] <- mean(heavy)
}
comparison <- data.frame(n,p.reject, p.heavy)
knitr::kable(comparison)

## -----------------------------------------------------------------------------
alpha1 <- 4 
alpha2 <- 10
n <- 300
m <- 1500
set.seed(1234)
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05)) 
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(0.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon 
  e <- epsilon[j]
sktests <- xbar <- m3 <- m2 <- sk <- numeric(m)
for (i in 1:m) { #for each replicate
alpha <- sample(c(alpha1, alpha2), replace = TRUE, 
                size = n, prob = c(1-e, e))
x <- rbeta(n, alpha, alpha)
xbar[i] <- mean(x)
m3[i] <- mean((x-xbar[i])^3)
m2[i] <- mean((x-xbar[i])^2)
sk[i] <- m3[i] / ((m2[i])^1.5)
sktests[i] <- as.integer(abs(sk[i]) >= cv) }
        pwr[j] <- mean(sktests)
}
#plot power vs epsilon 
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon)) 
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors 
lines(epsilon, pwr+se, lty = 3,col = "red")
lines(epsilon, pwr-se, lty = 3,col = "red")

## -----------------------------------------------------------------------------
n1 <- 4 
n2 <- 40
n <- 300
m <- 1500
set.seed(1234)
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05)) 
N <- length(epsilon)
pwr <- numeric(N)
#critical value for the skewness test
cv <- qnorm(0.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) { #for each epsilon 
  e <- epsilon[j]
sktests <- xbar <- m3 <- m2 <- sk <- numeric(m)
for (i in 1:m) { #for each replicate
nn <- sample(c(n1, n2), replace = TRUE, 
                size = n, prob = c(1-e, e))
x <- rt(n, nn)
xbar[i] <- mean(x)
m3[i] <- mean((x-xbar[i])^3)
m2[i] <- mean((x-xbar[i])^2)
sk[i] <- m3[i] / ((m2[i])^1.5)
sktests[i] <- as.integer(abs(sk[i]) >= cv) }
        pwr[j] <- mean(sktests)
}
#plot power vs epsilon 
plot(epsilon, pwr, type = "b",
xlab = bquote(epsilon)) 
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors 
lines(epsilon, pwr+se, lty = 3,col = "red")
lines(epsilon, pwr-se, lty = 3,col = "red")

## -----------------------------------------------------------------------------
count5test <- function(x, y) {
  X <- x - mean(x)
  Y <- y - mean(y)
  outx <- sum(X > max(Y)) + sum(X < min(Y))
  outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0)
  return(as.integer(max(c(outx, outy)) > 5))
}



# generate samples under H1 to estimate power
sigma1 <- 1
sigma2 <- 1.5
m<-1000
power1 <- mean(replicate(m, expr={
  x <- rnorm(20, 0, sigma1)
  y <- rnorm(20, 0, sigma2)
  count5test(x, y)
}))

print(power1)

x <- rnorm(20, 0, sigma1)
y <- rnorm(20, 0, sigma2)
var.test(x, y, ratio = 1,
         alternative = c("two.sided", "less", "greater"),
         conf.level = 0.95)


## -----------------------------------------------------------------------------
n <- c(10, 20, 30, 50, 100, 500) #sample sizes
cv <- qnorm(.975, 0, sqrt(6/n)) #crit. values for each n
sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
}

#n is a vector of sample sizes
#we are doing length(n) different simulations
p.reject <- numeric(length(n)) #to store sim. results
m <- 10000 #num. repl. each sim.
for (i in 1:length(n)) {
  sktests <- numeric(m) #test decisions
  for (j in 1:m) {
    x <- rnorm(n[i])
    #test decision is 1 (reject) or 0
    sktests[j] <- as.integer(abs(sk(x)) >= cv[i] )
  }
  p.reject[i] <- mean(sktests) #proportion rejected
}
cv <- qnorm(.975, 0, sqrt(6*(n-2) / ((n+1)*(n+3)))) 
round(cv, 4)

## -----------------------------------------------------------------------------
alpha <- .1
n <- 30
m <- 2500
epsilon <- c(seq(0, .15, .01), seq(.15, 1, .05))
N <- length(epsilon)
pwr <- numeric(N) #critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))
for (j in 1:N) {
  #for each epsilon
  e <- epsilon[j]
  sktests <- numeric(m)
  for (i in 1:m) {
    #for each replicate
    sigma <- sample(c(1, 10), replace = TRUE, size = n, prob = c(1-e, e))
    x <- rnorm(n, 0, sigma)
    sktests[i] <- as.integer(abs(sk(x)) >= cv)
  }
  pwr[j] <- mean(sktests)
} 

#plot power vs epsilon
plot(epsilon, pwr, type = "b", xlab = bquote(epsilon), ylim = c(0,1))
abline(h = .1, lty = 3)
se <- sqrt(pwr * (1-pwr) / m) #add standard errors
lines(epsilon, pwr+se, lty = 3)
lines(epsilon, pwr-se, lty = 3)

## -----------------------------------------------------------------------------
# initialize input and output
library(energy)
alpha <- .1
n <- 30
m <- 2500 #try smaller m for a trial run
epsilon <- .1
test1 <- test2 <- test3 <- numeric(m)

sk <- function(x) {
  #computes the sample skewness coeff.
  xbar <- mean(x)
  m3 <- mean((x - xbar)^3)
  m2 <- mean((x - xbar)^2)
  return( m3 / m2^1.5 )
} 

#critical value for the skewness test
cv <- qnorm(1-alpha/2, 0, sqrt(6*(n-2) / ((n+1)*(n+3))))

# estimate power
for (j in 1:m) {
  e <- epsilon
  sigma <- sample(c(1, 10), replace = TRUE, size = n, prob = c(1-e, e))
  x <- rnorm(n, 0, sigma)
  test1[j] <- as.integer(abs(sk(x)) >= cv)
  test2[j] <- as.integer(shapiro.test(x)$p.value <= alpha)
  test3[j] <- as.integer(mvnorm.etest(x, R=200)$p.value <= alpha)
}
print(c(epsilon, mean(test1), mean(test2), mean(test3)))
detach(package:energy)



## -----------------------------------------------------------------------------
data(law,package = "bootstrap")
n<-nrow(law)
LSAT<-law$LSAT
GPA<-law$GPA
theta.hat<-cor(LSAT,GPA)

# jackknife estimate of bias
theta.jack<-numeric(n)
for (i in 1:n) {
  theta.jack[i]<-cor(LSAT[-i],GPA[-i])
  bias<-(n-1)*(mean(theta.jack)-theta.hat)
}
print(bias)

# jackknife estimate of standard error
se<-sqrt((n-1)*mean((theta.jack-mean(theta.jack))^2))
print(se)

## -----------------------------------------------------------------------------
library(boot)
data(aircondit, package = "boot")
x<-aircondit$hours
qqnorm(x);
qqline(x);

boot.obj <- boot(x, R = 2000,
                 statistic = function(x, i){
                   mean(x[i])})
print(boot.ci(boot.obj, type=c("norm","basic","perc","bca")))

## -----------------------------------------------------------------------------
library(boot)
library(bootstrap)
data(scor, package="bootstrap")

n <- nrow(scor)
df <- as.data.frame(scor)
theta.hat <- eigen(cov(df))$value[1] / sum(eigen(cov(df))$value)

#Jackknife estimate of bias and standard error
theta.jack<-numeric(n)
for (i in 1:n) {
  x<-df[-i,]
  lambda <- eigen(cov(x))$values
  theta.jack[i] <- lambda[1] / sum(lambda)
}
bias <- (n - 1) * (mean(theta.jack) - theta.hat)
se <- (n - 1) * sqrt(var(theta.jack) / n)

rbind(bias,se)

## -----------------------------------------------------------------------------
library(DAAG);
attach(ironslag)

n <- length(magnetic) #in DAAG ironslag
N <- combn(n, 2) 
m<-dim(N)[2]
e1 <- e2 <- e3 <- e4 <- numeric(m)

# for n-fold cross validation
# fit models on leave-two-out samples
for (k in 1:m) {
  lto<-N[,k]
  y <- magnetic[-lto]
  x <- chemical[-lto] 
  
  J1 <- lm(y ~ x)
  yhat1 <- J1$coef[1] + J1$coef[2] * chemical[lto]
  e1[k] <- sum((magnetic[lto] - yhat1)^2)
  
  J2 <- lm(y ~ x + I(x^2))
  yhat2 <- J2$coef[1] + J2$coef[2] * chemical[lto]+J2$coef[3]           * chemical[lto]^2 
  e2[k] <- sum((magnetic[lto] - yhat2)^2)
  
  J3 <- lm(log(y) ~ x)
  logyhat3 <- J3$coef[1] + J3$coef[2] * chemical[lto]
  yhat3 <- exp(logyhat3)
  e3[k] <- sum((magnetic[lto] - yhat3)^2)
  
  J4 <- lm(log(y) ~ log(x))
  logyhat4 <- J4$coef[1] + J4$coef[2] * log(chemical[lto])
  yhat4 <- exp(logyhat4)
  e4[k] <- sum((magnetic[lto] - yhat4)^2)
}

c(mean(e1),mean(e2),mean(e3),mean(e4))

## -----------------------------------------------------------------------------
model<-lm(magnetic~chemical+I(chemical^2),data = ironslag)
summary(model)

## -----------------------------------------------------------------------------
n1 <- 100; n2 <- 150
set.seed(12345)
m <- 500
count5test <- function(x, y, s) {
X <- x - mean(x)
Y <- y - mean(y)
outx <- sum(X > max(Y)) + sum(X < min(Y))
outy <- sum(Y > max(X)) + sum(Y < min(X))
# return 1 (reject) or 0 (do not reject H0) 
return(as.integer(max(c(outx, outy)) > s))
}
x <- rnorm(n1)
y <- rnorm(n2)
s <- 5:15
R <- 100
q <- numeric(R)
alphahat <- pwr <- numeric(length(s))
for (j in 1:length(s)) {
  ss <- s[j]
  alphahat[j] <- count5test(x, y, ss) 
  z <- c(x, y)
  K <- 1:(n1+n2); n<-length(x)
  for (i in 1:R) {
  k <- sample(K, size = n, replace = FALSE)
  x1 <- z[k]; y1 <- z[-k] #complement of x1
  x1 <- x1 - mean(x1) 
  #centered by sample mean 
  y1 <- y1 - mean(y1)
  q[i] <- count5test(x1, y1, ss)
 }
 pwr[j] <- mean(c(alphahat[j], q))
}
plot(s, pwr, col = "red")

## -----------------------------------------------------------------------------
library(RANN) 
library(energy)
library(boot)
library(Ball)

## -----------------------------------------------------------------------------
Te <- function(z, x, sizes,k) {
  n1 <- sizes[1]
  n2 <- sizes[2]
  n <- n1 + n2
  if(is.vector(z)) z <- data.frame(z,0);
  z <- z[x, ]
  NN <- nn2(data=z, k=k+1) 
  block1 <- NN$nn.idx[1:n1,-1] 
  block2 <- NN$nn.idx[(n1+1):n,-1] 
  i1 <- sum(block1 < n1 + .5)
  i2 <- sum(block2 > n1+.5) 
  return((i1 + i2) / (k * n))
}

set.seed(12345)
m<-100
k<-3
n1<-n2<-50
n <- n1+n2
N = c(n1,n2)
eqdist.nn <- function(z,sizes,k){
   boot.obj <- boot(data=z,statistic=Te,R=999,
   sim = "permutation", sizes = sizes,k=k)
   ts <- c(boot.obj$t0,boot.obj$t)
   p.value <- mean(ts>=ts[1])
   list(statistic=ts[1],p.value=p.value)
}
 p.values <- matrix(0,m,3)
 power.comp<-function(mu1,mu2,sigma1,sigma2,alpha){
   x <- rnorm(n1,mu1,sigma1)
   y <- rnorm(n2,mu2,sigma2)
   z <- c(x,y)
   for(i in 1:m){

   p.values[i,1] <- eqdist.nn(z,N,k)$p.value
   p.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
   p.values[i,3] <-bd.test(x,y,R=999,seed=i*12345)$p.value
}
 pow<-apply(p.values<alpha,2,mean)
 names(pow)<-c("NN","energy","Ball")
 return(pow)
 }
 

power.comp(0,0,1,1.5,0.055) #Unequal variances and equal expectations


power.comp(0.5,0,1,1.5,0.02) #Unequal variances and unequal expectations


## -----------------------------------------------------------------------------
set.seed(12345)
m<-100
k<-3
n1<-n2<-50
n <- n1+n2
N = c(n1,n2)

eqdist.nn <- function(z,sizes,k){
   boot.obj <- boot(data=z,statistic=Te,R=999,
   sim = "permutation", sizes = sizes,k=k)
   ts <- c(boot.obj$t0,boot.obj$t)
   p.value <- mean(ts>=ts[1])
   list(statistic=ts[1],p.value=p.value)
}
 pt.values <- matrix(0,m,3)

 
 for(i in 1:m){
   x <- rt(n1,df=1)
   y <- rt(n2,df=1)
   z <- c(x,y)
   pt.values[i,1] <- eqdist.nn(z,N,k)$p.value
   pt.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
   pt.values[i,3]<-bd.test(x,y,R=999,seed=i*12345)$p.value
}
alpha.t <- 0.2
pow.t<-apply(pt.values<alpha.t,2,mean)
names(pow.t)<-c("NN","energy","Ball")
pow.t

## -----------------------------------------------------------------------------
 pb.values <- matrix(0,m,3)
 for(i in 1:m){
   x <- 0.5*rnorm(n1,0,1)+0.5*rnorm(n1,1,1)
   y <- 0.5*rnorm(n1,0,1.5)+0.5*rnorm(n1,1,1.5)
   z <- c(x,y)
   pb.values[i,1] <- eqdist.nn(z,N,k)$p.value
   pb.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
   pb.values[i,3] <-bd.test(x,y,R=999,seed=i*12345)$p.value
}
 alpha.b <- 0.2
 pow.b<-apply(pb.values<alpha.b,2,mean)
 names(pow.b)<-c("NN","energy","Ball")
 pow.b

## -----------------------------------------------------------------------------
set.seed(12345)
n1<-10
n2<-100
n <- n1+n2
N = c(n1,n2)
m<-50
power.comp<-function(mu1,mu2,s1,s2,a){
   x <- rnorm(n1,mu1,s1)
   y <- rnorm(n2,mu2,s2)
   z <- c(x,y)
   for(i in 1:m){
      
   pt.values[i,1] <- eqdist.nn(z,N,3)$p.value
   pt.values[i,2] <- eqdist.etest(z,sizes=N,R=999)$p.value
   pt.values[i,3] <-bd.test(x,y,R=999,seed=i*12345)$p.value
}
 alpha <- a
 pow<-apply(pt.values<a,2,mean)
 names(pow)<-c("NN","energy","Ball")
 return(pow)
}
power.comp(0,0,1,1.5,0.1)

## -----------------------------------------------------------------------------
lap.f <- function(x){
  # prop density function of Laplace distribution
  return(1/2*exp(-abs(x)))
}

## -----------------------------------------------------------------------------
rw.Me <- function(sigma, x0, N){ 
  # function to generate a random walk metropolis chain
  x <- numeric(N) 
  x[1] <- x0 
  u <- runif(N) 
  k <- 0 
  for (i in 2:N){ 
    y <- rnorm(1, x[i-1], sigma) 
    if (u[i] <= (lap.f(y)/lap.f(x[i-1]))){
      x[i] <- y 
    }else{ 
      x[i] <- x[i-1] 
      k <- k + 1 } 
    } 
  return(list(x=x, k=k)) 
} 

N <- 2000 
sigma <- c(.05, .5, 2, 16)

x0 <- 25 
rw1 <- rw.Me(sigma[1], x0, N) 
rw2 <- rw.Me(sigma[2], x0, N) 
rw3 <- rw.Me(sigma[3], x0, N) 
rw4 <- rw.Me(sigma[4], x0, N)

## ----eval=FALSE---------------------------------------------------------------
#  par(mfrow=c(2,2))
#  plot((1:N),rw1$x,type = "l",main = expression(paste(sigma,"=",0.05)))
#  plot((1:N),rw2$x,type = "l",main = expression(paste(sigma,"=",0.5)))
#  plot((1:N),rw3$x,type = "l",main = expression(paste(sigma,"=",2)))
#  plot((1:N),rw4$x,type = "l",main = expression(paste(sigma,"=",16)))

## -----------------------------------------------------------------------------
rw.k <- c(rw1$k, rw2$k, rw3$k, rw4$k)
rate.acceptance <- N/(rw.k+N)
rbind(sigma,rate.acceptance)

## -----------------------------------------------------------------------------
Gelman.Rubin <- function(psi) {
  # psi[i,j] is the statistic psi(X[i,1:j])
  # for chain in i-th row of X
  psi <- as.matrix(psi)
  n <- ncol(psi)
  k <- nrow(psi)
  psi.means <- rowMeans(psi)     #row means
  B <- n * var(psi.means)        #between variance est.
  psi.w <- apply(psi, 1, "var")  #within variances
  W <- mean(psi.w)               #within est.
  v.hat <- W*(n-1)/n + (B/n)     #upper variance est.
  r.hat <- v.hat / W             #G-R statistic
  return(r.hat)
}

## -----------------------------------------------------------------------------
normal.chain <- function(sigma, N, X1) {
  #generates a Metropolis chain for the standard Laplace distribution
  #with Normal(X[t], sigma) proposal distribution
  #and starting value X1
  x <- rep(0, N)
  x[1] <- X1
  set.seed(1122)
  u <- runif(N)
  for (i in 2:N) {
    xt <- x[i-1]
    y <- rnorm(1, xt, sigma) #candidate point
    r1 <- lap.f(y) * dnorm(xt, y, sigma)
    r2 <- lap.f(xt) * dnorm(y, xt, sigma)
    r <- r1 / r2
    if (u[i] <= r) 
      x[i] <- y 
    else
      x[i] <- xt
    }
  return(x)
}

## ----eval=FALSE---------------------------------------------------------------
#  sigma <- 1   #parameter of proposal distribution
#  k <- 4       #number of chains to generate
#  n <- 15000   #length of chains
#  b <- 1000    #burn-in length
#  
#  #choose overdispersed initial values
#  x0 <- c(-10, -5, 5, 10)
#  
#  #generate the chains
#  X <- matrix(0, nrow=k, ncol=n)
#  for (i in 1:k)
#    X[i, ] <- normal.chain(sigma, n, x0[i])
#  
#  #compute diagnostic statistics
#  psi <- t(apply(X, 1, cumsum))
#  for (i in 1:nrow(psi))
#    psi[i,] <- psi[i,] / (1:ncol(psi))
#  print(Gelman.Rubin(psi))
#  
#  #plot psi for the four chains
#  par(mfrow=c(2,2))
#  for (i in 1:k)
#    plot(psi[i, (b+1):n], type="l",
#         xlab=i, ylab=bquote(psi))
#  par(mfrow=c(1,1)) #restore default
#  
#  #plot the sequence of R-hat statistics
#  rhat <- rep(0, n)
#  for (j in (b+1):n)
#    rhat[j] <- Gelman.Rubin(psi[,1:j])
#  plot(rhat[(b+1):n], type="l", xlab="", ylab="R")
#  abline(h=1.2, lty=2)

## -----------------------------------------------------------------------------
cupper <- function(k,a){
  return(sqrt(a^2*k/(k+1-a^2)))
}
f1 <- function(u){
  (1+u^2/(k-1))^(-k/2)
}
f2 <- function(u){
  (1+u^2/k)^(-(k+1)/2)
}


kt <- c(4:25,100)
n <- length(kt)
A  <- numeric(n)


sol <- function(a){
   # the toot of sol2 is A(k)
   1-pt(cupper(k-1,a),k-1)-1+pt(cupper(k,a),k)
  }
for (i in 1:n) {
   k <- kt[i]
   A[i] <- uniroot(sol,c(1e-5,sqrt(k)-1e-5))$root
 }

cbind(df.t=kt,root.ex11.4=A
      
      )

## -----------------------------------------------------------------------------
nA. <- 444; nB. <- 132; nOO <- 361; nAB <- 63
p <- q <- r <- numeric(100)
p[1] <- 0.2; q[1] <- 0.2; r[1] <- (1- p[1]- q[1])
   #Given initial value of iteration
f <- function(a,b) {
  return((nB.*b/(2-b-2*a)+nB.+nAB)/(nA.*a/(2-a-2*b)+nA.+nAB))
}
g <- function(a,b) {
 return(((1-a/(2-a-2*b))*nA.+(1-b/(2-b-2*a))*nB.+2*nOO)/((nB.*b/(2-b-2*a)+
                          nB.+nAB)))
}
threshold <- 1e-5
#Given the threshold
for (k in 2:100) {
   p[k] <- 1/(1+f(p[k-1],q[k-1])*(1+g(p[k-1],q[k-1])))
   q[k] <- f(p[k-1],q[k-1])/(1+f(p[k-1],q[k-1])*(1+g(p[k-1],q[k-1])))
   r[k] <- 1- p[k] - q[k]
   #Through the theoretical steps of the EM algorithm, the relationship between the iteration value at each step and the previous iteration value is obtained.
   if((p[k]-p[k-1] <= threshold) & (q[k]-q[k-1] <= threshold) &
      (r[k]-r[k-1] <= threshold))
   #If the difference between two iterations of p, q, r is less than a given threshold, stop iteration
       {print(c(k, p[k], q[k],r[k]))
       break
    }
}

## -----------------------------------------------------------------------------
x <- seq(1,k,1)
plot(x, p[1:k], "b", col = "red",ylim=c(0,0.6), main = "The log-maximum likelihood values in M-steps" , xlab = "The number of iteration", ylab = "The value of iteration")
lines(x, q[1:k], "b", col = "blue")
lines(x, r[1:k], "b", col = "green")
legend("topright", legend = c("p", "q", "r"),lty = 1, col = c("red", "blue", "green"))

## ----warning=FALSE,message=FALSE----------------------------------------------
data(mtcars)
formulas = list(mpg ~ disp,mpg ~ I(1 / disp),mpg ~ disp + wt,mpg ~ I(1 / disp) + wt)
# for loop
result_1 = list()
for (i in formulas){
  result_1 = c(result_1,list(lm(data=mtcars,i)))
}
#lapply
result_2 = lapply(formulas,lm,data=mtcars)
print(result_1)
print(result_2)

## ----warning=FALSE, message=FALSE---------------------------------------------
trials <- replicate(
100,
t.test(rpois(10, 10), rpois(7, 10)),
simplify = FALSE
)
# Using anonymous function 
result_anonymous = sapply(trials, function(i) i$p.value)
# Don't use anonymous fnction
result_NO_anonymous = sapply(trials,'[[',i=3)
print(result_anonymous)
print(result_NO_anonymous)

## -----------------------------------------------------------------------------
library(parallel)
boot_df <- function(x) x[sample(nrow(x), rep = T), ]
rsquared <- function(mod) summary(mod)$r.square
boot_lm <- function(i) {
rsquared(lm(mpg ~ wt + disp, data = boot_df(mtcars)))
}

system.time(lapply(1:500, boot_lm))

## -----------------------------------------------------------------------------
set.seed(1201)
rw_MeR <- function(sigma, x0, N) {
  x <- numeric(N)
  x[1] <- x0
  u <- runif(N)
  k <- 0
  for (i in 2:N) {
    y <- rnorm(1, x[i-1], sigma)
    if (u[i] <= exp(-(abs(y) - abs(x[i-1]))))
      x[i] <- y else {
        x[i] <- x[i-1]
        k <- k + 1
      }
  }
  return(list(x=x, k=k))
}

## -----------------------------------------------------------------------------
library(Rcpp)
#include <Rcpp.h>
#using namespace Rcpp;
#// [[Rcpp::export]]
cppFunction('NumericVector rw_MeC(double sigma, double x0, int N) {
  NumericVector x(N);
  x[0] = x0;
  double u, y;
  int k = 0;
  for (int i = 1; i < N; i++) 
  {
    y = rnorm(1, x[i-1], sigma)[0];
    u = runif(1)[0];
    if (u <= exp(-(abs(y) - abs(x[i-1])))) 
    {
      x[i] = y; 
    }
    else 
    {
      x[i] = x[i-1];
      k++;
    }
  }
  return x;
}')

## ----eval=FALSE---------------------------------------------------------------
#  library(microbenchmark)
#  N = 2000
#  sigma <- c(0.5,1,10,100)
#  x0 = 25
#  for (i in 1:length(sigma)) {
#  ts = microbenchmark(rwR = rw_MeR(sigma[i], x0, N)$x,
#                      rwC = rw_MeC(sigma[i], x0, N))
#  print(summary(ts)[, c(1,3,5,6)])
#  
#  rwR = rw_MeR(sigma[i], x0, N)$x
#  rwC = rw_MeC(sigma[i], x0, N)
#  par(mfrow = c(2, 2))
#  b <- 1000
#  y <- (rwR)[b:N]
#  a <- ppoints(500)
#  QR <- ifelse(a <= 1/2, log(2*a), -log(2-2*a))
#  Q1 <- quantile(rwR, a)
#  qqplot(QR, Q1, main=paste("R  sigma=",sigma[i]))
#  abline(a=0, b=1)
#  
#  y <- (rwC)[b:N]
#  a <- ppoints(500)
#  QR <- ifelse(a <= 1/2, log(2*a), -log(2-2*a))
#  Q2 <- quantile(rwC, a)
#  qqplot(QR, Q2, main=paste("C  sigma=",sigma[i]),col="green")
#  abline(a=0, b=1)
#  
#  qqplot(Q1, Q2, main=paste("C-R  sigma=",sigma[i]),col="red")
#  abline(a=0, b=1)
#  }

