library('microbenchmark')


#______pour le calcul des efficacites relatives
coutDeSimulation <- function(f,g, n = 1){
  a <- microbenchmark(f(n),g(n))
  m <- summary(a)
  as.data.frame(list(m[1],m[4]))
}

t <- 8
lambda <- 2
n <- 1000


N <- 10000
n <- 10000
#pour le calcul de l'efficacité relative par la methode de monte carlo classique
f1 <- function(n){
  S <- rpois(n, 3.7)
  x <- numeric(n)
  for(i in 1:n){
    x[i] <- sum(rweibull(S[i], 0.5, 2))
  }
  y <- (x < 3)
  MC.estim(y)
}
#pour le calcul de l'efficacité relative par la methode strates proportionnelle 

strat_prop <- function(n, nk, level = 0.95 ){
  nk <- nk[nk>1]
  K <- length(nk)
  sigma2 <- numeric(K)
  delta <- numeric(K)
  
  for(k in 1:K){
    x <- numeric(nk[k])
    for(i in 1:nk[k]){
      x[i] <- sum(rweibull(k, 0.5, 2))
    }
    y <- x < 3
    delta[k] <- MC.estim(y)$delta
    sigma2[k] <- MC.estim(y)$Var #estimation des variances intra strat pour le calcul de la variance de l'estimateur
  }
  
  t <- 1:K
  p <- exp(-3.7) * (3.7)^(t) / factorial(t)
  q <- qnorm((1 + level) / 2)
  Var <- sum(p * p * sigma2 / nk)
  i1 <- sum(delta * p) + ppois(0,3.7) - q * sqrt(Var / n)
  i2 <- sum(delta * p) + ppois(0,3.7) + q * sqrt(Var / n)
  list(delta = sum(delta * p) + ppois(0,3.7), Var = Var, IC =c(i1,i2) )
}

#______________________
t <- 1:K
p <- exp(-3.7) * (3.7)^(t) / factorial(t)
nk <- round(p * rep(n, K), 0)
strat_prop(n, nk, 0.95)

strat_opti <- function(n, nk, level = 0.95, sig(K)){
  sigma2<-sig(K)
  sigma2 <- sigma2[nk>0]
  nk <- nk[nk>0]
  K <- length(nk)
  delta <- numeric(K)
  
  for(k in 1:K){
    x <- numeric(nk[k])
    for(i in 1:nk[k]){
      x[i] <- sum(rweibull(k, 0.5, 2))
    }
    y <- x < 3
    delta[k] <- MC.estim(y)$delta
  }
  
  t <- 1:K
  p <- exp(-3.7) * (3.7)^(t) / factorial(t)
  q <- qnorm((1 + level) / 2)
  Var <- sum(p * p * sigma2 / nk)
  i1 <- sum(delta * p) + ppois(0,3.7) - q * sqrt(Var / n)
  i2 <- sum(delta * p) + ppois(0,3.7) + q * sqrt(Var / n)
  list(delta = sum(delta * p) + ppois(0,3.7), Var = Var, IC =c(i1,i2) )
}

K <- 20 # K = +infini, 

#etape1:estimation des variances intra-strate
m <- 100
sig<-function(K){
  sigma2 <- numeric(K)
  for(k in 1:K){
    x <- numeric(m)
    for(i in 1:m){
      x[i] <- sum(rweibull(k, 0.5, 2))
    }
    y <- x < 3
    sigma2[k] <- MC.estim(y)$Var
  }
  return(sigma2)
}
#etape2:estimation de delta


t <- 1:K
p <- exp(-3.7) * (3.7)^(t) / factorial(t)
qk <- p * sqrt(sig(K))
qk <- qk/sum(qk)
nk <- round(qk*rep(n, K), 0)
strat_opti(n, nk, 0.95 ,sig(K))


microbenchmark(strat_opti(n, nk, 0.95 ,sig(K)), strat_prop(N,nk, 0.95), times = 10L)