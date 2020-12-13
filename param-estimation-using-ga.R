#install.packages("Metrics")
#install.packages("GA")
#install.packages("kSamples")
#install.packages("survival")

library("Metrics")
library("GA")
library("survival")
library("kSamples")

# Import Data Set
inp <- scan("AirCraft.data", list(0,0))
r <- 86
n <- 151
time<-inp[[1]]
cen<-inp[[2]]
y <- cbind(time,cen)
air <- data.frame(time,cen)
S <- survfit(Surv(time,cen, type="right")~1 , data=air)    #kaplan-meier
adtest <- 0
SEED <- 5001
type <- 1         # 1: GA-MSE ,  2: GA-MAE, 3: GA-MAPE   

# Cumulative Distribution Function
G <- function( t, p, alfa1, beta1, alfa2, beta2 )
{
  result <- ( 1-exp(-(t/alfa2)^beta2) ) + p * ( exp(-(t/alfa2)^beta2) - exp(-(t/alfa1)^beta1) ) 
  result
}

# Probability Density Function
g <- function( t, p, alfa1, beta1, alfa2, beta2 )
{
  result <- p * ( (beta1/alfa1) * ((t/alfa1)^(beta1-1)) * exp((t/alfa1)^beta1) ) + 
    (1-p) * ( (beta2/alfa2) * ((t/alfa2)^(beta2-1)) * exp((t/alfa2)^beta2) ) 
  result
}

thm2 <- function( p, alfa1, beta1, alfa2, beta2 )
{
  ys <- 0
  t1 <- 0
  for (i in 1:n) {
    t <- time[i]
    delta <- cen[i]
    if (delta > 0) { t1 <- t }
    ys[i] <- 1- G( t, p, alfa1, beta1, alfa2, beta2 )
  }
  ys
}

thm <- function( p, alfa1, beta1, alfa2, beta2 )
{
  ys <- 0
  for (i in 1:n) {
    t <- time[i]
    delta <- cen[i]
    ys[i] <- g( t, p, alfa1, beta1, alfa2, beta2 )^delta * 
      (1- G( t, p, alfa1, beta1, alfa2, beta2 ))^(1-delta)
  }
  ys
}

thm3 <- function( tm, cn, p, alfa1, beta1, alfa2, beta2 )
{
  ys <- 0
  for (i in 1:length(tm)) {
    t <- tm[i]
    ys[i] <- 1- G( t, p, alfa1, beta1, alfa2, beta2 )
  }
  ys
}

# Fitness Function
f <- function( p, alfa1, beta1, alfa2, beta2 )
{
  ys <- thm2( p, alfa1, beta1, alfa2, beta2 )
  uyg <- switch(type, mse( S$surv , ys ), mae( S$surv , ys ), mape( S$surv , ys ) ) 
  uyg
}


GA2 <- function()
{
  GA <- ga(type = "real-valued", fitness = function(x) -f(x[1], x[2], x[3], x[4], x[5]),
           lower = c(0, 5, 0, 0, 0), upper = c(0.5, 10, 10, 10, 10),
           #min = c(0, 5, 0, 0, 0), max = c(0.5, 10, 10, 10, 10),
           #min = c(0.5, 5, 0, 0, 0), max = c(1, 10, 10, 10, 10),
           popSize = 100, maxiter = 30, pmutation=0.5, seed=SEED)
  sm <- summary(GA)
  
  # Results of Genetic Algorithm
  w <- sm$sol[1,1]
  alfa1 <- sm$sol[1,2]
  beta1 <- sm$sol[1,3]
  alfa2 <- sm$sol[1,4]
  beta2 <- sm$sol[1,5]
  fit   <- sm$fitness
  
  ys <- thm2( w, alfa1, beta1, alfa2, beta2 )
  adt <- ad.test( S$surv, ys )
  sn <- array( c(fit, w, alfa1, beta1, alfa2, beta2,adt$ad[1,1], adt$ad[1,3]) , dim=c(8,1) )
  rownames(sn) <- c("Fitness","w","alfa1","beta1", "alfa2", "beta2","A2", "A2 p-value")
  colnames(sn) <- c("Results")
  sn
}

groups <- function()
{
  # groups for 10-fold cross validation
  g <- matrix(, nrow = 10, ncol = 0)
  gc <- ceiling(n/10)
  
  xn <- n
  xx <- 1:n
  set.seed(5001)
  gxc <- ceiling( xn/10 )
  xn <- xn - gxc
  c1 <<- sort(sample( xx,gxc )) 
  
  gxc <- ceiling( xn/9 )
  xn <- xn - gxc
  c2 <<- sort(sample( xx[-c1],gxc ) )
  
  gxc <- ceiling( xn/8 )
  xn <- xn - gxc
  c3 <<- sort(sample( xx[c(-c1,-c2)],gxc ))
  
  gxc <- ceiling( xn/7 )
  xn <- xn - gxc
  c4 <<- sort(sample( xx[c(-c1,-c2,-c3)],gxc ) )
  
  gxc <- ceiling( xn/6 )
  xn <- xn - gxc
  c5 <<- sort(sample( xx[c(-c1,-c2,-c3,-c4)],gxc ) )
  
  gxc <- ceiling( xn/5 )
  xn <- xn - gxc
  c6 <<- sort(sample( xx[c(-c1,-c2,-c3,-c4,-c5)],gxc ) )
  
  gxc <- ceiling( xn/4 )
  xn <- xn - gxc
  c7 <<- sort(sample( xx[c(-c1,-c2,-c3,-c4,-c5,-c6)],gxc ) )
  
  gxc <- ceiling( xn/3 )
  xn <- xn - gxc
  c8 <<- sort(sample( xx[c(-c1,-c2,-c3,-c4,-c5,-c6,-c7)],gxc ) )
  
  gxc <- ceiling( xn/2 )
  xn <- xn - gxc
  c9 <<- sort(sample( xx[c(-c1,-c2,-c3,-c4,-c5,-c6,-c7,-c8)],gxc )) 
  
  c10 <<- xx[c(-c1,-c2,-c3,-c4,-c5,-c6,-c7,-c8, -c9)]
}

newdata <- function( test_data )
{
  # Training Dataset for a fold
  time <<- times[-test_data]
  cen  <<- cens[-test_data]
  n    <<- length( time )
  
  y    <<- cbind(time,cen)
  air  <<- data.frame(time,cen)
  S    <<- survfit(Surv(time,cen, type="right")~1 , data=air)
  adtest<<- 0
}

CV <- function()
{
  # 10-fold Cross Validation algorithm for the Method
  SEED  <<- 5001
  
  times <<- time
  cens  <<- cen
  
  print("");print("---------- Cross Validation Fold 1 -------------")
  newdata( c1 )
  train1 <- GA2()
  colnames(train1) <- c("Training-1")
  
  ys  <- thm3( times[c1], cens[c1], train1[2], train1[3], train1[4], train1[5], train1[6] )
  air <- data.frame(times[c1], cens[c1])
  St  <- survfit(Surv(times[c1],cens[c1], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test1 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 2 -------------")
  newdata( c2 )
  train2 <- GA2()
  colnames(train2) <- c("Training-2")
  
  ys  <- thm3( times[c2], cens[c2], train2[2], train2[3], train2[4], train2[5], train2[6] )
  air <- data.frame(times[c2], cens[c2])
  St  <- survfit(Surv(times[c2],cens[c2], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test2 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 3 -------------")
  newdata( c3 )
  train3 <- GA2()
  colnames(train3) <- c("Training-3")
  
  ys  <- thm3( times[c3], cens[c3], train3[2], train3[3], train3[4], train3[5], train3[6] )
  air <- data.frame(times[c3], cens[c3])
  St  <- survfit(Surv(times[c3],cens[c3], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test3 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 4 -------------")
  newdata( c4 )
  train4 <- GA2()
  colnames(train4) <- c("Training-4")
  
  ys  <- thm3( times[c4], cens[c4], train4[2], train4[3], train4[4], train4[5], train4[6] )
  air <- data.frame(times[c4], cens[c4])
  St  <- survfit(Surv(times[c4],cens[c4], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test4 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 5 -------------")
  newdata( c5 )
  train5 <- GA2()
  colnames(train5) <- c("Training-5")
  
  ys  <- thm3( times[c5], cens[c5], train5[2], train5[3], train5[4], train5[5], train5[6] )
  air <- data.frame(times[c5], cens[c5])
  St  <- survfit(Surv(times[c5],cens[c5], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test5 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 6 -------------")
  newdata( c6 )
  train6 <- GA2()
  colnames(train6) <- c("Training-6")
  
  ys  <- thm3( times[c6], cens[c6], train6[2], train6[3], train6[4], train6[5], train6[6] )
  air <- data.frame(times[c6], cens[c6])
  St  <- survfit(Surv(times[c6],cens[c6], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test6 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 7 -------------")
  newdata( c7 )
  train7 <- GA2()
  colnames(train7) <- c("Training-7")
  
  ys  <- thm3( times[c7], cens[c7], train7[2], train7[3], train7[4], train7[5], train7[6] )
  air <- data.frame(times[c7], cens[c7])
  St  <- survfit(Surv(times[c7],cens[c7], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test7 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 8 -------------")
  newdata( c8 )
  train8 <- GA2()
  colnames(train8) <- c("Training-8")
  
  ys  <- thm3( times[c8], cens[c8], train8[2], train8[3], train8[4], train8[5], train8[6] )
  air <- data.frame(times[c8], cens[c8])
  St  <- survfit(Surv(times[c8],cens[c8], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test8 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 9 -------------")
  newdata( c9 )
  train9 <- GA2()
  colnames(train9) <- c("Training-9")
  
  ys  <- thm3( times[c9], cens[c9], train9[2], train9[3], train9[4], train9[5], train9[6] )
  air <- data.frame(times[c9], cens[c9])
  St  <- survfit(Surv(times[c9],cens[c9], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test9 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  print("");print("---------- Cross Validation Fold 10 -------------")
  newdata( c10 )
  train10 <- GA2()
  colnames(train10) <- c("Training-10")
  
  ys  <- thm3( times[c10], cens[c10], train10[2], train10[3], train10[4], train10[5], train10[6] )
  air <- data.frame(times[c10], cens[c10])
  St  <- survfit(Surv(times[c10],cens[c10], type="right")~1 , data=air)
  adt <- ad.test( St$surv, ys ) 
  test10 <- c( adt$ad[1,1], adt$ad[1,3] )
  
  # Output
  train <- cbind(train1, train2,train3,train4,train5,train6,train7,train8,train9,train10)
  test <- cbind(test1, test2,test3,test4,test5,test6,test7,test8,test9,test10)
  rownames(test) <- c("A2", "A2 p-value")
  
  list_data <- list(train, test)
  return(list_data)
}

####### Running 10-fold Cross-Validation  ################

type=1    # 1: GA-MSE ,     2: GA-MAE,     3: GA-MAPE
c1  <- 0
c2  <- 0
c3  <- 0
c4  <- 0
c5  <- 0
c6  <- 0
c7  <- 0
c8  <- 0
c9  <- 0
c10 <- 0
times <- 0
cens  <- 0
groups()
CV_result <- CV()

print("");print("Cross-Validation Folds for Training and Test")
CV_result

print("");print("Best Training Model according to Anderson-Darling")
CV_result[[1]][,which.min( CV_result[[1]][7,] )]

print("");print("Mean of A2 of Training")
apply(CV_result[[1]][7:8,], MARGIN = 1, mean)

print("");print("Mean of A2 of Testing")
apply(CV_result[[2]], MARGIN = 1, mean)


