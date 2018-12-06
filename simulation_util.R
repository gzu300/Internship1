library(ggplot2)
library(maSigPro)
library(tidyverse)
source('MaSigPro_util.R')

#################
#generate data  #
#################


sampling <- function(v, i, j, r){
  collector <- vector()
  for (each in v){
    base_trend <- rep(each, each = r)
    #print(base_trend)
    draw_sample <- diag(rnorm(length(base_trend), 0, 0))
    at_one_group <- c(draw_sample %*% base_trend)
    collector <- c(collector, at_one_group)
  }
  remaining <- c(rnorm(i*j*r-(length(collector)), 0, 0))
  combined <- c(collector, remaining)
  #print(combined)
  return(combined)
}

wrap_metabolites <- function(v, i, j, r){
  col.collector <- list()
  for (each in 1:length(v)){
    col.collector[[each]] <- sampling(v[[each]], i, j, r)
  }
  df <- data.frame(col.collector)
  return(df)
}

wrap_randoms <- function(v,i,j,r){
  #make 100 columns random
  nrow <- i*j*r
  ncol <- 20
  mx <- matrix(rnorm(nrow*ncol, 0, 0), nrow = nrow, ncol = ncol, dimnames = list(c(), noquote(c(paste('met', (length(v))+1:ncol, sep = '')))))
  return(mx)
}

##############
#create trend#
##############

create_points <- function(j){
  b <- c(runif(j, 1, 2))
  diff <- c(runif(j, 0.201, 0.3))#b*diff keeps the increament from 0.22 to 0.6 which is a siginificant change
  points.vs <- c(1)
  for (each in 1:j){
    a <- points.vs[each]+(b[each]*diff[each])
    points.vs <- c(points.vs, a)
  }
  return(points.vs)
}

create_pattern <- function(j){
  #choose 4 from 5 w/o replacement and shuffle created points to create pattern
  sign <- sample(c(1, -1), 1)
  points <- create_points(j)*sign
  trend <- c(1, sample(points, size = j-1))
  if (as.integer(trend[3]) == -1){#mimics challenge test
    trend[3] = 1
    trend[4] = 1
  }
  return(trend)
}

wrap_pattern <- function(n, j){
  #n is the total number of trends
  #multi is the number of multi pattern
  #return lists in the list
  wrap <- list()
  for (each in 1:n){
    ls <- list()
    multi <- sample(c(1,1,1,1,2,3), 1)
    for (every in 1:multi){
      a <- create_pattern(j)
      ls[[every]] = a
      wrap[[each]] = ls
    }
  }
  return(wrap)
}

###########
#pipeline #
###########

pipeline <- function(v, i, j, r){
  df <- design_exp(i, j, r)
  selected <- wrap_metabolites(v, i, j, r)
  col.length <- ncol(selected)
  colnames(selected) <- noquote(paste('met', 1:col.length, sep = ''))
  random <- wrap_randoms(v,i,j,r)
  df_ <- cbind(df, selected, random)
  df_final <- t(df_)
  return(df_final)
}

pipeline.raw <- function(v, i, j, r){
  df <- design_exp.raw(i, j, r)
  selected <- wrap_metabolites(v, i, j, r)
  col.length <- ncol(selected)
  colnames(selected) <- noquote(paste('met', 1:col.length, sep = ''))
  #random <- wrap_randoms(v,i,j,r)
  df_final <- cbind(df, selected)
  #df_final <- t(df_)
  return(df_final)
}

####################
#data generation#
####################
#### GENERATE TIME COURSE DATA
## generate n random gene expression profiles of a data set with 
## one control plus 3 treatments, 3 time points and r replicates per time point.
##r replicates
##
tc.GENE <- function(n, r,
                    var11 = 0.01, var12 = 0.01,var13 = 0.01,
                    var21 = 0.01, var22 = 0.01, var23 =0.01,
                    var31 = 0.01, var32 = 0.01, var33 = 0.01,
                    var41 = 0.01, var42 = 0.01, var43 = 0.01,
                    a1 = 0, a2 = 0, a3 = 0, a4 = 0,
                    b1 = 0, b2 = 0, b3 = 0, b4 = 0,
                    c1 = 0, c2 = 0, c3 = 0, c4 = 0)
{
  
  tc.dat <- NULL
  for (i in 1:n) {
    Ctl <- c(rnorm(r, a1, var11), rnorm(r, b1, var12), rnorm(r, c1, var13))  # Ctl group
    Tr1 <- c(rnorm(r, a2, var21), rnorm(r, b2, var22), rnorm(r, c2, var23))  # Tr1 group
    Tr2 <- c(rnorm(r, a3, var31), rnorm(r, b3, var32), rnorm(r, c3, var33))  # Tr2 group
    Tr3 <- c(rnorm(r, a4, var41), rnorm(r, b4, var42), rnorm(r, c4, var43))  # Tr3 group
    gene <- c(Ctl, Tr1, Tr2, Tr3)
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat
}

