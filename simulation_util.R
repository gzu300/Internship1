library(ggplot2)
library(maSigPro)
library(tidyverse)

#################
#generate data  #
#################

design_exp <- function(i, j, r){
  #i:number of groups
  #j:number of time points
  #r:number of replicates
  #v:list(list_i(c(j1,j2,j3), c(j1,j2,j3)), list(c(j1,j2,j3)))
  Group <- factor(paste('G', rep(c(1:i), each=j*r), sep = ''))
  Time <- rep(c(2^(1:j)), each=r)
  Replicate <- factor(rep(c(1:r), j))
  df <- data.frame(Group=Group, Time=Time, Replicate=Replicate) %>% 
    unite('index', c(Group, Time, Replicate), sep = '_') %>% 
    column_to_rownames(var = 'index')
  return(df)
}

design_exp.raw <- function(i, j, r){
  #i:number of groups
  #j:number of time points
  #r:number of replicates
  #v:list(list_i(c(j1,j2,j3), c(j1,j2,j3)), list(c(j1,j2,j3)))
  Group <- factor(paste('G', rep(c(1:i), each=j*r), sep = ''))
  Time <- rep(c(2^(1:j)), each=r)
  Replicate <- factor(rep(c(1:r), j))
  df <- data.frame(Group=Group, Time=Time, Replicate=Replicate)
  return(df)
}

sampling <- function(v, i, j, r){
  collector <- vector()
  for (each in v){
    base_trend <- rep(each, each = r)
    #print(base_trend)
    draw_sample <- diag(rnorm(length(base_trend), 1, 0.01))
    at_one_group <- c(draw_sample %*% base_trend)
    collector <- c(collector, at_one_group)
  }
  remaining <- c(rnorm(i*j*r-(length(collector)), 1, 0))
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
  mx <- matrix(rnorm(nrow*ncol, 1, 0), nrow = nrow, ncol = ncol, dimnames = list(c(), noquote(c(paste('x', (length(v))+1:ncol, sep = '')))))
  return(mx)
}


###########
#pipeline #
###########

pipeline <- function(v, i, j, r){
  df <- design_exp(i, j, r)
  selected <- wrap_metabolites(v, i, j, r)
  col.length <- ncol(selected)
  colnames(selected) <- noquote(paste('x', 1:col.length, sep = ''))
  random <- wrap_randoms(v,i,j,r)
  df_ <- cbind(df, selected, random)
  df_final <- t(df_)
  return(df_final)
}

pipeline.raw <- function(v, i, j, r){
  df <- design_exp.raw(i, j, r)
  selected <- wrap_metabolites(v, i, j, r)
  col.length <- ncol(selected)
  colnames(selected) <- noquote(paste('x', 1:col.length, sep = ''))
  #random <- wrap_randoms(v,i,j,r)
  df_final <- cbind(df, selected)
  #df_final <- t(df_)
  return(df_final)
}

####################
#asca design matrix#
####################

asca.design.matrix <- function(i,j,r){
  rows <- i*j*r
  blank <- vector(mode = 'list')
  for (col in c(i,j,i*j)){#iterate through different number of design matrix's columns
    mx <- matrix(0,nrow=rows, ncol=col)
    blank <- c(blank,list(mx)) 
  }
  #output <- list(blank[[1]]+c(rep(rep(1,j*r),rep(0,j*r*1),1),rep(0,j*r)))
  output <- list(blank[[1]]+c(rep(1,j*r),rep(0,j*r*1),rep(0,j*r)))
  output <- c(output,list(blank[[2]]+c(rep(c(rep(1,r),rep(0,r*3)),2),rep(0,r))))
  output <- c(output,list(blank[[3]]+c(rep(1,r),rep(0,r*7),rep(0,r))))
  output
}

a <- asca.design.matrix(2,4,4)
a[[3]]
a <- design_exp(2,4,4)
