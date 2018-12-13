#################
#data simulation#
#################

#a)single group continuous induction
#b)single group transitory repression
#c)differential multi-group induction
library(ggplot2)
library(maSigPro)
library(tidyverse)

###############
#design matrix#
###############
create_dummy <- function(i, j, r){
  dummy.groups <- list()
  rows <- i*j*r
  for (k in 1:i){
    n <- k
    dummy.groups[[n]] <- rep(0, rows)
    if (n == 1){
      dummy.groups[[n]][1:(j*r)] <- 1
    }
    else{
      x <- 1+(j*r)*(k-1)
      y <- (j*r)*k
      dummy.groups[[n]][x:y] <- 1
    }
  }
  return(dummy.groups)
}


design_exp.raw <- function(i, j, r,time){
  #i:number of groups
  #j:number of time points
  #r:number of replicates
  #v:list(list_i(c(j1,j2,j3), c(j1,j2,j3)), list(c(j1,j2,j3)))
  Group <- factor(rep(c(1:i), each=j*r))
  Time <- factor(rep(time, each=r))
  Replicate <- factor(rep(c(1:r), j))
  df <- data.frame(Group=Group, Time=Time, Replicate=Replicate)
  return(df)
}

design_exp <- function(i, j, r,time){
  #i:number of groups
  #j:number of time points
  #r:number of replicates
  #v:list(list_i(c(j1,j2,j3), c(j1,j2,j3)), list(c(j1,j2,j3)))
  df.raw <- design_exp.raw(i,j,r,time)
  df <- df.raw %>% 
    unite('index', c(Group, Time, Replicate), sep = '_') %>% 
    column_to_rownames(var = 'index')
  return(df)
}

design_matrix <- function(i, j, r, time){
  time.col <- rep(time, each=r)
  replicate.col <-  rep(1:(i*j), each = r)
  group <- paste('treatment', 1:i, sep = '')
  dummy.groups <- create_dummy(i, j, r)
  design.matrix <- data.frame(time=time.col, replicate=replicate.col, dummy.groups)
  colnames(design.matrix) <- noquote(c('Time', 'Replicate', group))
  rownames(design.matrix) <- rownames(design_exp(i,j,r,time))
  return(design.matrix)
}




