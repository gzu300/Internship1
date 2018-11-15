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

design_matrix <- function(i, j, r, df){
  time.col <- rep(c(2^c(1:j)), each=r)
  replicate.col <-  rep(1:(i*j), each = r)
  group <- c(paste('G', c(1:i), sep = ''))
  group.col <- rep(group, each = (j*r))
  dummy.groups <- create_dummy(i, j, r)
  design.matrix <- data.frame(time=time.col, replicate=replicate.col, dummy.groups)
  colnames(design.matrix) <- noquote(c('Time', 'Replicate', group))
  rownames(design.matrix) <- colnames(df)
  return(design.matrix)
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


