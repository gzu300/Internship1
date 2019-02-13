#########################data##################

library(tidyverse)
library(purrr)
#source('simulation_util.R')
source('MaSigPro_util.R')
source('ASCA-genes.1.2.1/sourceASCA.R')

######
#generate simulation data

#6 time points;2 treatments;4 replicates
# 
# rep <- 4
# sd <- 0
# ctrl <- 0.01*time+0.1
# tre1 <- -0.01*time-0.1
# fluctuate <- list(2,-2,2,-2,2,-2)
# notrend <- rep(list(0),6)
# ctrl <- notrend
# tre1 <- fluctuate
# tre2 <- fluctuate
# 
# trend1 <- append(ctrl,tre1)
# trend2 <- append(ctrl,tre2)
# trend3 <- append(ctrl,tre2)
# trend4 <- append(ctrl,ctrl)
#create number of features in each pattern and names of the patterns
ftrs_in_pat <- list(30,30,1,100)
pat_names <- list('one','two','three','flat') 

create.simulation <- function(pat_names.list,ftrs_in_pat.list,replicates,sd,...){
  dimnames <- pat_names.list %>% 
    map2(.,ftrs_in_pat.list,~rep(.x,.y)) %>% 
    flatten_chr() %>% 
    list(NULL,.)
  
  #create the simulation matrix
  wrap_tre <- function(dimnames,replicates,sd,...){
    trend_list <- list(...)
      trend_list %>% 
      map2(pat,.,~rep(.y,.x)) %>% 
      flatten() %>% 
      map(.,~rnorm(n = replicates,mean = .x,sd = sd)) %>% 
      flatten_dbl() %>% 
      matrix(.,nrow = length(trend_list[[1]])*replicates,dimnames = dimnames) %>% 
      t(.)
  }
  wrap_tre <- wrap_tre(dimnames = dimnames,replicates = rep,sd = sd,trend1,trend2,trend3,trend4)
  output <- list(df=wrap_tre,groups=dimnames)
  output
}
# a <- create.simulation(pat_names.list = pat_names,ftrs_in_pat,rep,sd,trend1,trend2,trend3,trend4)
# b <- a$df
