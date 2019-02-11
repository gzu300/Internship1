#########################data##################

library(tidyverse)
library(purrr)
#source('simulation_util.R')
source('MaSigPro_util.R')
source('ASCA-genes.1.2.1/sourceASCA.R')

######
#generate simulation data
######
#6 time points;2 treatments;4 replicates
#simulation generating pipelines:
#1.1)pmap for all the features,'ftrs_in_pat' is the number of featurs in each pattern.
#  example:if patterns <- c(3,2,1). the output 'ftrs_in_pat' is list(3,3,3,2,2,1)
ftrs_in_pat <- as.vector(patterns,mode='list')
#  2)n(n=replicates) is the number of independent drawings from rnorm. n must have the same length as mean
n <- as.vector(rep(n,length(mean)))
#  3)mean is the mean of rnorm.
#   example:as.vector(c(1,2,3),mode='list'). means for the 3 time points are 1,2 and 3
#  4)sd has the same structure as n
mean <- as.vector(time,mode = 'list')
sd <- as.vector(rep(sd,length(mean)))
#2.different treatments

