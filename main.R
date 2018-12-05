#########################data##################
######
#simulation data
######

i=2#treatment
j=5#time
r=4#replicate
time=c(0,2,6,10,24)#time points

#total of 100 metabolites(columns)
#total of 2*5*4=40 rows

######
#generate data
######
base.mx <- matrix(rnorm(4000,0,0),nrow = i*j*r, ncol = 100)
trend1 <- c(rep(c(1,2,3,4),each=j),rep(0,j*r))
trend2 <- c(rep(0,j*r),rep(c(1,2,3,4),each=j))
base.mx[,1:2] <- base.mx[,1:2]+cbind(trend1,trend2)

########
#real data
########
source('simulation_util.R')
library(readxl)
df <- read_xlsx('macrophage%20polarization.xlsx',range = 'PeakArea!I4:BL43',col_names = c(paste('met',1:56,sep = '')))
###pretreatment####


#####################masigpro#############

######
#design matrix
######
source('MaSigPro_util.R')
masigpro.design <- design_matrix(i,j,r,time)
design <- make.design.matrix(edesign = masigpro.design,degree = 2)
dis <- design$dis
edesign <- design$edesign

########
#fit model
########

#######
#plot
######

#######################asca-gene#############
######
#design matrix
#####
mx.i <- matrix(0,nrow = i*j*r, ncol = i, dimnames = list(c(),c('control','treatment')))+c(rep(1,j*r),rep(0,i*j*r))
mx.j <- matrix(0,nrow = i*j*r, ncol = j, dimnames = list(c(),c('T1','T2','T3','T4','T5')))+c(rep(c(rep(1,r),rep(0,j*r-r)),i),rep(0,r))
mx.ij <- matrix(0,nrow = i*j*r, ncol = i*j, dimnames = list(c(),c(paste('inter',1:(i*j)))))+c(rep(1,r),rep(0,i*j*r))

######
#fit model
######

######
#plot
######
