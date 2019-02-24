library(tidyverse)
library(readxl)
library(knitr)
#setwd('../../asca/')
##########data transformation and exploration####

colnames <- colnames(read_xlsx('data/macrophage%20polarization.xlsx',range = 'PeakArea!I1:BL1'))
df <- read_xlsx('data/macrophage%20polarization.xlsx',range = 'PeakArea!I4:BL43',col_names = colnames)

df$togroup <- rep(1:10,each=4)#add a dummy tag for calculation

df[20,] <- df[17:19,] %>% 
  summarise_all(.,funs(mean))#some of the measurements at time 24, replicate 4 are extremely low. suspected technical error.
                            #Other measurement looked normal are comparable with the other 3. 
                            #So replaced this with the mean of the other 3 replicates

df.clean <- df %>% 
  group_by(togroup) %>% 
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))#fill NA with mean of replicates. 

df.clean[1:4,"Fructose-1,6-bisphosphate"] <- df.clean$`Fructose-1,6-bisphosphate`[5]
# F16B is missing values at entire time 1 control. filled with values at time 1 treatment.
# the assumption is that value at time 0 for both control and treatment is the same. 
# plot of this metabolite seems ok with this replacement

plot.mean_sd <- function(df,...){
  mean_sd <- df %>% 
    group_by(togroup) %>% 
    summarise_all(.,funs(mean,sd)) %>% 
    select(.,-togroup) %>% 
    flatten()
  
  plot(mean_sd[1:560],mean_sd[561:1120],type = 'p',xlab = 'mean',ylab = 'sd',...)
  
}
df.log <- log(df.clean)#log transform

plot.mean_sd(df.log, main = 'log transformed', sub = 'log transformed shows homoscedastic')
plot.mean_sd(df.clean, main = 'original data', sub = 'heteroscedastic')#check the homoscedasticity

mean_sd.df <- df.log %>% 
  group_by(togroup) %>% 
  summarise_all(.,funs(mean,sd))

df.clean$togroup <- NULL#remove dummy tags
df.clean$Leucine <- NULL#leucine has lots of NA. so removed
df.log$Leucine <- NULL
df.log$togroup <- NULL
colnames <- colnames['Leucine' != colnames]

df.final <- t(scale(df.log))#for asca. autoscaling is done. data transposed for algorithm usage purpose

##############model#########
source('ASCA-genes.1.2.1/utils.R')#import all the functions in this folder
i=2 #treatment
j=5 #time
r=4 #replicates
time=c(0,2,6,10,24)

mx <- asca.design.matrix(i,j,r,time)#create design matrix. will prompt warning but is ok
Fac=c(1,1,2,2)#Number of components in each submodel. time, treatment, interaction and residual
              #if type is 1. the third value controls the number of components in submodel b.ab
type=1 #type 1 indicates submodel treatment(b) and interaction(ab) is analysed together.
       # type 2 means they will be analysed seperately

asca.fit <- ASCA.2f_leverage(X = t(df.final),Designa = mx$j, Designb = mx$i,type = type, Fac = Fac)
# ASCA.2F_leverage is the function for improved leverage. This function calls 
# ASCAfun12_leverage and PCA.GENE.unormed_loading function(in PCA-GENES script)
# ASCA.2f_leverage --> ASCAfun12_leverage --> PCA.GENE.unormed_loading

# Function for original leverage is called ASCA.2f. 
# ASCA.2f --> ASCAfun12 --> PCA.GENE

#Designa is time, Designb is treatment

permut_wrapped <- wrap.permutation.info(df.final = df.final,
                                        asca.fit = asca.fit,
                                        groups = colnames, #expressed patterns of features
                                        R = 100,#rounds of permutation
                                        which_leverage = ASCA.2f_leverage)#either ASCA.2f_leverage or
                                                                          #ASCA.2f
# this function wraps up the permutation for leverage(leverage.lims.R). 
# create leverage vs spe plots, loading plots and violin plots to view null distribution after permutation
# And some other data might useful.

#######visualize results########
plot(permut_wrapped$plots_for_submodel.leverage_spe)
plot(permut_wrapped$plots_for_submodel.loadings$PC1)
plot(permut_wrapped$Null_distribution_plot)
print(permut_wrapped$selected_features)