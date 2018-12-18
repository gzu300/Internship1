#########################data##################
########
#real data
########
setwd('../')
library(readxl)
colnames <- read_xlsx('macrophage%20polarization.xlsx',range = 'PeakArea!I1:BL1')
df <- read_xlsx('macrophage%20polarization.xlsx',range = 'PeakArea!I4:BL43',col_names = colnames(colnames))
sums <- apply(df,1,sum, na.rm=T)
df <- sweep(df*100,1,sums,"/")
i=2
j=5
r=4
time=c(0,2,6,10,24)
###exploration###
plot(apply(log(df),2,mean),apply(log(df),2,sd))
plot(apply(df,2,mean),apply(df,2,sd))
hist(df$Adenosine[1:19])
hist(log(df$Adenosine[1:19]))
###pretreatment####
library(tidyverse)
df$togroup <- rep(1:10,each=4)
df1 <- df %>% 
  group_by(togroup) %>% 
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
df1[1:4,"Fructose-1,6-bisphosphate"] <- df1$`Fructose-1,6-bisphosphate`[5]#fill NA with mean of replicates. F16B is missing values at entire time 1 control. filled with values at time 1 treatment
df1[20,] <- df1[19,]
df1$togroup <- NULL
df1$Leucine <- NULL#leucine has lots of NA. so removed
#df.log <- log(df1)#log transform
df.final <- t(scale(df1))

#####################masigpro#############

######
#design matrix
######
source('main.R')
masigpro.design <- design_matrix(i,j,r,time)
design <- make.design.matrix(edesign = masigpro.design,degree = 3)
dis <- design$dis
edesign <- design$edesign
colnames(df.final) <- rownames(edesign)

########
#fit model
########
masigpro.fit <- maSigPro(df.final,masigpro.design,degree=3)
print(masigpro.fit$summary)
toplot <- rownames(masigpro.fit$sig.genes$treatment2vstreatment1$sig.profiles)
#######
#plot
######
see.genes(masigpro.fit$sig.genes$treatment2vstreatment1,show.fit = T,dis = edesign,k = 3)
for (each in c(toplot)){
  PlotGroups(df.final[rownames(df.final)==each,],edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector,main = each)
}
AA <- df.final["Ascorbic Acid",]
PlotGroups(AA,edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector)
#######################asca-gene#############

######
#design matrix
#####
# mx.i <- matrix(0,nrow = i*j*r, ncol = i, dimnames = list(c(),c('control','treatment')))+c(rep(1,j*r),rep(0,i*j*r))
# mx.j <- matrix(0,nrow = i*j*r, ncol = j, dimnames = list(c(),c('T0','T2','T6','T10','T24')))+c(rep(c(rep(1,r),rep(0,j*r-r)),i),rep(0,r))
# mx.ij <- matrix(0,nrow = i*j*r, ncol = i*j, dimnames = list(c(),c(paste('inter',1:(i*j)))))+c(rep(1,r),rep(0,i*j*r))
mx <- asca.design.matrix(i,j,r,time)
######
#fit model
######
source('ASCA-genes.1.2.1/sourceASCA_mac.R')
asca.fit <- ASCA.2f(X = t(df.final),Designa = mx$j, Designb = mx$i,type = 1, Fac = c(1,1,1,2))
Fac=c(1,1,1,2)
type=1
######
#plot
######
#source('main.R')
# # lev.lim <- leverage.lims(df.final,R=10,FUN = ASCA.2f,Designa = mx.j, Designb = mx.i,Fac = c(1,1,3,2),type = 1,alpha = 0.05)$Cutoff[[2]]
# # spe.lim <- SPE.lims(my.asca = asca.fit,alpha = 0.05)[[2]]
# # leverage <- asca.fit$Model.bab$leverage
# # spe <- asca.fit$Model.bab$SPE
# # lev.spe.toplot <- data.frame(leverage=leverage,spe=spe)
# # lev.spe.toplot$metabolites <- rownames(df.final)
# 
# ggplot(data = lev.spe.toplot,aes(x=leverage,y=spe))+
#   geom_point()+
#   geom_hline(yintercept = spe.lim)+
#   geom_vline(xintercept = lev.lim)+
#   geom_text(aes(label=metabolites))

plot.leverage_spe(df.final,asca.fit,colnames)
plot.submodels_score(asca.fit,i,j)
plot.submodels_loading(asca.fit,colnames)
