#########################data##################
########
#real data
########
library(readxl)
colnames <- read_xlsx('macrophage%20polarization.xlsx',range = 'PeakArea!I1:BL1')
df <- read_xlsx('macrophage%20polarization.xlsx',range = 'PeakArea!I4:BL43',col_names = colnames(colnames))
rm(colnames)

###exploration###
plot(apply(log(df),2,mean),apply(log(df),2,sd))
hist(df$met1)
hist(log(df$met1))
###pretreatment####
library(tidyverse)
df$togroup <- rep(1:10,each=4)
df1 <- df %>% 
  group_by(togroup) %>% 
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
df1[1:4,"Fructose-1,6-bisphosphate"] <- df1$`Fructose-1,6-bisphosphate`[5]#fill NA with mean of replicates. F16B is missing values at entire time 1 control. filled with values at time 1 treatment
df1$togroup <- NULL
df1$Leucine <- NULL
df.log <- log(df1)#log transform
df.final <- t(scale(df.log))

#####################masigpro#############

######
#design matrix
######

masigpro.design <- design_matrix(i,j,r,time)
design <- make.design.matrix(edesign = masigpro.design,degree = 2)
dis <- design$dis
edesign <- design$edesign
colnames(df.final) <- rownames(edesign)

########
#fit model
########
masigpro.fit <- maSigPro(df.final,masigpro.design,degree=3)
toplot <- rownames(masigpro.fit$sig.genes$treatment2vscontrol$sig.profiles)
#######
#plot
######
see.genes(masigpro.fit$sig.genes$treatment2vscontrol,show.fit = T,dis = edesign,k = 3)
for (each in c(toplot)){
  PlotGroups(df.final[rownames(df.final)==each,],edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector)
}
AA <- df.final["Ascorbic Acid",]
PlotGroups(AA,edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector)
#######################asca-gene#############

######
#design matrix
#####
mx.i <- matrix(0,nrow = i*j*r, ncol = i, dimnames = list(c(),c('control','treatment')))+c(rep(1,j*r),rep(0,i*j*r))
mx.j <- matrix(0,nrow = i*j*r, ncol = j, dimnames = list(c(),c('T0','T2','T6','T10','T24')))+c(rep(c(rep(1,r),rep(0,j*r-r)),i),rep(0,r))
mx.ij <- matrix(0,nrow = i*j*r, ncol = i*j, dimnames = list(c(),c(paste('inter',1:(i*j)))))+c(rep(1,r),rep(0,i*j*r))

######
#fit model
######
asca.fit <- ASCA.2f(X = t(df.final),Designa = mx.j, Designb = mx.i,type = 1, Fac = c(3,2,2,2))
######
#plot
######
lev.lim <- leverage.lims(df.final,R=10,FUN = ASCA.2f,Designa = mx.j, Designb = mx.i,Fac = c(3,2,2,2),type = 1,alpha = 0.05)$Cutoff[[2]]
spe.lim <- SPE.lims(my.asca = asca.fit,alpha = 0.05)[[2]]
leverage <- asca.fit$Model.bab$leverage
spe <- asca.fit$Model.bab$SPE
lev.spe.toplot <- data.frame(leverage=leverage,spe=spe)
lev.spe.toplot$metabolites <- rownames(df.final)

ggplot(data = lev.spe.toplot,aes(x=leverage,y=spe))+
  geom_point()+
  geom_hline(yintercept = spe.lim)+
  geom_vline(xintercept = lev.lim)+
  geom_text(aes(label=metabolites))