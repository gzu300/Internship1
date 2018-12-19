#########################data##################
########
#real data
########
library(readxl)
colnames <- read_xlsx('macrophage%20polarization.xlsx',range = 'PeakArea!I1:BL1')
df <- read_xlsx('macrophage%20polarization.xlsx',range = 'PeakArea!I4:BL43',col_names = colnames(colnames))
rm(colnames)

# exploration----------------
plot(apply(log(df),2,mean),apply(log(df),2,sd))
plot(apply(df,2,mean),apply(df,2,sd))
plot(apply(sqrt(df),2,mean),apply(sqrt(df),2,sd))
hist(df$Adenosine[1:19])
hist(log(df$Adenosine[1:19]))
hist(sqrt(df$Adenosine[1:19]))

# pretreatment--------------------------------

library(tidyverse)
sums <- apply(df,MARGIN = 1,sum,na.rm=T)
df <- sweep(df,MARGIN = 1,STATS = sums,FUN = "/")
i=2
j=5
r=4
time=c(0,2,6,10,24)
df$togroup <- rep(1:10,each=4)
df1 <- df %>% 
  group_by(togroup) %>% 
  mutate_all(~ifelse(is.na(.x), mean(.x, na.rm = TRUE), .x))
df1[1:4,"Fructose-1,6-bisphosphate"] <- df1$`Fructose-1,6-bisphosphate`[5]#fill NA with mean of replicates. F16B is missing values at entire time 1 control. filled with values at time 1 treatment
df1[20,] <- df1[19,]
df1$togroup <- NULL
df1$Leucine <- NULL#leucine has lots of NA. so removed
df.log <- log(df1)#log transform
df.final.masigpro <- t(scale(df.log,center = F,scale = apply(df.log, 2, sd, na.rm = TRUE)))
df.final <- t(scale(df.log))
plot(apply(df.log,2,mean),apply(df.log,2,sd))
colnames <- rownames(df.final)
# masigpro ------------------------------

source('main.R')
masigpro.design <- design_matrix(i,j,r,time)
design <- make.design.matrix(edesign = masigpro.design,degree = 3)
dis <- design$dis
edesign <- design$edesign
colnames(df.final.masigpro) <- rownames(edesign)

masigpro.fit <- maSigPro(df.final.masigpro,masigpro.design,degree=3)
print(masigpro.fit$summary)
toplot <- rownames(masigpro.fit$sig.genes$treatment2vstreatment1$sig.profiles)
# plot--------------------------

see.genes(masigpro.fit$sig.genes$treatment2vstreatment1,show.fit = T,dis = edesign,k = 3)
for (each in c(toplot)){
  PlotGroups(df.final[rownames(df.final)==each,],edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector,main = each)
}
AA <- df.final.masigpro["Succinic Acid",]
PlotGroups(AA,edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector)
# asca-gene---------------------

mx <- asca.design.matrix(i,j,r,time)

source('ASCA-genes.1.2.1/sourceASCA_mac.R')
asca.fit <- ASCA.2f(X = t(df.final),Designa = mx$j, Designb = mx$i,type = 1, Fac = c(1,1,4,2))
Fac=c(1,1,4,2)
type=1

# asca section 2 ----------------------------------------------------------
asca.fit1 <- ASCA.2f(X = t(df.final),Designa = mx$j, Designb = mx$i,type = 2, Fac = c(1,1,4,2))
barplot(t(asca.fit1$Model.b$loadings))
# plot------------------------------------------
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

asca.selected <- plot.leverage_spe(df.final,asca.fit,colnames,20)
plot.submodels_score(asca.fit,i,j)
plot.submodels_loading(asca.fit,colnames)
matplot(asca.fit$Model.bab$data-asca.fit$Model.bab$scores%*%t(asca.fit$Model.bab$loadings),type = 'l')
