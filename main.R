#########################data##################
######
#simulation data
######
source('simulation_util.R')
source('MaSigPro_util.R')
source('ASCA-genes.1.2.1/sourceASCA.R')
setwd('../')
i=2#treatment
j=5#time
r=4#replicate
time=c(0,2,6,10,24)#time points
y=time+time^2

#total of 100 metabolites(columns)
#total of 2*5*4=40 rows

######
#generate data
######
tc.GENE <- function(n, r,
                    var11 = 0.01, var12 = 0.01,var13 = 0.01,var14 = 0.01,var15 = 0.01,
                    var21 = 0.01, var22 = 0.01, var23 =0.01,var24 = 0.01,var25 = 0.01,
                    #var31 = 0.01, var32 = 0.01, var33 = 0.01,var34 = 0.01,var35 = 0.01,
                    #var41 = 0.01, var42 = 0.01, var43 = 0.01,var44 = 0.01,var45 = 0.01,
                    #var51 = 0.01, var52 = 0.01, var53 = 0.01,var54 = 0.01,var55 = 0.01,
                    a1 = 0, a2 = 0, a3 = 0, a4 = 0,
                    b1 = 0, b2 = 0, b3 = 0, b4 = 0,
                    c1 = 0, c2 = 0, c3 = 0, c4 = 0,
                    d1 = 0, d2 = 0, d3 = 0, d4 = 0,
                    e1 = 0, e2 = 0, e3 = 0, e4 = 0)
{
  
  tc.dat <- NULL
  for (i in 1:n) {
    Ctl <- c(rnorm(r, a1, var11), rnorm(r, b1, var12), rnorm(r, c1, var13),rnorm(r, d1, var14),rnorm(r, e1, var15))  # Ctl group
    Tr1 <- c(rnorm(r, a2, var21), rnorm(r, b2, var22), rnorm(r, c2, var23),rnorm(r, d2, var24),rnorm(r, e2, var25))  # Tr1 group
    #Tr2 <- c(rnorm(r, a3, var31), rnorm(r, b3, var32), rnorm(r, c3, var33))  # Tr2 group
    #Tr3 <- c(rnorm(r, a4, var41), rnorm(r, b4, var42), rnorm(r, c4, var43))  # Tr3 group
    gene <- c(Ctl, Tr1)
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat
}
trend1 <- tc.GENE(2,4,a2=y[1],b2=y[2],c2=y[3],d2=y[4],e2=y[5])
flat <- tc.GENE(10,4)
df.final <- rbind(flat,trend1)
rownames(df.final) <- paste('gene',c(1:12))

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
AA <- df.final["gene 12",]
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
asca.fit <- ASCA.2f(X = t(df.final),Designa = mx.j, Designb = mx.i,type = 2, Fac = c(1,1,1,1))
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
  
