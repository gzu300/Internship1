#########################data##################
######
#simulation data
######
source('simulation_util.R')
source('MaSigPro_util.R')
source('ASCA-genes.1.2.1/sourceASCA.R')
setwd('../')
i=4#treatment
j=5#time
r=4#replicate
time=c(0,2,6,10,24)#time points
y1 <- 0.05*time-0.005*time^2
y1.t <- y1+4
y1.t1 <- y1+3
y2 <- -0.05*time+0.005*time^2
y2.t <- y2-1
y3 <- 0.05*time
y3.t <- y3-2

#total of 100 metabolites(columns)
#total of 2*5*4=40 rows

######
#generate data
######
tc.GENE <- function(n, r,
                    var11 = 0.2, var12 = 0.2, var13 = 0.2,var14 = 0.2,var15 = 0.2,
                    var21 = 0.2, var22 = 0.2, var23 = 0.2,var24 = 0.2,var25 = 0.2,
                    var31 = 0.2, var32 = 0.2, var33 = 0.2,var34 = 0.2,var35 = 0.2,
                    var41 = 0.2, var42 = 0.2, var43 = 0.2,var44 = 0.2,var45 = 0.2,
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
    Tr2 <- c(rnorm(r, a3, var31), rnorm(r, b3, var32), rnorm(r, c3, var33),rnorm(r, d3, var34),rnorm(r, e3, var35))  # Tr2 group
    Tr3 <- c(rnorm(r, a4, var41), rnorm(r, b4, var42), rnorm(r, c4, var43),rnorm(r, d4, var44),rnorm(r, e4, var45))  # Tr3 group
    gene <- c(Ctl, Tr1, Tr2,Tr3)
    tc.dat <- rbind(tc.dat, gene)
  }
  tc.dat
}
onediff <- tc.GENE(10,4,
                   #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                   a2=y1.t[1],b2=y1.t[2],c2=y1.t[3],d2=y1.t[4],e2=y1.t[5],
                   a3=2,b3=2,c3=2,d3=2,e3=2)
#a4=y2[1],b4=y2[2],c4=y2[3],d4=y2[4],e4=y2[5]
twodiff.2way <- tc.GENE(10,4,
                        #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                        a2=y1.t[1],b2=y1.t[2],c2=y1.t[3],d2=y1.t[4],e2=y1.t[5],
                        #a3=y2[1],b3=y2[2],c3=y2[3],d3=y2[4],e3=y2[5],
                        a4=y2[1],b4=y2[2],c4=y2[3],d4=y2[4],e4=y2[5])
difftime <- tc.GENE(10,4,
                    #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                    a2=0.2+y2[1],b2=0.2+y2[2],c2=0.2+y2[3],d2=0.2+y2[4],e2=0.2+y2[5],
                    a3=0.3+y2[1],b3=0.3+y2[2],c3=0.3+y2[3],d3=0.3+y2[4],e3=0.3+y2[5],
                    a4=0.4+y2[1],b4=0.4+y2[2],c4=0.4+y2[3],d4=0.4+y2[4],e4=0.4+y2[5],
                    var41 = 0.2, var42 = 0.2, var43 = 0.2,var44 = 0.2,var45 = 0.2)
small <- tc.GENE(1,4,a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                 a2=y1[1],b2=y1[2],c2=y1[3],d2=y1[4],e2=y1[5],
                 a3=y2[1],b3=y2[2],c3=y2[3],d3=y2[4],e3=y2[5],
                 a4=y2[1],b4=y2[2],c4=y2[3],d4=y2[4],e4=y2[5])
flat <- tc.GENE(100,4)
df.final <- rbind(onediff,twodiff.2way,difftime,small,flat)
rownames(df.final) <- c(paste('onediff',c(1:10)),paste('twodiff1',c(1:10)),paste('difftime',c(1:10)),'small',paste('flat',1:100))

#####################masigpro#############

######
#design matrix
######

masigpro.design <- design_matrix(i,j,r,time)
design <- make.design.matrix(edesign = masigpro.design,degree = 3)
dis <- design$dis
edesign <- design$edesign
colnames(df.final) <- rownames(edesign)

########
#fit model
########
masigpro.fit <- maSigPro(df.final,masigpro.design,degree=3,step.method = 'forward')
toplot <- rownames(masigpro.fit$sig.genes$treatment2vscontrol$sig.profiles)
#######
#plot
######
see.genes(masigpro.fit$sig.genes$treatment2vscontrol,show.fit = T,dis = edesign,k = 4)
for (each in c(toplot)){
  PlotGroups(df.final[rownames(df.final)==each,],edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector)
}
AA <- df.final["twodiff1 1",]
PlotGroups(AA,edesign = edesign, show.fit = T, dis = dis, groups.vector = design$groups.vector)
#######################asca-gene#############

######
#design matrix
#####
mx.i <- matrix(0,nrow = i*j*r, ncol = i, dimnames = list(c(),c('control','treatment1','treatment2','treatment3')))+c(rep(1,j*r),rep(0,i*j*r))
mx.j <- matrix(0,nrow = i*j*r, ncol = j, dimnames = list(c(),c('T0','T2','T6','T10','T24')))+c(rep(c(rep(1,r),rep(0,j*r-r)),i),rep(0,r))
mx.ij <- matrix(0,nrow = i*j*r, ncol = i*j, dimnames = list(c(),c(paste('inter',1:(i*j)))))+c(rep(1,r),rep(0,i*j*r))

######
#fit model
######
asca.fit <- ASCA.2f(X = t(df.final),Designa = mx.j, Designb = mx.i,type = 1, Fac = c(1,1,2,2))
######
#plot
######
lev.lim <- leverage.lims(df.final,R=10,FUN = ASCA.2f,Designa = mx.j, Designb = mx.i,Fac = c(2,2,2,2),type = 1,alpha = 0.05)$Cutoff[[2]]
spe.lim <- SPE.lims(my.asca = asca.fit,alpha = 0.05)[[2]]
leverage <- asca.fit$Model.bab$leverage
spe <- asca.fit$Model.bab$SPE
lev.spe.toplot <- data.frame(leverage=leverage,spe=spe)
lev.spe.toplot$metabolites <- rownames(df.final)
lev.spe.toplot$label <- c(rep('onediff',10),rep('twodiff',20),rep('small',1),rep('flat',100))

ggplot(data = lev.spe.toplot,aes(x=leverage,y=spe,color=label))+
  geom_point()+
  geom_hline(yintercept = spe.lim)+
  geom_vline(xintercept = lev.lim)

bab.toplot <- data.frame(score=asca.fit$Model.bab$scores,time=rep(time,4),treatments=rep(c('ctrl','treat1','treat2','treat3'),each=5))
ggplot(data = bab.toplot,aes(x=time,y=score.2,color=treatments))+
  geom_line()
