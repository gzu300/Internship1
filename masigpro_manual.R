################
#This is the main spread created on 16/10/2018#
################

source('Google Drive/BDA_internship/data/simulation_util.R')
source('Google Drive/BDA_internship/data/MaSigPro_util.R')
set.seed(999)


i=5
j=4
r=4
v <- list(list(c(2,4,6,8)))
df_manual <- pipeline(v, i,j,r)
df <- df_manual

#design matrix
design.matrix <- design_matrix(i, j, r, df)
design <- make.design.matrix(design.matrix, degree = 3)
a <- maSigPro(df, design.matrix, degree = 3, k = 4)

see.genes(a$sig.genes$G3vsG1, design.matrix, k = 4)


####################
#step wise masigpro#
####################
design <- make.design.matrix(design.matrix, degree = 3)

fit <- p.vector(df, design)
fit$dat
tstep <- T.fit(fit, step.method = 'forward')
tstep$coefficients

sigs <- get.siggenes(tstep, vars = 'groups', rsq = 0.7)
sigs$sig.genes$G1$sig.profiles

c <- see.genes(sigs$sig.genes$G2vsG1, dis = design$dis, k=4)
c$cut
x1 <- df[rownames(df)=='x1']
PlotGroups(x1, edesign = design.matrix, dis = design$dis)

###############
#plot raw data#
###############
df.raw <- pipeline.raw(v,i,j,r)
for (each in 3:ncol(df.raw)){
  name <- colnames(df.raw[each])
  plt <- ggplot(data = data.frame(df.raw), aes(x = Time, y = df.raw[, each], colour = Group, group = Group))+
    geom_point()+
    stat_summary(fun.y = 'median',geom = 'line')+
    ylab(name)
  print(plt)
  #dev.copy(jpeg, paste(name, '.jpeg'))
  #dev.off()
}
