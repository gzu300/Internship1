#########################data##################
######
#simulation data
######
library(tidyverse)
library(purrr)
#source('simulation_util.R')
source('MaSigPro_util.R')
source('ASCA-genes.1.2.1/sourceASCA.R')
setwd('../')



######
#generate data
######
generate_data.nocorrelation <- function(p1,p2,p3,neg){
  #6 time points;2 treatments;4 replicates
  #no correlation between metabolites
  tc.GENE <- function(n, r,
                      var11 = 0.2, var21 = 0.2, var31 = 0.2,var41 = 0.2,var51 = 0.2,
                      var12 = 0.2, var22 = 0.2, var32 = 0.2,var42 = 0.2,var52 = 0.2,
                      var13 = 0.2, var23 = 0.2, var33 = 0.2,var43 = 0.2,var53 = 0.2,
                      var14 = 0.2, var24 = 0.2, var34 = 0.2,var44 = 0.2,var54 = 0.2,
                      var15 = 0.2, var25 = 0.2, var35 = 0.2,var45 = 0.2,var55 = 0.2,
                      var16 = 0.2, var26 = 0.2, var36 = 0.2,var46 = 0.2,var56 = 0.2,
                      a1 = 0,      a2 = 0,      a3 = 0,     a4 = 0,
                      b1 = 0,      b2 = 0,      b3 = 0,     b4 = 0,
                      c1 = 0,      c2 = 0,      c3 = 0,     c4 = 0,
                      d1 = 0,      d2 = 0,      d3 = 0,     d4 = 0,
                      e1 = 0,      e2 = 0,      e3 = 0,     e4 = 0,
                      f1 = 0,      f2 = 0,      f3 = 0,     f4 = 0)
  {
    
    tc.dat <- NULL
    for (i in 1:n) {
      #Ctl <- c(rnorm(r, a1, var11), rnorm(r, b1, var12), rnorm(r, c1, var13),rnorm(r, d1, var14),rnorm(r, e1, var15))  # Ctl group
      #Tr1 <- c(rnorm(r, a2, var21), rnorm(r, b2, var22), rnorm(r, c2, var23),rnorm(r, d2, var24),rnorm(r, e2, var25),rnorm(r,f2,var26))  # Tr1 group
      Tr2 <- c(rnorm(r, a3, var31), rnorm(r, b3, var32), rnorm(r, c3, var33),rnorm(r, d3, var34),rnorm(r, e3, var35),rnorm(r,f3,var36))  # Tr2 group
      Tr3 <- c(rnorm(r, a4, var41), rnorm(r, b4, var42), rnorm(r, c4, var43),rnorm(r, d4, var44),rnorm(r, e4, var45),rnorm(r, f4,var46))  # Tr3 group
      gene <- c(Tr2, Tr3)
      tc.dat <- rbind(tc.dat, gene)
    }
    tc.dat
  }
  one <- tc.GENE(p1,4,
                 #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                 #a2=y1.t[1],b2=y1.t[2],c2=y1.t[3],d2=y1.t[4],e2=y1.t[5],
                 a3=y1[1],b3=y1[2],c3=y1[3],d3=y1[4],e3=y1[5],f3=y1[6],
                 a4=y12[1],b4=y12[2],c4=y12[3],d4=y12[4],e4=y12[5],f4=y12[6])
  two <- tc.GENE(p2,4,
                 #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                 #a2=y1.t[1],b2=y1.t[2],c2=y1.t[3],d2=y1.t[4],e2=y1.t[5],
                 a3=y2[1],b3=y2[2],c3=y2[3],d3=y2[4],e3=y2[5],
                 a4=y22[1],b4=y22[2],c4=y22[3],d4=y22[4],e4=y22[5],f4=y22[6])
  three <- tc.GENE(p3,4,
                   #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                   #a2=0.2+y2[1],b2=0.2+y2[2],c2=0.2+y2[3],d2=0.2+y2[4],e2=0.2+y2[5],
                   a3=y3[1],b3=y3[2],c3=y3[3],d3=y3[4],e3=y3[5],f3=y3[6],
                   a4=y32[1],b4=y32[2],c4=y32[3],d4=y32[4],e4=y32[5],f4=y32[6])
  #var41 = 0.2, var42 = 0.2, var43 = 0.2,var44 = 0.2,var45 = 0.2)
  # four <- tc.GENE(1,4,a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
  #                  a2=y1[1],b2=y1[2],c2=y1[3],d2=y1[4],e2=y1[5],
  #                  a3=y2[1],b3=y2[2],c3=y2[3],d3=y2[4],e3=y2[5],
  #                  a4=y2[1],b4=y2[2],c4=y2[3],d4=y2[4],e4=y2[5])
  flat <- tc.GENE(neg,4)
  df.final <- rbind(one,two,three,flat)
  groups <- c(rep('one',nrow(one)),rep('two',nrow(two)),rep('three',nrow(three))
              ,rep('flat',nrow(flat)))
  rownames(df.final) <- c(paste('one',1:nrow(one)),paste('two',1:nrow(two)),paste('three',1:nrow(three))
                          ,paste('flat',1:nrow(flat)))
  output <- vector(mode='list',length = 2)
  names(output) <- c('df','groups')
  output[[1]] <- df.final
  output[[2]] <- groups
  output
}

generate_data.same_var <- function(p1,p2,p3,neg,var){
  #6 time points;2 treatments;4 replicates
  #no correlation between metabolites
  tc.GENE <- function(n, r,
                      var,
                      a1 = 0,      a2 = 0,      a3 = 0,     a4 = 0,
                      b1 = 0,      b2 = 0,      b3 = 0,     b4 = 0,
                      c1 = 0,      c2 = 0,      c3 = 0,     c4 = 0,
                      d1 = 0,      d2 = 0,      d3 = 0,     d4 = 0,
                      e1 = 0,      e2 = 0,      e3 = 0,     e4 = 0,
                      f1 = 0,      f2 = 0,      f3 = 0,     f4 = 0)
  {
    
    tc.dat <- NULL
    for (i in 1:n) {
      #Ctl <- c(rnorm(r, a1, var11), rnorm(r, b1, var12), rnorm(r, c1, var13),rnorm(r, d1, var14),rnorm(r, e1, var15))  # Ctl group
      #Tr1 <- c(rnorm(r, a2, var21), rnorm(r, b2, var22), rnorm(r, c2, var23),rnorm(r, d2, var24),rnorm(r, e2, var25),rnorm(r,f2,var26))  # Tr1 group
      Tr2 <- c(rnorm(r, a3, var), rnorm(r, b3, var), rnorm(r, c3, var),rnorm(r, d3, var),rnorm(r, e3, var),rnorm(r,f3,var))  # Tr2 group
      Tr3 <- c(rnorm(r, a4, var), rnorm(r, b4, var), rnorm(r, c4, var),rnorm(r, d4, var),rnorm(r, e4, var),rnorm(r, f4,var))  # Tr3 group
      gene <- c(Tr2, Tr3)
      tc.dat <- rbind(tc.dat, gene)
    }
    tc.dat
  }
  one <- tc.GENE(p1,4,
                 #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                 #a2=y1.t[1],b2=y1.t[2],c2=y1.t[3],d2=y1.t[4],e2=y1.t[5],
                 a3=y1[1],b3=y1[2],c3=y1[3],d3=y1[4],e3=y1[5],f3=y1[6],
                 a4=y12[1],b4=y12[2],c4=y12[3],d4=y12[4],e4=y12[5],f4=y12[6],
                 var=var)
  two <- tc.GENE(p2,4,
                 #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                 #a2=y1.t[1],b2=y1.t[2],c2=y1.t[3],d2=y1.t[4],e2=y1.t[5],
                 a3=y2[1],b3=y2[2],c3=y2[3],d3=y2[4],e3=y2[5],
                 a4=y22[1],b4=y22[2],c4=y22[3],d4=y22[4],e4=y22[5],f4=y22[6],
                 var=var)
  three <- tc.GENE(p3,4,
                   #a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
                   #a2=0.2+y2[1],b2=0.2+y2[2],c2=0.2+y2[3],d2=0.2+y2[4],e2=0.2+y2[5],
                   a3=y3[1],b3=y3[2],c3=y3[3],d3=y3[4],e3=y3[5],f3=y3[6],
                   a4=y32[1],b4=y32[2],c4=y32[3],d4=y32[4],e4=y32[5],f4=y32[6],
                   var=var)
  #var41 = 0.2, var42 = 0.2, var43 = 0.2,var44 = 0.2,var45 = 0.2)
  # four <- tc.GENE(1,4,a1=y1[1],b1=y1[2],c1=y1[3],d1=y1[4],e1=y1[5],
  #                  a2=y1[1],b2=y1[2],c2=y1[3],d2=y1[4],e2=y1[5],
  #                  a3=y2[1],b3=y2[2],c3=y2[3],d3=y2[4],e3=y2[5],
  #                  a4=y2[1],b4=y2[2],c4=y2[3],d4=y2[4],e4=y2[5])
  flat <- tc.GENE(neg,4,var=var)
  df.final <- rbind(one,two,three,flat)
  groups <- c(rep('one',nrow(one)),rep('two',nrow(two)),rep('three',nrow(three))
              ,rep('flat',nrow(flat)))
  rownames(df.final) <- c(paste('one',1:nrow(one)),paste('two',1:nrow(two)),paste('three',1:nrow(three))
                          ,paste('flat',1:nrow(flat)))
  output <- vector(mode='list',length = 2)
  names(output) <- c('df','groups')
  output[[1]] <- df.final
  output[[2]] <- groups
  output
}

create.simulation <- function(pat_names.list,ftrs_in_pat.list,replicates,sd,...){
  dimnames <- pat_names.list %>% 
    map2(.,ftrs_in_pat.list,~rep(.x,.y)) %>% 
    flatten_chr() %>% 
    list(NULL,.)
  
  #create the simulation matrix
  wrap_tre <- function(dimnames,replicates,sd,...){
    trend_list <- list(...)
    trend_list %>% 
      map2(ftrs_in_pat.list,.,~rep(.y,.x)) %>% 
      flatten() %>% 
      map(.,~rnorm(n = replicates,mean = .x,sd = sd)) %>% 
      flatten_dbl() %>% 
      matrix(.,nrow = length(trend_list[[1]])*replicates,dimnames = dimnames) %>% 
      t(.)
  }
  wrap_tre <- wrap_tre(dimnames = dimnames,replicates = rep,sd = sd,trend1,trend2,trend3,trend4)
  output <- list(df=wrap_tre,groups=dimnames[[2]])
  output
}

#######################asca-gene#############

######
#design matrix
#####
##functions
asca.design.matrix <- function(i,j,r,time){
  #time: a vector consists of time points in experiment design.
  mx <- vector(mode = 'list')
  mx[[1]] <- matrix(0,nrow = i*j*r, ncol = i, dimnames = list(c(),paste('treatment',1:i)))+c(rep(1,j*r),rep(0,i*j*r))
  mx[[2]] <- matrix(0,nrow = i*j*r, ncol = j, dimnames = list(c(),paste('T',time)))+c(rep(c(rep(1,r),rep(0,j*r-r)),i),rep(0,r))
  mx[[3]] <- matrix(0,nrow = i*j*r, ncol = i*j, dimnames = list(c(),c(paste('inter',1:(i*j)))))+c(rep(1,r),rep(0,i*j*r))
  names(mx) <- c('i','j','ij')
  mx
}

######
#plot
######
#leveragevsspe
plot.leverage_spe_original <- function(df.final,asca.fit,groups, R=1,No.sample=NULL){
  ###########
  #OBLOSETED#
  ###########
  #original leverage normalised to 1
  lev.lim <- leverage.lims(df.final,R=R,FUN = ASCA.2f_leverage,Designa = mx$j, Designb = mx$i,Fac = Fac,type = type,alpha = 0.05,showvar = F, showscree = F)$Cutoff[[2]]
  spe.lim <- SPE.lims(my.asca = asca.fit,alpha = 0.01)[[2]]
  leverage <- asca.fit$Model.bab$leverage
  spe <- asca.fit$Model.bab$SPE
  lev.spe.toplot <- data.frame(leverage=leverage,spe=spe)
  lev.spe.toplot$metabolites <- rownames(df.final)
  lev.spe.toplot$groups <- groups
  plot <- ggplot(data = lev.spe.toplot,aes(x=leverage,y=spe,color=groups))+
    geom_point()+
    geom_hline(yintercept = spe.lim)+
    geom_vline(xintercept = lev.lim)+
    #geom_text(aes(label=c(1,2,3,rep('',100)),hjust=-1.2))+
    labs(title = paste('original leverage and SPE with',ncol(asca.fit$Model.bab$scores),'PCs'))+
    theme(legend.text = element_text(size=5),legend.key.size = unit(0.1,'cm'))
  print(plot)
  output <- groups[leverage>lev.lim]
  output
  #lev.spe.toplot
}


######fitted asca_gene data frame######
fitted <- function(df.final,asca.fit,groups,R){
  #attention: df.final. rows are features(eg.metabolites), columns are samples.
  #groups is a vector of strings or numbers specifies patterns of variables. 
  #it needs to be the same length as variables. same pattern gives the same tag.
  #R is the number of permutations
  #which_leverage argument: improved leverage is 'ASCA.2f_leverage'. original leverage is 'ASCA.2f'
  
  ##original leverage loading normalised to 1
  ##improved leverage scores normalised to 1
  lev.lim <- leverage.lims(df.final,R=R,FUN = which_leverage,Designa = mx$j, Designb = mx$i,Fac = Fac,type = type,alpha = 0.05,showvar = F, showscree = F)$Cutoff[[2]]
  spe.lim <- SPE.lims(my.asca = asca.fit,alpha = 0.01)[[2]]
  asca.fit_leverage <- ASCA.2f_leverage(t(df.final),Designa = mx$j, Designb = mx$i,Fac = Fac,type = type,showvar = F, showscree = F)
  
  leverage <- asca.fit_leverage$Model.bab$leverage
  leverage.original <- asca.fit$Model.bab$leverage
  spe <- asca.fit$Model.bab$SPE
  #assemble datafram for ggplot2
  lev.spe.toplot <- data.frame(leverage=leverage,
                               leverage.original=leverage.original,
                               spe=spe,
                               metabolites=1:nrow(df.final),
                               groups=groups,
                               predicted=(leverage>lev.lim),
                               truth=(groups != 'flat'),
                               score=asca.fit$Model.bab$scores
                               )
  output <- list(lev.lim,spe.lim,lev.spe.toplot)
  names(output) <- c('lev_limit','spe_lim','stats_for_plot')
  output
}
########stats###########
stats <- function(fitted.data){
  FP <- sum(fitted.data$predicted&!fitted.data$truth)
  FN <- sum(!fitted.data$predicted&fitted.data$truth)
  TP <- sum(fitted.data$predicted&fitted.data$truth)
  TN <- sum(!fitted.data$predicted&!fitted.data$truth)
  output <- list(TP,TN,FP,FN)
  names(output) <- c('TP','TN','FP','FN')
  output
}

#####plots########

#plot.leverage_spe(df.final,asca.fit, groups)
plot.leverage_spe <- function(asca.fitted, sd=sd,title = paste('improved leverage and SPE with',Fac[3],'PCs')){
  plot <- ggplot(data = asca.fitted$stats_for_plot,aes(x=leverage,y=spe,color=groups))+
    geom_point()+
    geom_hline(yintercept = asca.fitted$spe_lim)+
    geom_vline(xintercept = asca.fitted$lev_limit)+
    #geom_text(aes(label=c(1,2,3,rep('',100)),hjust=-1.2))+
    labs(title = title)+
    theme(legend.text = element_text(size=5),legend.key.size = unit(0.1,'cm'))
  plot
}
#score plot
plot.submodels_score <- function(asca.fit,i,j,title = paste('score plot for submodel b.ab')){
  #plot score vs time of all the PCs for submodel b.ab
  scores <- asca.fit$Model.bab$scores
  PCs <- ncol(scores)
  bab.toplot <- data.frame(scores=scores,time=rep(time,i),treatments=rep(paste('treatment',1:i),each=j))
  output <- as.vector(1:PCs,mode = 'list')
  for (each in 1:PCs){
    plot <- ggplot(data = bab.toplot,aes(x=time,y=bab.toplot[[each]],color=treatments))+
      geom_line()+
      ylab(label = paste('PC',each))+
      labs(title = title)
    output[[each]] <- plot
  }
  output
}
#plot.submodels_score(asca.fit,i)
  #loading plot
plot.submodels_loading <- function(asca.fit,groups=NULL,title = paste('loading plot for submodel b.ab')){
  #plot loadings for all the PCs
  #groups is a vector of strings or numbers specifies patterns of variables. 
  #it needs to be the same length as variables. same pattern gives the same tag.
  bab.loadings <- data.frame(loading=asca.fit$Model.bab$loadings)
  PCs <- ncol(bab.loadings)
  bab.loadings$metabolites <- 1:nrow(bab.loadings)
  bab.loadings$groups <- groups
  output <- as.vector(1:PCs,mode = 'list')
  for (each in 1:PCs){
    plot <- ggplot(bab.loadings,aes(x=metabolites,y=bab.loadings[[each]],fill=groups))+
      geom_bar(stat = 'identity',position = 'dodge')+
      ylab(paste('PC',each))+
      labs(title = title)+
      theme(legend.text = element_text(size=5),legend.key.size = unit(0.1,'cm'))
    output[[each]] <- plot
  }
  output
}

plot.submodels <- function(asca.fitted,asca.fit, Fac=Fac...){
  output <- list(leverage_spe=plot.leverage_spe(asca.fitted = asca.fitted),
                 scores=plot.submodels_score(asca.fit,i,j),
                 loadings=plot.submodels_loading(asca.fit,groups=groups))
  output
}

