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
  wrap_tre <- wrap_tre(dimnames = dimnames,replicates = rep,sd = sd,...)
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

######fitted asca_gene data frame######
wrap.permutation.info <- function(df.final,asca.fit,groups,R,which_leverage,alpha=0.05,...){
  #attention: df.final. rows are features(eg.metabolites), columns are samples.
  #groups is a vector of strings or numbers specifies patterns of variables. for real data, feature names could be filled in
  #it needs to be the same length as variables. same pattern gives the same tag.
  #R is the number of permutations
  #which_leverage argument: improved leverage is 'ASCA.2f_leverage'. original leverage is 'ASCA.2f'. without the quote in argument
  
  ##original leverage loading normalised to 1
  ##improved leverage scores normalised to 1
  
  #...are the arguments to be passed into 'labels' layer of ggplot in plot.submodel function and
  #its downstream functions. such as: 'title', 'tag', etc
  
  permutated.data <- leverage.lims(df.final,R=R,FUN = which_leverage,Designa = mx$j, Designb = mx$i,Fac = Fac,type = type,alpha = alpha,showvar = F, showscree = F)
  permu.data <- permutated.data$NullDistribution$Model.bab
  lev.lim <- permutated.data$Cutoff[[2]]
  spe.lim <- SPE.lims(my.asca = asca.fit,alpha = alpha)[[2]]
  
  leverage <- asca.fit$Model.bab$leverage
  spe <- asca.fit$Model.bab$SPE
  #assemble dataframe for leverage vs spe as well as stats for prediction accuracy. FP, FN... can be calculated in stats function below
  lev.spe.toplot <- data.frame(leverage=leverage,
                               spe=spe,
                               metabolites=1:nrow(df.final),
                               patterns=groups,
                               predicted=(leverage>lev.lim),
                               truth=(groups != 'flat')
                               )
  selected_features <- sort(groups[leverage>lev.lim])
  Nulldist_plot <- plot.NullDistribution(permu.data,asca.fit$Model.bab$leverage,groups,lev.lim,R,alpha,size=6)
  output <- list(lev_limit=lev.lim, spe_lim=spe.lim, stats_for_plot=lev.spe.toplot, selected_features=selected_features, Null_distribution_plot=Nulldist_plot)
  submodel_plots <- plot.submodels(output,asca.fit, groups=groups, Fac=Fac, size=6,...)
  output <- c(output, plots_for_submodel=submodel_plots)
  output
}

########stats###########
stats <- function(permut_wrapped){
  FP <- sum(permut_wrapped$predicted&!permut_wrapped$truth)
  FN <- sum(!permut_wrapped$predicted&permut_wrapped$truth)
  TP <- sum(permut_wrapped$predicted&permut_wrapped$truth)
  TN <- sum(!permut_wrapped$predicted&!permut_wrapped$truth)
  output <- list(TP,TN,FP,FN)
  names(output) <- c('TP','TN','FP','FN')
  output
}

#####plots########

#plot.leverage_spe(df.final,asca.fit, groups)
plot.leverage_spe <- function(permut_wrapped, sd=sd,title = paste('improved leverage and SPE with',Fac[3],'PCs'),size,...){
  plot <- ggplot(data = permut_wrapped$stats_for_plot,aes(x=leverage,y=spe,color=patterns))+
    geom_point()+
    geom_hline(yintercept = permut_wrapped$spe_lim)+
    geom_vline(xintercept = permut_wrapped$lev_limit)+
    labs(title = title,...)+
    theme(axis.title = element_text(size=size),panel.background = element_blank(),legend.title = element_text(size=size),legend.text = element_text(size=3),legend.key.size = unit(0.08,'cm'), plot.title = element_text(size=size))
  plot
}

#score plot
plot.submodels_score <- function(asca.fit,i,j,title = paste('score plot for submodel b.ab')){
  #plot score vs time of all the PCs for submodel b.ab
  scores <- asca.fit$Model.bab$scores
  PCs <- 1:ncol(scores)
  bab.toplot <- data.frame(scores=scores,
                           time=rep(time,i),
                           treatments=rep(paste('treatment',1:i),each=j))
  output <- as.vector(PCs,mode = 'list') %>% 
    set_names(paste('PC',PCs,sep = ''))
  for (each in PCs){
    plot <- ggplot(data = bab.toplot,aes(x=time,y=bab.toplot[[each]],color=treatments))+
      geom_line()+
      ylab(label = paste('PC',each))+
      labs(title = title)
    output[[each]] <- plot
  }
  output
}

#loading plot
plot.submodels_loading <- function(asca.fit,groups=NULL,title = paste('loading plot for submodel b.ab'),size=10,...){
  #plot loadings for all the PCs
  #groups is a vector of strings or numbers specifies patterns of variables. 
  #it needs to be the same length as variables. same pattern gives the same tag.
  bab.loadings <- data.frame(loading=asca.fit$Model.bab$loadings,
                             metabolites=1:length(groups),
                             groups = groups)
  PCs <- 1:ncol(asca.fit$Model.bab$loadings)
  output <- as.vector(PCs,mode = 'list') %>% 
    set_names(paste('PC',PCs,sep = ''))
  for (each in PCs){
    plot <- ggplot(bab.loadings,aes(x=metabolites,fill=groups))+
      geom_col(aes_string(y = colnames(bab.loadings)[each]))+
      ylab(paste('PC',each))+
      labs(title = title,...)+
      theme(axis.title = element_text(size=size),panel.background = element_blank(),legend.title = element_text(size=size),legend.text = element_text(size=3),legend.key.size = unit(0.08,'cm'), plot.title = element_text(size=size))
    output[[each]] <- plot
  }
  output
}

plot.submodels <- function(permut_wrapped,asca.fit, Fac=Fac, groups,size=10,...){
  output <- list(leverage_spe=plot.leverage_spe(permut_wrapped = permut_wrapped,size = size,...),
                 scores=plot.submodels_score(asca.fit,i,j),
                 loadings=plot.submodels_loading(asca.fit,groups=groups,size = size,...))
  output
}

plot_metabolites <- function(df,range,...){
  df.final <- df
  df.toplot <- data.frame(t(df.final[range,]))
  df.toplot$time <- rep(rep(time,each=r),i)
  df.toplot$treatment <- factor(rep(1:2,each=r*j))
  a <- df.toplot %>% 
    gather(key = metabolites,value = value,1:length(range))
  
  output <- ggplot(a,aes(x=time,y=value,color=treatment))+
    geom_point(size=0.5)+
    stat_summary(fun.y = mean,geom = 'line')+
    theme(legend.title = element_text(size = 5),panel.background = element_blank())+
    facet_wrap(metabolites~., scales = 'free')+
    geom_smooth(method = 'lm', formula = y~poly(x,2),se = F,linetype = '3313')+
    labs(...)
  output
}
plot_a_metabolite <- function(df,FUN,which,size=6,...){
  df.final <- df
  trend.toplot <- data.frame(replicate=rep(1:(i*j),each=r),time=rep(time,each=r),treatment=factor(rep(1:i,each=j*r)),metabolite=df.final[which,])
  ggplot(trend.toplot,aes(x=time,y=metabolite,color=treatment))+
    geom_point()+
    stat_summary(fun.y = FUN,geom = 'line')+
    geom_smooth(method = 'lm', formula = y~poly(x,2),se = F,linetype = '3313')+
    theme(legend.position=c(0.9,0.1),axis.title = element_text(size=size),legend.title = element_text(size=size),legend.text = element_text(size=3),legend.key.size = unit(0.08,'cm'), plot.title = element_text(size=size),panel.background = element_blank())+
    labs(...)
}

plot.NullDistribution <- function(permu.data,model.data,colnames,cutoff,R,alpha,size){
  permu.data <- permu.data
  model.data <- model.data
  model.leverage <- data.frame(metabolites=factor(colnames[1:length(model.data)]),leverage=model.data)
  Nulldist <- data.frame(permu.data) %>% 
    gather(.,key = metabolites,value = leverage) %>% 
    mutate(.,metabolites=rep(colnames[1:nrow(permu.data)],each=ncol(permu.data))) %>% 
    ggplot(.,aes(metabolites,leverage,color=metabolites))+
    geom_violin(draw_quantiles = c(1-alpha))+
    geom_point(data = model.leverage,aes(x = metabolites,y = leverage))+
    geom_hline(yintercept = cutoff)+
    theme(axis.text.x.bottom = element_text(angle = 90,size = 7,hjust = 1,vjust = 0.5),panel.background = element_blank(),legend.position=c(0.5,0.8),axis.title = element_text(size=size),legend.title = element_text(size=size),legend.text = element_text(size=3),legend.key.size = unit(0.08,'cm'), plot.title = element_text(size=4))+
    labs(title = paste('Null distribution by',R,'rounds of permutation. alpha:',alpha))
  Nulldist
}
