######
# 13/11/18
# permutate, calcualte_levergeSPE, and permutate_x_times modules created
###

library(tidyverse)

permutate <- function(mx_raw){
  #rows are samples;columns are variables
  m <- nrow(mx_raw)
  places <- sample(1:m)
  m.output <- mx_raw[places,]
  return(m.output)
}

permutate_x_times <- function(mx,R,matrix.a,matrix.b){
  #16/11 change from rbind to cbind to make quantile calculation easier. this is bullshit. add tags for each
  #gene is easier
  #2007 anova-gene calculates overall 99%quantile of a gene's 99% quantile leverage
  #permutate->split(xa,xb,xab)->pca->leverage_spe
  #mx:NxM. M sample, N genes
  for (r in 1:R){
    data <- permutate(mx)#NxM
    X.. <- mean(data)
    Xa <- ((data%*%matrix.a)/(nrow(matrix.a)/ncol(matrix.a)))%*%t(matrix.a)-X..#NxM
    Xb <- ((data%*%matrix.b)/(nrow(matrix.b)/ncol(matrix.b)))%*%t(matrix.b)-X..
    Xab <- data-Xa-Xb+X..
    each_distri <- calculate_leverageSPE(t(Xb+Xab))#MxN
    if (r==1){
      sum_distri <- each_distri
      #print(r)
    }
    else{
      sum_distri <- rbind(sum_distri, each_distri)
    }
  }
  return(sum_distri)
}

calculate_leverageSPE <- function(mx, pcs=1){
  #16/11/18 naming of the columns are not important anymore due to rbind to cbind
  #mx is MxN. M sample. N genes
  pca.model <- prcomp(mx, center = T, scale. = F)
  loadings <- pca.model$rotation[,1:pcs]
  scores <- pca.model$x[,1:pcs]
  df <- as.data.frame(diag(loadings%*%t(loadings)))
  colnames(df) <- 'leverage'
  resi <- mx-(scores%*%t(loadings))
  SPE <- diag(t(resi)%*%resi)
  df$SPE <- SPE
  df$gene <- c(1:ncol(mx))
  return(df)
}

leverage_lim <- function(permutated.data){
  #16/11/18 rows are different genes. columns are different permutations
  quant_per_gene <- permutated.data %>% 
    group_by(gene) %>% 
    summarise(cutoff.per_gene = quantile(leverage, probs = 0.95))
  limit <- quantile(quant_per_gene$cutoff.per_gene,probs = 0.95)
  return(limit)
}

spe_lim <- function(permutated.data){
  m <- mean(permutated.data$SPE)
  v <- var(permutated.data$SPE)
  h <- 2*m*m/v
  g <- v/(2*m)
  limit <- g*qchisq(p = 0.95,df = h)
  return(limit)
}

design.a <- function(Designa){
  design.a <- as.matrix(Designa[1:32,])
  return(design.a)
}

design.b <- function(Designb){
  design.b <- as.matrix(Designb[1:32,])
  return(design.b)
}

design.ab <- function(design.a){
  chunk <- design.a[1:16,]
  blank <- matrix(0,32,8)
  blank[1:16,1:4] <- blank[1:16,1:4]+chunk
  blank[17:32,5:8] <- blank[17:32,5:8]+chunk
  return(blank)
}

create_Xa <- function(Designa){
  #14/11/18
  #rows: time factors. 4 levels. 2 treatments, 4 replicates
  #columns: Number of genes:100
  design.a <- design.a(Designa)
  mx <- matrix(c(rnorm(100, -3, 0), rnorm(100,-1,0),rnorm(100,1,0),rnorm(100, 3, 0)), nrow = 4, ncol = 100, byrow = T)
  #mean center the data
  mx <- as.matrix(scale(mx, scale = F))
  #broadcast to 32X100
  xa <- design.a%*%mx
  return(xa)
}

create_Xb <- function(Designb){
  ##14/11/18
  #rows:treatments
  #columns: genes
  #3 different trends of genes.
  #gene 1-25 increase
  #gene 26-50 decrease
  #gene 1 extra increase
  #base noise
  mx <- matrix(rnorm(200, 0, 0),nrow = 2, ncol = 100)
  #patterns
  abn1 <- matrix(rnorm(25, 2, 0), nrow = 1, ncol = 25)#induction pattern 1:25
  abn2 <- matrix(rnorm(25,-1,0),nrow = 1,ncol = 25)#reduction pattern 26-50
  abn_extra1 <- 2#add extra noise
  mx[2,1:25] <- mx[2,1:25]+abn1
  mx[2,26:50] <- mx[2,26:50]+abn2
  mx[2,1] <- mx[2,1]+abn_extra1
  #mean center
  mx <- scale(mx,scale = F)
  #mx[,2] <- mx[,2]+abn_extra2
  #broadcast to 32x100
  design.b <- design.b(Designb)
  xb <- design.b%*%mx
  return(xb)
}

create_Xab <- function(Designa){
  #design matrix for interaction term
  #rows: interaction Xij. 8 rows
  #columns 100 columns
  design.ab <- matrix(0,nrow = 32, ncol = 8)
  design.ab[1:16,1:4] <- design.ab[1:16,1:4]+Designa[1:16]
  design.ab[17:32,5:8] <- design.ab[17:32,5:8]+Designa[1:16]
  #8x100
  mx <- matrix(rnowm(800,0,1),nrow = 8,ncol = 100)
}

generate_dataset <- function(Designa, Designb){
  a <- create_Xa(Designa)
  b <- create_Xb(Designb)
  mx <- a+b
  return(mx)
}