######
# 13/11/18
# permutate, calcualte_levergeSPE, and permutate_x_times modules created
###

library(tidyverse)

permutate <- function(mx_raw){
  #perterb columns
  # 1,2,3     1,3,2
  # 1,2,3  -> 1,3,2 
  m <- ncol(mx_raw)
  places <- sample(1:m)
  m.output <- mx_raw[,places]
  return(m.output)
}

permutate_x_times <- function(mx,R,matrix.a,matrix.b,matrix.ab){
  #16/11 change from rbind to cbind to make quantile calculation easier. this is bullshit. add tags for each
  #gene is easier
  #2007 anova-gene calculates overall 99%quantile of a gene's 99% quantile leverage
  #permutate->split(xa,xb,xab)->pca->leverage_spe
  #mx:NxM. M sample, N genes
  for (r in 1:R){
    data <- permutate(mx)#NxM
    X.. <- apply(data,1,mean)
    Xa <- ((data%*%matrix.a)/(nrow(matrix.a)/ncol(matrix.a)))%*%t(matrix.a)-X..#NxM
    Xb <- ((data%*%matrix.b)/(nrow(matrix.b)/ncol(matrix.b)))%*%t(matrix.b)-X..
    Xab <- ((data%*%matrix.ab)/(nrow(matrix.ab)/ncol(matrix.ab)))%*%t(matrix.ab)-Xa-Xb+X..
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

calculate_leverageSPE <- function(mx, pcs=2){
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
######################################################
design.a <- function(Designa){
  design.a <- as.matrix(Designa[1:32,])
  return(design.a)
}

design.b <- function(Designb){
  design.b <- as.matrix(Designb[1:32,])
  return(design.b)
}

design.ab <- function(Designa){
  design.a <- as.matrix(Designa[1:32,])
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
  a <- design.a%*%mx
  return(a)
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
  mx <- matrix(rnorm(200, 0, 1),nrow = 2, ncol = 100)
  #patterns
  abn1 <- matrix(rnorm(25, 2, 0), nrow = 1, ncol = 25)#induction pattern 1:25
  abn2 <- matrix(rnorm(25,-2,0),nrow = 1,ncol = 25)#reduction pattern 26-50
  abn_extra1 <- 2#add extra noise
  mx[2,1:25] <- mx[2,1:25]+abn1
  mx[2,26:50] <- mx[2,26:50]+abn2
  #mx[2,1] <- mx[2,1]+abn_extra1
  #mean center
  mx <- scale(mx,scale = F)
  #mx[,2] <- mx[,2]+abn_extra2
  #broadcast to 32x100
  design.b <- design.b(Designb)
  b <- design.b%*%mx
  return(b)
}

create_Xab <- function(Designa){
  #design matrix for interaction term
  #rows: interaction Xij. 8 rows
  #columns 100
  #8x100
  mx <- matrix(rnorm(8*77,0,0),nrow = 8,ncol = 77)
  abn1 <- matrix(c(-2,-1,1,2,2,1,-1,-2),8,10)
  abn2 <- matrix(c(0,0,0,0,0,0,0,0),8,10)
  abn3 <- matrix(c(-1,3,-1,-1,1,-3,1,1),8,3)
  mx <- cbind(abn1,abn2,abn3,mx)
  mx <- scale(mx,scale = F)
  design.ab <- design.ab(Designa)
  ab <- design.ab%*%mx
  return(ab)
}

generate_dataset <- function(Designa, Designb, X..=0){
  a <- create_Xa(Designa)
  b <- create_Xb(Designb)
  ab <- create_Xab(Designa)
  mx <- a+b+ab
  #mx <- scale(mx,scale = F)
  return(mx)
}