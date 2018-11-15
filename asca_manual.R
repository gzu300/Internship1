######
# 13/11/18
# permutate, calcualte_levergeSPE, and permutate_x_times modules created
###


permutate <- function(mx_raw){
  #rows are samples;columns are variables
  m <- nrow(mx_raw)
  places <- sample(1:m)
  m.output <- mx_raw[places,]
  return(m.output)
}

calculate_leverageSPE <- function(mx, pcs=2){
  pca.model <- prcomp(mx, center = T, scale. = F)
  loadings <- pca.model$rotation[,1:pcs]
  scores <- pca.model$x[,1:pcs]
  df <- as.data.frame(diag(loadings%*%t(loadings)))
  colnames(df) <- 'leverage'
  resi <- mx-(scores%*%t(loadings))
  SPE <- diag(t(resi)%*%resi)
  df$SPE <- SPE
  return(df)
}

permutate_x_times <- function(mx,R){
  for (r in 1:R){
    m.output <- permutate(mx)
    each_distri <- calculate_leverageSPE(m.output)
    if (r==1){
      sum_distri <- each_distri
      print(r)
    }
    else{
      sum_distri <- rbind(sum_distri, each_distri)
    }
  }
  return(sum_distri)
}

design.a <- function(){
  load('Example.RData')
  rm(data.example)
  design.a <- as.matrix(Designa[1:32,])
  return(design.a)
}

design.b <- function(){
  load('Example.RData')
  rm(data.example)
  design.b <- as.matrix(Designb[1:32,])
  return(design.b)
}

create_Xa <- function(){
  #14/11/18
  #rows: time factors. 4 levels. 2 treatments, 4 replicates
  #columns: Number of genes:100
  design.a <- design.a()
  mx <- matrix(c(rnorm(100, -3, 1), rnorm(100,-1,1),rnorm(100,1,1),rnorm(100, 3, 1)), nrow = 4, ncol = 100, byrow = T)
  #mean center the data
  mx <- as.matrix(scale(mx, scale = F))
  #broadcast to 32X100
  xa <- design.a%*%mx
  return(xa)
}

create_Xb <- function(){
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
  abn1 <- matrix(rnorm(25, 2, 1), nrow = 1, ncol = 25)#induction pattern 1:25
  abn2 <- matrix(rnorm(25,-1,1),nrow = 1,ncol = 25)#reduction pattern 26-50
  abn_extra1 <- 2#add extra noise
  mx[2,1:25] <- mx[2,1:25]+abn1
  mx[2,26:50] <- mx[2,26:50]+abn2
  mx[2,1] <- mx[2,1]+abn_extra1
  #mean center
  mx <- scale(mx,scale = F)
  #mx[,2] <- mx[,2]+abn_extra2
  load('Example.RData')
  rm(data.example)
  #broadcast to 32x100
  design.b <- design.b()
  xb <- design.b%*%mx
  return(xb)
  # abn2 <- matrix(c(rnorm(300, 5, 1), rnorm(200, -2, 1)), nrow = 50, ncol = 10, byrow = TRUE)
  # abn3 <- matrix(rnorm(1000, 5, 1), nrow = 50, ncol = 20)
  # 
  # group1 <- mx[1:50,1:20]
  # group2 <- mx[51:100, 21:30]
  # group3 <- mx[101:150, 31:50]
  # groups <- c(rep('g4',2),
  #             rep('g1', ncol(group1)-2),
  #             rep('g2', ncol(group2)),
  #             rep('g3', ncol(group3)),
  #             rep('ref', 40))
  # 
  # treatments <- c(rep('t1', nrow(group1)),
  #                 rep('t2', nrow(group2)),
  #                 rep('t3', nrow(group3)),
  #                 rep('ref', 50))
  # 
  # mx[1:50,1:20] <- group1+abn1
  # mx[51:100, 21:30] <- group2+abn2
  # mx[101:150, 31:50] <- group3+abn3
  # mx[,1] <- mx[,1]+abn_extra1
  # mx[,2] <- abn_extra2
}

create_Xab <- function(){
  #design matrix for interaction term
  #rows: interaction Xij. 8 rows
  #columns 100 columns
  design.ab <- matrix(0,nrow = 32, ncol = 8)
  load('Example.RData')
  rm(data.example)
  design.ab[1:16,1:4] <- design.ab[1:16,1:4]+Designa[1:16]
  design.ab[17:32,5:8] <- design.ab[17:32,5:8]+Designa[1:16]
  #8x100
  mx <- matrix(rnowm(800,0,1),nrow = 8,ncol = 100)
}

generate_dataset <- function(){
  a <- create_Xa()
  b <- create_Xb()
  mx <- a+b
  return(mx)
}