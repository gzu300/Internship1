RMSECV = function(X,LV,split){
  X = as.matrix(X)
  size = dim(X)
  I = size[1]
  J = size[2]
  NTest = ceiling(I/split);
  TESTSET = 1:(NTest*split)*NA
  TESTSET[1:I] = (1:I);
  TEST = matrix(TESTSET,nrow = split)
  PRESS = 1:10 * 0;
  RMSECV = PRESS
  Xhat = X*0
  XX = X;
  for(lv in 1:LV){
    for(s in 1:split){
      Train = 1:I;
      Test = TEST[s,];
      Test = Test[!is.na(Test)]
      Train = Train [! Train %in% Test]
      Xtrain = XX[Train,]
      Xtest = XX[Test,]
      
      # PCA on Train
      svd = svd(Xtrain);
      v = svd$v[,1]
      
      for (j in 1:J){
        # for each variable j; predict score of testset without using this variable
        # and then predict Xhat test set for this variable
        # This guarantees independence between Xhat and X.
        vj = v[j]
        v_j = v[-j]
        t_test_j = as.matrix(Xtest[,-j]) %*% as.vector(v_j) / sum(v_j^2)
        Xhat[Test,j] = Xhat[Test,j] + t_test_j * vj
      } # end J
      
    } # end split
    
    E = X - Xhat;
    PRESS[lv] = sum(as.vector(E^2))
    RMSECV[lv] = sqrt(PRESS[lv]/(I*J));
    
    # Deflate XX with 1 component
    svd_all = svd(XX);
    XX = XX - svd_all$d[1] * (svd_all$u[,1]%*%t(svd_all$v[,1]))
  } # end LV
 
  return(RMSECV)
} # end function
