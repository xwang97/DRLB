ASSO<-function(MAT,DIM=1000,Thres){
  
  if(min(ncol(MAT),nrow(MAT))<DIM){
    DIM<-min(ncol(MAT),nrow(MAT))-1
  }
  
  
  
  X<-MAT %*% t(MAT)
  VEC<-diag(X)
  X<-sweep(X,2,VEC,FUN = "/")
  X<-apply(X>Thres, 2, as.numeric)
  X[is.na(X)]<-0
  
  A<-MAT
  B<-NULL
  C<-NULL
  D<-A*0
  Thres_use<-0.05*sum(A)
  COL<-0
  NOW<-0
  for (i in 1:ncol(X)) {
    B_TEMP<-X[,i]
    C_TEMP<-rep(1,ncol(A))
    A_TEMP<-B_TEMP%o%C_TEMP
    TEMP<-colSums(A_TEMP*A*2-A_TEMP)
    C_TEMP<-as.numeric(TEMP>0)
    if(NOW<sum(TEMP[TEMP>0])){
      B<-B_TEMP
      C<-C_TEMP
      NOW<-sum(TEMP[TEMP>0])
      COL<-i
    }
  }
  X<-X[,-COL]
  A_stop<-sum(apply((A-B%o%C)>0,2,as.numeric))
  
  
  for (m in 2:DIM) {
    
    if(A_stop>Thres_use){
      B_use<-NULL
      C_use<-NULL
      COL<-0
      for (i in 1:ncol(X)) {
        B_TEMP<-cbind(B,X[,i])
        C_TEMP<-rbind(C,rep(1,ncol(A)))
        A_TEMP<-apply((B_TEMP%*%C_TEMP)>0,2,as.numeric)
        TEMP<-colSums(A_TEMP*A*2-A_TEMP)
        if(NOW<sum(TEMP[TEMP>0])){
          NOW<-sum(TEMP[TEMP>0])
          B_use<-X[,i]
          C_use<-as.numeric(TEMP>0)
          COL<-i
        }
      }
      
      
      if(COL==0){
        break
      }
      B<-cbind(B,B_use)
      C<-rbind(C,C_use)
      X<-X[,-COL]
      A_stop<-sum(apply((A-B%*%C)>0,2,as.numeric))
    }
  }
  
  result<-list()
  result[[1]]<-B
  result[[2]]<-C
  return(result)
}
