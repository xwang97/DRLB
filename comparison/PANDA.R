GAMMA<-function(MAT_A,MAT_AC,Core){
  TEMP<-Core[[1]]%o%Core[[2]]
  TEMP<-TEMP+MAT_AC
  TEMP<-apply(TEMP>0, 2, as.numeric)
  TEMP<-abs(MAT_A-TEMP)
  return(sum(TEMP))
}

Find_Core<-function(MAT_A,MAT_AR,MAT_AC){
  Core<-list()
  E<-NULL
  S<-order(rowSums(MAT_AR),decreasing = T)
  Cl<-rep(0,nrow(MAT_AR))
  Cr<-rep(0,ncol(MAT_AR))
  Cl[S[1]]<-1
  Cr<-MAT_AR[S[1],]
  Core[[1]]<-Cl
  Core[[2]]<-Cr
  G<-GAMMA(MAT_A,MAT_AC,Core)+sum(Core[[1]])+sum(Core[[2]])
  for (i in 2:length(S)) {
    Core1<-list()
    Cl1<-Cl
    Cr1<-Cr
    Cl1[S[i]]<-1
    Cr1[which(MAT_AR[S[i],]==0)]<-0
    Core1[[1]]<-Cl1
    Core1[[2]]<-Cr1
    G1<-GAMMA(MAT_A,MAT_AC,Core1)+sum(Core1[[1]])+sum(Core1[[2]])
    if(G1<G){
      Core<-Core1
      G<-G1
    }else{
      E<-append(E,S[i])
    }
  }
  
  Core[[3]]<-E
  names(Core)<-c("left","right","Extension")
  return(Core)
}


Extend_Core<-function(Core,MAT_A,MAT_AC){
  
  Cl<-Core[[1]]
  Cr<-Core[[2]]
  
  G<-GAMMA(MAT_A,MAT_AC,Core)+sum(Core[[1]])+sum(Core[[2]])
  for (i in Core[[3]]) {
    Core1<-Core[1:2]
    Core1[[1]]<-Cl
    Core1[[1]][i]<-1
    G1<-GAMMA(MAT_A,MAT_AC,Core1)+sum(Core1[[1]])+sum(Core1[[2]])
    if(G1<G){
      Core[1:2]<-Core1
      G<-G1
    }
  }
  for (i in which(Core[[2]]==0)) {
    Core1<-Core[1:2]
    Core1[[2]]<-Cr
    Core1[[2]][i]<-1
    G1<-GAMMA(MAT_A,MAT_AC,Core1)+sum(Core1[[1]])+sum(Core1[[2]])
    if(G1<G){
      Core[1:2]<-Core1
      G<-G1
    }
  }
  return(Core)
} 


###PANDA
PANDA<-function(MAT,DIM=1000){
  
  if(min(ncol(MAT),nrow(MAT))<DIM){
    DIM<-min(ncol(MAT),nrow(MAT))-1
  }
  Thres_use<-0.05*sum(MAT)
  MAT_A<-MAT
  MAT_B<-NULL
  MAT_C<-NULL
  MAT_AR<-MAT_A
  MAT_AC<-0*MAT
  G<-sum(MAT_A)
  G1<-0
  for (i in 1:DIM) {
    if(sum(MAT_AR)<=Thres_use){
      result<-list()
      result[[1]]<-MAT_B
      result[[2]]<-MAT_C
      return(result)
      break()
    }else{
      Core<-Find_Core(MAT_A,MAT_AR,MAT_AC)
      Core<-Extend_Core(Core,MAT_A,MAT_AC)
      MAT_B<-cbind(MAT_B,Core[[1]])
      MAT_C<-rbind(MAT_C,Core[[2]])
      
      TEMP<-Core[[1]] %o% Core[[2]]
      MAT_AR<-MAT_AR-TEMP
      MAT_AR<-apply(MAT_AR>0,2,as.numeric)
      MAT_AC<-MAT_AC+TEMP
      MAT_AC<-apply(MAT_AC>0,2,as.numeric)
      
      G1<-sum(abs(MAT_A-MAT_AC))+sum(MAT_B)+sum(MAT_C)
      if(G1<G){
        G<-G1
      }else{
        result<-list()
        result[[1]]<-MAT_B
        result[[2]]<-MAT_C
        return(result)
        break()
      }
    }
  }
  result<-list()
  result[[1]]<-MAT_B
  result[[2]]<-MAT_C
  return(result)
}
