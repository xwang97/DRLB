
#BIND functions
#############
Block_filter<-function(MAT,Thres_l=0.01,Thres_h=0.8)
{
  ROW<-apply(MAT,1,mean)
  COL<-apply(MAT,2,mean)
  tg_r<-which((ROW>Thres_l)&(ROW<Thres_h))
  tg_c<-which((COL>Thres_l)&(COL<Thres_h))
  MAT<-MAT[tg_r,tg_c]
  return(MAT)
}

########
filter_size<-function(MAT,rr,cr)
{
	tr<-min(floor(nrow(MAT)*rr)+1,nrow(MAT))
	tc<-min(floor(ncol(MAT)*rr)+1,ncol(MAT))
	for(i in 1:nrow(MAT))
	{
		MAT[i,order(MAT[i,])[1:(nrow(MAT)-tr)]]<-0
	}
	for(i in 1:ncol(MAT))
	{
		MAT[order(MAT[,i])[1:(ncol(MAT)-tr)],i]<-0
	}
	return(MAT)
}


Block_filter2<-function(MAT,Thres_lc=0.01,Thres_hc=0.8,Thres_lr=0.01,Thres_hr=0.8)
{
  ROW<-apply(MAT,1,mean)
  COL<-apply(MAT,2,mean)
  tg_r<-which((ROW>Thres_lr)&(ROW<Thres_hr))
  tg_c<-which((COL>Thres_lc)&(COL<Thres_hc))
  MAT<-MAT[tg_r,tg_c]
  return(MAT)
}

BinBlock<-function(MAT){
  
  #MAT_orig<-MAT
   
  ppc<-apply(MAT,1,mean)
  ppr<-apply(MAT,2,mean)
  
  xc0<-sort(ppc)
  xr0<-sort(ppr)
  
  EM_xc0<-Generate_Edistr_base(xc0)
  EM_xr0<-Generate_Edistr_base(xr0)

  KS_c<-c()
  for(i in 1:ncol(MAT))
  {
    xc1<-sort(ppc[which(MAT[,i]==1)])
    KS_c<-c(KS_c,KS_dist_sum4(xc1,EM_xc0))
	if(i%%100==1)
	{
	print(i)
	}
  }
  KS_r<-c()
  for(i in 1:nrow(MAT))
  {
    xr1<-sort(ppr[which(MAT[i,]==1)])
    KS_r<-c(KS_r,KS_dist_sum4(xr1,EM_xr0))
	if(i%%100==1)
	{
	print(i)
	}
  }
  names(KS_r)<-rownames(MAT)
  names(KS_c)<-colnames(MAT)
  SS<-list(KS_r,KS_c)
  return(SS)
}

Generate_Edistr_base<-function(distr_base,kk=100)
{
xc000<-c()
for(i in 1:length(distr_base))
{
xc000<-c(xc000,rep(distr_base[i],floor(distr_base[i]*kk)+1))
}
rr0<-sort(xc000)
return(rr0)
}

######

BinBlock2_BS<-function(MAT)
{
  #MAT_orig<-MAT
   
  ppc<-apply(MAT,1,mean)
  ppr<-apply(MAT,2,mean)
  
  xc0<-sort(ppc)
  xr0<-sort(ppr)
  
  EM_xc0<-Generate_Edistr_base(xc0)
  EM_xr0<-Generate_Edistr_base(xr0)
  
  EM_xc0<-sort(EM_xc0)	
  KS_c<-c()
  BS_col<-c()
  temp_c<-ppc*0
  for(i in 1:ncol(MAT))
  {
      xc1<-sort(ppc[which(MAT[,i]==1)])
	ww<-0
	if(length(xc1)>1)
	{
		diff_c<-KS_dist_sum5(xc1,EM_xc0)
		ww<--sum(diff_c[which(diff_c<0)])
 		temp_c<-ppc*0
		temp_c[names(diff_c)]<-diff_c
		ccc<-temp_c[order(ppc)]
		ccc[length(ccc)]<-1e-10
		tgs<-which(ccc!=0)
		tgs0<-c(0,tgs[-length(tgs)])
		index<-tgs-tgs0
		ids<-1:length(index)
		new_ids<-unlist(mapply(rep,ids,index))
		diff_c1<-c(diff_c,1e-10)
		BS_c<-diff_c1[new_ids]
		names(BS_c)<-names(temp_c[order(ppc)])
		BS_col<-cbind(BS_col,BS_c[rownames(MAT)])
	}
      KS_c<-c(KS_c,ww)
	if(i%%500==1)
	{
		print(i)
	}
  }
  colnames(BS_col)<-colnames(MAT)

  EM_xr0<-sort(EM_xr0)
  KS_r<-c()
  BS_row<-c()
  temp_c<-ppr*0
  for(i in 1:nrow(MAT))
  {
      xr1<-sort(ppr[which(MAT[i,]==1)])
	ww<-0
	if(length(xr1)>1)
	{
		diff_c<-KS_dist_sum5(xr1,EM_xr0)
		ww<--sum(diff_c[which(diff_c<0)])
		temp_c<-ppr*0
		temp_c[names(diff_c)]<-diff_c
		ccc<-temp_c[order(ppr)]
		ccc[length(ccc)]<-1e-10
		tgs<-which(ccc!=0)
		tgs0<-c(0,tgs[-length(tgs)])
		index<-tgs-tgs0
		ids<-1:length(index)
		new_ids<-unlist(mapply(rep,ids,index))
		diff_c1<-c(diff_c,1e-10)
		BS_c<-diff_c1[new_ids]
		names(BS_c)<-names(temp_c[order(ppr)])
		BS_row<-rbind(BS_row,BS_c[colnames(MAT)])
	}
      KS_r<-c(KS_r,ww)
	if(i%%1000==1)
	{
		print(i)
	}
  }
  rownames(BS_row)<-rownames(MAT)
  names(KS_r)<-rownames(MAT)
  names(KS_c)<-colnames(MAT)
  SS<-list(KS_r,KS_c, BS_row, BS_col)
  return(SS)
}

BinBlock2_PLOT<-function(MAT)
{
  #MAT_orig<-MAT
   
  ppc<-apply(MAT,1,mean)
  ppr<-apply(MAT,2,mean)
  
  xc0<-sort(ppc)
  xr0<-sort(ppr)
  
  EM_xc0<-Generate_Edistr_base(xc0)
  EM_xr0<-Generate_Edistr_base(xr0)
  
  EM_xc0<-sort(EM_xc0)
  KS_c<-c()
  list_cc<-list()
  for(i in 1:ncol(MAT))
  {
      xc1<-sort(ppc[which(MAT[,i]==1)])
	ww<-0
	list_cc[[i]]<-0
	if(length(xc1)>1)
	{
		diff_c<-KS_dist_sum5(xc1,EM_xc0)
		list_cc[[i]]<-diff_c
		ww<--sum(diff_c[which(diff_c<0)])
	}
      KS_c<-c(KS_c,ww)
	if(i%%1000==1)
	{
	print(i)
	}
  }

  EM_xr0<-sort(EM_xr0)
  list_rr<-list()
  KS_r<-c()
  for(i in 1:nrow(MAT))
  {
      xr1<-sort(ppr[which(MAT[i,]==1)])
	ww<-0
	list_rr[[i]]<-0
	if(length(xr1)>1)
	{
		diff_c<-KS_dist_sum5(xr1,EM_xr0)
		list_rr[[i]]<-diff_c
		ww<--sum(diff_c[which(diff_c<0)])
	}
      KS_r<-c(KS_r,ww)
	if(i%%1000==1)
	{
		print(i)
	}
  }
  names(KS_r)<-rownames(MAT)
  names(KS_c)<-colnames(MAT)
  names(list_rr)<-rownames(MAT)
  names(list_cc)<-colnames(MAT)
  SS<-list(KS_r,KS_c, list_rr, list_cc)
  return(SS)
}


BinBlock2<-function(MAT){
  
  #MAT_orig<-MAT
   
  ppc<-apply(MAT,1,mean)
  ppr<-apply(MAT,2,mean)
  
  xc0<-sort(ppc)
  xr0<-sort(ppr)
  
  EM_xc0<-Generate_Edistr_base(xc0)
  EM_xr0<-Generate_Edistr_base(xr0)
  
  EM_xc0<-sort(EM_xc0)
  KS_c<-c()
  for(i in 1:ncol(MAT))
  {
      xc1<-sort(ppc[which(MAT[,i]==1)])
	ww<-0
	if(length(xc1)>1)
	{
		diff_c<-KS_dist_sum5(xc1,EM_xc0)
		ww<--sum(diff_c[which(diff_c<0)])
	}
      KS_c<-c(KS_c,ww)
	if(i%%1000==1)
	{
	print(i)
	}
  }

  EM_xr0<-sort(EM_xr0)
  #pp_r<-1:length(EM_xr0)/length(EM_xr0)
  KS_r<-c()
  for(i in 1:nrow(MAT))
  {
      xr1<-sort(ppr[which(MAT[i,]==1)])
	ww<-0
	if(length(xr1)>1)
	{
		diff_c<-KS_dist_sum5(xr1,EM_xr0)
		ww<--sum(diff_c[which(diff_c<0)])
	}
      KS_r<-c(KS_r,ww)
	if(i%%1000==1)
	{
		print(i)
	}
  }
  names(KS_r)<-rownames(MAT)
  names(KS_c)<-colnames(MAT)
  SS<-list(KS_r,KS_c)
  return(SS)
}

KS_dist_sum5_past<-function(tg_distr,rr_bg)
{
	rr0<-rr_bg
	rr<-tg_distr
	tg_aa<-floor((1:(length(rr)-1)/length(rr))*length(rr0))
	tg_aa1<-tg_aa+1
	tt1<-(rr0[tg_aa]+rr0[tg_aa1])/2-rr[1:(length(rr)-1)]
	tt2<-(rr0[tg_aa]+rr0[tg_aa1])/2
	tt<-tt1/(1-tt2)
	names(tt)<-names(tg_distr)
	return(tt)
}


KS_dist_sum5<-function(tg_distr,rr_bg)
{
	  rr0<-rr_bg
        rr<-tg_distr
        tg_aa<-floor((1:(length(rr))/(length(rr)+1))*length(rr0))
        tg_aa1<-tg_aa
        tt1<-(rr0[tg_aa]+rr0[tg_aa1])/2-rr[1:(length(rr))]
        tt2<-(rr0[tg_aa]+rr0[tg_aa1])/2
        tt<-tt1/(1-tt2)
        names(tt)<-names(tg_distr)
        return(tt)
}



KS_dist_sum4<-function(tg_distr,rr0)
{
  rr0<-sort(rr0)
  pp<-1:length(rr0)/length(rr0)
  rr<-sort(tg_distr)
  tt_dist<-0
  if(length(rr)>1){
    for(i in 1:(length(rr)-1))
    {
      rs<-i/length(rr)
      ii<-sum(rs>=pp)
      tt<-mean(rr0[ii],rr0[ii+1])-rr[i]
      tt<-tt/(1-mean(rr0[ii],rr0[ii+1]))
      if(tt<0)
      {
        tt_dist<-tt_dist-tt
      }
    }
  }
  return(tt_dist)
}


denoise_with_bind <- function(data) {
  xx <- data
  rownames(xx) <- paste("R", 1:nrow(xx), sep = "")
  colnames(xx) <- paste("C", 1:nrow(xx), sep = "")
  BSxx <- BinBlock2_BS(xx)
  BS_BIND <- t(t(BSxx[[1]])) %*% BSxx[[2]]
  xx2 <- xx * (BS_BIND > quantile(BS_BIND, 0.8))
  return(xx2)
}

denoise_with_new <- function(data) {
  xx <- data
  rownames(xx) <- paste("R", 1:nrow(xx), sep = "")
  colnames(xx) <- paste("C", 1:nrow(xx), sep = "")
  BSxx <- BinBlock2_BS(xx)
  
  aaa1 <- BSxx[[3]]
  aaa2 <- BSxx[[4]]
  
  aaa1 <- -aaa1
  aaa2 <- -aaa2
  
  ccc <- aaa1 + aaa2
  ccc <- abs(aaa1 + aaa2) * (sign(aaa1) + sign(aaa2))
  # boxplot(ccc[which(xx0 == 1)], ccc[which(xx0 == 0)])
  
  ccc_new <- ccc
  cut_q <- 0.8
  ccc_new[which(ccc_new < quantile(ccc_new, cut_q))] <- 0
  # boxplot(ccc_new[which(xx0 == 1)], ccc_new[which(xx0 == 0)])
  
  BS_new <- ccc_new
  xx2 <- xx * (BS_new > quantile(BS_new, 0.8))
  
  ccc_new1 <- filter_size(ccc_new, 0.25, 0.25)
  BS_new1 <- ccc_new1
  xx3 <- xx * (BS_new1 > quantile(BS_new1, 0.2))
  
  ############
  aaa1 <- BSxx[[3]]
  aaa2 <- BSxx[[4]]
  
  aaa1 <- -aaa1
  aaa2 <- -aaa2
  
  ccc <- aaa1 + aaa2
  ppbg <- t(t(apply(xx, 1, mean))) %*% apply(xx, 2, mean)
  ccc <- abs(aaa1 + aaa2) * (sign(aaa1) + sign(aaa2)) * ppbg
  
  ccc_new <- ccc
  ccc_new[which(ccc_new < 0)] <- 0
  ccc_new2 <- filter_size(ccc_new, 0.25, 0.25)
  ccc_new2[which(ccc_new < quantile(ccc_new, 0.8))] <- 0
  
  BS_new2 <- ccc_new2
  xx4 <- xx * (BS_new2 > quantile(BS_new2, 0.2))
  return(list(xx2, xx3, xx4))
}