VTM <-function(vc, dm) {
  matrix(vc, ncol=length(vc), nrow=dm, byrow=T)
}

##################### Splus function : cumsum2 #####################
## cumsum2 gets cumsum of each column, that is, apply(X, 2, cumsum)

cumsum2 <- function(mydat)     #cumsum by row, col remains the same
{
  if(is.null(dim(mydat))) return(cumsum(mydat))
  else{
    out <- matrix(cumsum(mydat), nrow=nrow(mydat))
    out <- out - VTM(c(0, out[nrow(mydat), -ncol(mydat)]), nrow(mydat))
    return(out)
  }
}

#################### End of function : cumsum2 #####################
sum.I <- function(yy,FUN,Yi,Vi=NULL)
{
  if (FUN=="<"|FUN==">=") { yy <- -yy; Yi <- -Yi}
  pos <- rank(c(yy,Yi),ties.method='f')[1:length(yy)]-rank(yy,ties.method='f')
  if (substring(FUN,2,2)=="=") pos <- length(Yi)-pos
  if (!is.null(Vi)) {
    if(substring(FUN,2,2)=="=") tmpind <- order(-Yi) else  tmpind <- order(Yi)
    ##Vi <- cumsum2(as.matrix(Vi)[tmpind,])
    Vi <- apply(as.matrix(Vi)[tmpind,,drop=F],2,cumsum)
    return(rbind(0,Vi)[pos+1,])
  } else return(pos)
}

## ======================================================= ##
## Calculating the Sampling Weights with Possible Matching ##
## ======================================================= ##  
FNCC.WGT.FUN <- function(data,V.IND, Iik0)
{
  xi = data[,1];  di = data[,2]; vi = data[,3]; nv = sum(vi); n1 = sum(di)
  xk = xi[vi==1]; dk = di[vi==1]; phatk = dk; ind0 = dk == 0; nv0 = sum(ind0)
  ## tl: event time for the cases; nl: (matched) risk set size for tl
  tl = V.IND[,1]; n.Rl = V.IND[,2]; 
  #if(is.null(Zmatch.ind))
  #  {
  #    pl = c(1,cumprod(1-m.match/(n.Rl-1))); phatk[ind0]=1-pl[sum.I(xk[ind0],">",tl)+1]
  #  }else{
  junk = Iik0*log(1-m.match/pmax(n.Rl-1,m.match)); junk[is.na(junk)] = 0
  phatk[ind0]=1-exp(apply(junk,2,sum))
  #  }
  phatk
}

## ====================================================== ##
## Perturbing the Sampling Weights with Possible Matching ##
## ====================================================== ##  
PtbNCC.WGT.FUN <- function(data, V.IND, Iik0, Zmatch=NULL,B0=500,wgtk = NULL,B.dist=myB.dist)
{ 
  ## V.IND: 1st col = ti; 2nd col = n(ti);  3:(nv0+2): indicate xk is the control of ti ##
  ## Iik0: indicates the matched risk set of ti 
  ti = V.IND[,1]; n.Ri = V.IND[,2]; id.ii = V.IND[,3]; V.IND = V.IND[,-(1:3)]; n1 = length(ti)
  
  dj.i = xj.i = matrix(NA,ncol=m.match,nrow=n1)
  V.IND.vec = na.omit(c(V.IND))
  xj.i[!is.na(c(V.IND))] = data[V.IND.vec,1] 
  dj.i[!is.na(c(V.IND))] = data[V.IND.vec,2] 
  ind.control = c(dj.i!=1) ## try to identify those dj.i = 1 and sampled as control 
  id.0j.all = sort(unique(na.omit(V.IND[ind.control]))); id.all = sort(unique(c(id.ii,id.0j.all)))
  x0j = c(xj.i)[match(id.0j.all,c(V.IND))]; nsub = length(id.all)
  
  wptb.mat = matrix(0,ncol=B0,nrow=nsub); 
  for(bb in 1:B0)
  {
    Vjptb = rep(0,nsub); Bii = rexp(n1); 
    if(B.dist=="exp"){B0mat = matrix(rexp(n1*m.match),ncol=m.match)}else{B0mat = matrix(rchisq(n1*m.match,2)*0.5,ncol=m.match)} 
    V0jptb = 1 - tapply(c(1-B0mat[ind.control]), c(V.IND[ind.control]), prod,na.rm=T)
    Vjptb[match(id.0j.all,id.all)] = V0jptb	
    Vjptb[match(id.ii,    id.all)] = Bii	
    dLam.ti = apply(B0mat*!is.na(V.IND), 1, sum)/(n.Ri); Lam0j.ptb = apply(dLam.ti*Iik0,2,sum)
    pjptb = rep(1,nsub); pjptb[match(id.0j.all,id.all)] = 1-exp(-Lam0j.ptb)
    wptb.mat[,bb] = Vjptb/pjptb;  
  }
  if(is.null(wgtk))
  {
    return(cbind(id.all,wptb.mat)[order(id.all),])
  }else{
    return(list("wptb"=cbind(id.all,wptb.mat)[order(id.all),],"wptb.naive"=wgtk*matrix(rexp(nsub*B0),nrow=nsub)))
  }
}

coxph.ALASSO = function(data,Gi,regularize,const)
{
  data = as.matrix(data); nn = nrow(data); y = data[,1]; d = data[,2]; x = data[,-(1:2)]; pp = ncol(x)
  bini = coxph(Surv(y,d)~x,weight=Gi,robust=T); I.info = solve(bini$naive.var); 
  if(regularize)
  {
    Db = diag(abs(bini$coef)); I.half = svd(I.info); I.half = I.half$u%*%diag(sqrt(I.half$d))%*%I.half$v
    pseudoY = c(I.half%*%bini$coef); pseudoX = I.half%*%Db 
    tmpfit = lars(pseudoX, pseudoY,type="lasso",normalize=F,intercept=F)
    lam.all = c(seq(0,max(tmpfit$lambda),length=200)); m0 = length(lam.all);
    b.all = predict(tmpfit,s=lam.all,type="coefficients",mode="lambda")$coef; 
    df.all = apply(b.all!=0,1,sum); 
    BIC.lam = apply((pseudoY - pseudoX%*%t(b.all))^2,2,sum) + log(sum(d))*df.all
    m.opt = (1:m0)[BIC.lam==min(BIC.lam)]; bhat = b.all[m.opt,]*abs(bini$coef); lamhat = lam.all[m.opt]; ##print(lamhat)
  }else{
    bhat = bini$coef	  	
  }
  bhat
}

PTB.Gamma.FUN <- function(wgt.ptb,data,gammahat,regularize,newy,betaonly=F,const)
{
  n.t0 = length(t0); xk = data[,1]; dk = data[,2]; wgtk = data[,3]; yk = data[,-(1:3),drop=F]; pp = ncol(yk)
  if(is.null(newy)){newy=rep(0,pp)}; newy = as.matrix(newy); n.y0 = ncol(newy) ## pp x n.y0
  cox.ptb = function(wgtp,xk,dk,yk,wgtk,regularize,newy,betaonly,const)
  {
    ##tmpind = 1:length(xk) 
    tmpind = wgtp>0; datanew = cbind(xk,dk,yk)[tmpind,]; 
    out=betap=coxph.ALASSO(datanew,Gi=wgtp[tmpind],regularize=regularize,const=const);##print(length(betap)); 
    if(!betaonly){
      ebyk = c(exp(yk%*%betap)); tj = xk[dk==1]; pi.tj = PI.k.FUN(tj,ebyk,xk,yk,wgtp)/sum(wgtp)
      Lamt0p = sum.I(t0, ">=", tj, wgtp[dk==1]/pi.tj)/sum(wgtp); 
      out=c(log(Lamt0p)+rep(c(t(newy)%*%betap),rep(n.t0,n.y0)),out) ## n.t0 x n.y0
    }
    out
  }
  apply(wgt.ptb,2,cox.ptb,xk=xk,dk=dk,yk=yk,wgtk=wgtk,regularize=regularize,newy=newy,betaonly=betaonly,const=const)
}

NCC.Cox.FUN.Regularize <- function(data,V.IND,Iik0,wgtk.ptb=NULL,B0=500,y0.cut=NULL,Zmatch.ind=NULL,rtn="ALL",regularize=F,betaonly=F,pos.only=F)
{
  xi = data[,1]; di = data[,2]; vi = data[,3]; yi = data[,-c(1:3,Zmatch.ind),drop=F]; n.t0 = length(t0); py = ncol(yi)
  NN = length(xi); nv = sum(vi); tj = V.IND[,1]; pi.tj = V.IND[,2]/NN
  wgtk = 1/FNCC.WGT.FUN(data,V.IND,Iik0); ##print(sum(wgtk))
  xk = xi[vi==1]; dk = di[vi==1]; yk = yi[vi==1,,drop=F]; idk = (1:NN)[vi==1];   
  gammahat=betahat=coxph.ALASSO(data=cbind(xk,dk,yk),Gi=wgtk,regularize=regularize,const=0); 
  if(!betaonly){
    ebyk = c(exp(yk%*%betahat)); pi.tj = PI.k.FUN(tj,ebyk,xk,yk,wgtk)/sum(wgtk)
    LamCox.t0 = sum.I(t0, ">=", tj, 1/pi.tj)/sum(wgtk); 
    if(is.null(y0.cut)){y0.cut=rep(0,py)}; y0.cut = as.matrix(y0.cut); n.y0 = ncol(y0.cut) ## py x n.y0
    out = gammahat = c(log(LamCox.t0)+rep(t(y0.cut)%*%betahat,rep(n.t0,n.y0)),gammahat)
  }
  if(rtn=="EST"){return(gammahat)}else{
    wgtk.ptb = PtbNCC.WGT.FUN(data,V.IND,Iik0,B0=B0); wgtk.ptb = wgtk.ptb[match(idk,wgtk.ptb[,1]),-1]
    ## ==================================================================================================== ##
    ## to handle negative weights, perturb w/ wgtk.ptb + const*wgtk --> sd(gammahat)/(1+tau), hence rescale ##
    ## ==================================================================================================== ##
    ind.pos = apply(wgtk.ptb,2,min)>0; wgtk.ptb.neg = wgtk.ptb[,!ind.pos]; wgtk.ptb.pos = wgtk.ptb[,ind.pos];
    ##scale0 = mean((apply(wgtk.ptb,1,sd)/sqrt(apply(wgtk.ptb^2*(wgtk.ptb>0),1,mean)-apply(wgtk.ptb*(wgtk.ptb>0),1,mean)^2))) 
    scale0 = mean((apply(abs(wgtk.ptb-wgtk),1,mean)/apply(abs(wgtk.ptb-wgtk)*(wgtk.ptb>0),1,mean)))
    if(pos.only & mean(ind.pos)>0.5){
      gammaptb = PTB.Gamma.FUN(wgtk.ptb.pos,cbind(xk,dk,wgtk,yk),
                               gammahat=gammahat,regularize=regularize,newy=y0.cut,betaonly=betaonly,const=0) 
    }else{
      if(sum(!ind.pos)>0){const = max((1-wgtk.ptb.neg)/wgtk)}else{const=0}; 
      gammaptb = PTB.Gamma.FUN(wgtk.ptb,cbind(xk,dk,wgtk,yk),gammahat=gammahat,
                               regularize=regularize,newy=y0.cut,betaonly=betaonly,const=const) 
      tmpind = -(1:(n.y0*n.t0))						 
      gammaptb[tmpind,] = (gammaptb[tmpind,] - gammahat[tmpind])*scale0+gammahat[tmpind]
    }
    gamma.sd = sqrt(apply( (gammaptb-gammahat)^2*(apply(gammaptb+gammahat==0,1,mean)<0.5),1,mean, na.rm=T)); 
    if(!betaonly){
      cutoff.sup = NULL
      for(kk in 1:n.y0){
        tmpind = 1:n.t0+n.t0*(kk-1);
        z95.suplogLam = apply(na.omit(abs(gammaptb[tmpind,,drop=F]-gammahat[tmpind])/gamma.sd[tmpind]),2,max)
        z95.suplogLam = quantile(z95.suplogLam,0.95,na.rm=T)
        cutoff.sup = c(cutoff.sup,z95.suplogLam)
      }
    }else{cutoff.sup=NULL}
    wgtk.ptb.naive = wgtk*diag(nv)
    gammaptb.naive = PTB.Gamma.FUN.explicit(wgtk.ptb.naive,cbind(xk,dk,wgtk,yk),logLamhat=log(LamCox.t0),betahat=betahat,
                                            regularize=regularize,newy=y0.cut,betaonly=betaonly)
    gamma.sd.naive = sqrt(apply(gammaptb.naive^2,1,sum))
    out = c(cutoff.sup,gammahat,gamma.sd,gamma.sd.naive)
  }
  out	
}

PI.k.FUN <- function(tt,ebyi,xi,yi,wgt.ptb=NULL,k0=0)
{
  out = ebyi; yi=as.matrix(yi); py = ncol(yi); 
  if(!is.null(wgt.ptb))
  { 
    wgt.ptb=as.matrix(wgt.ptb); pv=ncol(wgt.ptb); 
    if(k0==1){yi=yi[,rep(1:py,pv),drop=F]; out=out*wgt.ptb[,rep(1:pv,rep(py,pv))] # v1*z1, v1*z2, ...
    }else{ out=out*wgt.ptb }
  }
  if(k0==1){out=out*yi}
  if(k0==2){out=c(out)*yi[,rep(1:py,py)]*yi[,rep(1:py,rep(py,py))]}
  as.matrix(sum.I(tt,"<=",xi,out))
}

PTB.Gamma.FUN.explicit <- function(wgt.ptb,data,logLamhat,betahat,regularize,newy,betaonly)
{
  xi = data[,1]; di = data[,2]; wgti = data[,3]; yi = data[,-(1:3),drop=F]; ebyi = c(exp(yi%*%betahat)); py = ncol(yi)
  tmpind = di==1; tj = xi[tmpind]; wgt.ptb.j = wgt.ptb[tmpind,]; pv = ncol(wgt.ptb)
  pi0.tj   = c(PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=0))/sum(wgti) 
  pi1.tj   =   PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=1)/sum(wgti)
  Ahat = PI.k.FUN(tj,ebyi,xi,yi,wgti,k0=2)/sum(wgti)*pi0.tj - pi1.tj[,rep(1:py,py)]*pi1.tj[,rep(1:py,rep(py,py))]
  Ahat = matrix(apply(Ahat/pi0.tj^2,2,sum),ncol=py)
  term1 = t(wgt.ptb.j)%*%yi[tmpind,]
  term2 = t(wgt.ptb.j)%*%(pi1.tj/pi0.tj)
  pi1.tj.ptb = PI.k.FUN(tj,ebyi,xi,yi,wgt.ptb=wgt.ptb,k0=1)/sum(wgti) ## py x nptb
  pi0.tj.ptb = PI.k.FUN(tj,ebyi,xi,yi,wgt.ptb=wgt.ptb,k0=0)/sum(wgti);## nptb x 1 
  term3 = t(rep(1,length(tj)))%*%((pi1.tj.ptb*pi0.tj-pi0.tj.ptb[,rep(1:pv,rep(py,pv))]*pi1.tj[,rep(1:py,pv)])/pi0.tj^2)
  term3 = t(matrix(term3,nrow=py))
  if(regularize){
    betastar = matrix(0,nrow=py,ncol=pv); tmpind = betahat!=0
    if(sum(tmpind)>0){
      betastar[tmpind,] = solve(Ahat[tmpind,tmpind])%*%(t(term1-term2-term3)[tmpind,])  # py x nptb
    }
  }else{
    betastar = solve(Ahat)%*%t(term1 - term2 - term3)  # py x nptb
  }
  out=betastar
  if(!betaonly)
  {
    term1 = sum.I(t0,">=",tj,wgt.ptb.j/pi0.tj)     # n.t0 x nptb
    term2 = sum.I(t0,">=",tj,pi0.tj.ptb/pi0.tj^2)  # n.t0 x nptb
    term3 = sum.I(t0,">=",tj,pi1.tj/pi0.tj^2);term3 = term3%*%betastar      # n.t0 x py
    logLamstar = (term1-term2-term3)/sum(wgti)/exp(logLamhat)  # n.t0 x nptb
    if(is.null(newy)){newy=rep(0,py)}; newy = as.matrix(newy); n.y0 = ncol(newy) ## py x n.y0
    logLamstar = logLamstar[rep(1:n.t0,n.y0),] + (t(newy)%*%betastar)[rep(1:n.y0,rep(n.t0,n.y0)),]
    out=rbind(logLamstar,out)
  }
  out
}

NCC.Cox.FUN.clogit <- function(data,V.IND,y0.cut=NULL,Zmatch.ind=NULL)
{
  xi = data[,1]; di = data[,2]; vi = data[,3]; yi = data[,-c(1:3,Zmatch.ind),drop=F]; pp=ncol(yi)
  if(is.null(y0.cut)){y0.cut=matrix(0,nrow=pp,ncol=1)}; n.y0 = ncol(y0.cut); 
  V.IND = V.IND[order(V.IND[,1]),]
  nv = sum(vi); tj = V.IND[,1]; n.tj = V.IND[,2]; ncase=length(tj); ind.Rji = V.IND[,-(1:2),drop=F] 
  data.cluster = data.frame("Xj" = xi[c(t(ind.Rji))],"status"=rep(c(1,rep(0,m.match)),ncase),
                            "stratum"=rep(1:ncase,rep(m.match+1,ncase)),"yy"=yi[c(t(ind.Rji)),])
  tmpfm = paste0("Surv(Xj,status)~",paste(names(data.cluster)[-(1:3)],collapse="+"),"+strata(stratum)")						  
  prob.Stilde.tj = 1/pmax(choose(n.tj-1,m.match),1); prob.Pi.hat.tj = n.tj/(m.match+1)
  data.cluster = na.omit(data.cluster) ## remove NA observations due to clusters with less than m.match controls
  clogit.fit = coxph(as.formula(tmpfm),data=data.cluster)    
  betahat = clogit.fit$coef; 
  Pi.tj.kk = function(dat0,kk,bhat)
  {
    yy = as.matrix(dat0[,-(1:3)]); ind.clus = dat0$stratum; eby = c(exp(yy%*%bhat)); pp=ncol(yy)
    ind.select = cumsum(table(ind.clus)) ## index to select cumsum of yy ##
    switch(as.character(kk),"0"={tmpout=eby},"1"={tmpout=eby*yy},
           "2"={tmpout=eby*yy[,rep(1:pp,pp)]*yy[,rep(1:pp,rep(pp,pp))]})
    tmpout = as.matrix(cumsum2(tmpout))[ind.select,,drop=F]
    tmpout - rbind(0,tmpout[-nrow(tmpout),,drop=F])
  }
  Stilde.tj =  Pi.hat.tj = as.list(1:3); names(Stilde.tj)=names(Pi.hat.tj)=paste0("K",0:2)
  for(kk in 1:3){
    tmp.pi.tj = Pi.tj.kk(data.cluster,kk-1,bhat=betahat)
    Stilde.tj[[kk]] = tmp.pi.tj*prob.Stilde.tj; Pi.hat.tj[[kk]] = tmp.pi.tj*prob.Pi.hat.tj }
  ################################################
  ## Calculating Information Matrix for Betahat ##
  ################################################
  V.Rj.tj = Stilde.tj$K2/c(Stilde.tj$K0) - Stilde.tj$K1[,rep(1:pp,pp)]*Stilde.tj$K1[,rep(1:pp,rep(pp,pp))]/c(Stilde.tj$K0^2)
  Ibb.inv = solve(matrix(apply(V.Rj.tj,2,sum),ncol=pp)); sd.beta = sqrt(diag(Ibb.inv))
  ###############################################################
  ## Calculating Cumulative Hazard at y0.cut (pp x n.y0 matrix)##
  ###############################################################
  ind.t0 = pmax(sum.I(t0,">=",tj),1); 
  hhat.tj.y0 = VTM(exp(c(t(betahat)%*%y0.cut)),ncase)/c(Pi.hat.tj$K0)            ## n.tj x n.y0 
  Lamhat.t0.y0 = cumsum2(hhat.tj.y0)[ind.t0,]									   ## n.t0 x n.y0 
  
  bhat.tj.y0 = VTM(c(y0.cut),ncase)-(Pi.hat.tj$K1/c(Pi.hat.tj$K0))[,rep(1:pp,n.y0)] ## n.tj x (pp x n.y0)
  bhat.tj.y0 = hhat.tj.y0[,rep(1:n.y0,rep(pp,n.y0))]*bhat.tj.y0
  bhat.t0.y0 = cumsum2(bhat.tj.y0)[ind.t0,]                                      ## n.t0 x (pp x n.y0)
  bhat.t0.y0 = matrix(t(bhat.t0.y0),nrow=pp)                                     ## pp x (n.y0 x n.t0)
  var.Lamhat.y0 = cumsum2(hhat.tj.y0^2)[ind.t0,]+ matrix(diag(t(bhat.t0.y0)%*%Ibb.inv%*%bhat.t0.y0),byrow=T,ncol=n.y0) ## n.t0 x n.y0
  
  out = cbind("Est"=c(log(Lamhat.t0.y0),betahat), "SE"=c(sqrt(var.Lamhat.y0)/Lamhat.t0.y0,sd.beta))
}
