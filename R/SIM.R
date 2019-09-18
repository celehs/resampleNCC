SIM.FUN <- function(nn=5000,DepCen=F, Zmatch=F,a0=NULL,model="cox",p.add.noise=0)
{
  ## ===================================================== ##
  ## Z is the matching vector; Y is the marker of interest ##
  ## ===================================================== ##
  y0i = rnorm(nn)/2; sd.y = 1; yi = (y0i+matrix(rnorm(nn*length(gam.true),sd=sd.y),nrow=nn))
  if(yes.Z){
    Zi = 1*(pnorm(yi[,1]+rnorm(nn))>0.5); 
    Zi = cbind(Zi,pnorm(yi[,2]+rnorm(nn))*5); Zi = round(Zi) 
    if(is.null(a0)){a0 = rep(0,ncol(Zi))}
  }else{Zi = NULL}    
  if(!DepCen){
    ci = pmin(rgamma(nn,2,2)+0.2,runif(nn,0.5,2))
    # ci = rexp(nn,1)
  }else{
    ##ci = pmin(exp(apply(as.matrix(yi),1,sum)-2),runif(nn,0.5,2))
    ci = apply(yi,1,sum); ci = pmin(rgamma(nn,2,5)*exp(ci/2)+0.1,runif(nn,0.5,2))
  }
  switch(model,
         "cox"={ti = -c(yi%*%gam.true)+log(-log(runif(nn))); ti = exp(ti/2+1.5)},
         "aft"={ti = exp(yi%*%gam.true+rnorm(nn)+2.6)})
  xi = pmin(ti,ci); di = 1*(ti<=ci); bi = vi = di; ind.case = (1:nn)[di==1]; ##print(mean(di))
  tivec = nivec = rep(0,sum(di)); 
  
  IND.ik = matrix(NA,nrow=sum(di),ncol=m.match+1) ## id of cases and the corresponding controls
  Iik0 = matrix(0,nrow=sum(di),ncol=nn-sum(di))  ## indicates whether xk is in the matched risk set for ti
  for(l in 1:sum(di)){
    tmpind = ind.case[l]; risksetind = xi>xi[tmpind];  IND.ik[l,1]=tmpind;  
    ## =========================================================================== ##
    ## if matching, additional constraint of |Zi - Zl| <= a0 needs to be satisfied ##
    ## =========================================================================== ##
    if(Zmatch){risksetind = risksetind&(apply(abs(Zi-VTM(Zi[tmpind,],nn))<=VTM(a0,nn),1,prod)==1)}  
    Iik0[l,] = 1*risksetind[di==0]
    risksetind = (1:nn)[risksetind]; nl = length(risksetind);nivec[l]=nl;tivec[l]=ti[tmpind];  
    ## =========================================================================== ##
    ## if riskset is empty, no control would be selected                           ##
    ## if riskset size < m.match, only select # of available for Finite Population ##
    ## =========================================================================== ##
    if(length(risksetind)>0){
      controlind = as.numeric(sample(as.character(risksetind),min(m.match,nl))); ##Finite Population
      IND.ik[l,-1] = c(controlind,rep(NA,m.match-length(controlind)))
      vi[controlind]=1
    }
    
  }
  Iik0 = Iik0[,vi[di==0]==1];  if(p.add.noise>0){yi = cbind(yi,y0i+matrix(rnorm(nn*p.add.noise,sd=sd.y),nrow=nn))}
  list("data"=cbind(xi,di,vi,yi,Zi), "Vi.k0.IND"=cbind(tivec,nivec,IND.ik),"Iik0"=Iik0) 
}    
