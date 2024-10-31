
library(ismev)
library(lmomco)
library(FAdist)
library(Rsolnp)
#library(numDeriv)
library(EnvStats) 
library(MASS)
# --------------------------------------

test = welmet.rgev(xdat=maha10, numr=5, alpha.qm=0.5, w.Hlme1=0.5, w.Hqm= 0.5,
                  numB= 200, qqt= c(.95,.99,.995))

test2 = welmet.rgev(maha10, numr=3)
#---------------------------------------------------------------------------------------

# --- main program -------------------------------------------               

welmet.rgev = function(xdat, numr=NULL, alpha.qm=0.5, w.Hlme1=0.5, w.Hqm= 0.5,
                       numB= 200, qqt= c(.95,.99,.995))
{

  z=list(); mle.rgev=list()
  numq=length(qqt)
  if( is.null(numr) ) numr = dim(xdat)[2]
  
  nsample=nrow(xdat)  # nrow(xdat) = sample size n; ncol(xdat) = r
                         
# ++++++ lme1 +++++++++++++
  
  lme1 = lme1.gev(datr=xdat, qqt=qqt)
  
  z$lme1.rl = lme1$lme1.rl
  z$lme1.theta = lme1$lme1.theta
  theta.BM = lme1$lme1.theta
  
# --- Parametric bootstrap for S inverse ------
  
    theta =matrix(NA, numB,3); lme.rl =matrix(NA, numB,numq)

    for (ib in 1:numB){

      Bsam = gen.rgev.hosking(par= z$lme1.theta, sim_r=numr, 
                              sim_n=nsample, sim_k=1)

      theta[ib,1:3] = lme1.gev(datr=Bsam[[1]], qqt=qqt)$lme1.theta
      lme.rl[ib,1:numq] = quagev(qqt, vec2par(theta[ib,1:3],'gev'))
    }
  
    Hth = cov(theta)
    if( det(Hth) <= 0 ){
      cat("trouble in cov of theta.lme1","\n")
    }
    Sinv = solve(Hth)
    z$PBse.lme1.theta = sqrt(c(Hth[1,1], Hth[2,2],Hth[3,3]))
    
    Hrl = cov(lme.rl)
    for (kq in 1:numq){
      z$PBse.lme1.rl[kq] = sqrt(Hrl[kq,kq])
    }

# +++++ welmet_QM ++++++++++
    
  wel.qm = QM.trsf(xdat, numr=numr, alpha.qm=alpha.qm, 
                   quant=qqt, theta.BM=theta.BM, Sinv=Sinv)

  z$welmet = wel.qm

# +++++++  rmle +++++++++
    
    mle.rgev = rgevmle.park(xdat, numr=numr, ntry=10, lowb= -1.0, 
                            const=T, qpro=qqt)
    
    z$rmle.rl = mle.rgev$rmle.rl 
    z$rmle.theta = mle.rgev$rmle.theta
#   if(z$rmle.theta[3] < -1.0) z$rmle.theta[3] = -1.0

# ++++  mle1 +++++++++++++  

    mle1 = gev.max.consT(xdat[,1], ntry=10, lowb=-1.0, const=F)$mle   # hosking style xi

   if( mle1[3] <= -1.0) {
     mle1= gev.max.consT(xdat[,1], ntry=10, lowb=-1.0, const=T)$mle   # hosking style xi
   }
  
  z$mle1.rl = quagev(qqt, vec2par(mle1,'gev'))
  z$mle1.theta = mle1
  
# +++++ Hybrid of rmle and lme +++++++

    z$rl.Hlme1 = w.Hlme1*z$rmle.rl + (1-w.Hlme1)*z$lme1.rl
    z$theta.Hlme1 = w.Hlme1*z$rmle.theta + (1-w.Hlme1)*z$lme1.th
  
    z$rl.Hqm = w.Hqm*z$rmle.rl + (1-w.Hqm)*wel.qm$welmet.rl
    z$theta.Hqm = w.Hqm*z$rmle.theta + (1-w.Hqm)*wel.qm$welmet.th

  return(z)
}

#--------------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

rlarg.fit.consT.stnry = function (xdat, r = dim(xdat)[2], init=NULL, ydat = NULL, mul = NULL, sigl = NULL, 
            shl = NULL, mulink = identity, siglink = identity, shlink = identity, 
            muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
            method = "Nelder-Mead", maxit = 10000, lowb=lowb, const=NULL, ...) 
  {
    
    lowb = -1.0

    z <- list()                             # coles style para
    npmu <- length(mul) + 1
    npsc <- length(sigl) + 1
    npsh <- length(shl) + 1
    z$trans <- FALSE
    in2 <- sqrt(6 * var(xdat[, 1]))/pi
    in1 <- mean(xdat[, 1]) - 0.57722 * in2
    if (is.null(mul)) {
      mumat <- as.matrix(rep(1, dim(xdat)[1]))
      if (is.null(muinit)) 
        muinit <- in1
    }
    else {
      z$trans <- TRUE
      mumat <- cbind(rep(1, dim(xdat)[1]), ydat[, mul])
      if (is.null(muinit)) 
        muinit <- c(in1, rep(0, length(mul)))
    }
    if (is.null(sigl)) {
      sigmat <- as.matrix(rep(1, dim(xdat)[1]))
      if (is.null(siginit)) 
        siginit <- in2
    }
    else {
      z$trans <- TRUE
      sigmat <- cbind(rep(1, dim(xdat)[1]), ydat[, sigl])
      if (is.null(siginit)) 
        siginit <- c(in2, rep(0, length(sigl)))
    }
    if (is.null(shl)) {
      shmat <- as.matrix(rep(1, dim(xdat)[1]))
      if (is.null(shinit)) 
        shinit <- 0.1
    }
    else {
      z$trans <- TRUE
      shmat <- cbind(rep(1, dim(xdat)[1]), ydat[, shl])
      if (is.null(shinit)) 
        shinit <- c(0.1, rep(0, length(shl)))
    }
    xdatu <- xdat[, 1:r, drop = FALSE]
    
#    init <- c(muinit, siginit, shinit)    # park modified these 2 lines
    init = init
    
    z$model <- list(mul, sigl, shl)
    z$link <- deparse(substitute(c(mulink, siglink, shlink)))
    u <- apply(xdatu, 1, min, na.rm = TRUE)
    
    rlarg.lik <- function(a) {
      mu <- mulink(drop(mumat %*% (a[1:npmu])))
      sc <- siglink(drop(sigmat %*% (a[seq(npmu + 1, length = npsc)])))
      xi <- shlink(drop(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)])))
      if (any(sc <= 0)) 
        return(10^6)
      y <- 1 + xi * (xdatu - mu)/sc
      if (min(y, na.rm = TRUE) <= 0) 
        l <- 10^6
      else {
        y <- (1/xi + 1) * log(y) + log(sc)
        y <- rowSums(y, na.rm = TRUE)
        l <- sum((1 + xi * (u - mu)/sc)^(-1/xi) + y)
      }
      l
    }
    
    #  x <- optim(init, rlarg.lik, hessian = TRUE, method = method, 
    #             control = list(maxit = maxit, ...))
    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    x <- optim(init, rlarg.lik, hessian = TRUE, method = c("L-BFGS-B"), 
               lower= c(-Inf, 0, -1.0), upper=c(Inf, Inf, -lowb),      # coles style para
               control = list(maxit = maxit, ...))
    #++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    mu <- mulink(drop(mumat %*% (x$par[1:npmu])))
    sc <- siglink(drop(sigmat %*% (x$par[seq(npmu + 1, length = npsc)])))
    xi <- shlink(drop(shmat %*% (x$par[seq(npmu + npsc + 1, length = npsh)])))
    z$conv <- x$convergence
    z$nllh <- x$value
    #  z$data <- xdat
    if (z$trans) {
      for (i in 1:r) z$data[, i] <- -log((1 + (as.vector(xi) * 
                                                 (xdat[, i] - as.vector(mu)))/as.vector(sc))^(-1/as.vector(xi)))
    }
    z$mle <- x$par
    z$cov <- solve(x$hessian)
    z$se <- sqrt(diag(z$cov))
    #  z$vals <- cbind(mu, sc, xi)
    z$r <- r
    if (show) {
      if (z$trans) 
        print(z[c(2, 3)])
      print(z[4])
      if (!z$conv) 
        print(z[c(5, 7, 9)])
    }
    class(z) <- "rlarg.fit"
    invisible(z)
    
    return(z)
  }
# --------------------------------------------------------------------------
#-------------------------------------------------------------  
rgevmle.park = function(xdat, numr=NULL, ntry=20, lowb=-1.0, const=T, qpro=NULL){
  
  zz=list(); k=list(); z=list()
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat[,1],ntry)
  nllh= rep(NA, ntry)
  
  tryCatch( 
    for(i in 1:nrow(init)){
      
     value= try( rlarg.fit.consT.stnry(xdat[,1:numr],r=numr, init=init[i,1:3], 
                                       show=F, lowb=lowb, const=F) 
                , silent=T)                                        # coles style para
     
     if(is(value)[1]=="try-error"){
       k[[i]] <- list(value=10^6)
#       cat("i try error= ", i,"\n" )
     }else{
       k[[i]] <- value
       nllh[i]= k[[i]]$nllh
     }
     
    } #for  
  ) #tryCatch
  
  selc_num = which.min(nllh)

  x  <-k[[selc_num]]

  z$conv <- x$conv
  z$nllh <- x$nllh
  z$mle <- x$mle             # coles style parameter
  z$mle[3] = - x$mle[3]       # Hosking style para
  

  if( z$mle[3] <= lowb & const==T ) {
    
    nllh= rep(NA, ntry)
    tryCatch( 
      for (i in 1:nrow(init)) {
        
        value= try( rlarg.fit.consT.stnry(xdat[,1:numr],r=numr, init=init[i,1:3], 
                                          show=F, lowb=lowb, const=T) 
                    , silent=T)                                        # coles style para
        
        if(is(value)[1]=="try-error"){
          k[[i]] <- list(value=10^6)
        }else{
          k[[i]] <- value
          nllh[i]= k[[i]]$nllh
        }
        
      } #for
    ) #tryCatch

    selc_num = which.min(nllh)
    x  <-k[[selc_num]]
    
    z$conv <- x$conv
    z$nllh <- x$nllh
    z$mle <- x$mle              # coles style parameter
    z$mle[3] = - x$mle[3]       # Hosking style para
    
  }
  
  zz$rmle.rl = quagev(qpro, vec2par(z$mle, type='gev'))
  zz$rmle.theta = z$mle                                # hosking style parameter
  
  return(zz)
}
#------------------------------------------------------------
gev.max.consT=function (xdat, ntry=20, lowb= -1.0, const=T) 
{
  z <- list();  k =list()           # hosking style para
  n=ntry

  nsample=length(xdat)
  z$nsample=nsample
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat,ntry)
  
  #-------------------------------------------------  
  gev.lik.max <- function(a) {
    
    mu <- a[1]      #mulink(mumat %*% (a[1:npmu]))
    sc <- a[2]      #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- a[3]      #shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    y <- (xdat - mu)/sc
    y <- 1 - xi * y        # park modify to negative, for xi in hosking
    
    for (i in 1:nsample){
      y[i] = max(0, y[i], na.rm=T) }
    
    if (any(y <= 0) || any(sc <= 0)) 
      return(10^6)
    
    if( abs(xi) >= 10^(-5) ) {ooxi= 1/xi       # park modify to xi in hosking
    }  else  {ooxi=sign(xi)*10^5}
    
    zz=nsample*(log(sc)) + sum( exp(ooxi *log(y)) ) + sum(log(y) * (1-(ooxi)) ) 
    
    return(zz)
  }
  #-------------------------------------------------------------
  tryCatch(
    for(i in 1:nrow(init)){
      
      value <- try(solnp(init[i,], fun=gev.lik.max, 
                         LB =c(-Inf,0,lowb),UB =c(Inf,Inf,1.0),     # hosking style para
                         control=list(trace=0, outer.iter=10,
                                      delta=1.e-7, inner.iter=40, tol=1.e-5) ))
      
      if(is(value)[1]=="try-error"){
        k[[i]] <- list(value=10^6)
      }else{
        k[[i]] <- value
      }
      
    } #for
  ) #tryCatch
  
  optim_value  <-data.frame(num=1:n,value=sapply(k, function(x) x$value[which.min(x$value)]))
  
  optim_table1 <-optim_value[order(optim_value$value),]
  selc_num  <- optim_table1[1,"num"]
  
  x  <-k[[selc_num]]
  
  #  mu <- x$par[1];  sc <- x$par[2];  xi <- x$par[3] 
  
  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$mle <- x$par                            # hosking style parameter
  
  return(z)
}
#------------------------------------------------------
#---------------------------------------------------------------------------------------
#----------------------------------------------------------------

 com.rl = function(xdat, td.hap, td.cbd, numr2=NULL, quant=NULL,
                  theta.BM=NULL, Sinv=NULL){

  z=list()

  numq=length(quant)
  cbd.rl = matrix(NA, numq, numr2);
  ma.rl = matrix(NA, numq, numr2)
  rk.gev = matrix(NA, nrow=numr2, ncol=3)
  ma.theta = rep(NA, 3)
  ld=matrix(NA, numr2, 3); gld= rep(NA, numr2)

#  td.hap.new = c(xdat[,1], td.hap)
  td.cbd.new = cbind(xdat[,1], td.cbd)
  
  for (kw in 1:(numr2) ) {
         kung = lmoms(td.cbd.new[,kw], nmom=3)
         
         if( are.lmom.valid(kung, checkt3t4=TRUE) == F){
           cbd.rl[,kw] = NA
           cat(" Invalid L-moms", "\n")
         }else{

           rk.gev[kw,1:3] = pargev(kung, checklmom=F)$para     # Lme

         } # end if mom.valid
         
       cbd.rl[,kw] = quagev(quant, vec2par(rk.gev[kw,],type="gev")  )

       ld[kw,1:3] = theta.BM[1:3]-rk.gev[kw,1:3]
       gld[kw] = exp( -( t(ld[kw,1:3]) %*% Sinv %*% ld[kw,1:3] )/2 )
  } # end for kw
  
  gld[1] =1.0

  id =seq(1,numr2)
  numid=length(id)

  if(numid==0){
    cat("check LME for BM data= ",  rk.gev[1,1:3],"\n")
    stop
    
  }else{
    
    wlme=rep(0, numr2)
    wlme[id] = gld[id]/sum(gld[id])
    
    z$gld = -2*log(gld)
    z$di = gld
    z$wlme = wlme
    z$welmet.rl = cbd.rl[,id] %*% wlme[id]     # rl with wlemet
    z$welmet.th = wlme[id] %*% rk.gev[id,]     # theta with welmet

    wt2= rep(1/numid, numr2)
    z$ma.rl = cbd.rl[,id] %*% wt2[id]     # rl with simple average
    z$ma.theta = wt2[id] %*% rk.gev[id,]  # theta with simple average
    
    z$each.rl = cbd.rl      # rl for each component
    z$each.theta = rk.gev

  }
  return(z)
 }
#------------------------------------------------------------
#--------------------------------------------------------------
#-------------------------------------------------------------------
QM.trsf = function(xdat, numr=NULL, alpha.qm=0.5, quant=NULL,
                  theta.BM=NULL, Sinv=NULL) { 
  
  z=list()
  alpha = alpha.qm
  nsample = nrow(xdat)
  ndim = dim(xdat)[2]
  count=0
  
  if(ndim ==2){ 
    cat("Qm is not available, ndim= ", ndim,"\n")
    quit
  }
  
  td.fin=matrix(NA, nrow=nsample, ncol=ndim)
  numq=length(quant)
  
  aver1= mean(xdat[,1]); std1= sqrt(var(xdat[,1]))
  td.hap=NULL
  td.cbd=NULL
  r1= xdat[,1]
  
  r1.gev = pargev(lmoms(r1, nmom=3),checklmom=F)

  numr2= numr-1
  if(numr >= ndim) {
 #   cat("numr should be < ndim. So changed to numr=ndim-1","\n")
    numr2 = ndim-1
  }
  
  for (kq in 2:numr2){

    rpre= xdat[,kq-1]
    r2= xdat[,kq]
    r3= xdat[,kq+1]

# using empirical cdf 'pemp' function from EnvStat package ------------------   
    
     r2c=  pemp(r2, r2)*alpha + pemp(r2, r3)*(1-alpha)

    r2c[ which(r2c >= 0.999999) ] = 0.999
    r2c[ which(r2c <= 0.000001) ] = 0.001
    
    td.fin[,kq] = quagev(r2c, r1.gev)
    
    td.hap= c(td.hap, td.fin[,kq])
    td.cbd= cbind(td.cbd, td.fin[,kq] )
    
  } # end for kq
  
  if(numr2 == numr-1){
    
    z= com.rl(xdat, td.hap, td.cbd, numr2=numr2, quant=quant,
              theta.BM=theta.BM, Sinv=Sinv)
    
  }else{
    cat("wrong in qm.trsf, numr2 is not same to numr-1 =", numr2,"\n")
  }

  return(z)
}
#-------------------------------------------------------------  
#-------------------------------------------------------------
  gen.rgev.hosking = function(par, sim_r, sim_n, sim_k ){
    
    sim_par_umat  <-array(runif(sim_n *sim_r *sim_k),c(sim_n, sim_r, sim_k))
    sim_par_umat1 <-lapply(seq(sim_k), function(x) sim_par_umat[ , , x]) # array to list
    sim_par_umat2 <-lapply(seq(sim_k), function(x) t(apply(sim_par_umat1[[x]],1,cumprod)))   # r-largest

    sim_par_sample <-lapply(seq(sim_k), function(x) quagev(f=sim_par_umat2[[x]][,1:sim_r],
                                                           vec2par(par, 'gev') ) )
                                        #     location = par[1],scale = par[2],shape = par[3]))
    
    return(sim_par_sample)
  }
#--------------------------------------------------------------
#------------------------------------------------------
  ginit.max <-function(data,ntry){
    
    n=ntry
    init <-matrix(rep(0,n*3),ncol=3)
    
    lmom_init = lmoms(data,nmom=5)
    lmom_est <- pargev(lmom_init)
    
    init[1,1]    <-lmom_est$para[1]
    init[1,2]    <-lmom_est$para[2]
    init[1,3]    <-lmom_est$para[3]
    
    maxm1=ntry; maxm2=maxm1-1
    init[2:maxm1,1] <- init[1,1]+ rnorm(n=maxm2,mean=0,sd = 5)
    init[2:maxm1,2] <- abs( init[1,2]+ rnorm(n=maxm2,mean=5,sd = 5)) +1
    init[2:maxm1,3] <- runif(n=maxm2,min= -0.5,max=0.5)
#    init[2:maxm1,2] = max(0.1, init[2:maxm1,2])
    
    return(init)
  }
#----------------------------------------
lme1.gev = function (datr=NULL, qqt=NULL){
  
  z=list()
  kung= lmoms(datr[,1], nmom=3)
  
  if( are.lmom.valid(kung, checkt3t4=TRUE) == F){
    cat("----trouble in lme1.gev----","\n")
    z$lme1.rl = NA
    
  }else{
    hos= pargev(kung, checklmom=F)
    
    z$lme1.rl= quagev(qqt, hos)
    z$lme1.theta = hos$para
    
  }
  return(z)
}
