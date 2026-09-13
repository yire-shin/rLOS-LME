# library(lmomco)
# library(EnvStats)
# library(ismev)

# xdat <-readRDS("sancheong.RDS")
# xdat <-na.omit(xdat[,-1])
# s1=Sys.time()
#   z= welmet.rgev(xdat, opt.weight.hybrid=FALSE) 
# Sys.time()-s1
# 
# s1=Sys.time()
#   z= welmet.rgev(xdat, opt.weight.hybrid=TRUE) 
# Sys.time()-s1

# --- main program ------------------------------------------------              

welmet.rgev = function(xdat, numr=NULL, alpha.qm=0.5, 
                       w.Hlme1=0.5, w.Hqm= 0.5,
                       numB= 200, qqt= c(.95,.99,.995,.999),
                       ntry=5, opt.weight.hybrid=TRUE)
{

  z=list()
  mle.rgev=list()
  numq=length(qqt)
  xdat <- as.matrix(xdat)
  
  if(is.null(numr)) numr <- ncol(xdat)
  
  stopifnot(
    numr >= 1L,
    numr <= ncol(xdat),
    numB >= 2L
  )

  nsample=nrow(xdat)  # nrow(xdat) = sample size n; ncol(xdat) = r
                         
# ++++++ lme1 +++++++++++++
  
  lme1 = lme1.gev(datr=xdat, qqt=qqt)
  z$lme1.rl = lme1$lme1.rl
  z$lme1.theta = lme1$lme1.theta

# --- Parametric bootstrap for S inverse ------
  
  theta =matrix(NA, numB,3)
  lme.rl =matrix(NA, numB,numq)
  
  theta.boot=z$lme1.theta
      
  Bsam = gen.rgev.hosking(par= theta.boot, sim_r=1,
                          sim_n=nsample, sim_k=numB)
      
    theta.BM = theta.boot

    for (ib in 1:numB){

      fit= lme1.gev(datr=as.matrix(Bsam[[ib]]), qqt=qqt)
      theta[ib,1:3] = fit$lme1.theta
      lme.rl[ib,1:numq] = fit$lme1.rl 
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
    
  wel.qm=list()  
  if(numr==1){
    wel.qm$welmet.rl = z$lme1.rl 
    wel.qm$welmet.th = z$lme1.theta
    
  }else if(numr >= 2){

    wel.qm = QM.trsf(xdat, numr=numr, alpha.qm=alpha.qm, 
                   quant=qqt, theta.BM=theta.BM, Sinv=Sinv)
  }
  z$welmet = wel.qm
  
# ++++  mle1 +++++++++++++  
  mle1 = gev.max.consT(xdat[,1], ntry=ntry, lowb=-1.0)$mle   # hosking style xi
  
  z$mle1.rl = quagev(qqt, vec2par(mle1,'gev'))
  z$mle1.theta = mle1

# +++++++  rmle +++++++++

  if(numr >= 2){
    mle.rgev = rgevmle.park(xdat, numr=numr, ntry=ntry, lowb= -1.0, 
                            const=TRUE, qpro=qqt, 
                            start.para=z$mle1.theta)

    z$rmle.rl = mle.rgev$rmle.rl 
    z$rmle.theta = mle.rgev$rmle.theta
    
  }else if(numr ==1){
    z$rmle.rl = z$mle1.rl 
    z$rmle.theta = z$mle1.theta
  }

# +++++ Hybrid of rmle and lme +++++++
    
    if(opt.weight.hybrid == TRUE){
      
      w.Hqm = ( 0.5359 +0.0053*nsample + 0.4098*z$rmle.theta[3] 
                - 0.5161*(z$rmle.theta[3])^2 )
      w.Hlme1 = (0.3826 + 0.0025*nsample 
                 + 1.8770*(z$rmle.theta[3])^2 )
      
      w.Hqm   <- pmin(1, pmax(0, w.Hqm))
      w.Hlme1 <- pmin(1, pmax(0, w.Hlme1))
    }
    
    z$rl.Hlme1 = w.Hlme1*z$rmle.rl + (1-w.Hlme1)*z$lme1.rl
    z$theta.Hlme1 = w.Hlme1*z$rmle.theta + (1-w.Hlme1)*z$lme1.theta
  
    z$rl.Hqm = w.Hqm*z$rmle.rl + (1-w.Hqm)*wel.qm$welmet.rl
    z$theta.Hqm = w.Hqm*z$rmle.theta + (1-w.Hqm)*wel.qm$welmet.theta
    
    z$w.Hlme1= w.Hlme1
    z$w.Hqm = w.Hqm

  return(z)
}
#--------------------------------------------------------------
#  The following function is just a simple modification of 
#  'rlarg.fit' function in the 'ismev' package.
#--------------------------------------------------------------

rlarg.fit.consT.stnry = function (xdat, r = dim(xdat)[2], init=NULL, ydat = NULL, mul = NULL, sigl = NULL, 
            shl = NULL, mulink = identity, siglink = identity, shlink = identity, 
            muinit = NULL, siginit = NULL, shinit = NULL, show = TRUE, 
            method = "Nelder-Mead", maxit = 1000, lowb=lowb, const=NULL, ...) 
  {

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

    # +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++  
    
    x <- optim(init, rlarg.lik, hessian = FALSE, method = c("L-BFGS-B"), 
               lower= c(-Inf, 1e-8, lowb), upper=c(Inf, Inf, 1),      # coles style para
               control = list(maxit = maxit, 
                              factr = 1e7, pgtol = 0))
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
    # z$cov <- solve(x$hessian)
    # z$se <- sqrt(diag(z$cov))
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
# ----------------------------------------------------------------
#-------------------------------------------------------------  
rgevmle.park = function(xdat, numr=NULL, ntry=5, lowb=lowb, 
                        const=TRUE, qpro=NULL, start.para=NULL){
  
  zz=list(); k=list(); z=list()
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat[,1],ntry)
  if(!is.null(start.para)) init[ntry,1:3] = start.para
  nllh= rep(NA, ntry)
  
  tryCatch( 
    for(i in 1:nrow(init)){
      
     value= try( rlarg.fit.consT.stnry(xdat[,1:numr],r=numr, init=init[i,1:3], 
                                       show=F, lowb=lowb, const=FALSE) 
                , silent=TRUE)                                        # coles style para
     
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
  z$mle <- x$mle             # coles style parameter
  z$mle[3] = - x$mle[3]       # Hosking style para
  

  if( z$mle[3] <= lowb & const==TRUE ) {
    
    nllh= rep(NA, ntry)
    tryCatch( 
      for (i in 1:nrow(init)) {
        
        value= try( rlarg.fit.consT.stnry(xdat[,1:numr],r=numr, init=init[i,1:3], 
                                          show=F, lowb=lowb, const=TRUE) 
                    , silent=TRUE)                                        # coles style para
        
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
  
  zz$conv = z$conv
  zz$nllh = z$nllh
  zz$rmle.rl = quagev(qpro, vec2par(z$mle, type='gev'))
  zz$rmle.theta = z$mle                                # hosking style parameter
  
  return(zz)
}
#------------------------------------------------------------
gev.max.consT=function (xdat, ntry=5, lowb= lowb, 
                        const=TRUE) 
{
  z <- list();  k =list()           # hosking style para
  n=ntry

  nsample=length(xdat)
  z$nsample=nsample
  
  init= matrix(0, nrow=ntry, ncol=3)
  init <- ginit.max(xdat,ntry)
  
  #--------------------------------------------------------- 
  # The following function is a simple modification of 
  # 'gev.lik' function in the 'ismev' package.
  # --------------------------------------------------------
  gev.lik.max <- function(a) {
    
    mu <- a[1]      #mulink(mumat %*% (a[1:npmu]))
    sc <- a[2]      #siglink(sigmat %*% (a[seq(npmu + 1, length = npsc)]))
    xi <- a[3]      #shlink(shmat %*% (a[seq(npmu + npsc + 1, length = npsh)]))
    
    y <- 1 - xi * (xdat - mu)/sc       # park modify to negative, for xi in hosking
    
    # for (i in 1:nsample){
    #   y[i] = max(0, y[i], na.rm=T) }
    
    if (any(y <= 0) || sc <= 0) 
      return(10^6)
    
    if( abs(xi) >= 10^(-5) ) {ooxi= 1/xi       # park modify to xi in hosking
    }  else  {ooxi=sign(xi)*10^5}
    
    zz=nsample*(log(sc)) + sum( exp(ooxi *log(y)) ) + sum(log(y) * (1-(ooxi)) ) 
    
    return(zz)
  }
  #-------------------------------------------------------------
  tryCatch(
    for(i in 1:nrow(init)){
      
    value= try( optim(
        init[i, ], gev.lik.max, method = "L-BFGS-B",
        lower = c(-Inf, 1e-8, lowb),
        upper = c( Inf, Inf, 1),
        control = list( maxit = 1000, factr = 1e7)
      ) )
      
      # value <- try(solnp(init[i,], fun=gev.lik.max, 
      #                    LB =c(-Inf,0,lowb),UB =c(Inf,Inf,1.0),     # hosking style para
      #                    control=list(trace=0, outer.iter=10,
      #                            delta=1.e-7, inner.iter=40, 
      #                            tol=1.e-5) ))
      
      if(is(value)[1]=="try-error"){
        k[[i]] <- list(value=10^6)
      }else{
        k[[i]] <- value
      }
      
    } #for
  ) #tryCatch
  
  optim_value  <-data.frame(num=1:n,value=sapply(k, 
                            function(x) x$value[which.min(x$value)]))
  
  optim_table1 <-optim_value[order(optim_value$value),]
  selc_num  <- optim_table1[1,"num"]
  
  x  <- k[[selc_num]]

  z$conv <- x$convergence
  z$nllh <- x$value[which.min(x$value)]
  z$mle <- x$par                            # hosking style parameter
  
  return(z)
}
#------------------------------------------------------
#----------------------------------------------------------------
#----------------------------------------------------------------

 com.rl = function(xdat, td.hap, td.cbd, numr2=NULL, quant=NULL,
                  theta.BM=NULL, Sinv=Sinv){

  zc=list()

  numq=length(quant)
  cbd.rl = matrix(NA, numq, numr2)
  rk.gev = matrix(NA, nrow=numr2, ncol=3)
  gld= rep(NA, numr2)

  td.cbd.new = cbind(xdat[,1], td.cbd)
  
  for (kw in 1:(numr2) ) {
       rkw=td.cbd.new[,kw]
         kung = lmoms(rkw, nmom=3)
         
         if( are.lmom.valid(kung, checkt3t4=TRUE) == F){
           cbd.rl[,kw] = NA
           cat(" Invalid L-moms", "\n")
         }else{

            rk.gev[kw,1:3] = pargev(kung, checklmom=F)$para     # Lme

         } # end if lmom.valid
         
       cbd.rl[,kw] = quagev(quant, vec2par(rk.gev[kw,],type="gev")  )

       d= theta.BM - rk.gev[kw,]
       gld[kw] = exp( -0.5 * drop(crossprod( d, Sinv %*% d ) ) )
                   
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

    zc$welmet.rl = cbd.rl[,id] %*% wlme[id]     # rl with wlemet
    zc$welmet.th = wlme[id] %*% rk.gev[id,]     # theta with welmet

    zc$each.rl = cbd.rl      # rl for each component
    zc$each.theta = rk.gev

  }
  return(zc)
 }
#------------------------------------------------------------
#--------------------------------------------------------------
#-------------------------------------------------------------------
QM.trsf = function(xdat, numr=NULL, alpha.qm=NULL, quant=NULL,
                  theta.BM=NULL, Sinv=Sinv) { 
  
  z=list()
  alpha = alpha.qm
  ndim = dim(xdat)[2]
  
  if(ndim ==2){ 
    cat("Qm is not available, ndim= ", ndim,"\n")
    quit
  }
  
  r1= xdat[,1]
  
  r1.gev=list()
    r1.gev = pargev(lmoms(r1, nmom=3),checklmom=F)
    
  numr2= numr-1
  if(numr >= ndim) {
    numr2 = ndim-1
  }
  
  if(is.null(alpha)){
    stop("alpha should be specified when alpha_optim==FALSE")
    
  }else if(!is.null(alpha)){  
    
    td= make.td(xdat, numr2, r1.gev, alpha)
 
  # }else if(alpha_optim==TRUE){
  #   
  #   alpha= alfin = alpha_find(xdat, numr2, r1.gev)
  #   td= make.td(xdat, numr2, r1.gev, alpha=alfin)

  } # end if
    
   if(numr2 == numr-1){
    
    z= com.rl(xdat, td$td.hap, td$td.cbd, numr2=numr2, quant=quant,
              theta.BM=theta.BM, Sinv=Sinv)
    
  }else{
    cat("wrong in qm.trsf, numr2 is not same to numr-1 =", numr2,"\n")
  }
  z$alpha= alpha
  return(z)
}
#-------------------------------------------------------------
make.td = function(xdat, numr2, r1.gev, alpha=NULL){  
  
  td.fin=matrix(NA, nrow=nrow(xdat), ncol=ncol(xdat))
  td.hap=NULL
  td.cbd=NULL
  
  for (kq in 2:numr2){
    
    r2= xdat[,kq]
    r3= xdat[,kq+1]
    
    # using empirical cdf 'pemp' function from EnvStats package 
    
    r2c= suppressWarnings( EnvStats::pemp(r2, r2)*alpha
                         + EnvStats::pemp(r2, r3)*(1-alpha) )
    
    r2c[ which(r2c >= 0.999999) ] = 0.999
    r2c[ which(r2c <= 0.000001) ] = 0.001
    
    td.fin[,kq] = quagev(r2c, r1.gev)
    
    td.hap= c(td.hap, td.fin[,kq])
    td.cbd= cbind(td.cbd, td.fin[,kq] )
    
  } # end for kq
  return(list(td.hap=td.hap,td.cbd=td.cbd))
}
  #-------------------------------------------------------------
  # gen.rgev.hosking = function(par, sim_r, sim_n, sim_k ){
  #   
  #   sim_par_umat  <-array(runif(sim_n *sim_r *sim_k),c(sim_n, sim_r, sim_k))
  #   sim_par_umat1 <-lapply(seq(sim_k), function(x) sim_par_umat[ , , x]) # array to list
  #   sim_par_umat2 <-lapply(seq(sim_k), function(x) t(apply(sim_par_umat1[[x]],1,cumprod)))   # r-largest
  #   
  #   sim_par_sample <-lapply(seq(sim_k), function(x) quagev(f=sim_par_umat2[[x]][,1:sim_r],
  #                                                          vec2par(par, 'gev') ) )
  #   # location = par[1], scale = par[2], shape = par[3]))
  #   
  #   return(sim_par_sample)
  # }
  
#--------------------------------------------------------------  
gen.rgev.hosking <- function(par, sim_r, sim_n, sim_k) {
    
    sim_par_umat <- array(
      runif(sim_n * sim_r * sim_k),
      c(sim_n, sim_r, sim_k)
    )
    
    sim_par_sample <- lapply(seq(sim_k), function(x) {
      
      umat <- sim_par_umat[, , x]
      
      if (sim_r > 1) {
        
        umat_cum <- t(apply(umat, 1, cumprod))
        
        quagev(
          f = umat_cum[, 1:sim_r],
          para = vec2par(par, "gev")
        )
        
      } else {
        
        matrix(
          quagev(
            f = umat,
            para = vec2par(par, "gev")
          ),
          ncol = 1
        )
      }
    })
    
    return(sim_par_sample)
  }
  
#--------------------------------------------------------------
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

    return(init)
  }
#----------------------------------------------------------
lme1.gev = function (datr=NULL, qqt=c(.99,.995)){
  
  z=list()
  kung= lmomco::lmoms(datr[,1], nmom=3)
  
  if( lmomco::are.lmom.valid(kung, checkt3t4=TRUE) == F){
    warning("Invalid L-moments in lme1.gev")
    z$lme1.rl = rep(NA, length(qqt))
    z$lme1.theta = rep(NA, 3)
    
  }else{
    hos= lmomco::pargev(kung, checklmom=F)
    
    z$lme1.rl= lmomco::quagev(qqt, hos)
    z$lme1.theta = hos$para
    
  }
  return(z)
}
#-------------------------------------------------------------
#-------------------------------------------------------------  
# alpha_find = function(xdat, numr2, r1.gev){
#   
#   alstar= rep(NA, numr2-2+1)
#   for (kq in 2:numr2){
#     
#     rpre= xdat[,kq-1]
#     r2= xdat[,kq]
#     r3= xdat[,kq+1]
#     
#     r2c1 = EnvStats::pemp(r2, r2)
#     r2c2 = EnvStats::pemp(r2, r3)
#     
#     ialp=0
#     KL=rep(NA,7)
#     
#     for (alpha in seq(0.2,0.8,0.1)){
#       
#       ialp=ialp+1
#       r2c = r2c1*alpha + r2c2*(1-alpha)
#       
#       r2c[ which(r2c >= 0.999999) ] = 0.999
#       r2c[ which(r2c <= 0.000001) ] = 0.001
#       
#       ztilda= quagev(r2c, r1.gev)
#       
#       KL[ialp]= emp_sym_kl(ztilda, xdat[,1])
#       
#       #        cat("alpha,KL=",alpha, KL[ialp],"\n")
#     } # end for
#     
#     al.id= which.min(KL)
#     alstar[kq-1]= (al.id+1)/10
#     #      cat("kq, alstar=",kq,alstar[kq-1],"\n")
#   } # end for kq
#   
#   mean(alstar)
# }
# #-------------------------------------------------------------
# emp_sym_kl <- function(x, y, breaks = "FD", eps = 1e-8) {
#   0.5 * emp_kl(x, y, breaks, eps) +
#     0.5 * emp_kl(y, x, breaks, eps)
# }
# 
# emp_kl <- function(x, y, breaks = "FD", eps = 1e-8) {  # 
#   # x: transformed pseudo-BM, e.g. F1^{-1}(T_{i alpha})
#   # y: original BM, e.g. Z_1
#   
#   br <- hist(c(x, y), breaks = breaks, plot = FALSE)$breaks
#   
#   px <- hist(x, breaks = br, plot = FALSE)$counts
#   py <- hist(y, breaks = br, plot = FALSE)$counts
#   
#   px <- px / sum(px)
#   py <- py / sum(py)
#   
#   
#   # avoid log(0)
#   px <- px + eps
#   py <- py + eps
#   px <- px / sum(px)
#   py <- py / sum(py)
#   
#   sum(px * log(px / py))
# }
