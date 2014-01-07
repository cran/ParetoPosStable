coef.PPSfit <-
  function (object, ...) 
  {
    if (!class(object) == "PPSfit") stop("Object must belong to class PPS")            
    print.default(format(object$estimate), print.gap = 2, quote = FALSE)
    invisible(object)
  }

dPPS <-
  function(x,lam,sc,v){
    if (lam<=0 | sc<=0 | v<=0) salida<-0
    else{
      salida<-rep(NA,length(x))
      salida[(x<sc)==1]<-rep(0,sum(x<sc))
      salida[(x<sc)==0]<-lam*v*(log(x[(x<sc)==0]/sc))^(v-1)*x[(x<sc)==0]^(-1)*exp(-lam*(log(x[(x<sc)==0]/sc))^v)
    }
    return(salida)
  }

GoF.PPSfit <-
  function (PPSfit, k = 2000, show.iters = TRUE) {
    if (is.null(PPSfit$Pareto)) PPSfit$Pareto <- FALSE
    if (PPSfit$Pareto==TRUE & !is.null(PPSfit$sigma)) pars <- c(as.numeric(PPSfit$estimate), PPSfit$sigma, 1)
    if (PPSfit$Pareto==TRUE & is.null(PPSfit$sigma)) pars <- c(as.numeric(PPSfit$estimate), 1)
    if (PPSfit$Pareto==FALSE & is.null(PPSfit$sigma)) pars <- as.numeric(PPSfit$estimate)
    if (PPSfit$Pareto==FALSE & !is.null(PPSfit$sigma)) pars <- c(PPSfit$estimate[[1]], PPSfit$sigma, PPSfit$estimate[[2]])
    datos <- PPSfit$obs
    estadistico.ks <- ks.test(datos, pPPS, lam = pars[1], sc = pars[2], v = pars[3])$statistic
    estadistico.ad <- ad.test(datos, pPPS, lam = pars[1], sc = pars[2], v = pars[3])$statistic
    rango <- function(y) -(length(y)*ecdf(y)(y)-(length(y)+1))
    #z <- log(log(datos/pars[2]))
    #y <- log(-log(rango(datos)/(PPSfit$n+1)))
    z <- log(log(datos[datos!=pars[2]]/pars[2]))
    y <- log(-log(rango(datos[datos!=pars[2]])/(PPSfit$n+1)))
    estadistico.pps <- sum((y-(log(pars[1])+pars[3]*z))^2)
    metodo <- PPSfit$estim.method
    n <- PPSfit$n
    ks.sample <- c()
    ad.sample <- c()
    pps.sample <- c()
    for (i in 1:k){
      if (show.iters == TRUE) cat("Step",i,"out of", k, "\n")
      datos.sim <- rPPS(n, pars[1], pars[2], pars[3])#Simulo datos con los parámetros estimados (H0)
      if (PPSfit$Pareto == FALSE & is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, Pareto = PPSfit$Pareto)
        pars.sim <- as.numeric(ajuste.sim$estimate)
      }
      if (PPSfit$Pareto == FALSE & !is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, sigma = PPSfit$sigma, Pareto = PPSfit$Pareto)
        pars.sim <- c(ajuste.sim$estimate[[1]], PPSfit$sigma, ajuste.sim$estimate[[2]])   
      }
      if (PPSfit$Pareto == TRUE & is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, Pareto = PPSfit$Pareto)
        pars.sim <- c(as.numeric(ajuste.sim$estimate),1)
      }
      if (PPSfit$Pareto == TRUE & !is.null(PPSfit$sigma)){
        ajuste.sim <- PPS.fit(datos.sim, estim.method = metodo, sigma = PPSfit$sigma, Pareto = PPSfit$Pareto)
        pars.sim <- c(as.numeric(ajuste.sim$estimate), PPSfit$sigma, 1) 
      }
      test.ks<-ks.test(datos.sim, pPPS, lam = pars.sim[1], sc = pars.sim[2], v = pars.sim[3])
      ks.sample <- c(ks.sample, test.ks$statistic)
      test.ad<-ad.test(datos.sim, pPPS, lam = pars.sim[1], sc = pars.sim[2], v = pars.sim[3])
      ad.sample <- c(ad.sample, test.ad$statistic)
      #z <- log(log(datos.sim/pars.sim[2]))
      #y <- log(-log(rango(datos.sim)/(PPSfit$n+1)))
      z <- log(log(datos.sim[datos.sim!=pars.sim[2]]/pars.sim[2]))
      y <- log(-log(rango(datos.sim[datos.sim!=pars.sim[2]])/(PPSfit$n+1)))
      pps.sample <- c(pps.sample, sum((y-(log(pars.sim[1])+pars.sim[3]*z))^2))
    }
    return(list(ks.statistic = estadistico.ks, ad.statistic = estadistico.ad, pps.statistic = estadistico.pps, ks.p.value = sum(estadistico.ks<ks.sample)/(k+1), ad.p.value = sum(estadistico.ad<ad.sample)/(k+1), pps.p.value = sum(estadistico.pps<pps.sample)/(k+1), PPSfit))
  }


hPPS <-
  function(x,lam,sc,v){
    if (lam<=0 | sc<=0 | v<=0) salida<-0
    else{
      salida<-rep(NA,length(x))
      salida[(x<sc)==1]<-rep(0,sum(x<sc))
      salida[(x<sc)==0]<-lam*v*(log(x[(x<sc)==0]/sc))^(v-1)*x[(x<sc)==0]^(-1)
    }
    return(salida)
  }

logLik.PPSfit <-
  function (object, ...) 
  {
    val <- object$loglik
    attr(val, "nobs") <- object$n
    attr(val, "df") <- length(object$estimate)
    class(val) <- "logLik"
    val
  }

pareto.fit <-
  function(x, estim.method = "MLE", sigma = NULL, start, ...){
    if (!is.element(estim.method,c("MLE", "OLS"))) stop("The estimation method in method.estim is not OLS nor MLE")
    
    if (missing(x) || length(x) == 0L || mode(x) != "numeric" || any(x<=0)) stop("'x' must be a non-empty numeric vector of positive values")
    if (any(!is.finite(x))) stop("'x' contains missing or infinite values")
    
    n <- length(x)
    Fn<-ecdf(x)
    y<-log(1-n*Fn(x)/(n+1))
    xName <- paste(deparse(substitute(x), 500), collapse = "\n")
    if (missing(start)) start <- NULL
    dots <- names(list(...))
    
    if (estim.method == "OLS"){
      if (is.null(sigma)){
        sce<-function(sigma){
          xx<-log(sigma)-log(x)
          ajuste <- lm(y ~ -1 + xx)
          sce<-sum(ajuste$residuals^2)
        }
        sigma <- optimize(f=sce,interval=c(0,min(x)))$minimum
        xx<-log(sigma)-log(x)
        ajuste <- lm(y ~ -1 + xx)
        lambda <- ajuste$coefficients[[1]]
        loglik <- n*(log(lambda)+lambda*log(sigma))-(lambda+1)*sum(log(x))
        return(structure(list(estimate = list(lambda = lambda, sigma = sigma), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS", Pareto = TRUE), class = "PPSfit"))
      }
      if (!is.null(sigma)){
        sc<-sigma
        xx<-log(sigma)-log(x)
        ajuste <- lm(y ~ -1 + xx)
        lambda <- ajuste$coefficients[[1]]
        loglik <- n*(log(lambda)+lambda*log(sigma))-(lambda+1)*sum(log(x))
        return(structure(list(estimate = list(lambda = lambda), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS", sigma = sigma, Pareto = TRUE), class = "PPSfit"))
      }
    }
    
    if (estim.method == "MLE"){
      missing(start)
      sigma <- min(x)
      lambda <- n/sum(log(x/sigma))
      loglik <- n*(log(lambda)+lambda*log(sigma))-(lambda+1)*sum(log(x))
      return(structure(list(estimate = list(lambda = lambda, sigma = sigma), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "MLE", Pareto = TRUE), class = "PPSfit"))
    }
  }


plot.PPSfit <-
  function (x, which = 1:4, ask = prod(par("mfcol")) < length(which) && dev.interactive(), ylim, breaks, ...)
  {
    if (!class(x) == "PPSfit") 
      stop("Object must belong to class PPSfit")
    if (ask) {
      op <- par(ask = TRUE)
      on.exit(par(op))
    }
    par(mar = c(6, 4, 4, 2) + 0.1)
    show <- rep(FALSE, 4)
    show[which] <- TRUE
    lam <- x$estimate[[1]]
    if (is.null(x$Pareto)) x$Pareto <- FALSE
    if (x$Pareto==FALSE & !is.null(x$sigma)){
      sc <- x$sigma
      v <- x$estimate[[2]]
    }
    if (x$Pareto==FALSE & is.null(x$sigma)){   
      sc <- x$estimate[[2]]
      v <- x$estimate[[3]] 
    }  
    if (x$Pareto==TRUE & !is.null(x$sigma)){
      sc <- x$sigma
      v <- 1
    }
    if (x$Pareto==TRUE & is.null(x$sigma)){   
      sc <- x$estimate[[2]]
      v <- 1
    }
    obs <- x$obs
    obsName <- x$obsName
    if (missing(breaks)) breaks <- hist(obs, plot = FALSE, right = FALSE)$breaks
    else breaks <- hist(obs, plot = FALSE, right = FALSE, breaks)$breaks
    PPSDens <- function(y) dPPS(y, lam, sc, v)
    PPSProbs <- function(y) pPPS(sort(y), lam, sc, v)
    if (missing(ylim)) {
      ymax <- 1.06 * max(PPSDens(obs), na.rm = TRUE)
      ylim <- c(0, ymax)
    }
    if (show[1]) {
      hist.default(obs, right = FALSE, freq = FALSE, 
                   ylim = ylim, 
                   main=NULL,
                   xlab="x",
                   breaks, 
                   ...)       
      curve(PPSDens, min(breaks) - 1, max(breaks) + 1, add = TRUE, 
            ylab = NULL,col=2,lty=2,lwd=2)
      #legend("topright", legend = substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    if (show[2]) {
      ecdf2<-function(y) length(y)*ecdf(y)(y)/(length(y)+1)
      plot(obs, ecdf2(obs), 
           xlab="x",
           ylab="Cumulative probability", 
           ...)
      lines(sort(obs), PPSProbs(sort(obs)),col=2,lty=2,lwd=2)
      #legend("bottomright",legend=substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    if (show[3]) {
      rango<-function(y) -(length(y)*ecdf(y)(y)-(length(y)+1))
      plot(log(obs), log(rango(obs)), 
           xlab="log(x)",
           ylab="log(rank(x))",
           ...)
      lines(log(sort(obs)), log(x$n*(1-PPSProbs(sort(obs)))),col=2,lty=2,lwd=2)
      #legend("bottomleft", legend = substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    if (show[4]) {
      rango<-function(y) -(length(y)*ecdf(y)(y)-(length(y)+1))
      if (min(obs) < sc+0.000001) obs <- obs[which(obs >= sc+0.000001)]
      plot(log(log(obs/sc)), log(-log(rango(obs)/(x$n+1))), 
           xlab="log(log(x/scale))",
           ylab="log(-log(rank(x)/(n+1)))",
           ...)
      lines(log(log(sort(obs)/sc)), log(-log(1-PPSProbs(sort(obs)))),col=2,lty=2,lwd=2)
      #legend("topleft",legend=substitute(paste("PPS(",lambda, " = ", nn1 ,", ", sigma, " = ", nn2, ", ", nu, " = ", nn3,")"), list(nn1=round(lam,3),nn2=round(sc,3),nn3=round(v,3))), bty = "n", col = 2, lty = 2, lwd = 2)
    }
    invisible()
  }


pPPS <-
  function(x,lam,sc,v){
    if (lam<=0 | sc<=0 | v<=0) salida<-0
    else{
      salida<-rep(NA,length(x))
      salida[(x<sc)==1]<-rep(0,sum(x<sc))
      salida[(x<sc)==0]<-1-exp(-lam*(log(x[(x<sc)==0]/sc))^v)
    }
    return(salida)
  }

PPS.fit <-
  function(x, estim.method = "MLE", sigma = NULL, start, Pareto=FALSE, ...){
    
    if (Pareto==TRUE) pareto.fit(x, estim.method, sigma, start, ...)
    else{
      
#      if (!is.element(estim.method,c("MLE", "OLS", "iMLE", "LMOM", "KSMIN"))) stop("The estimation method in method.estim is not OLS, MLE, iMLE, LMOM nor KSMIN")
      if (!is.element(estim.method,c("MLE", "OLS", "iMLE", "LMOM"))) stop("The estimation method in method.estim is not OLS, MLE, iMLE nor LMOM")
      
      if (missing(x) || length(x) == 0L || mode(x) != "numeric" || any(x<=0)) stop("'x' must be a non-empty numeric vector of positive values")
      if (any(!is.finite(x))) stop("'x' contains missing or infinite values")
      
      if (missing(start)) start <- NULL
      xName <- paste(deparse(substitute(x), 500), collapse = "\n")
      dots <- names(list(...))
      n <- length(x)
      Fn <- ecdf(x)
      y <- log(-log(1-n*Fn(x)/(n+1)))
      
      
      # OLS ---------------------------------------------------------------------
      OLS.estimate<-function(x, sigma){
        
        if (is.null(sigma)){
          sce <- function(sigma){
            z <- log(x/sigma) 
            v.ols <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
            lam.ols <- exp(mean(y)-v.ols*mean(log(z)))
            sce <- sum(lm(y~log(z))$residuals^2)
          }
          sigma <- optimize(f=sce,interval=c(0,min(x)))$minimum
          z <- log(x/sigma) 
          v <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
          lam <- exp(mean(y) - v*mean(log(z)))
          loglik <- n*(log(lam)+log(v))+(v-1)*sum(log(z))-lam*sum(z^v)-sum(log(x))
          return(structure(list(estimate = list(lambda = lam, sigma = sigma, nu = v), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS"), class = "PPSfit"))
        }
        
        if (!is.null(sigma)){
          sc <- sigma
          z <- log(x/sc)
          v <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
          lam <- exp(mean(y)-v*mean(log(z)))
          loglik <- n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sc)))-lam*sum(log(x/sc)^v)-sum(log(x))
          return(structure(list(estimate = list(lambda = lam, nu = v), sigma = sigma, loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS"), class = "PPSfit"))
        }
      }
      
      
      # MLE ---------------------------------------------------------------------
      MLE.estimate <- function(x, sigma, start){
        
        if (!is.null(sigma)){
          sc <- sigma
          z <- log(x/sc)
          v.ols <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
          lam.ols <- exp(mean(y)-v.ols*mean(log(z)))
          obj <- function(v) (1/v+(1/n)*sum(log(z))-sum(z^v*log(z))/sum(z^v))^2
          optimum <- nlm(p = v.ols, f = obj)
          v<-optimum$estimate
          lam <- 1/(sum(z^v)/n)
          loglik <- n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sc)))-lam*sum(log(x/sc)^v)-sum(log(x))
          return(structure(list(estimate = list(lambda = lam, nu = v), sigma = sigma, loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "MLE", convergence = optimum$code), class = "PPSfit"))
        }
        
        if (is.null(sigma)){
          if (is.null(start)){
            sc <- min(x)/2
            Fn <- ecdf(x)
            y <- log(-log(1-n*Fn(x)/(n+1)))
            z <- log(x/sc)
            v <- sum((log(z)-mean(log(z)))*(y-mean(y)))/sum((log(z)-mean(log(z)))^2)
            lam <- exp(mean(y)-v*mean(log(z)))
            start <- list(lambda = lam, sigma = sc, nu = v)
            start <- start[!is.element(names(start), dots)]
          }
          nm <- names(start)
          if (sum(is.element(nm,c("lambda","sigma","nu")))<3) stop("'start' doesn't specify the names of the PPS parameters as 'lambda', 'sigma' and 'nu'")
          
          logL <- function(p){
            lam <- exp(p[1])
            sc <- (1/(1+exp(-p[2])))*min(x)
            v <- exp(p[3])
            -(n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sc)))-lam*sum(log(x/sc)^v)-sum(log(x)))
          }
          p0 <- c(log(start[[1]]), -log((min(x)-start[[2]])/start[[2]]), log(start[[3]]))
          res <- optim(p0, logL, hessian=FALSE)
          l.est <- exp(res$par[1])
          s.est <- (1/(1+exp(-res$par[2])))*min(x)
          v.est <- exp(res$par[3])
          estimate <- list(lambda = l.est, sigma = s.est, nu = v.est)
          return(structure(list(estimate = estimate, loglik = -res$value, n = n, obs = x, obsName = xName, estim.method = "MLE", convergence = res$convergence), class = "PPSfit"))
        }
      }
      
      
      # iMLE --------------------------------------------------------------------
      iMLE.estimate<-function(x){
        logver <- function(sigma){
          z <- log(x/sigma)
          obj <- function(v) 1/v+(1/n)*sum(log(z))-sum(z^v*log(z))/sum(z^v)
          cota.inf <- 1/(max(log(z))-mean(log(z)))
          cota.sup <- 1/(mean(log(z))-min(log(z)))
          if (cota.inf*cota.sup>=0){
            cota.inf <- 0
            cota.sup <- 100
          }
          v.mle <- uniroot(obj, lower = cota.inf, upper = cota.sup)$root
          lam.mle <- 1/(sum(z^v.mle)/n)
          return(- (n * (log(lam.mle)+log(v.mle)) + (v.mle-1)*sum(log(log(x/sigma))) - lam.mle*sum(log(x/sigma)^v.mle) - sum(log(x))))
        }
        sigma.mle.lv <- optimize(f=logver,interval=c(0,min(x)))$minimum
        z <- log(x/sigma.mle.lv)
        obj <- function(v) 1/v+(1/n)*sum(log(z))-sum(z^v*log(z))/sum(z^v)
        v.mle.lv <- uniroot(obj, lower = 0, upper = 100)$root
        lam.mle.lv <- 1/(sum(z^v.mle.lv)/n)
        loglik.optimum <- n*(log(lam.mle.lv)+log(v.mle.lv))+(v.mle.lv-1)*sum(log(log(x/sigma.mle.lv)))-lam.mle.lv*sum(log(x/sigma.mle.lv)^v.mle.lv)-sum(log(x))
        return(structure(list(estimate = list(lambda = lam.mle.lv, sigma = sigma.mle.lv, nu = v.mle.lv), loglik = loglik.optimum, n = n, obs = x, obsName = xName, estim.method = "iMLE"), class = "PPSfit"))
      }
      
      
      # LMOM --------------------------------------------------------------------
      LMOM.estimate <- function(x, start){
        code1<-1
        code2<-1
        code3<-1
        l1 <- function(l,s,v) {
          obj1 <- function(p) s*exp((-(1/l)*log(1-p))^(1/v))
          if((integrate(obj1,0,1,stop.on.error = FALSE)$message=="the integral is probably divergent")==TRUE){
            print("Warning: L1 moment evaluated where it is probably divergent")
            return(code1 <- 0)
          }
          else integrate(obj1,0,1)$value
        }
        l2 <- function(l,s,v) {
          obj2 <- function(p) (s*exp((-(1/l)*log(1-p))^(1/v)))*(2*p-1)
          if((integrate(obj2,0,1,stop.on.error = FALSE)$message=="the integral is probably divergent")==TRUE){
            print("Warning: L2 moment evaluated where it is probably divergent")
            return(code2 <- 0)
          }
          else integrate(obj2,0,1)$value
        }
        l3 <- function(l,s,v) {
          obj3 <- function(p) (s*exp((-(1/l)*log(1-p))^(1/v)))*(6*p^2-6*p+1)
          if((integrate(obj3,0,1,stop.on.error = FALSE)$message=="the integral is probably divergent")==TRUE){
            print("Warning: L3 moment evaluated where it is probably divergent")
            return(code3 <- 0)
          }
          else integrate(obj3,0,1)$value
        }
        if ((code1==0)|(code2==0)|(code3==0)) print("The LMOM method has no solution")
        else{
          obj <- function(p) {
            l<-exp(p[1])
            s<-min(x)/(1+exp(-p[2]))
            v<-exp(p[3])
            (samlmu(x)[[1]]-l1(l,s,v))^2 + (samlmu(x)[[2]]-l2(l,s,v))^2 + (samlmu(x)[[3]]-l3(l,s,v))^2
          }
          if (is.null(start)) pars<-iMLE.estimate(x)$estimate
          else pars <- start
          minimo<-optim(par = c(log(pars[[1]]), -log(min(x)/pars[[2]]-1), log(pars[[3]])), fn = obj)
          lam<-exp(minimo$par[1])
          sigma<-min(x)/(1+exp(-minimo$par[2]))
          v<-exp(minimo$par[3])
          minimum <- minimo$value
          loglik <- (n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sigma)))-lam*sum(log(x/sigma)^v)-sum(log(x)))
          return(structure(list(estimate = list(lambda = lam, sigma = sigma, nu = v), minimo = minimum, loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "LMOM",convergence = minimo$convergence), class = "PPSfit"))
        }
      }
      
      
      # KSMIN -------------------------------------------------------------------
#       KSMIN.estimate <- function(x, start){
#         if (is.null(start)) pars<-iMLE.estimate(x)$estimate
#         else pars <- start
#         obj <- function(p) {
#           l<-exp(p[1])
#           s<-min(x)/(1+exp(-p[2]))
#           v<-exp(p[3])
#           ks.test(x,pPPS,l,s,v)$statistic
#         }
#         minimo<-optim(par = c(log(pars[[1]]), -log(min(x)/pars[[2]]-1), log(pars[[3]])), fn = obj)
#         lam<-exp(minimo$par[1])
#         sigma<-min(x)/(1+exp(-minimo$par[2]))
#         v<-exp(minimo$par[3])
#         minimum <- minimo$value
#         loglik <- (n*(log(lam)+log(v))+(v-1)*sum(log(log(x/sigma)))-lam*sum(log(x/sigma)^v)-sum(log(x)))
#         return(structure(list(estimate = c(lam, sigma, v), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "KSMIN", convergence = minimo$convergence), class = "PPSfit"))
#       }
      
      
      # Return ------------------------------------------------------------------
      if (estim.method == "MLE") return(MLE.estimate(x, sigma, start))
      if (estim.method == "iMLE") return(iMLE.estimate(x))
      if (estim.method == "OLS") return(OLS.estimate(x, sigma))
      if (estim.method == "LMOM") return(LMOM.estimate(x, start))
#       if (estim.method == "KSMIN") return(KSMIN.estimate(x, start))
    }
  }


print.PPSfit <-
  function (x, digits = max(3, getOption("digits") - 3), ...) 
  {
    if (!class(x) == "PPSfit") {
      stop("Object must belong to class PPS")
    }
    cat("\nData:     ", x$obsName, "\n")
    cat("\n")
    if (!is.null(x$sigma)){
      cat("Sigma:\n")
      print.default(format(x$sigma, digits = digits), print.gap = 2, quote = FALSE)
    }
    cat("\n")
    cat("Parameter estimates:\n")
    print.default(format(x$estimate, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    if (is.null(x$Pareto)) x$Pareto <- FALSE
    if (x$Pareto == TRUE){
      cat("Pareto:\n")
      print.default(TRUE)
    }
    cat("\n")
    cat("Log-likelihood:\n")
    print.default(format(x$loglik, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    cat("Sample size:\n")
    print.default(format(x$n, digits = digits), print.gap = 2, quote = FALSE)
    cat("\n")
    cat("\nEstimation method:     ", x$estim.method, "\n") 
    invisible(x)
  }


qPPS <-
  function(p,lam,sc,v){
    if (lam<=0 | sc<=0 | v<=0) salida<-0
    else{
      salida<-rep(NA,length(p))
      salida[((p<0)|(p>1))==1]<-rep(0,sum(((p<0)|(p>1))))
      salida[((p<0)|(p>1))==0]<-sc*exp((-(1/lam)*log(1-p[((p<0)|(p>1))==0]))^(1/v))
    }
    return(salida)
  }

rPPS <-
  function(n,lam,sc,v){
    weib<-rweibull(n,shape=v,scale=1)
    sc*exp(lam^(-1/v)*weib)
  }

se.PPSfit <-
  function(PPSfit, k=2000, show.iters = TRUE){
    if (is.null(PPSfit$Pareto)) PPSfit$Pareto <- FALSE
    if (PPSfit$Pareto == FALSE){
      if (is.null(PPSfit$sigma)){
        lambdas <- c()
        sigmas <- c()
        nus <- c()
        i <- 1
        while (i < k){
          if (show.iters == TRUE) cat("Step",i,"out of", k, "\n")
          datos <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
          ajuste <- PPS.fit(datos, estim.method = PPSfit$estim.method, start = PPSfit$estimate)
          lambdas <- c(lambdas, ajuste$estimate[[1]])
          sigmas <- c(sigmas, ajuste$estimate[[2]])
          nus <- c(nus, ajuste$estimate[[3]])
          i <- i+1
        }
        se.lambda <- sqrt(sum((lambdas - PPSfit$estimate[[1]])^2)/k)
        se.sigma <- sqrt(sum((sigmas - PPSfit$estimate[[2]])^2)/k)
        se.nu <- sqrt(sum((nus - PPSfit$estimate[[3]])^2)/k)
        return(list(se.lambda = se.lambda, se.sigma = se.sigma, se.nu = se.nu))
      }
      else {
        lambdas <- c()
        nus <- c()
        i <- 1
        while (i < k){
          if (show.iters == TRUE) cat(c("Step ",i,"out of ",k, "\n"))
          datos <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
          ajuste <- PPS.fit(datos, estim.method = PPSfit$estim.method, sigma = PPSfit$sigma, start = PPSfit$estimate)
          lambdas <- c(lambdas, ajuste$estimate[[1]])
          nus <- c(nus, ajuste$estimate[[2]])
          i <- i+1
        }
        se.lambda <- sqrt(sum((lambdas - PPSfit$estimate[[1]])^2)/k)
        se.nu <- sqrt(sum((nus - PPSfit$estimate[[2]])^2)/k)
        return(list(se.lambda = se.lambda, se.nu = se.nu))
      }
    }
    if (PPSfit$Pareto == TRUE){
      if (is.null(PPSfit$sigma)){
        lambdas <- c()
        sigmas <- c()
        i <- 1
        while (i < k){
          if (show.iters == TRUE) cat("Step",i,"out of", k, "\n")
          datos <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
          ajuste <- PPS.fit(datos, estim.method = PPSfit$estim.method, Pareto = TRUE)
          lambdas <- c(lambdas, ajuste$estimate[[1]])
          sigmas <- c(sigmas, ajuste$estimate[[2]])
          i <- i+1
        }
        se.lambda <- sqrt(sum((lambdas - PPSfit$estimate[[1]])^2)/k)
        se.sigma <- sqrt(sum((sigmas - PPSfit$estimate[[2]])^2)/k)
        return(list(se.lambda = se.lambda, se.sigma = se.sigma))
      }
      else {
        lambdas <- c()
        i <- 1
        while (i < k){
          if (show.iters == TRUE) cat(c("Step ",i,"out of ",k, "\n"))
          datos <- sample(x = PPSfit$obs, size = PPSfit$n, replace = TRUE)
          ajuste <- PPS.fit(datos, estim.method = PPSfit$estim.method, sigma = PPSfit$sigma, Pareto = TRUE)
          lambdas <- c(lambdas, ajuste$estimate[[1]])
          i <- i+1
        }
        se.lambda <- sqrt(sum((lambdas - PPSfit$estimate[[1]])^2)/k)
        return(list(se.lambda = se.lambda))
      }
    }
  }
