outExpGP <- function(x, y, uy, out) {

  fit      = out$fit
  method   = out$method
  xGP      = out$xGP
  prior_PD = out$prior_PD
  # lasso    = out$lasso
  if(prior_PD != 0)
    return(NULL)

  if(method == 'sample') {

    theta   = rstan::extract(fit,'theta')[[1]]
    yGP     = rstan::extract(fit,'yGP')[[1]]
    sigma   = mean(rstan::extract(fit,'sigma')[[1]])
    resid   = rstan::extract(fit,'resid')[[1]]
    mod     = rstan::extract(fit,'m')[[1]]
    dL      = rstan::extract(fit,'dL')[[1]]
    lp      = rstan::extract(fit,'lp__')[[1]]

    map     = which.max(lp)
    mod0 = mod[map,]
    mExp = expDecayModel(x,theta[map,1:3],dataType)
    dL0  = dL[map,]

    sum  = rstan::summary(fit,
                          use_cache=FALSE,
                          probs=c()
    )$summary[1:3,c(1,3)]


  } else {

    theta = fit$par$theta
    mod   = fit$par$m
    resid = fit$par$resid
    dL    = fit$par$dL
    yGP   = fit$par$yGP
    sigma = fit$par$sigma

    mExp = expDecayModel(x,theta,dataType)
    mod0 = mod
    dL0 = dL

    pars = c('theta')
    opt = list()
    for (par in pars)
      opt[[par]] = fit$par[[par]]
    opt = unlist(opt,use.names = TRUE)

    se = rep(NA, length(opt))
    if(!is.null(fit$hessian)) {
      H = fit$hessian
      tags = colnames(H)
      tags = gsub('\\.','',tags)
      colnames(H) = rownames(H) = tags
      se = list()
      for (par in names(opt))
        se[[par]]  = sqrt(-1/H[par,par])
      se = unlist(se)
    }
    sum = data.frame(mean = opt, sd = se)

  }

  return(
    list(
      sum = sum,
      mExp = mExp,
      mod = mod0,
      dL = dL0
    )
  )
}
