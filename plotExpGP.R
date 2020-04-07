# Ctrl output

if(fitOut$prior_PD == 0) {
  sink(file=paste0('Results/',tag,'_ctrl.txt'),append = TRUE)

  cat('\n ExpGP parameters:\n')
  fit = fitGP$fit
  if(fitGP$method == 'sample' | fitGP$method == 'vb') {
    pars = c('theta','yGP','lambda','sigma','br')
    print(fit,pars=pars)
    br = mean(rstan::extract(fit,'br')[[1]])
  } else {
    cat('theta  :',fit$par$theta,'\n')
    cat('yGP    :',fit$par$yGP,'\n')
    cat('lambda :',fit$par$lambda,'\n')
    cat('sigma  :',fit$par$sigma,'\n')
    cat('br     :',br<-fit$par$br,'\n')
  }
  cat('\n\n')

  # Probability Interval for Birge's ratio
  br = FitOCTLib::printBr(fit)
  sink()
}

# Plots
nMC = 100  # Nb spaghetti lines to plot (if any)

fit      = fitOut$fit
method   = fitOut$method
xGP      = fitOut$xGP
prior_PD = fitOut$prior_PD
# lasso    = fitOut$lasso

# Build tag from options
tag1 =''
if(prior_PD != 0)
  tag1=paste0(tag1,'_priPD')
tagOut = paste0('Results/',tag,'_ExpGP',tag1)

pars = c('theta','yGP','lambda','sigma','br','lp__')
if(prior_PD != 0)
  pars = c('theta','yGP','lambda','sigma','lp__')
S = as.matrix(fit,pars=pars)
pdf(file = paste0(tagOut,'_ctrl.pdf'))
print(rstan::traceplot(fit, inc_warmup=TRUE, pars = pars))
pairs(fit, pars = pars, gap=0)
par(mfrow=c(1,1),pty='m')
rgumlib::SAPlot(S)
dev.off()

png(filename = paste0(tagOut,'_results.png'),
    width=2000, height=1000*(2-prior_PD))
gPars$plot_title = 'Modulated exp. fit'
sum = FitOCTLib::plotExpGP(
  x, y, uy, ySpl, out=fitOut, modScale=modRange, gPars=gPars,
  dataType = dataType, br = br
)
gPars$plot_title = ' '
dev.off()
