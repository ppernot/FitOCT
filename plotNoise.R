# Ctrl output
sink(file=paste0('Results/',tag,'_ctrl.txt'))
cat('\n Noise fit parameters:\n')
a = fits$theta
for(i in 1:length(a))
  cat(paste0('a_',i,' : '),signif(a[i],3),'\n')
cat('\n\n')
sink()

# Plots
tagOut  = paste0('Results/',tag,'_EstimateNoise')
png(filename = paste0(tagOut,'_results.png'),
    width=2000, height=1000)
gPars$plot_title = 'Noise estimation'
SNR = FitOCTLib::plotNoise(x, y, uy, ySpl, gPars, dataType)
gPars$plot_title = ''
dev.off()
