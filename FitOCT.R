rm(list = ls()); gc() # Clean environment

# Libraries ####
libs =c('parallel','rstan','inlmisc','FNN')
for (lib in libs ) {
  if(!require(lib,character.only = TRUE))
    install.packages(lib,dependencies=TRUE)
  library(lib,character.only = TRUE)
}
lib ='FitOCTLib'
if(!require(lib,character.only = TRUE))
  devtools::install_github("ppernot/FitOCTlib")
library(lib,character.only = TRUE)

# Options ####
options(
  mc.cores = parallel::detectCores(),
  width = 90,
  stringsAsFactors = FALSE
)
rstan_options(auto_write = TRUE)
set.seed(1234) # Initialise la graine du RNG

# Graphical parameters ####
gPars = list(
  cols    = inlmisc::GetColors(8),
  col_tr  = inlmisc::GetColors(8,alpha=0.1), # Light, for spaghetti curves
  col_tr2 = inlmisc::GetColors(8,alpha=0.4), # Darker, for legends...
  pty     = 's',
  mar     = c(3,3,1.6,.2),
  mgp     = c(2,.75,0),
  tcl     = -0.5,
  lwd     = 6,
  cex     = 3.5,
  cex.leg = 0.8,
  xlabel  = 'stromal depth (µm)',
  plot_title = '',
  graphTable = TRUE
)

# Control parameters ####

### Default values / Set values here
ctrlPars = list(
  depthSel    = NULL, # Otherwise c(xmin,xmax)
  dataType    = 2,    # Intensity
  subSample   = 1,
  smooth_df   = 15,
  method      = c('sample','optim','vb')[1],
  nb_warmup   = 500,
  nb_sample   = 1000,
  modRange    = 0.5,
  ru_theta    = 0.05,
  lambda_rate = 0.1,
  gridType    = 'internal',
  Nn          = 10,
  rho_scale   = 0.1,
  priPost     = TRUE, # Compare prior and posterior pdf ?
  priorType   = 'abc'
)

# Override parameters with control file
ctrlFile = 'ctrlParams.yaml'
if (file.exists(ctrlFile)) {
  ## Get from file
  lPars = rlist::list.load(ctrlFile)
  ## Expose parameters
  for (n in names(lPars))
    ctrlPars[[n]]= lPars[[n]]
}
cat('Configuration Parameters\n')
cat('------------------------\n')
str(ctrlPars, give.head=FALSE, give.length=FALSE)

# Expose parameters
for (n in names(ctrlPars))
  assign(n,rlist::list.extract(ctrlPars,n))



# RUN ####
dataDirs = c("DataWl","Data1","DataSynth")[3]
for (dataDir in dataDirs) {

  dataSets = list.dirs(path=dataDir,full.names = FALSE)[-1]

  # Results table: one table per dataDir
  resultsTable = data.frame(
    date = NA,
    tag = NA,
    SNR = NA,
    Mono_br = NA,
    Mono_brCI95 = NA,
    Mono_alert = NA,
    GP_br = NA,
    GP_brCI95 = NA,
    GP_alert = NA,
    C0 = NA,
    uC0 = NA,
    A0 = NA,
    uA0 = NA,
    Ls = NA,
    uLs = NA,
    eta = NA
  )

  for(dataSet in dataSets) {

    tag = paste0(dataDir,'_',dataSet)
    cat(tag,'------------------------','\n')

    ### Get ans select Data
    D = read.csv(paste0(dataDir,'/',dataSet,'/Courbe.csv'))
    C = FitOCTLib::selX(D[,1],D[,2],depthSel,subSample)
    x = C$x; y = C$y
    depth = max(x)-min(x)

    ### Estimate data uncertainty
    fits = FitOCTLib::estimateNoise(x, y, df = smooth_df)
    uy   = fits$uy           # Used by next stages
    ySpl = fits$ySmooth      # Used by plotMonoExp
    source ("./plotNoise.R") # Side effect: provides SNR

    ### MAP Inference of exponential decay parameters
    fitm   = FitOCTLib::fitMonoExp(x, y, uy, dataType = dataType)
    theta0    = fitm$best.theta   # Used by next stage
    cor.theta = fitm$cor.theta
    unc.theta = fitm$unc.theta
    source ("./plotMonoExp.R")
    mono_br = br

    if(is.null(br$alert)) {
      # Monoexp fit OK
      resultsTable = rbind(
        resultsTable,
        data.frame(
          date = date(),
          tag = tag,
          SNR = signif(SNR$SNR,3),
          Mono_br = signif(br$br,3),
          Mono_brCI95 = paste0(signif(br$CI95[1],3),'-',
                          signif(br$CI95[2],3)),
          Mono_alert = "Mono OK",
          GP_br = NA,
          GP_brCI95 = NA,
          GP_alert = NA,
          C0  = signif(theta0[1],4),
          uC0 = signif(unc.theta[1],2),
          A0  = signif(theta0[2],4),
          uA0 = signif(unc.theta[2],2),
          Ls  = signif(theta0[3],4),
          uLs = signif(unc.theta[3],2),
          eta = signif(theta0[3]/depth,4)
        )
      )
      df = cbind(
        x = x,
        y = y,
        uy = uy,
        mExp = fitm$fit$par$m
      )
      write.csv(
        df,
        file = paste0('Results/',tag,'_Curves.csv'),
        row.names = FALSE
      )
      next # No need to try fitExpGP
    }

    ### Prior
    priExp = FitOCTLib::estimateExpPrior(
      x, uy, dataType, priorType,
      out = fitm, ru_theta = ru_theta,
      eps = 1e-3
    )

    # - Posterior Distribution
    fitGP = FitOCTLib::fitExpGP(
      x, y, uy,
      dataType      = dataType,
      Nn            = Nn,
      gridType      = gridType,
      method        = method,
      theta0        = priExp$theta0,
      Sigma0        = priExp$Sigma0,
      lambda_rate   = lambda_rate,
      rho_scale     = ifelse(rho_scale==0, 1./Nn, rho_scale),
      nb_warmup     = nb_warmup,
      nb_iter       = nb_warmup + nb_sample,
      prior_PD      = 0,
      open_progress = FALSE
    )
    fitOut = fitGP
    source ("./plotExpGP.R")

    # - Compare prior/posterior marginal pdfs
    if(priPost & fitGP$method != 'optim')
       source("./priPost.R")

    source("./outExpGP.R")
    params = outExpGP(x, y, uy, out = fitGP)

    resultsTable = rbind(
      resultsTable,
      data.frame(
        date = date(),
        tag = tag,
        SNR = signif(SNR$SNR,3),
        Mono_br = signif(mono_br$br,3),
        Mono_brCI95 = paste0(signif(mono_br$CI95[1],3),'-',
                        signif(mono_br$CI95[2],3)),
        Mono_alert = "Mono PB !",
        GP_br = signif(br$br,3),
        GP_brCI95 = paste0(signif(br$CI95[1],3),'-',
                             signif(br$CI95[2],3)),
        GP_alert = ifelse(is.null(br$alert),"GP OK","GP PB !"),
        C0  = signif(params$sum[1,1],4),
        uC0 = signif(params$sum[1,2],2),
        A0  = signif(params$sum[2,1],4),
        uA0 = signif(params$sum[2,2],2),
        Ls  = signif(params$sum[3,1],4),
        uLs = signif(params$sum[3,2],2),
        eta = signif(params$sum[3,1]/depth,4)
      )
    )
    df = cbind(
      x = x,
      y = y,
      uy = uy,
      mExp = params$mExp,
      mGP  = params$mod,
      dL   = params$dL
    )
    write.csv(
      df,
      file = paste0('Results/',tag,'_Curves.csv'),
      row.names = FALSE
    )
  }

  # Save results table
  write.csv(
    resultsTable[-1,], # Remove first empty line
    file = paste0('Results/',dataDir,'_resultsTable.csv'),
    append = TRUE
  )
}

