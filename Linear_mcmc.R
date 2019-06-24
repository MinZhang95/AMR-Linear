##------------------------------------------------
## This code is for Model 1: i--index year; j--index observation within each year
## Likelihood: 
##  Y_{ij}     ~ I(c_{ij} = 1) N(\beta_{1i}, \sigma_1^2) + I(c_{ij} = 0)N(\beta_{0i}, \sigma_0^2)
##  \beta_{1i} ~ N(\mu_1, \tau_1^2), \beta_{0i} ~ N(\mu_0, \tau_0^2)
##   logit(p_i) =  alpha_i
##   alpha_i ~ N(\theta, \nu^2)
## Parameters: \sigma_1^2, \sigma_0^2, \mu_1, \mu_0, \tau_1^2, \tau_0^2, \theta,\nu^2 
## priors: independet conjugate prior 
##------------------------------------------------



##---------------------------------------
## Normal with known variance σ2
## data:  y[n]|theta, sigma^2 ~  iid N(theta, sigma^2)
## prior:  theta|a, b   ~ N(a, b)
## posterior: theta|y[], sigma^2 ~ N(m, v)
## where v = 1 / (n/sigma^2 + 1/b), m = v*(sum(y)/sigma^2 + a/b) 
##---------------------------------------
updateNormalMean <- function(y, sigma, a, b) {
  ## @y[n]
  ## @sigma >0, b >0, a  
  n <- length(y)
  v <- 1 / (n / (sigma * sigma) + 1 / b)
  m <- v * (sum(y) / (sigma * sigma) + a / b)
  rnorm(1, mean = m, sd = sqrt(v))
}

##---------------------------------------
## Normal with known mean μ
## data:  y[n]|mu, sigma^2 ~  iid N(mu, sigma^2)
## prior:  sigma^2|a, b   ~ IG(shape = a, rate = b), i.e. E(1/sigma^2) = a/b, var(1/sigma2) = a / b^2
## posterior: sigma^2|y[], mu  ~ IG(a + 0.5 * n, b + 0.5 * sum((y-mu)^2))
##---------------------------------------
updateNormalVarianceIG <- function(y, mu , a, b) {
  n <- length(y)
  sse <- sum((y - mu) ^ 2)
  gamma.shape <- a + 0.5 * n
  gamma.rate  <- b + 0.5 * sse
  1 / rgamma(1, shape = gamma.shape, rate = gamma.rate)
}


## Update unobserved continuous observation for each censored observation
## reference: r package: `mixAK`
rTNorm <- function(y, y1, mu, sigma, censoredType) {
  ## @y: (numeric)observed value 
  ## @y1: (numeric)upper bound for interval censored data, it's zero for other case 
  ## @mu: (numeric)normal mean
  ## @sigma: (numeric)normal sd 
  ## @censoredType: (int)0 -- right censored; 1 --- left-censored; 2 --- interval censored 
  ## Use inverse cdf sampling method to get a sample from censored normal 
  N_limit <- 8 
  N_prob0 <- 1e-15 
  N_prob1 <- 1 - 1e-15
  if (censoredType == 0) {
    ## right censored 
    ## centerd observation y 
    Zy <- (y - mu) / sigma
    Phiy <- pnorm(Zy, mean = 0, sd = 1)  
    U <- Phiy + (1 - Phiy) * runif(1)   ## U ~ Unif(Phiy, 1) ## 
    if (U > N_prob1) {
      ans <- mu + sigma * N_limit
    }  else {
      if (U < N_prob0) {
        ans <- mu - sigma * N_limit
      } else {
        ans <- mu + sigma * qnorm(U, 0, 1) 
      }
    }
  }
  
  if (censoredType == 1) {
    ## left-censored 
    Zy <- (y - mu) / sigma
    Phiy <- pnorm(Zy, 0, 1)
    U <- Phiy * runif(1)  ## U ~ Unif(0, Phiy)
    if (U < N_prob0) {
      ans <- mu - sigma * N_limit
    } else {
      if (U > N_prob1) {
        ans <- mu + sigma * N_limit
      } else {
        ans <- mu + sigma * qnorm(U, 0, 1)
      }
    }
  }
  
  if (censoredType == 2) {
    ## interval censored 
    Zy0 <- (y - mu) / sigma
    Zy1 <- (y1 - mu) / sigma
    PhiZy0 <- pnorm(Zy0, 0, 1)
    PhiZy1 <- pnorm(Zy1, 0, 1)
    U <- PhiZy0 + (PhiZy1 - PhiZy0) * runif(1)    ## U ~ Unif(PhiZy0, PhiZy1)
    if (U < N_prob0) {
      ans <- mu - sigma * N_limit 
    } else {
      if (U > N_prob1) {
        ans <- mu + sigma * N_limit 
      } else {
        ans <- mu + sigma * qnorm(U, 0, 1)
      }
    }
  }
  return(ans)
}

## Step 1: 
## updating the latent (censored) observations 
updateCensObs <- function(y, y1, beta, sigma, c, censor, yearLabel) {
  ## @y: [] observed values 
  ## @y1: [] upper limit bound for interval censored 
  ## @beta: [numUniqueYear, 2] matrix: normal mean, first column is non-resistant group, second column is resistant group  
  ## @sigma:[2] vecotr: normal sd for two compoenent, [non-resistant, resistant] 
  ## @c:[] componenent allocation 1 or 0 
  ## @censor:[] censor type 
  ## @yearLabel:[] group/year label for each observation 
  
  n <- length(y)
  yUpdated <- rep(0, n)
  ## sample ## 
  for (i in 1:n) {
    if (c[i] == 0) {
      ## first component 
      yUpdated[i] <- rTNorm(y[i], y1[i], beta[yearLabel[i], 1], sigma[1], censor[i])
    } else {## second component 
      yUpdated[i] <- rTNorm(y[i], y1[i], beta[yearLabel[i], 2], sigma[2], censor[i])
    }
  }
  return(yUpdated)
}


## update allocation: c_{ij}, P(c_{ij} = 1) = p_i --> resistant group using log trick 
updateAlloc <- function(y, p, beta, sigma, yearLabel) {
  ## @y: [] observation
  ## @p: [yearNum] mixture weight for resistant component(second one) in each year
  ## @beta[, 2]: normal mean for continous observation y
  ## @sigma[2]: sd 
  ## @yearLabel: [] 
  ## return updated allocation for each observation 
  
  cUpdated <- rep(0, length(y))
  for (i in 1:length(y)) {
    ##f1_log <- 1 - log(p[yearLabel[i]]) + dnorm(y[i], mean = beta[yearLabel[i], 1], sd = sigma[1], log = TRUE)
    f1_log <- log(1 - p[yearLabel[i]]) + dnorm(y[i], mean = beta[yearLabel[i], 1], 
                                               sd = sigma[1], log = TRUE)
    f2_log <- (log(p[yearLabel[i]])) + dnorm(y[i], mean = beta[yearLabel[i], 2], 
                                             sd = sigma[2], log = TRUE)
    ## normalize 
    mx = max(c(f1_log, f2_log))
    f1_log_new <- f1_log - mx 
    f2_log_new <- f2_log - mx 
    h <- exp(f2_log_new) / (exp(f1_log_new) + exp(f2_log_new))
    cUpdated[i] <- rbinom(1, 1, h)
  }
  return(cUpdated)
}


## Hierical structure for p_i 
##     logit(p_i) = alpha_i 
##      alpha_i ~ N(\theta, \nu^2)
## Prior: 
##      \theta ~ N(0, 10^6)
##      \nu^2 ~ IG(a, b)

## update alpha mean theta 
updateMuAlpha <- function(alpha, sigmaAlpha, a = 0, b = 10^6) {
  ## @alpha[numYear]
  updateNormalMean(alpha, sigmaAlpha, a, b)
}

## update alpha sd nu  
updateSigmaAlpha <- function(alpha, muAlpha, a = 0.0001, b = 0.0001) {
  ## alpha[numYear]
  sqrt(updateNormalVarianceIG(alpha, muAlpha, a, b))
}


logitInv <- function(x) {
  ##x[]:
  ## calculate logit inverse: logit(a) = x --> a = exp(x) / (1 + exp(x))
  n <- length(x)
  res <- sapply(seq(1, n), function(s){exp(x[s]) / (1 + exp(x[s]))})
  return(res)
}

## update mixture weight p (prob for second resistant group) 
updateWeight <- function(alpha) {
  ## alpha[numYear]
  ## 
  logitInv(alpha)
}

## Random walk MH for updating alpha 
log_prior <- function(alpha, muAlpha, sigmaAlpha) {
  ## alpha: numeric 
  dnorm(alpha, mean = muAlpha, sd = sigmaAlpha, log = T)
  
}

log_likilihood <- function(c, propIntercept) {
  
  numSecondComp <- sum(c == 1)
  totobs <- length(c)
  res <- propIntercept * numSecondComp - totobs * log(1 + exp(propIntercept))
  return(res)
}

updateAlpha <- function(c, alpha, muAlpha, sigmaAlpha, prop.sd = 1) {
  ## use random walk 
  ## @alpha: numeric  
  ## return alphaUpdate
  
  ## draw candidate from symmetric MV t-dist distribution 
  ##Alpha.star <- c(alpha + rt(1, 2, 3))
  ##Alpha.star <- rnorm(1, mean = alpha, sd = prop.sd)
  Alpha.star <- alpha + rt(1, df = 3)
  ## log prior 
  # lpold <- dnorm(alpha, mean = 0, sd = 100, log = T)
  # lpnew <- dnorm(Alpha.star, mean = 0, sd = 100, log = T)
  lpold <- log_prior(alpha, muAlpha, sigmaAlpha)
  lpnew <- log_prior(Alpha.star, muAlpha, sigmaAlpha)
  
  lold <- log_likilihood(c, alpha) 
  lnew <- log_likilihood(c, Alpha.star)
  ratio <- (lnew + lpnew) - (lold + lpold)
  acc <- ifelse(log(runif(1)) <= ratio, 1, 0)
  
  return(ifelse(acc == 1, Alpha.star, alpha))
}

updateAlphaVec <- function(c,  uniqueYearLength, yearLabel, alpha, 
                           muAlpha, sigmaAlpha,
                           prop.sd = 1) {
  ## use random walk 
  ## return alphaUpdate[]
  alphaUpdate <- rep(NA, uniqueYearLength)
  for (i in 1:uniqueYearLength) {
    cSubset <- c[yearLabel == i]
    alphaUpdate[i] <- updateAlpha(cSubset, alpha[i], muAlpha, sigmaAlpha, prop.sd) 
  }
  return(alphaUpdate)
}

## update mu: hierarchical normal mean for beta
##        \beta_{li} ~ N(\mu_l, \tau_l^2) i = 1 .. I
## prior:  \mu_l ~ N(a, b)
## posterior: \mu_l|. ~ N(m, v), where v = (I/tau^2 + 1 / b)^{-1}, m = v * ((sum beta)/tau^2 + a / b)
## revise for linear model
updateHierarchicalMu <- function(mu, beta, tau, a, b, yearLabel) {
  ## beta[numUniqueYear, 2]: 
  ## tau[2]: sd, 
  ## muUpdata[[1]]: mu.0i = gamma0 + gamma1*time.i
  ## muUpdata[[2]]: mu.1
  ## muUpdata[[3]]: gamma0, gamma1
  time <- c(1:length(unique(yearLabel))) # c(1, 2, ..., 16) always start with 1
  muUpdate <- list(rep(0, nrow(beta)), 0, rep(0, 2))
  gamma1 <- mu[[3]][2]
  # update gamma0
  gamma0Update <- updateNormalMean(beta[, 1]-gamma1*time, tau[1], a=0, b=10^6)
  gamma1Update <- updateNormalMean((beta[1, 1]-gamma0Update)/1, tau[1]/1, a=0, b=10^2)
  for (t in time[-1]) {
    gamma1Update <- updateNormalMean((beta[t, 1]-gamma0Update)/t, tau[1]/t, a=0, b=10^2)
  }
  muUpdate[[1]] <- gamma0Update + gamma1Update*time
  muUpdate[[2]] <- updateNormalMean(beta[, 2], tau[2], a=0, b=10^6)
  muUpdate[[3]] <- c(gamma0Update, gamma1Update)
  return(muUpdate)
}

## update beta: data normal mean 
## revise for linear model 
updateBeta <- function(y, c, sigma, mu, tau, uniqueYearLength, yearLabel) {
  ## y[]:
  ## c[]: allocation
  ## beta[numUniqueYear, 2]:
  ## mu[2]: beta hierarchical mean  
  ## sigma[2]: beta hierarchical sd 
  ## tau[2]
  betaUpdated <- matrix(0, nrow = uniqueYearLength, ncol = 2)
  
  for (i in 1:uniqueYearLength) {
    ## 
    ySubset <- y[yearLabel == i]
    cSubset <- c[yearLabel == i]
    ##  
    first.comp.obs <- cSubset == 0
    second.comp.obs <- cSubset == 1
    ## number of first component 
    num.first.comp <- sum(first.comp.obs)
    ## number of second component 
    num.second.comp <- sum(second.comp.obs)
    if (num.first.comp == 0) {
      ## no observations in the first component 
      ## draw from prior 
      betaUpdated[i, 1] <- rnorm(1, mean = mu[[1]][i], sd = tau[1])
      betaUpdated[i, 2] <- updateNormalMean(ySubset, sigma[2], mu[[2]], tau[2]) 
      
    } else if (num.second.comp == 0) {
      ## no observations in the second component 
      ## draw from prior 
      betaUpdated[i, 2] <- rnorm(1, mean = mu[[2]], sd = tau[2])
      betaUpdated[i, 1] <- updateNormalMean(ySubset, sigma[1], mu[[1]][i], tau[1])
    } else {
      betaUpdated[i, 1] <- updateNormalMean(ySubset[cSubset == 0], sigma[1], 
                                            mu[[1]][i], tau[1])
      betaUpdated[i, 2] <- updateNormalMean(ySubset[cSubset == 1], sigma[2], 
                                            mu[[2]], tau[2])
    }
  }
  return(betaUpdated)
}

## update mixture sigma: 
## sigma_l|. ~ IG(sigmaPriorA + 0.5 * \sum_i \sum_j I(c_ij = l), sigmaPriorB + 0.5 * \sum_i \sum_j I(c_ij = l)(y_ij - beta_il)^2) 
updateSigma <- function(y, c, yearLabel, uniqueYearLength, beta, 
                        sigmaPriorA = 0.0001, sigmaPriorB = 0.0001){
  ##y: []
  ##c: [], label componenet 
  ##yearLabel: []
  ##uniqueYearLength: int
  ##beta[uniqueYear, 2]: mix mean matrix 
  
  sigmaUpdated <- c(0, 0)
  ## get Mean for each obs 
  yMean <- sapply(seq(1, length(y)), function(s){ifelse(c[s] == 0, beta[yearLabel[s], 1], 
                                                        beta[yearLabel[s], 2])})
  
  sigmaUpdated[1] <- sqrt(updateNormalVarianceIG(y[c == 0], yMean[c == 0], sigmaPriorA,
                                                 sigmaPriorB))
  sigmaUpdated[2] <- sqrt(updateNormalVarianceIG(y[c == 1], yMean[c == 1], sigmaPriorA,
                                                 sigmaPriorB))
  return(sigmaUpdated)
}


## update hierarchical sd tau|. ~ IG(a + I/2, b + 0.5 * sum(beta - mu)^2)
## revise for linear model
updateTau <- function(beta, mu, a = 0.0001, b = 0.0001, yearLabel) {
  ## beta[uniqueYear, 2]: mix mean matrix 
  ## mu[2]: hierarchical mean 
  tauUpdated <- c(0, 0)
  tauUpdated[1] <- sqrt(updateNormalVarianceIG(beta[, 1], mu[[1]], a, b))
  tauUpdated[2] <- sqrt(updateNormalVarianceIG(beta[, 2], mu[[2]], a, b))
  return(tauUpdated)
}
#-----------------------
# for trace plot
#-----------------------
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

model2_mcmc <- function(y0, y1, c, p, beta, sigma, mu, tau, yearLabel, censor,  
                        muAlpha, sigmaAlpha, alpha, ## initial values 
                        iterMax, output, prop.sd, seed = 2019) {
  set.seed(seed)
  uniqueYearLength <- length(unique(yearLabel)) 
  
  p.keep <- list()
  beta.keep <- list()
  sigma.keep <- list()
  mu.keep <- list()
  tau.keep <- list()
  muAlpha.keep <- list()
  sigmaAlpha.keep <- list() 
  
  iter <- 0
  while (iter < iterMax) {
    ## update observations 
    yupdate <- updateCensObs(y0, y1, beta, sigma, c, censor, yearLabel)
    ## update allocation 
    cupdate <- updateAlloc(yupdate, p, beta, sigma, yearLabel) 
    ## update mean of alpha 
    muAlphaUpdate <- updateMuAlpha(alpha, sigmaAlpha)
    ## update sd of alpha 
    sigmaAlphaUpdate <- updateSigmaAlpha(alpha, muAlphaUpdate)
    ## update alpha 
    alphaUpdate <- updateAlphaVec(cupdate,  uniqueYearLength, yearLabel, alpha, 
                                  muAlphaUpdate, sigmaAlphaUpdate, prop.sd) 
    ## update mixture weight 
    pupdate <- updateWeight(alphaUpdate) 
    ## update hierarchical sd tau 
    tauUpdate <- updateTau(beta, mu, a = 0.0001, b = 0.0001, yearLabel)
    ## update hierarchical mean mu
    muUpdate <- updateHierarchicalMu(mu, beta, tauUpdate,  a = 0, b = 10^6, yearLabel) 
    ## update mixture normal mean beta 
    betaUpdate <- updateBeta(yupdate, cupdate, sigma, muUpdate, tauUpdate, 
                             uniqueYearLength, yearLabel)
    ## update variance first 
    sigmaUpdate <- updateSigma(yupdate, cupdate, yearLabel, uniqueYearLength, betaUpdate)
    
    ## for next iteration
    beta <- betaUpdate
    sigma <- sigmaUpdate
    mu <- muUpdate
    tau <- tauUpdate
    alpha <- alphaUpdate
    muAlpha <- muAlphaUpdate
    sigmaAlpha <- sigmaAlphaUpdate
    p <- pupdate
    c <- cupdate
    
    iter <- iter + 1
    
    ## keep 
    muAlpha.keep[[iter]] <- muAlphaUpdate
    sigmaAlpha.keep[[iter]] <- sigmaAlphaUpdate
    p.keep[[iter]] <- pupdate
    beta.keep[[iter]] <- betaUpdate
    sigma.keep[[iter]] <- sigmaUpdate
    mu.keep[[iter]] <- muUpdate
    tau.keep[[iter]] <- tauUpdate
    
    # if (iter %% 1 == 0) {
    #   cat("iter is ", iter, "\n")
    #   print(pupdate)
    #   print(betaUpdate)
    #   print(sigmaUpdate)
    #   print(muUpdate)
    #   print(tauUpdate)
    #   cat("muAlpha is ", muAlphaUpdate, ", and sigmaAlpha is ", sigmaAlphaUpdate, "\n")
    # }
  }
  
  ## only save results  after burnin
  # muAlpha.res.removeBurnin <- muAlpha.keep[seq(burnin + 1, iterMax)]
  # sigmaAlpha.res.removeBurnin <- sigmaAlpha.keep[seq(burnin + 1, iterMax)]
  # p.res.removeBurnin <- p.keep[seq(burnin + 1, iterMax)]
  # beta.res.removeBurnin <- beta.keep[seq(burnin + 1, iterMax)]
  # sigma.res.removeBurnin <- sigma.keep[seq(burnin + 1, iterMax)]
  # mu.res.removeBurnin <- mu.keep[seq(burnin + 1, iterMax)]
  # tau.res.removeBurnin <- tau.keep[seq(burnin + 1, iterMax)]
  save(muAlpha.keep,     file = paste0(output, "muAlpha.keep.RData"))
  save(sigmaAlpha.keep,  file = paste0(output, "sigmaAlpha.keep.RData"))
  save(p.keep,           file = paste0(output, "p.keep.RData"))
  save(beta.keep,        file = paste0(output, "beta.keep.RData"))
  save(sigma.keep,       file = paste0(output, "sigma.keep.RData"))
  save(mu.keep,          file = paste0(output, "mu.keep.RData"))
  save(tau.keep,         file = paste0(output, "tau.keep.RData"))
  
}



