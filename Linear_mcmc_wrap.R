library("dplyr")
# setwd("~/Documents/ISU/RA_AMR/SwineSalData")
anti_Name <- c("AMI", "AMP", "ATM", "AUG", "AXO", "AZM", "CAZ", "CCV", "CEP", "CEQ", "CHL", "CIP", "COT", "CTC",
               "CTX", "FEP", "FIS", "FOX", "GEN", "IMI", "KAN", "NAL", "PTZ", "SMX", "STR", "TET", "TIO")

for (i in 1:2) { 
  for (j in c(11, 27)) { #c(1:27)[-c(1, 6, 8, 10, 14, 16, 20, 26)] 
    print(paste0(c("SalT", "Sal4")[i], "/ ", anti_Name[j]))
    serotype <- c("SalT", "Sal4")[i]
    antibiotic <-  anti_Name[j]
    ##----------------------------------------------
    ## Loading data set 
    ##----------------------------------------------
    datFileName <- paste0("DataBySeroAnti/", serotype, "_", antibiotic, ".csv")
    dat <- read.csv(datFileName, stringsAsFactors = F, header = T)
    dat <- dat %>% filter(Year >= 2002)
    y0 <- dat$l_vec
    y1 <- dat$u_vec
    censor <- dat$censored
    minYear <- min(dat$Year)
    yearLabel <- as.numeric(as.factor(as.numeric(dat$Year)))
    uniqueYearLength <- length(unique(yearLabel))
    ##----------------------------------------------
    ## initial values : estimated from data set 
    ##----------------------------------------------
    dat_group <- dat %>% mutate(Rslt_log2 = log(Rslt, base = 2),
                                cGroup = ifelse(Concl=="R", 1, 0)) %>% 
      group_by(Year, cGroup) %>% 
      summarise(meanMIC = mean(Rslt_log2), 
                num = n(),
                meany0MIC = mean(l_vec)) %>%  
      data.frame() 
    
    dat_totObs <- dat %>% group_by(Year) %>% summarise(totObs = n()) %>% data.frame() 
    dat_full <- left_join(dat_group, dat_totObs) %>% mutate(prop = num / totObs)
    ## proportion of resistant 
    p_tmp <-  left_join(data.frame(Year = sort(unique(dat$Year))), 
                        dat_full %>% filter(cGroup == 1) %>% dplyr::select(Year, prop))
    if (any(is.na(p_tmp$prop))) {
      p_tmp$prop[is.na(p_tmp$prop)] <- 0
    }
    p <- p_tmp$prop
    if (any(p < 1e-15)) {
      p[p < 1e-15] <- 0.000001
    }
    if (any(p > 1-(1e-15))) {
      p[p > 1-(1e-15)] <- 1-0.000001
    }
    ## beta 
    beta_tmp <- data.frame(Year = sort(unique(dat$Year)))
    beta_1 <- left_join(beta_tmp, dat_full %>% filter(cGroup == 0) %>% 
                          dplyr:: select(Year, meany0MIC) %>% rename(beta1 = meany0MIC)) 
    beta_1_2 <- left_join(beta_1, dat_full %>% filter(cGroup == 1) %>% 
                            dplyr:: select(Year, meany0MIC) %>% rename(beta2 = meany0MIC)) 
    if (any(is.na(beta_1_2$beta1))) {
      beta_1_2$beta1[which(is.na(beta_1_2$beta1))] <- mean(beta_1_2$beta1, na.rm = T)
    }
    if (any(is.na(beta_1_2$beta2))) {
      beta_1_2$beta2[which(is.na(beta_1_2$beta2))] <- mean(beta_1_2$beta2, na.rm = T)
    }
    beta <- beta_1_2[, c(2, 3)]
    ## initial value for allocation 
    c <- ifelse(dat$Concl=="R", 1, 0)
    ## sigma
    sigma <- dat %>% mutate(cGroup = ifelse(Concl=="R", 1, 0)) %>% 
      group_by(cGroup) %>% summarise(se = sd(l_vec)) %>% data.frame () %>% .[, "se"] 
    SMALLVAL <- 0.0001
    if (any(sigma == 0)) {
      sigma <- sigma + SMALLVAL 
    }
    if (any(is.na(sigma))) {
      sigma[is.na(sigma)] <- 0.1
    }
    ## mu
    lmod <- lm(beta$beta1~c(1:nrow(beta)))
    mu[[3]] <- lmod$coefficients
    mu[[2]] <- beta$beta2 %>% mean()
    mu[[1]] <- lmod$coefficients[1] + lmod$coefficients[2]*c(1:nrow(beta))
    ## tau
    tau <- c(summary(lmod)$sigma, apply(beta, 2, sd)[2])
    ## calculate logit(p)
    logit <- function(x) {
      n <- length(x) 
      if (any(x < 1e-15)) {
        x[x < 1e-15] <- 1e-15
      }
      sapply(seq(1, n), function(s){log(x[s] / (1 - x[s]))})
    }
    alpha <- logit(p)
    muAlpha <- mean(alpha)
    sigmaAlpha <- sd(alpha)
  
    ##--------------------------------------
    source("Linear_mcmc.R")
    if (!dir.exists(paste0("LinearModelOutput/From2002/", serotype, "_", antibiotic, "_Res2"))) {
      dir.create(paste0("LinearModelOutput/From2002/", serotype, "_", antibiotic, "_Res2"))
    }
    
    output <- paste0("LinearModelOutput/From2002/", serotype, "_", antibiotic, "_Res2/")
    
    iterMax <- 10000
    
    model2_mcmc(y0, y1, c, p, beta, sigma, mu, tau, yearLabel, censor, 
                muAlpha, sigmaAlpha, alpha, ## initial values
                iterMax, output, prop.sd = 5, seed = 2019)
  }
}
system("say The process if finished")
