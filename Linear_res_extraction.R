##---------------------------------------------------
## This Code is used to extract the output for hierarchical logit 
## set working directory as ./Model1_Paper1
##---------------------------------------------------
##rm(list = ls())

##---------------------------------------------------
library("dplyr")
library("ggplot2")
anti_Name <- c("AMI", "AMP", "ATM", "AUG", "AXO", "AZM", "CAZ", "CCV", "CEP", "CEQ", "CHL", "CIP", "COT", "CTC",
               "CTX", "FEP", "FIS", "FOX", "GEN", "IMI", "KAN", "NAL", "PTZ", "SMX", "STR", "TET", "TIO")
All_Sal <- c("SalT", "Sal4")

##---------------------------------------------------

##---------------------------------------------------
## Extract Hierarchical logit model 
##---------------------------------------------------
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
                     ncol = cols, nrow = ceiling(numPlots/cols), byrow = T)
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

extractHierarchicalBeta <- function(yearVec, path) {
  load(paste0(path, "beta.keep.RData")) 
  load(paste0(path, "mu.keep.RData")) 
  load(paste0(path, "p.keep.RData")) 
  
  beta.res.removeBurnin <- tail(beta.keep, 6000)
  mu.res.removeBurnin <- tail(mu.keep, 6000)
  p.res.removeBurnin <- tail(p.keep, 6000)
  ## extract 
  I <- length(yearVec)
  beta <- data.frame(beta0X25 = apply(matrix(unlist(lapply(beta.res.removeBurnin, 
                                                           function(s){s[,1]})), ncol = I, byrow = T), 2, 
                                      function(s){quantile(s, probs = 0.025)}),
                     beta0X975 = apply(matrix(unlist(lapply(beta.res.removeBurnin, 
                                                            function(s){s[,1]})), ncol = I, byrow = T), 2, 
                                       function(s){quantile(s, probs = 0.975)}),
                     beta0mean = colMeans(matrix(unlist(lapply(beta.res.removeBurnin, 
                                                               function(s){s[,1]})), ncol = I, byrow = T)),
                     beta1X25 = apply(matrix(unlist(lapply(beta.res.removeBurnin, 
                                                           function(s){s[,2]})), ncol = I, byrow = T), 2, 
                                      function(s){quantile(s, probs = 0.025)}),
                     beta1X975 = apply(matrix(unlist(lapply(beta.res.removeBurnin, 
                                                            function(s){s[,2]})), ncol = I, byrow = T), 2, 
                                       function(s){quantile(s, probs = 0.975)}),
                     beta1mean = colMeans(matrix(unlist(lapply(beta.res.removeBurnin, 
                                                               function(s){s[,2]})), ncol = I, byrow = T)), 
                     
                     mu0X25 = apply(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T), 2, 
                                    function(s) quantile(s, 0.025))[1:I],

                     mu0X975 = apply(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T), 2, 
                                    function(s) quantile(s, 0.975))[1:I],
                     
                     mu0mean = colMeans(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T))[1:I],
                     
                     ga0X25 = apply(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T), 2, 
                                    function(s) quantile(s, 0.025))[I+2],
                     
                     ga0X975 = apply(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T), 2, 
                                     function(s) quantile(s, 0.975))[I+2],
                     
                     ga0mean = colMeans(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T))[I+2],
                     
                     ga1X25 = apply(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T), 2, 
                                    function(s) quantile(s, 0.025))[I+3],
                     
                     ga1X975 = apply(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T), 2, 
                                     function(s) quantile(s, 0.975))[I+3],
                     
                     ga1mean = colMeans(matrix(unlist(mu.res.removeBurnin), nrow = length(mu.res.removeBurnin), byrow = T))[I+3],
                     
                     pX25 = apply(matrix(unlist(p.res.removeBurnin), nrow = length(p.res.removeBurnin), byrow = T), 2, 
                                  function(s) quantile(s, 0.025)),
                     pX975 = apply(matrix(unlist(p.res.removeBurnin), nrow = length(p.res.removeBurnin), byrow = T), 2, 
                                   function(s) quantile(s, 0.975)),
                     pmean = apply(matrix(unlist(p.res.removeBurnin), nrow = length(p.res.removeBurnin), byrow = T), 2, mean)
                     
  )
  
  beta <- beta[, c(3, 1, 2, 6, 4, 5, 9, 7, 8, 12, 10, 11, 15, 13, 14, 18, 16, 17)]        
  beta <- cbind(yearVec, beta)
  rm(list = ls(pattern = "*RData"))
  return(beta)
}

calculateDiff <- function(yearVec, path) {
  load(paste0(path, "beta.keep.RData")) 
  beta.res.removeBurnin <- tail(beta.keep, 6000)
  diff.res.removeBurnin <- lapply(beta.res.removeBurnin, function(t) t[-1,]-t[-nrow(beta.res.removeBurnin[[1]]),])
  I <- length(yearVec)
  if (I>2) {
    diff <- data.frame(diff1X25 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                             function(s){s[,1]})), ncol = I-1, byrow = T), 2, 
                                        function(s){quantile(s, probs = 0.025)}),
                       diff1X975 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                              function(s){s[,1]})), ncol = I-1, byrow = T), 2, 
                                         function(s){quantile(s, probs = 0.975)}),
                       diff1mean = colMeans(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                                 function(s){s[,1]})), ncol = I-1, byrow = T)),
                       diff2X25 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                             function(s){s[,2]})), ncol = I-1, byrow = T), 2, 
                                        function(s){quantile(s, probs = 0.025)}),
                       diff2X975 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                              function(s){s[,2]})), ncol = I-1, byrow = T), 2, 
                                         function(s){quantile(s, probs = 0.975)}),
                       diff2mean = colMeans(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                                 function(s){s[,2]})), ncol = I-1, byrow = T)))
    
  } else {
    diff <- data.frame(diff1X25 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                             function(s){s[1]})), ncol = I-1, byrow = T), 2, 
                                        function(s){quantile(s, probs = 0.025)}),
                       diff1X975 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                              function(s){s[1]})), ncol = I-1, byrow = T), 2, 
                                         function(s){quantile(s, probs = 0.975)}),
                       diff1mean = colMeans(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                                 function(s){s[1]})), ncol = I-1, byrow = T)),
                       diff2X25 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                             function(s){s[2]})), ncol = I-1, byrow = T), 2, 
                                        function(s){quantile(s, probs = 0.025)}),
                       diff2X975 = apply(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                              function(s){s[2]})), ncol = I-1, byrow = T), 2, 
                                         function(s){quantile(s, probs = 0.975)}),
                       diff2mean = colMeans(matrix(unlist(lapply(diff.res.removeBurnin, 
                                                                 function(s){s[2]})), ncol = I-1, byrow = T)))
    
  }
  
  diff <- diff[, c(3,1,2,6,4,5)]
  diff <- cbind(Year=yearVec[-1], diff)
  rm(list = ls(pattern = "*RData"))
  return(diff)
}



for (i in 1:2) {
  for (j in c(11, 27)) { 
    path <- paste0("./LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/")
    yearSource <- read.csv(paste0("./DataBySeroAnti/", All_Sal[i], "_", anti_Name[j], ".csv")) %>% 
      filter(Year >= 2002)
    yearVec <- yearSource %>% 
      dplyr::select("Year") %>% 
      unique() %>% .[,"Year"]
    ParamEst <- extractHierarchicalBeta(yearVec, path)
    # sig.ind <- character()
    # if (nrow(ParamEst) > 1) {
    #   DiffEst <- calculateDiff(yearVec, path)
    #   sig.ind <- ifelse(sum(DiffEst$diff1X25 * DiffEst$diff1X975 > 0)>0, "*", "")
    #   sig.year <- DiffEst$Year[which(DiffEst$diff1X25 * DiffEst$diff1X975 > 0)]
    #   write.csv(DiffEst, paste0(path, "DiffEst", sig.ind, ".csv"))
    # }
    ## draw the means of non-resistance group
    # pdf(paste0(path, "BetaFirstEst.pdf"))
    # plot(ParamEst$yearVec, 2^(ParamEst$beta0mean), type = "l", col="red", 
    #      main = paste0(All_Sal[i], "_", anti_Name[j],"_","Beta First Comp Est"),
    #      xlab = "Year", ylab = "Estimated mean of 1st component")
    # points(ParamEst$yearVec, 2^(ParamEst$beta0mean), col="red", pch=20)
    # if (sig.ind == "*") {
    #   points(ParamEst$yearVec[which(DiffEst$diff1X25 * DiffEst$diff1X975 > 0)+1],
    #          2^(ParamEst$beta0mean[which(DiffEst$diff1X25 * DiffEst$diff1X975 > 0)+1]),
    #          col = "black", pch = "O")
    # }
    # lines(x=ParamEst$yearVec, y=2^(ParamEst$mu0mean), col="blue")
    # dev.off()
    ## write out the estimates
    slope.sig <- character()
    slope.sig <- ifelse(ParamEst$ga1X25[1] * ParamEst$ga1X975[1] > 0, "*", "")
    write.csv(ParamEst, paste0(path, "ParamEst.csv"))
    write.csv(ParamEst, paste0(path, "ParamEst", slope.sig, ".csv"))
    ## trace plot of beta0's (whole chain)
    # load(paste0("./Swine", All_Sal[i], "/", All_Sal[i], "_", anti_Name[j], "_Res/beta.keep.RData"))
    # trplt.beta0 <- list()
    # for (k in 1:length(yearVec)) {
    #   y <- lapply(beta.keep, function(t) t[k, 1]) %>% unlist()
    #   tracedf <- data.frame("x"=c(1:length(y)), "y"=y)
    #   trplt.beta0[[k]] <- ggplot(tracedf, aes(x=x, y=y)) +
    #     geom_line() + theme_bw() + labs(x = paste0("year", k), y ="beta_0")
    # }
    # if (length(yearVec) < 16) {
    #   for (k in (length(yearVec)+1):16) {
    #     trplt.beta0[[k]] <- ggplot(tracedf, aes(x=x, y=y)) +
    #       geom_blank()
    #   }
    # }
    # pdf(paste0("./Swine", All_Sal[i],"/", All_Sal[i], "_", anti_Name[j], "_Res/", All_Sal[i],"_", anti_Name[j],"_trplt_beta0.pdf"))
    # multiplot(trplt.beta0[[1]], trplt.beta0[[2]], trplt.beta0[[3]], trplt.beta0[[4]],
    #           trplt.beta0[[5]], trplt.beta0[[6]], trplt.beta0[[7]], trplt.beta0[[8]],
    #           trplt.beta0[[9]], trplt.beta0[[10]],trplt.beta0[[11]],trplt.beta0[[12]],
    #           trplt.beta0[[13]], trplt.beta0[[14]],trplt.beta0[[15]],trplt.beta0[[16]], cols = 4)
    # dev.off()
    # ## trace plot of gamma's
    # load(paste0("./Swine", All_Sal[i], "/", All_Sal[i], "_", anti_Name[j], "_Res/mu.keep.RData"))
    # trplt.gamma <- list() # list of 2
    # for (k in 1:2) {
    #   y <- lapply(mu.keep, function(t) t[[3]][k]) %>% unlist
    #   tracedf <- data.frame("x"=c(1:length(y)), "y"=y)
    #   trplt.gamma[[k]] <- ggplot(tracedf, aes(x=x, y=y)) +
    #     geom_line() + theme_bw() + labs(x = paste0("gamma", k-1), y ="")
    # }
    # pdf(paste0("./Swine", All_Sal[i],"/", All_Sal[i], "_", anti_Name[j], "_Res/", All_Sal[i],"_", anti_Name[j],"_trplt_gamma.pdf"))
    # multiplot(trplt.gamma[[1]], trplt.gamma[[2]], cols = 1)
    # dev.off()
  }
}

#------------------
# comb plot log2 scale
#------------------
for (i in 1:2) { 
  for (j in c(11, 27)) { #c(1:27)[-c(1, 6, 8, 10, 14, 16, 20, 26)]
    serotypeNames <- All_Sal[i]
    antibioticsDrugs <- anti_Name[j]
    path <- paste0("./LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/")
    datFileName <- paste0("DataBySeroAnti/", serotypeNames, "_", antibioticsDrugs, ".csv") 
    yearSource <- read.csv(datFileName) %>% filter(Year >= 2002)
    yearVec <- yearSource %>% 
      dplyr::select("Year") %>% 
      unique() %>% .[,"Year"]
    plotDatSource <- data.frame(Year = as.character(yearVec))
    dat <- read.csv(datFileName, stringsAsFactors = F, header = T) %>% filter(Year >= 2002)
    dat_group <- dat %>% mutate(Rslt_log2 = log(Rslt, base = 2),
                                cGroup = ifelse(Concl=="R", 1, 0), 
                                Year = as.factor(Year)) %>% 
      group_by(Year, cGroup) %>% 
      summarise(meanMIC = mean(Rslt_log2), 
                num = n(),
                meany0MIC = mean(l_vec)) %>%  data.frame() 
    dat_totObs <- dat %>% mutate(Year = as.factor(Year)) %>% 
      group_by(Year) %>% summarise(totObs = n()) %>% data.frame() 
    dat_full <- left_join(dat_group, dat_totObs) %>% mutate(prop = num / totObs)
    plotDatSource <- left_join(plotDatSource, dat_full[dat_full$cGroup==1, c("Year", "prop")])
    plotDatSource <- left_join(plotDatSource, dat_group[dat_group$cGroup==0, c("Year", "meanMIC")])
    beta0est <- read.csv(paste0("LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/", "ParamEst.csv"))[, 2:5]
    beta0est$yearVec <- as.factor(beta0est$yearVec)
    beta0est.org <- beta0est
    #beta0est.org[, 2:4] <- 2^(beta0est[, 2:4])
    line <- read.csv(paste0("LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/", "ParamEst.csv"))[, "mu0mean"]
    ga0mean <- read.csv(paste0("LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/", "ParamEst.csv"))[, "ga0mean"][1]
    ga1mean <- read.csv(paste0("LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/", "ParamEst.csv"))[, "ga1mean"][1]
    line.org <- line
    #line.org <- 2^line
    plotDatSource <- left_join(plotDatSource, beta0est.org, by = c("Year"="yearVec"))
    
    load(paste0("./LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/mu.keep.RData"))
    mu.res.removeBurnin <- tail(mu.keep, 6000)
    
    covOfGamma <- cov(matrix(unlist(mu.res.removeBurnin), nrow = 6000, byrow = T)[, 16], matrix(unlist(mu.res.removeBurnin), nrow = 6000, byrow = T)[, 17])
    varGamma0 <- matrix(unlist(mu.res.removeBurnin), nrow = 6000, byrow = T)[, 16] %>% var()
    varGamma1 <- matrix(unlist(mu.res.removeBurnin), nrow = 6000, byrow = T)[, 17] %>% var()
    
    t <- 1:14
    varOfLinear <- varGamma0 + varGamma1 * t^2 + 2 * covOfGamma * t
    centerOfLinear <- ga0mean + ga1mean * t
    upperOfLinear <- centerOfLinear + 1.96 * sqrt(varOfLinear)
    lowerOfLinear <- centerOfLinear - 1.96 * sqrt(varOfLinear)
    
    # ymin <- plotDatSource[, c("beta0mean", "beta0X25", "beta0X975", "meanMIC")] %>% min(na.rm=T)
    # ymax <- plotDatSource[, c("beta0mean", "beta0X25", "beta0X975", "meanMIC")] %>% max(na.rm=T)
    ymin <- min(lowerOfLinear)
    ymax <- max(upperOfLinear)
    pmax <- plotDatSource[, "prop"] %>% max(na.rm=T)
    
    # sig.ind <- character()
    # if (nrow(beta0est) > 1) {
    #   DiffEst <- calculateDiff(yearVec, path)
    #   sig.ind <- ifelse(sum(DiffEst$diff1X25 * DiffEst$diff1X975 > 0)>0, "*", "")
    #   sig.year <- DiffEst$Year[which(DiffEst$diff1X25 * DiffEst$diff1X975 > 0)]
    #   sig.year.nrow <- sig.year-min(yearVec)+1
    #   sig.df <- data.frame("sig.year" = sig.year.nrow, "sig.val" = (plotDatSource$beta0mean[sig.year.nrow]-ymin)*(pmax/(ymax-ymin)))
    # }
    plotDatSource <- plotDatSource %>% 
      mutate(prop = ifelse(is.na(prop), 0, prop))
    plotres <- ggplot(plotDatSource, aes(x = Year, y = prop)) +
      geom_bar(stat = "identity", fill = "grey85") +
      labs(y = "Observed proportion of resistant isolates", x = "Year") +
      geom_line(aes(y = (meanMIC-ymin)*(pmax/(ymax-ymin)), group = 2), col = "grey50")+
      geom_point(aes(y = (meanMIC-ymin)*(pmax/(ymax-ymin)), group = 2), col = "grey50", size = 1.5, shape = 4)+
      geom_line(aes(y = (beta0mean-ymin)*(pmax/(ymax-ymin)), group = 2), col = "red")+
      geom_point(aes(y = (beta0mean-ymin)*(pmax/(ymax-ymin)), group = 2), col = "red", size = 1.5)+
      # geom_line(aes(y = (beta0X25-ymin)*(pmax/(ymax-ymin)), group = 2), col = "red", linetype="dotted")+
      # geom_line(aes(y = (beta0X975-ymin)*(pmax/(ymax-ymin)), group = 2), col = "red", linetype="dotted")+
      geom_line(aes(y = (line.org-ymin)*(pmax/(ymax-ymin)), group = 2), col = "blue") +
      scale_y_continuous(limits=c(0,0.1),
                         sec.axis = sec_axis(trans = ~./(pmax/(ymax-ymin))+ymin, name = "Estimated mean log2(MIC) in non-resistant population")) +
      scale_color_manual(labels = c("Naive mean", "Est beta0"), values = c("black", "red")) +
      theme_light() +
      theme(axis.text.x = element_text(angle = 270, hjust = 1), 
            axis.title.y.right=element_text(color="red"),
            axis.text.y.right=element_text(color="red")) +
      labs(title = paste0(c("S. enterica serovar Typhimurium", "S. enterica serovar I,4,[5],12:i:-")[i],
                          ", ", antibioticsDrugs, ", log2(MIC) scale"))
    
    # plotres <- plotres + geom_point(aes(x = sig.year, y = sig.val), data = sig.df, shape = 1, size = 3, col = "black")
    
    plotres <- plotres +
      geom_ribbon(aes(x = t, 
                      ymin=(lowerOfLinear-ymin)*(pmax/(ymax-ymin)),
                      ymax=(upperOfLinear-ymin)*(pmax/(ymax-ymin))), 
                  fill="blue", alpha="0.1")  
    
    pdf(paste0("./LinearModelOutput/From2002/", All_Sal[i], "_", anti_Name[j], "_Res/", All_Sal[i], "_", anti_Name[j], "_LogComb_newRangeCI.pdf"))
    plotres
    dev.off()
  }
}

