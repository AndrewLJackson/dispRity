## Selects the items from the dispRity object needed for the model testings
select.model.list <- function(data, observed, cent.tend = median, rarefaction) {

    if(observed) {
        ## If observed is required
        central_tendency <- unlist(extract.dispRity(data, observed = TRUE))

        ## If disparity is a single value
        if(unique(unlist(lapply(data$disparity, lapply, lapply, length))) != 1) {
            ## Calculate the variance from the disparity data
            variance <- unlist(lapply(extract.dispRity(data, observed = FALSE), lapply, var))
        } else {
            ## Extract directly the variance from the data
            variance <- sapply(data[[3]], function(x) var(data[[1]][x[[1]]]))
        }

    } else {

        ## Getting the disparity
        if(!missing(rarefaction)) {
            disparity_tmp <- extract.dispRity(data, observed = FALSE, rarefaction = rarefaction)    
        } else {
            disparity_tmp <- extract.dispRity(data, observed = FALSE)
        }
        ## Calculating the central tendency
        central_tendency <- unlist(lapply(disparity_tmp, lapply, cent.tend))
        ## Calculating the variance
        variance <- unlist(lapply(disparity_tmp, lapply, var))
    }

    ## Getting the length of the samples
    summary_table <- summary(data)
    if(!missing(rarefaction)) {
        sample_length <- rep(rarefaction, length(central_tendency))
    } else {
        sample_length <- summary_table$n[which(!is.na(summary_table[,3]))]
    }

    ## Samples
    if(data$call$subsamples[1] == "continuous") {
        subsamples <- sort(as.numeric(names(data$subsamples)))
    } else {
        subsamples <- seq(1:length(data$subsamples))
    }

    ## Returns the data
    return(list("central_tendency" = central_tendency,
                "variance" = variance,
                "sample_size" = sample_length,
                "subsamples" = rev(subsamples)))
}

# Code up-dated from Gene Hunt's 'paleoTS' package. The code here

# These functions fit a model in trait evolution tracks a covariate (z) over time. There are two parameterizations: "Joint" and "AD"; all functions without "joint" in their names use the "AD" parameterization. The "joint" parameterization assumes the trait values are a linear function of the covariate, whereas the "AD" parameterization assumes that changes in the traits are a linear function of changes in the covariate.

# library(paleoTS)
# library(mnormt)


# pooled variance of all data in time series (modified from PaleoTS::pool.var)
pooled.variance <- function(data_list, rescale.variance = FALSE)  {

    ## Calculate the pooled variance
    pooled_variance <- sum(data_list$variance * (data_list$sample_size - 1))/sum(data_list$sample_size - 1)

    if (rescale.variance) {
        ## Rescale the variance
        data_out <- data_list
        data_out$variance <- rep(pooled_variance, length(data_list$central_tendency))
        return(data_out)
    } else {
        return(pooled_variance)
    }
}



    ############### OU MODEL    ###############
 
    MLrw <- function (data) {
    
        N <- length(data$central_tendency) - 1
        t <- (data$subsamples[N + 1] - data$subsamples[1])/N
        epsilon <- 2 * pooledVariance(data) / round(median(data$sample_size))
        xDifference <- diff(data$central_tendency)
        xMeanDifference <- mean(xDifference)
        muStep <- xMeanDifference / t
        sigmaSquaredStep <- (1/t) * ((1/N) * sum(xDifference ^ 2) - xMeanDifference^2 - epsilon)
        outputParams <- c(muStep, sigmaSquaredStep)
        names(outputParams) <- c("muStep", "sigmaSquaredStep")
        return(outputParams)
    
    }


    OULikelihoodJointOptim <- function (p, data, timeSplit=NULL, nOptima=1)  {
    
        timeSplit <- c(0, timeSplit, length(data$subsamples))
        nTheta <- nOptima
        ancState <- p[1]
        sigmaSquared <- p[2]
        alpha <- p[3]
        optima <- p[-c(1:3)]
        NN <- length(data$central_tendency)
        FF <- function(a, b) abs(a - b)
        VCV <- outer(data$subsamples, data$subsamples, FUN = FF)
        VCV <- exp(-alpha * VCV)
        ou.V <- function (sigmaSquared, alpha, time) (sigmaSquared / (2 * alpha)) * (1 - exp(-2 * alpha * time))  
        VCVd <- ou.V(sigmaSquared, alpha, data$subsamples)
        VCV2 <- outer(VCVd, VCVd, pmin)
        VCV <- VCV * VCV2
        diag(VCV) <- VCVd + data$variance / data$sample_size
        ou.M <- function (ancState, optima, alpha, time) optima * (1 - exp(-alpha * time)) + ancState * exp(-alpha * time)
        ouMean <- sapply(1:nTheta, function(x) ou.M(ancState, optima[x], alpha, data$subsamples))
        meanOU <- c()
        for(x in 1:nOptima) meanOU <- c(meanOU, ouMean[(timeSplit[x] + 1) : timeSplit[x + 1] , x])  
        S <- dmnorm(t(data$central_tendency), mean = meanOU, varcov = VCV, log = TRUE)
        return(S)
    
    }


    disparityOU <- function (data, nOptima=1, timeSplit=NULL, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {
    
        if (poolV)  data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputOptimisation <- MLrw(data)
        halfLife <- (data$subsamples[length(data$subsamples)] - data$subsamples[1])/4
        inputParams <- c(data$central_tendency[1], inputOptimisation[2]/10, log(2)/halfLife, rep(data$central_tendency[length(data$central_tendency)], nOptima))
        names(inputParams) <- c("ancState", "sigmaSquared", "alpha", paste0("theta", 1:nOptima))
        if (is.null(cl$ndeps)) cl$ndeps <- abs(inputParams/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") optimisedPrm <- optim(inputParams, fn =  OULikelihoodJointOptim, control = cl, method = meth, lower = c(NA, 1e-10, 1e-08, rep(NA, nOptima)), hessian = hess, data=data, timeSplit=timeSplit, nOptima=nOptima)
        else optimisedPrm <- optim(inputParams, fn =  OULikelihoodJointOptim, control = cl, method = meth, lower = c(NA, 1e-10, NA, 1e-08), hessian = hess, data=data)
        if (hess)   optimisedPrm$se <- sqrt(diag(-1 * solve(optimisedPrm$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "OU", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))

    }


    ############### BM MODEL    ###############

## Select the parameters for a BM model
BM.parameters <- function (data_list) {

    ## Create the output model
    output_model <- array(dim = 2)
    output_model[1] <- data_list$central_tendency[1]

    ## Number of subsamples
    n_subsamples <- length(data_list$central_tendency) - 1

    ## Setting up the simulation time
    time <- (data_list$subsamples[n_subsamples + 1] - data_list$subsamples[1])/n_subsamples
    
    ## parameters
    epsilon <- 2 * pooled.variance(data_list) / round(median(data_list$sample_size))
    x_difference <- diff(data_list$central_tendency)
    x_mean_difference <- mean(x_difference)
    sigma_squared_step <- (1/time) * ((1/n_subsamples) * sum(x_difference ^ 2) - epsilon)
    output_param <- c("sigma_squared_step" = sigma_squared_step)

    ## Output
    output_model[2] <- output_param
    names(output_model) <- c("anc_state", "sigma_squared")
    return(output_model)
}

    BMLikelihoodJointOptim <- function (p, data)  {
    
        ancState <- p[1]
        sigmaSquared <- p[2]
        NN <- length(data$central_tendency)
        VCV <- sigmaSquared * outer(data$subsamples, data$subsamples, FUN = pmin)
        diag(VCV) <- diag(VCV) + data$variance / data$sample_size
        M <- rep(ancState, NN)
        return(dmnorm(data$central_tendency, mean = M, varcov = VCV, log = TRUE))

    }

    disparityBM <- function (data, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {

        if (poolV) data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputOptimisation <- array(dim = 2)
        inputOptimisation[1] <- data$central_tendency[1]
        inputOptimisation[2] <- MLbm(data) #min(c(MLbm(data), 1e-07))
        names(inputOptimisation) <- c("ancState", "sigmaSquared")
        if (is.null(cl$ndeps))  cl$ndeps <- abs(inputOptimisation/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") {
            optimisedPrm <- optim(inputOptimisation, fn = BMLikelihoodJointOptim, control = cl, method = meth, lower = c(NA, 1e-6), hessian = hess, data = data)
            } else {
            optimisedPrm <- optim(inputOptimisation, fn = BMLikelihoodJointOptim, control = cl, method = meth, hessian = hess, data = data)
            }
        if (hess) 
            optimisedPrm$se <- sqrt(diag(-1 * solve(w$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "BM", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))

}

    ############### STASIS MODEL    ###############

    MLStasis <- function (data) {
    
        N <- length(data$central_tendency)
        varPooled <- pooledVariance(data)
        theta <- mean(data$central_tendency[2:N])
        omega <- var(data$central_tendency[2:N]) - varPooled/median(data$variance)
        ww <- c(omega, theta)
        names(ww) <- c("omega", "theta")
        return(ww)

    }


    StasisLikelihoodJointOptim <- function (p, data, timeSplit=NULL, nOptima=1) {
    
        timeSplit <- c(0, timeSplit, length(data$subsamples))
        theta <- p[-1]
        omega <- p[1]
        N <- length(data$central_tendency)
        VCV <- diag(omega + data$variance/data$sample_size)
        meanStasis <- c()
        for(x in 1:nOptima) meanStasis <- c(meanStasis, rep(theta[x], length(timeSplit[x] : timeSplit[x + 1]) - 1))    
        return(dmnorm(data$central_tendency, mean = meanStasis, varcov = VCV, log = TRUE))

    }

    disparityStasis <- function (data, nOptima=1, timeSplit=NULL, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {

        if (poolV) data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputOptimisation <- MLStasis(data)
        inputOptimisation <- c(inputOptimisation, rep(inputOptimisation[2], nOptima - 1))
        if (inputOptimisation[2] <= 0 || is.na(inputOptimisation[2]))   inputOptimisation[2] <- 1e-07
        if (is.null(cl$ndeps))  cl$ndeps <- abs(inputOptimisation/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-09
        if (meth == "L-BFGS-B") optimisedPrm <- optim(inputOptimisation, fn = StasisLikelihoodJointOptim, control = cl, method = meth, lower = c(0, rep(NA, nOptima)), hessian = hess, data = data, timeSplit=timeSplit, nOptima=nOptima)
        else optimisedPrm <- optim(inputOptimisation, fn =StasisLikelihoodJointOptim, control = cl, method = meth, hessian = hess, data = data, timeSplit=timeSplit)
        if (hess) optimisedPrm$se <- sqrt(diag(-1 * solve(optimisedPrm$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "Stasis", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))
        
}



    ############### EARLY BURST MODEL   ###############

    MLEarlyBurst <- function (data) {
    
        N <- length(data$central_tendency) - 1
        t <- (data$subsamples[N + 1] - data$subsamples[1])/N
        epsilon <- 2 * pooledVariance(data) / round(median(data$sample_size))
        xDifference <- diff(data$central_tendency)
        xMeanDifference <- mean(xDifference)
        sigmaSquaredStep <- (1/t) * ((1/N) * sum(xDifference ^ 2) - epsilon)
        a <- log(1e-5) / max(data$subsamples) * (1/2)
        ww <- c(sigmaSquaredStep, a)
        names(ww) <- c("sigmaSquaredStep", "a")
        return(ww)

    }

    EBLikelihoodJointOptim <- function (p, data)  {
    
        ancState <- p[1]
        sigmaSquared <- p[2]
        rRate <- p[3]   
        NN <- length(data$central_tendency)
        timeOut <-  data$subsamples
        M <- rep(ancState, NN) 
        VCV <- outer(sigmaSquared * ((exp(rRate * timeOut) - 1) / rRate), sigmaSquared * ((exp(rRate * timeOut) - 1) / rRate), FUN=pmin)    
        diag(VCV) <- diag(VCV) + data$variance / data$sample_size
        return(dmnorm(data$central_tendency, mean = M, varcov = VCV, log = TRUE))

    }


    disparityEB <- function (data, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {

        if (poolV) data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputOptimisation <- array(dim = 3)
        inputOptimisation[1] <- data$central_tendency[1]
        inputOptimisation[2:3] <- MLEarlyBurst(data) 
        names(inputOptimisation) <- c("ancState", "sigmaSquared", "a")
        if (    is.null(cl$ndeps))  cl$ndeps <- abs(inputOptimisation/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") {
            optimisedPrm <- optim(inputOptimisation, fn = EBLikelihoodJointOptim, control = cl, method = meth, lower = c(NA, 0, -20), upper=c(Inf, Inf, -1e-6) , hessian = hess, data = data)
            } else {
            optimisedPrm <- optim(inputOptimisation, fn = EBLikelihoodJointOptim, control = cl, method = meth, hessian = hess, data = data)
            }
    if (hess)   optimisedPrm$se <- sqrt(diag(-1 * solve(w$hessian)))    else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "EB", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))
    
   }
    
    
    
    
    
    
    shiftMode <- function(dataset, shiftTime, modeOne, modeTwo, fullDetails=FALSE) {
        
        startSeq <- which(abs(dataset[[4]] - shiftTime) == min(abs(dataset[[4]] - shiftTime)))[1]
        startSeq2 <- startSeq  + 1

        datasetOne <- lapply(dataset, function(x) x[1:startSeq])
        datasetTwo <- lapply(dataset, function(x) x[-c(1:startSeq)])
        datasetTwo[[4]] <- datasetTwo[[4]] - datasetTwo[[4]][1]
    
        if(modeOne == "OU") modelOne <- disparityOU(datasetOne)
        if(modeOne == "BM") modelOne <- disparityBM(datasetOne)
        if(modeOne == "EB") modelOne <- disparityEB(datasetOne)
        if(modeOne == "Trend")  modelOne <- disparityTrend(datasetOne)
        
        if(modeTwo == "OU") modelTwo <- disparityOU(datasetTwo)
        if(modeTwo == "BM") modelTwo <- disparityBM(datasetTwo)
        if(modeTwo == "EB") modelTwo <- disparityEB(datasetTwo)
        if(modeTwo == "Trend")  modelTwo <- disparityTrend(datasetTwo)
        
        
        k <- sum(modelOne$K, modelTwo$K)
        n <- length(dataset[[1]])
        AICModel <- (-2 * sum(modelOne$logL, modelTwo$logL)) + (2 * k)
        AICcModel <- (-2 * sum(modelOne$logL, modelTwo$logL)) + (2 * k) * (n / ( n - k - 1))
        
        details <- list()
        details$AIC <- AICModel
        details$AICc <- AICcModel
        if(fullDetails) { details$modelOne <- modelOne ; details$modelTwo <- modelTwo }
        return(details)
                
        }
        
        

    
    ############### TREND MODEL ###############

    MLtrend <- function (data)  {

        N <- length(data$central_tendency) - 1
        t <- (data$subsamples[N + 1] - data$subsamples[1])/N
        epsilon <- 2 * pooledVariance(data) / round(median(data$sample_size))
        xDifference <- diff(data$central_tendency)
        xMeanDifference <- mean(xDifference)
        trendParam <- xMeanDifference / t
        sigmaSquaredStep <- (1/t) * ((1/N) * sum(xDifference ^ 2) - epsilon)
        outputParams <- c(sigmaSquaredStep, trendParam)
        names(outputParams) <- c("sigmaSquaredStep", "trendParam")
        return(outputParams)
        
    }

    trendLikelihoodJointOptim <- function (p, data)  {
    
        ancState <- p[1]
        sigmaSquared <- p[2]
        trendParam <- p[3]
        NN <- length(data$central_tendency)
        VCV <- sigmaSquared * outer(data$subsamples, data$subsamples, FUN = pmin)
        diag(VCV) <- diag(VCV) + data$variance / data$sample_size
        M <- rep(ancState, NN) + trendParam * data$subsamples
        return(dmnorm(data$central_tendency, mean = M, varcov = VCV, log = TRUE))

    }

    disparityTrend <- function (data, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {

        if (poolV) data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputOptimisation <- array(dim = 3)
        inputOptimisation[1] <- data$central_tendency[1]
        inputOptimisation[2:3] <- MLtrend(data) #min(c(MLbm(data), 1e-07))
        names(inputOptimisation) <- c("ancState", "sigmaSquared", "trendParam")
        if (is.null(cl$ndeps))  cl$ndeps <- abs(inputOptimisation/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") {
            optimisedPrm <- optim(inputOptimisation, fn = trendLikelihoodJointOptim, control = cl, method = meth, lower = c(NA, 0, -100), hessian = hess, data = data)
            } else {
            optimisedPrm <- optim(inputOptimisation, fn = trendLikelihoodJointOptim, control = cl, method = meth, hessian = hess, data = data)
            }
        if (hess) 
            optimisedPrm$se <- sqrt(diag(-1 * solve(w$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "Trend", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))
    

        
}




########## multi-mode models

########## OU two Trend

    OU2Trend <- function (p, data, timeMin)  {
    
        timeMin <- which.min(abs(data$subsamples - timeMin))
    
        ancState <- p[1]
        sigmaSquared <- p[2]
        alpha <- p[3]
        optima <- p[4]
        trendParam <- p[5]
                        
        NN <- length(data$central_tendency)
        FF <- function(a, b) abs(a - b)
        VCV <- outer(data$subsamples, data$subsamples, FUN = FF)
        VCV <- exp(-alpha * VCV)
        ou.V <- function (sigmaSquared, alpha, time) (sigmaSquared / (2 * alpha)) * (1 - exp(-2 * alpha * time))  
        VCVd <- ou.V(sigmaSquared, alpha, data$subsamples)
        VCV2 <- outer(VCVd, VCVd, pmin)
        VCV <- VCV * VCV2
        diag(VCV) <- VCVd + data$variance / data$sample_size
        ou.M <- function (ancState, optima, alpha, time) optima * (1 - exp(-alpha * time)) + ancState * exp(-alpha * time)
        ouMean <- ou.M(ancState, optima, alpha, data$subsamples)
        
        ########
        ########
                    
        VCV[c((timeMin): dim(VCV)[1]), ] <- 0
        VCV[ ,c((timeMin): dim(VCV)[1])] <- 0   
        
        timeOut <-  data$subsamples[-c(1 : (timeMin - 1))]
        timeOutDiff <- diff(timeOut[1:2])
        timeOut2 <- timeOut - (min(timeOut) - timeOutDiff)
        
        
        NN <- length(timeOut)
        M2 <- rep(ancState, NN) 
         
        VCV2 <- sigmaSquared * outer(timeOut, timeOut, FUN = pmin)
        diag(VCV2) <- diag(VCV2) + data$variance[-c(1 : (timeMin - 1))] / data$sample_size[-c(1 : (timeMin - 1))]
        M2 <- M2 + trendParam * timeOut2
        
        VCV[timeMin:dim(VCV)[1], timeMin:dim(VCV)[1]] <- VCV2   
        
        MM <- c(ouMean[1:(timeMin - 1)], M2)    
        S <- dmnorm(t(data$central_tendency), mean = MM, varcov = VCV, log = TRUE)
        return(S)
    
    }
    
    
        disparityOU2Trend <- function (data, timeMin=NULL, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {
    
        if (poolV)  data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputOptimisation <- MLrw(data)
        halfLife <- (data$subsamples[length(data$subsamples)] - data$subsamples[1])/4
        inputParams <- c(data$central_tendency[1], inputOptimisation[2]/10, log(2)/halfLife, data$central_tendency[length(data$central_tendency)])
        inputParams[5] <- MLtrend(data)[2]
        names(inputParams) <- c("ancState", "sigmaSquared", "alpha", "theta", "trend")
        if (is.null(cl$ndeps)) cl$ndeps <- abs(inputParams/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") optimisedPrm <- optim(inputParams, fn =  OU2Trend, control = cl, method = meth, lower = c(NA, 1e-10, 1e-08, NA, -100), hessian = hess, data=data, timeMin=timeMin)
        else optimisedPrm <- optim(inputParams, fn =  OU2Trend, control = cl, method = meth, lower = c(NA, 1e-10, NA, 1e-08, -100), hessian = hess, data=data)
        if (hess)   optimisedPrm$se <- sqrt(diag(-1 * solve(optimisedPrm$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "disparityOU2Trend", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))

    }

    
########## OU 2 EB


    OU2EB <- function (p, data, timeMin)  {
    
        timeMin <- which.min(abs(data$subsamples - timeMin))
    
        ancState <- p[1]
        sigmaSquared <- p[2]
        alpha <- p[3]
        optima <- p[4]
        rRate <- p[5]
        
        NN <- length(data$central_tendency)
        FF <- function(a, b) abs(a - b)
        VCV <- outer(data$subsamples, data$subsamples, FUN = FF)
        VCV <- exp(-alpha * VCV)
        ou.V <- function (sigmaSquared, alpha, time) (sigmaSquared / (2 * alpha)) * (1 - exp(-2 * alpha * time))  
        VCVd <- ou.V(sigmaSquared, alpha, data$subsamples)
        VCV2 <- outer(VCVd, VCVd, pmin)
        VCV <- VCV * VCV2
        diag(VCV) <- VCVd + data$variance / data$sample_size
        ou.M <- function (ancState, optima, alpha, time) optima * (1 - exp(-alpha * time)) + ancState * exp(-alpha * time)
        ouMean <- ou.M(ancState, optima, alpha, data$subsamples)
                
        VCV[c(timeMin: dim(VCV)[1]), ] <- 0
        VCV[ ,c(timeMin: dim(VCV)[1])] <- 0
                
        
        ########
        ########
        
        
        timeOut <-  data$subsamples[-c(1 : (timeMin - 1))]
        timeOutDiff <- diff(timeOut[1:2])
        timeOut2 <- timeOut - (min(timeOut) - timeOutDiff)
        
        NN <- length(timeOut)
        M2 <- rep(ancState, NN) 
        VCV2 <- outer(sigmaSquared * ((exp(rRate * timeOut2) - 1) / rRate), sigmaSquared * ((exp(rRate * timeOut2) - 1) / rRate), FUN=pmin) 
        diag(VCV2) <- diag(VCV2) + data$variance[-c(1 : (timeMin - 1))] / data$sample_size[-c(1 : (timeMin - 1))]
        
        VCV[timeMin:dim(VCV)[1], timeMin:dim(VCV)[1]] <- VCV2
        
        MM <- c(ouMean[1:(timeMin - 1)], M2)        
        
        S <- dmnorm(t(data$central_tendency), mean = MM, varcov = VCV, log = TRUE)
        return(S)
   
    }
    
    
        disparityOU2EB <- function (data, timeMin=NULL, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {
    
        if (poolV)  data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputOptimisation <- MLrw(data)
        halfLife <- (data$subsamples[length(data$subsamples)] - data$subsamples[1])/4
        inputParams <- c(data$central_tendency[1], inputOptimisation[2]/10, log(2)/halfLife, data$central_tendency[length(data$central_tendency)])
        inputParams[5] <- MLEarlyBurst(data)[2]
        names(inputParams) <- c("ancState", "sigmaSquared", "alpha", "theta", "a")
        if (is.null(cl$ndeps)) cl$ndeps <- abs(inputParams/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") optimisedPrm <- optim(inputParams, fn =  OU2EB, control = cl, method = meth, lower = c(NA, 1e-10, 1e-08, NA, -20), upper=c(NA,NA,NA,NA, -1e-6), hessian = hess, data=data, timeMin=timeMin)
        else optimisedPrm <- optim(inputParams, fn =  OU2EB, control = cl, method = meth, lower = c(NA, 1e-10, 1e-08, NA, -20), upper=c(NA,NA,NA,NA, -1e-6), hessian = hess, data=data)
        if (hess)   optimisedPrm$se <- sqrt(diag(-1 * solve(optimisedPrm$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "disparityOU2EB", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))

}
    
    
    ##################
    
    ########## BM 2 EB


    BM2EB <- function (p, data, timeMin)  {
    
    
        ancState <- p[1]
        sigmaSquared <- p[2]
        rRate <- p[3]
        
        ancState <- p[1]
        sigmaSquared <- p[2]
        NN <- length(data$central_tendency)
        VCV <- sigmaSquared * outer(data$subsamples, data$subsamples, FUN = pmin)
        diag(VCV) <- diag(VCV) + data$variance / data$sample_size
        MOne <- rep(ancState, NN)
                
        VCV[c(timeMin: dim(VCV)[1]), ] <- 0
        VCV[ ,c(timeMin: dim(VCV)[1])] <- 0
                
        
        ########
        ########
        
        
        timeOut <-  data$subsamples[-c(1 : (timeMin - 1))]
        timeOutDiff <- diff(timeOut[1:2])
        timeOut2 <- timeOut - (min(timeOut) - timeOutDiff)
        
        NN <- length(timeOut)
        M2 <- rep(ancState, NN) 
        VCV2 <- outer(sigmaSquared * ((exp(rRate * timeOut2) - 1) / rRate), sigmaSquared * ((exp(rRate * timeOut2) - 1) / rRate), FUN=pmin) 
        diag(VCV2) <- diag(VCV2) + data$variance[-c(1 : (timeMin - 1))] / data$sample_size[-c(1 : (timeMin - 1))]
        
        VCV[timeMin:dim(VCV)[1], timeMin:dim(VCV)[1]] <- VCV2
        
        MM <- c(MOne[1:(timeMin - 1)], M2)      
        
        S <- dmnorm(t(data$central_tendency), mean = MM, varcov = VCV, log = TRUE)
        return(S)
   
    }
    
    
        disparityBM2EB <- function (data, timeMin=NULL, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {
    
        if (poolV)  data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputParams <- array(dim = 3)
        inputParams[1] <- data$central_tendency[1]
        inputParams[2] <- MLbm(data) #min(c(MLbm(data), 1e-07))
        inputParams[3] <- MLEarlyBurst(data)[2]
        names(inputParams) <- c("ancState", "sigmaSquared", "a")
        timeMin <- which.min(abs(data$subsamples - timeMin))
        if (is.null(cl$ndeps)) cl$ndeps <- abs(inputParams/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") optimisedPrm <- optim(inputParams, fn =  BM2EB, control = cl, method = meth, lower = c(NA, 0, -20), upper=c(NA, NA, -1e-6), hessian = hess, data=data, timeMin=timeMin)
        else optimisedPrm <- optim(inputParams, fn =  OU2EB, control = cl, method = meth, lower = c(NA, 0, -20), upper=c(NA,NA,-1e-6), hessian = hess, data=data)
        if (hess)   optimisedPrm$se <- sqrt(diag(-1 * solve(optimisedPrm$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "disparityBM2EB", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))

    }

########## BM 2 Trend


    BM2Trend <- function (p, data, timeMin)  {
    
    
        ancState <- p[1]
        sigmaSquared <- p[2]
        trendParam <- p[3]
    
        NN <- length(data$central_tendency)
        VCV <- sigmaSquared * outer(data$subsamples, data$subsamples, FUN = pmin)
        diag(VCV) <- diag(VCV) + data$variance / data$sample_size
        MOne <- rep(ancState, NN)
                
        VCV[c(timeMin: dim(VCV)[1]), ] <- 0
        VCV[ ,c(timeMin: dim(VCV)[1])] <- 0
                
        
        ########
        ########
        
        timeOut <-  data$subsamples[-c(1 : (timeMin - 1))]
        timeOutDiff <- diff(timeOut[1:2])
        timeOut2 <- timeOut - (min(timeOut) - timeOutDiff)
        
        NN <- length(timeOut)
        M2 <- rep(ancState, NN) 
         
        VCV2 <- sigmaSquared * outer(timeOut, timeOut, FUN = pmin)
        diag(VCV2) <- diag(VCV2) + data$variance[-c(1 : (timeMin - 1))] / data$sample_size[-c(1 : (timeMin - 1))]
        M2 <- M2 + trendParam * timeOut2
        
        VCV[timeMin:dim(VCV)[1], timeMin:dim(VCV)[1]] <- VCV2
        
        MM <- c(MOne[1:(timeMin - 1)], M2)      
        
        S <- dmnorm(t(data$central_tendency), mean = MM, varcov = VCV, log = TRUE)
        return(S)
   
    }
    
    
        disparityBM2Trend <- function (data, timeMin=NULL, poolV = TRUE, cl = list(fnscale = -1), meth = "L-BFGS-B", hess = FALSE) {
    
        if (poolV)  data <- pooledVariance(data, rescaleVarObject = TRUE)
        if (data$subsamples[1] != 0)  stop("use '0' for first age")
        inputParams <- array(dim = 3)
        inputParams[1] <- data$central_tendency[1]
        inputParams[2] <- MLbm(data) #min(c(MLbm(data), 1e-07))
        inputParams[3] <-MLtrend(data)[2]
        timeMin <- which.min(abs(data$subsamples - timeMin))
        names(inputParams) <- c("ancState", "sigmaSquared", "trend")
        if (is.null(cl$ndeps)) cl$ndeps <- abs(inputParams/10000)
        cl$ndeps[cl$ndeps == 0] <- 1e-08
        if (meth == "L-BFGS-B") optimisedPrm <- optim(inputParams, fn = BM2Trend, control = cl, method = meth, lower = c(NA, 0, -100), upper=c(NA,NA, 100), hessian = hess, data=data, timeMin=timeMin)
        else optimisedPrm <- optim(inputParams, fn =  OU2EB, control = cl, method = meth, lower = c(NA, 0, -100), upper=c(NA,NA, 100), hessian = hess, data=data)
        if (hess)   optimisedPrm$se <- sqrt(diag(-1 * solve(optimisedPrm$hessian)))
        else optimisedPrm$se <- NULL
        k <- length(optimisedPrm$par)
        n <- length(data[[1]])
        aic <- (-2 * optimisedPrm$value) + (2 * k)
        aicc <- (-2 * optimisedPrm$value) + (2 * k) * (n / ( n - k - 1))
        return(list(logL = optimisedPrm$value, parameters = optimisedPrm$par, modelName = "disparityBM2Trend", method = "Joint", K = k, n = length(data$central_tendency), se = optimisedPrm$se, AIC=aic, AICc=aicc))



    }




##########################################################################


    simulateOU <-  function (nSims=1, sigmaSquared=1,  ancState=0.15, alpha=1, optima=0.15,  nOptima=1, timeSplit=NULL, timeSpan=100, variance=1e-1, sample_size=100)  {
        
        NN <- length(1:timeSpan)
        nTheta <- nOptima

        if(length(optima) != NN) optima <- rep(optima, NN)
        if(length(variance) != NN) variance <- rep(variance, NN)
        if(length(sample_size) != NN) sample_size <- rep(sample_size, NN)      
        timeSplit   <- c(0, timeSplit, timeSpan)        
        timeSpan <- 0:(timeSpan-1)
        
        FF <- function(a, b) abs(a - b)
        VCV <- outer(timeSpan, timeSpan, FUN = FF)
        VCV <- exp(-alpha * VCV)
        ou.V <- function (sigmaSquared, alpha, time) (sigmaSquared / (2 * alpha)) * (1 - exp(-2 * alpha * time))  
        VCVd <- ou.V(sigmaSquared, alpha, timeSpan)
        VCV2 <- outer(VCVd, VCVd, pmin)
        VCV <- VCV * VCV2
        diag(VCV) <- VCVd + variance / sample_size
        ou.M <- function (ancState, optima, alpha, time) optima * (1 - exp(-alpha * time)) + ancState * exp(-alpha * time)
        ouMean <- sapply(1:nTheta, function(x) ou.M(ancState, optima[x], alpha, timeSpan))
        meanOU <- c()
        for(x in 1:nOptima) meanOU <- c(meanOU, ouMean[(timeSplit[x] + 1) : timeSplit[x + 1] , x])  
        central_tendency <- t(rmnorm(nSims, mean=meanOU, varcov = VCV))[1,]
        variance <- variance
        nn <- sample_size
        tt <- timeSpan
        data <- list(central_tendency, variance, nn, tt)
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
    
    }
    
    
    simulateBM <-  function (nSims=1, sigmaSquared=1, timeSpan=100, variance=1e-1, sample_size=100, ancState=0)  {
    
        NN <- length(1:timeSpan)
        if(length(variance) != NN) variance <- rep(variance, NN)
        if(length(sample_size) != NN) sample_size <- rep(sample_size, NN)
        timeSpan <- 0:(timeSpan-1)
        VCV <- sigmaSquared * outer(timeSpan, timeSpan, FUN = pmin)
        diag(VCV) <- diag(VCV) + variance / sample_size
        M <- rep(ancState, NN)
        central_tendency <- t(rmnorm(nSims, mean=M, varcov = VCV))[1,]
        variance <- variance
        nn <- sample_size
        tt <- timeSpan
        data <- list(central_tendency, variance, nn, tt)
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
    
    }
    
    
    simulateStasis <- function (nSims=1, theta=1, omega=1, timeSpan=100, timeSplit=NULL, nOptima=1, variance=1e-1, sample_size=100) {
                
        
        NN <- length(1:timeSpan)  
        if(length(variance) != NN) variance <- rep(variance, NN)
        if(length(sample_size) != NN) sample_size <- rep(sample_size, NN)
        timeSplit   <- c(0, timeSplit, timeSpan)
        VCV <- diag(omega + variance/sample_size)
        meanStasis <- c()
        for(x in 1:nOptima) meanStasis <- c(meanStasis, rep(theta[x], length(timeSplit[x] : timeSplit[x + 1]) - 1))    
        central_tendency <- t(rmnorm(nSims, mean=meanStasis, varcov = VCV))[1,]
        variance <- variance
        nn <- sample_size
        tt <- 0:(timeSpan-1)
        data <- list(central_tendency, variance, nn, tt)
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
    
   
    }
    
  
   simulateEB <- function (nSims=1, sigmaSquared=1, rRate=-0.1, timeSpan=100, variance=1e-1, sample_size=10, ancState=0)  {
        
        NN <- length(1:timeSpan)  
        if(length(variance) != NN) variance <- rep(variance, NN)
        if(length(sample_size) != NN) sample_size <- rep(sample_size, NN)       
        timeSpan <- 0:(timeSpan-1)
        VCV <- outer(sigmaSquared * ((exp(rRate * timeSpan) - 1) / rRate), sigmaSquared * ((exp(rRate * timeSpan) - 1) / rRate), FUN=pmin)  
        diag(VCV) <- diag(VCV) + variance / sample_size
        M <- rep(ancState, NN)
        central_tendency <- t(rmnorm(nSims, M, varcov = VCV))[1,]
        variance <- variance
        nn <- sample_size
        tt <- timeSpan
        data <- list(central_tendency, variance, nn, tt)
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
    
}


    simulateTrend <-  function (nSims=1, sigmaSquared=1, ancState=0.1, trendParam=0.1, timeSpan=100, variance=1e-1, sample_size=100)  {
    
        NN <- length(1:timeSpan)
        if(length(variance) != NN) variance <- rep(variance, NN)
        if(length(sample_size) != NN) sample_size <- rep(sample_size, NN)
        timeSpan <- 0:(timeSpan-1)
        VCV <- sigmaSquared * outer(timeSpan, timeSpan, FUN = pmin)
        diag(VCV) <- diag(VCV) + variance / sample_size
        M <- rep(ancState, NN) + trendParam * timeSpan
        central_tendency <- t(rmnorm(nSims, mean=M, varcov = VCV))[1,]
        variance <- variance
        nn <- sample_size
        tt <- timeSpan
        data <- list(central_tendency, variance, nn, tt)
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
    
    }
    
    
        simulateShiftMode <- function(shiftTime, modeOne, modeTwo, sigmaSquared=c(1, 1), ancState=0.1, trendParam=c(0.1, 0.1), variance=c(1e-1, 1e-1), rRate=c(-0.1, -0.1), alpha=c(1, 1), optima=c(0.15, 0.15), timeSpan=100, sample_size=c(100, 100)) {
                        
        nSimsOne <-     timeSpan - shiftTime
        nSimsTwo <-     timeSpan - nSimsOne

        timeSpanOne <- 0:(nSimsOne)
        timeSpanTwo <- 0:(nSimsTwo)
        
        nSimsOne <- length(timeSpanOne ) -1 
        nSimsTwo <- length(timeSpanTwo) -1
        
        if(modeOne == "OU") modelOne <- simulateOU(sigmaSquared=sigmaSquared[1],  ancState=ancState[1], alpha=alpha[1], optima=optima[1],  nOptima=1, timeSpan=nSimsOne[1], variance=variance[1], sample_size=sample_size[1])
        if(modeOne == "BM") modelOne <- simulateBM(sigmaSquared=sigmaSquared[1], nSimsOne=nSimsOne[1], variance=variance[1], sample_size=sample_size[1])
        if(modeOne == "EB") modelOne <- simulateEB(sigmaSquared=sigmaSquared[1], rRate=rRate[1], nSimsOne=nSimsOne[1], variance=variance[1], sample_size=sample_size[1], timeSpan=nSimsOne[1])
        if(modeOne == "Trend")  modelOne <- simulateTrend(sigmaSquared=sigmaSquared[1], ancState=ancState[1], trendParam=trendParam[1], variance=variance[1], sample_size=sample_size[1], timeSpan=nSimsOne[1])
        
        ancState2 <- tail(modelOne[[1]], 1)
        
        if(modeTwo == "OU") modelTwo <- simulateOU(sigmaSquared=sigmaSquared[2],  ancState=ancState2, alpha=alpha[2], optima=optima[2],  nOptima=1, timeSpan=nSimsTwo, variance=variance[2], sample_size=sample_size[2])
        if(modeTwo == "EB") modelTwo <- disparityEB(sigmaSquared=sigmaSquared[2], rRate=rRate[2], nSimsOne=nSimsTwo, variance=variance[2], sample_size=sample_size[2], ancState=ancState2)
        if(modeTwo == "Trend")  modelTwo <- simulateTrend(sigmaSquared=sigmaSquared[2], ancState=ancState2, trendParam=trendParam[2], variance=variance[2], sample_size=sample_size[2], timeSpan=nSimsTwo)
        
        datasetOne <- lapply(modelOne, function(x) x)
        datasetTwo <- lapply(modelTwo, function(x) x)
        datasetTwo[[4]] <- datasetTwo[[4]] - datasetTwo[[4]][1]
        
        timeSpanOne <- datasetOne[[4]]
        timeSpanTwo <- datasetTwo[[4]]
        
        central_tendency <- unlist(c(modelOne[[1]], modelTwo[[1]]))
        variance <- unlist(c(modelOne[[2]], modelTwo[[2]]))
        sample_size <- unlist(c(modelOne[[3]], modelTwo[[3]]))
        time <- c(timeSpanOne, tail(timeSpanOne, 1) + 1 + timeSpanTwo)
    
        data <- list(central_tendency, variance, sample_size, time)
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
                
        }       
        
            simulateOU2Trend <- function (nSims=1, sigmaSquared=1, ancState=0.15, alpha=1, optima=0.15, trendParam=0.1, timeMin = 50, timeSpan=100, variance=1e-1, sample_size=100)  {
        
        timeSpanAll <- 0:(timeSpan-1)
        
        if(length(optima) != length(timeSpanAll)) optima <- rep(optima, length(timeSpanAll))
        if(length(variance) != length(timeSpanAll)) variance <- rep(variance, length(timeSpanAll))
        if(length(sample_size) != length(timeSpanAll)) sample_size <- rep(sample_size, length(timeSpanAll))    
    
        splitOne <- timeMin - 1
        timeSpan <- timeSpanAll
        
        FF <- function(a, b) abs(a - b)
        VCV <- outer(timeSpan, timeSpan, FUN = FF)
        VCV <- exp(-alpha * VCV)
        ou.V <- function (sigmaSquared, alpha, time) (sigmaSquared / (2 * alpha)) * (1 - exp(-2 * alpha * time))  
        VCVd <- ou.V(sigmaSquared, alpha, timeSpan)
        VCV2 <- outer(VCVd, VCVd, pmin)
        VCV <- VCV * VCV2
        diag(VCV) <- VCVd + variance / sample_size
        ou.M <- function (ancState, optima, alpha, time) optima * (1 - exp(-alpha * time)) + ancState * exp(-alpha * time)
        ouMean <- ou.M(ancState, optima, alpha, timeSpan)
        
                    
        VCV[c((timeMin): dim(VCV)[1]), ] <- 0
        VCV[ ,c((timeMin): dim(VCV)[1])] <- 0   
        
        ########
        ########

        timeOut <-  timeSpanAll[-c(1 : (timeMin - 1))]
        timeOutDiff <- diff(timeOut[1:2])
        timeOut2 <- timeOut - (min(timeOut) - timeOutDiff)
        
        NN <- length(timeOut)
        M2 <- rep(ancState, NN) 
         
        VCV2 <- sigmaSquared * outer(timeOut, timeOut, FUN = pmin)
        diag(VCV2) <- diag(VCV2) + variance[-c(1 : (timeMin - 1))] / sample_size[-c(1 : (timeMin - 1))]
        M2 <- M2 + trendParam * timeOut2
        
        VCV[timeMin:dim(VCV)[1], timeMin:dim(VCV)[1]] <- VCV2   
                
        MM <- c(ouMean[1:(timeMin - 1)], M2)            
                
        central_tendency <- t(rmnorm(nSims, mean=MM, varcov = VCV))[1,]
        variance <- variance
        nn <- sample_size
        tt <- timeSpan
        data <- list(central_tendency, variance, nn, tt)
        
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
    
    }
    
    
    
    
        simulateOU2EB <- function (nSims=1, sigmaSquared=1, ancState=0.15, alpha=1, optima=0.15, rRate=-0.1, timeMin = 50, timeSpan=100, variance=1e-1, sample_size=100)  {
    
        timeSpanAll <- 0:(timeSpan-1)
        
        if(length(optima) != length(timeSpanAll)) optima <- rep(optima, length(timeSpanAll))
        if(length(variance) != length(timeSpanAll)) variance <- rep(variance, length(timeSpanAll))
        if(length(sample_size) != length(timeSpanAll)) sample_size <- rep(sample_size, length(timeSpanAll))    
    
        splitOne <- timeMin - 1
        timeSpan <- timeSpanAll
        
        FF <- function(a, b) abs(a - b)
        VCV <- outer(timeSpan, timeSpan, FUN = FF)
        VCV <- exp(-alpha * VCV)
        ou.V <- function (sigmaSquared, alpha, time) (sigmaSquared / (2 * alpha)) * (1 - exp(-2 * alpha * time))  
        VCVd <- ou.V(sigmaSquared, alpha, timeSpan)
        VCV2 <- outer(VCVd, VCVd, pmin)
        VCV <- VCV * VCV2
        diag(VCV) <- VCVd + variance / sample_size
        ou.M <- function (ancState, optima, alpha, time) optima * (1 - exp(-alpha * time)) + ancState * exp(-alpha * time)
        ouMean <- ou.M(ancState, optima, alpha, timeSpan)
        
                    
        VCV[c((timeMin): dim(VCV)[1]), ] <- 0
        VCV[ ,c((timeMin): dim(VCV)[1])] <- 0   
        
        ########
        ########

        timeOut <-  timeSpanAll[-c(1 : (timeMin - 1))]
        timeOutDiff <- diff(timeOut[1:2])
        timeOut2 <- timeOut - (min(timeOut) - timeOutDiff)
        
        NN <- length(timeOut)
        M2 <- rep(ancState, NN) 
        
        VCV2 <- outer(sigmaSquared * ((exp(rRate * timeOut2) - 1) / rRate), sigmaSquared * ((exp(rRate * timeOut2) - 1) / rRate), FUN=pmin) 
        diag(VCV2) <- diag(VCV2) + variance[-c(1 : (timeMin - 1))] / sample_size[-c(1 : (timeMin - 1))]
         
   
        MM <- c(ouMean[1:(timeMin - 1)], M2)    
        
        VCV[timeMin:dim(VCV)[1], timeMin:dim(VCV)[1]] <- VCV2   
                
                
        central_tendency <- t(rmnorm(nSims, mean=MM, varcov = VCV))[1,]
        variance <- variance
        nn <- sample_size
        tt <- timeSpan
        data <- list(central_tendency, variance, nn, tt)
        
        names(data) <- c("central_tendency", "variance", "sample_size", "time")
        return(data)
    
    }
    
    
