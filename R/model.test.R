#' @title Disparity Model Testing
#'
#' @description Test one or more models of disparity evolution
#'
#' @param data a \code{dispRity} object.
#' @param models one or more model names, can be \code{"EB"}, \code{"BM"}, \code{"OU"}, \code{"Trend"} and/or \code{"Shift"}.
#' @param parameters a \code{list} of model parameters to be passed to each model.
#' 
#' @details
#' \itemize{
#'      \item{\code{"EB"}} an Early-Burst model (see \link{\code{model.test.EB}} for more details)
#'      \item{\code{"BM"}} a Brownian-Motion model (see \link{\code{model.test.BM}} for more details)
#'      \item{\code{"OU"}} an Ornstein-Uhlenbeck model (see \link{\code{model.test.OU}} for more details)
#'      \item{\code{"Trend"}} a Trend model (see \link{\code{model.test.Trend}} for more details)
#'      \item{\code{"Shift"}} a model with a shift between any of the above models (see \link{\code{model.test.Shift}} for more details)
#' }
#' 
#' @examples
#'
#' @seealso \link{\code{test.dispRity}}, \link{\code{model.test.EB}}, \link{\code{model.test.BM}}, \link{\code{model.test.OU}}, \link{\code{model.test.Trend}}, \link{\code{model.test.Shift}}
#' 
#' @author Mark Puttick and Thomas Guillerme
#' @export


# ebModel <- disparityEB(listFull)$AICc
# bmModel <- disparityBM(listFull)$AICc
# ouModel <- disparityOU(listFull)$AICc
# trendModel <- disparityTrend(listFull)$AICc
# multiouModel <- disparityOU(listFull, nOptima=2, timeSplit=54)$AICc
# ou2trendModel <- disparityOU2Trend(listFull, timeMin=54)$AICc
# ou2ebModel <- disparityOU2EB(listFull, timeMin=54)$AICc
# bm2trendModel <- disparityBM2Trend(listFull, timeMin=54)$AICc
# bm2ebModel <- disparityBM2EB(listFull, timeMin=54)$AICc




model.test <- function(data, models, parameters) {

    return(0)
}


#' @title BM Model Testing
#'
#' @description Applies a Brownian Motion model to disparity data
#'
#' @param data a \code{dispRity} object or a \code{list} containing the mean, variance, samplesize and subsamples names on which to test the model.
#' @param pool.variance \code{logical} whether to pool the variance or not (default = \code{TRUE}).
#' @param control ???
#' @param ... the model parameters
#' 
#' @examples
#'
#' @seealso \link{\code{model.test}}
#' 
#' @author Mark Puttick and Thomas Guillerme
#' @export

model.test.BM <- function(data, pool.variance = TRUE, control = list(fnscale = -1), method = "L-BFGS-B", hessian = FALSE) {

        ## data is just to be consistent if the input is dispRity
        data_list <- data

        if(class(data_list) == "dispRity") {
            ## If the input is of class dispRity, first create the model list
            data_list <- select.model.list(data_list)
        }

        if (pool.variance) { #TG: why do we want that to be a default?
            ## Rescale the object's variance
            data_list <- pooled.variance(data_list, rescale.variance = TRUE)
        }

        #if (dataFile$time[1] != 0)  stop("use '0' for first age")
        #TG: I think this should not be a requirement select.model.list returns list of subsamples from smaller value to bigger

        #TG: I've stored the other inputOptimisation steps in BM.parameters for optimisations
        input_optimisation <- BM.parameters(data_list) #min(c(MLbm(dataFile), 1e-07))


        ## Setting the precision
        if (is.null(control$ndeps)) {
            control$ndeps <- abs(input_optimisation/10000)
        }

        control$ndeps[control$ndeps == 0] <- 1e-08
    
        if (method == "L-BFGS-B") {
            optimised_parameters <- optim(input_optimisation, fn = BMLikelihoodJointOptim, control = control, method = method, lower = c(NA, 1e-6), hessian = hessian, dataFile = dataFile)
            } else {
            optimised_parameters <- optim(input_optimisation, fn = BMLikelihoodJointOptim, control = control, method = method, hessian = hessian, dataFile = dataFile)
            }
        if (hessian) 
            optimised_parameters$se <- sqrt(diag(-1 * solve(w$hessian)))
        else optimised_parameters$se <- NULL
        k <- length(optimised_parameters$par)
        n <- length(dataFile[[1]])
        aic <- (-2 * optimised_parameters$value) + (2 * k)
        aicc <- (-2 * optimised_parameters$value) + (2 * k) * (n / ( n - k - 1))


        return(list(logL = optimised_parameters$value, parameters = optimised_parameters$par, modelName = "BM", method = "Joint", K = k, n = length(dataFile$meanVector), se = optimised_parameters$se, AIC=aic, AICc=aicc))
}


#' @title EB Model Testing
#'
#' @description Applies an Early Burst model to disparity data
#'
#' @param data a \code{dispRity} object.
#' @param ... the model parameters
#' 
#' @examples
#'
#' @seealso \link{\code{model.test}}
#' 
#' @author Mark Puttick and Thomas Guillerme
#' @export

model.test.EB <- function(data, ...) {

    return()
}

#' @title OU Model Testing
#'
#' @description Applies an Ornstein-Uhlenbeck model to disparity data
#'
#' @param data a \code{dispRity} object.
#' @param ... the model parameters
#' 
#' @examples
#'
#' @seealso \link{\code{model.test}}
#' 
#' @author Mark Puttick and Thomas Guillerme
#' @export

model.test.OU <- function(data, ...) {

    return()
}



#' @title Trend Model Testing
#'
#' @description Applies a Trend model to disparity data
#'
#' @param data a \code{dispRity} object.
#' @param ... the model parameters
#' 
#' @examples
#'
#' @seealso \link{\code{model.test}}
#' 
#' @author Mark Puttick and Thomas Guillerme
#' @export

model.test.Trend <- function(data, ...) {

    return()
}

#' @title Shift Model Testing
#'
#' @description Applies a Shift model to disparity data
#'
#' @param data a \code{dispRity} object.
#' @param shifts where to apply the shifts.
#' @param models the name of the models to shift (\code{"EB"}, \code{"BM"}, \code{"OU"} or \code{"Trend"}).
#' @param parameters a \code{list} of model parameters to be passed to each model.
#' 
#' @examples
#'
#' @seealso \link{\code{model.test}}
#' 
#' @author Mark Puttick and Thomas Guillerme
#' @export

model.test.Shift <- function(data, models, parameters, ...) {

    return()
}
