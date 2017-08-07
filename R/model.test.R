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

#' @title BM Model Testing
#'
#' @description Applies a Brownian Motion model to disparity data
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

model.test.BM <- function(data, ...) {

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
