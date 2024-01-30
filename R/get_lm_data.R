#' Build a landmark dataset
#'
#' @param data Data frame from which to construct landmark super dataset
#' @param outcome A list with items time and status, containing character
#'   strings identifying the names of time and status variables, respectively,
#'   of the survival outcome
#' @param lm The value of the landmark time point at which to construct the
#'   landmark dataset.
#' @param horizon Scalar, the value of the prediction window (ie predict risk
#'   within time w landmark points)
#' @param covs A list with items fixed and varying, containing character strings
#'   specifying column names in the data containing time-fixed and time-varying
#'   covariates, respectively.
#' @param format Character string specifying whether the original data are in
#'   wide (default) or in long format.
#' @param id Character string specifying the column name in data containing the
#'   subject id.
#' @param rtime Character string specifying the column name in data containing
#'   the (running) time variable associated with the time-varying variables;
#'   only needed if format = "long".
#' @param left.open Boolean (default = FALSE), indicating if the intervals for the
#'   time-varying covariates are open on the left (and closed on the right) or
#'   vice-versa.
#' @param split.data List of data split according to ID. Allows for faster
#'   computation.
#'
#' @details This function is based from [dynpred::cutLM()] with minor changes.
#'   The original function was authored by Hein Putter.
#' @references van Houwelingen HC, Putter H (2012). Dynamic Prediction in
#'   Clinical Survival Analysis. Chapman & Hall.
#' @return A landmark dataset.
#' @export
#'
# TODO: add examples
# TODO: add references


get_lm_data <- function (data, outcome, lm, horizon, covs, static_covs,
                         format = c("wide", "long"), id, rtime, right = TRUE) {
  format <- match.arg(format)
  if (format == "wide") {
    lmdata <- data
    if (!is.null(covs$varying)){
      for (col in covs$covarying)
        lmdata[[col]] <- 1 - as.numeric(lmdata[[col]] > lm)
    }


  }
  else {
    if (missing(id))
      stop("argument 'id' should be specified for long format data")
    if (missing(rtime))
      stop("argument 'rtime' should be specified for long format data")
    ord <- order(data[[id]], data[[rtime]])
    data <- data[ord, ] # order the data per id and rtime = time_index, not cut to lms
    ids <- unique(data[[id]]) # get all of the ids
    n <- length(ids) # number of different people
    lmdata <- data[which(!duplicated(data[[id]])), ] # take the first line for each id


    # add an empty time_to_ct column
    time_to_ct <- "time_to_ct"
    lmdata[time_to_ct] <- NA
    for (i in 1:n) {
      wh <- which(data[[id]] == ids[i])
      di <- data[wh, ] # get the data for each id


      # find the closest time_to_event to the landmark time
      closest <- which.min(abs(di[[rtime]] - lm))
      # get the time distance to the landmark
      time_to_lm <- di[[rtime]][closest] - lm

      # check if the closest is less than 3 years 
      if (abs(time_to_lm) <= 3){
        lmdata[i, ] <- di[closest, ]
        lmdata[i, time_to_ct] <- time_to_lm
      }
      else{
        lmdata[i, ] <- NA
        lmdata[i, time_to_ct] <- NA
        # static variables should be copied over
        lmdata[i, static_covs] <- di[closest, static_covs]
      }

    }




  }
  lmdata <- lmdata[lmdata[[outcome$time]] > lm, ]
  if (format == "long")
    lmdata <- lmdata[!is.na(lmdata[[id]]), ]
  lmdata[outcome$status] <- lmdata[[outcome$status]] *
    as.numeric(lmdata[[outcome$time]] <= horizon)
  lmdata[outcome$time] <- pmin(as.vector(lmdata[[outcome$time]]), horizon)
  lmdata$LM <- lm
  if (format == "long")
    cols <- match(c(id, outcome$time, outcome$status, covs$fixed,
                    covs$varying, rtime, "time_to_ct", "LM"), names(lmdata))
  else cols <- match(c(outcome$time, outcome$status, covs$fixed,
                       covs$varying, "time_to_ct", "LM"), names(lmdata))

  # replace NAs with the most frequent value of the column, this is done at the stack level
  lmdata <- lmdata %>% mutate_if(is.factor, function(x) replace(x, is.na(x), names(sort(table(x), decreasing = TRUE)[1])))
  return(lmdata[, cols])
}

