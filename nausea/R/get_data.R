## get fredmd data, from "fbi"

#' Get FRED-MD dataset
#'
#' @param file string address to FREDMD data
#' @param date_start Date denoting start of period, of format YYYY-MM-01
#' @param date_end Date denoting end of period, of format YYYY-MM-01
#' @param transform Boolean, apply stationarity transforms
#'
#' @return data.frame of FRED-MD data
#' @export
#'
#' @examples fred_data <- fredmd()
fredmd <- function(file = "https://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv",
                   date_start = NULL, date_end = NULL, transform = TRUE) {
  # Error checking
  if (!is.logical(transform))
    stop("'transform' must be logical.")
  if ((class(date_start) != "Date") && (!is.null(date_start)))
    stop("'date_start' must be Date or NULL.")
  if ((class(date_end) != "Date") && (!is.null(date_end)))
    stop("'date_end' must be Date or NULL.")

  if (class(date_start) == "Date") {
    if (as.numeric(format(date_start, "%d")) != 1)
      stop("'date_start' must be Date whose day is 1.")
    if (date_start < as.Date("1959-01-01"))
      stop("'date_start' must be later than 1959-01-01.")
  }

  if (class(date_end) == "Date") {
    if (as.numeric(format(date_end, "%d")) != 1)
      stop("'date_end' must be Date whose day is 1.")
  }


  # Prepare raw data
  #rawdata <- readr::read_csv(file, col_names = FALSE, col_types = cols(X1 = col_date(format = "%m/%d/%Y")),
  #                           skip = 2)
  rawdata <- read.csv(file, header = TRUE)
           #skip = 2)
  rawdata <- rawdata[2:(nrow(rawdata) - 1), ] # remove NA rows
  rawdata <- as.data.frame(rawdata)
  colnames(rawdata) <- c("date", colnames(rawdata))[1:ncol(rawdata)]

  # Import tcode tcodes is an internal data of the R package
  tcode <- c(5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 2, 2, 2, 5, 5, 2, 2, 5, 5, 5, 5, 5, 5, 5,
             5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 1, 2, 1, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 2, 6,
             6, 5, 6, 6, 7, 6, 6, 6, 2, 5, 5, 2, 5, 2, 2, 2, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 5, 5,
             5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 2, 6, 6, 6, 6, 1)


  # Subfunction transxf: data transformation based on tcodes
  transxf <- function(x, tcode) {
    # Number of observations (including missing values)
    n <- length(x)

    # Value close to zero
    small <- 1e-06

    # Allocate output variable
    y <- rep(NA, n)
    y1 <- rep(NA, n)

    # TRANSFORMATION: Determine case 1-7 by transformation code
    if (tcode == 1) {
      # Case 1 Level (i.e. no transformation): x(t)
      y <- x

    } else if (tcode == 2) {
      # Case 2 First difference: x(t)-x(t-1)
      y[2:n] <- x[2:n] - x[1:(n - 1)]

    } else if (tcode == 3) {
      # case 3 Second difference: (x(t)-x(t-1))-(x(t-1)-x(t-2))
      y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]

    } else if (tcode == 4) {
      # case 4 Natural log: ln(x)
      if (min(x, na.rm = TRUE) > small)
        y <- log(x)

    } else if (tcode == 5) {
      # case 5 First difference of natural log: ln(x)-ln(x-1)
      if (min(x, na.rm = TRUE) > small) {
        x <- log(x)
        y[2:n] <- x[2:n] - x[1:(n - 1)]
      }

    } else if (tcode == 6) {
      # case 6 Second difference of natural log:
      # (ln(x)-ln(x-1))-(ln(x-1)-ln(x-2))
      if (min(x, na.rm = TRUE) > small) {
        x <- log(x)
        y[3:n] <- x[3:n] - 2 * x[2:(n - 1)] + x[1:(n - 2)]
      }

    } else if (tcode == 7) {
      # case 7 First difference of percent change:
      # (x(t)/x(t-1)-1)-(x(t-1)/x(t-2)-1)
      y1[2:n] <- (x[2:n] - x[1:(n - 1)])/x[1:(n - 1)]
      y[3:n] <- y1[3:n] - y1[2:(n - 1)]
    }

    return(y)
  }


  # Transform data
  if (transform) {
    # Apply transformations
    N <- ncol(rawdata)
    data <- rawdata
    data[, 2:N] <- NA

    # Perform transformation using subfunction transxf (see below for
    # details)
    for (i in 2:N) {
      temp <- transxf(rawdata[, i], tcode[i - 1])
      data[, i] <- temp
    }

  } else {
    data <- rawdata
  }


  # Null case of date_start and date_end
  if (is.null(date_start))
    date_start <- as.Date("1959-01-01")
  if (is.null(date_end))
    date_end <- data[, 1][nrow(data)]


  # Subset data
  index_start <- which.max(data[, 1] == date_start)
  index_end <- which.max(data[, 1] == date_end)

  outdata <- data[index_start:index_end, ]
  class(outdata) <- c("data.frame", "fredmd")
  return(outdata)

}

# transmap <- function(x){
#   out1 <- x[-1,-1]
#   trans_vec <- x[1,-1]
#   for (trans in 1:length(trans_vec)) {
#     if( trans_vec[trans] == 2) out1[,trans] <- c(NA,diff(out1[,trans]))
#     if( trans_vec[trans] == 3) out1[,trans] <- c(NA,NA,diff(out1[,trans], differences =  2))
#     if( trans_vec[trans] == 4) out1[,trans] <- log(out1[,trans])
#     if( trans_vec[trans] == 5) out1[,trans] <- c(NA,diff( log(out1[,trans]) ))
#     if( trans_vec[trans] == 6) out1[,trans] <- c(NA,NA,diff( log(out1[,trans]), differences = 2 ))
#     if( trans_vec[trans] == 7) out1[,trans] <- c(NA, NA,diff(out1[-1,trans]/(out1[-nrow(out1),trans]) ))
#   }
#   out <-  data.frame(date = x[-1,1] ,out1)
#   return(out)
# }


## get data

#' Get entire nowcasting dataset
#'
#' @param observation_end Date denoting end of period, of format YYYY-MM-01
#'
#' @return Balanced panel of nowcasting dataset
#' @export
#'
#' @examples nowcasting_data <- get_data()
get_data <- function(observation_end = "2021-06-01"){
  #FRED-MD monthly data
  current <- fredmd("https://files.stlouisfed.org/files/htdocs/fred-md/monthly/current.csv",
                         date_start = NULL, date_end = NULL)
  #get NYFED series from ALFRED
  try_series <- alfred::get_fred_series("PAYEMS", series_name = "PAYEMS", observation_start = "2002-01-01", observation_end = observation_end)

  for (name in colnames(nowcasting::NYFED$base)[-1]) {
    new_series <- alfred::get_fred_series(name, series_name = name, observation_start = "2002-01-01", observation_end = observation_end)
    try_series <- dplyr::full_join(try_series, new_series, by = "date")
  }

  nyfed_trans <- nowcasting::NYFED$legend$Transformation
  nyfed_trans[25] <- 6 #smoothing transform
  alfred_ts <- ts(try_series[,-1], start = c(2002,1), end = c(2021, 6), frequency = 12)
  data_alfred <- nowcasting::Bpanel(base = alfred_ts,
                        trans = nyfed_trans,#NYFED$legend$Transformation,
                        aggregate = FALSE, na.prop = .75)
  #add surveys
  survey_series <- dplyr::full_join(alfred::get_fred_series("CFNAIDIFF", series_name = "CFNAIDIFF", observation_start = "2004-06-01", observation_end = observation_end),
                             alfred::get_fred_series("BACTSAMFRBDAL", series_name = "BACTSAMFRBDAL", observation_start = "2004-06-01", observation_end = observation_end))
  survey_ts <- ts(survey_series[,-1], start = c(2004,6), end = c(2021, 6), frequency = 12)
  survey_trans <- c(nyfed_trans, 0,0)
  data_survey <-  ts.intersect(data_alfred, survey_ts)

  # fredmd_0 <- current#read.csv("current.csv")
 #head(fredmd)
  # fredmd_1 <- fredmd_0[-1,]
  # ts.plot(fredmd_1[,2:5])
  # fredmd_2 <- transmap(fredmd_0[-(2:500),])
  # ts.plot(fredmd_2[-1,-1])
  #
  # fredmd_3 <- fredmd_2[-(1:17),]
  # fredmd_3[,1] <- as.Date(fredmd_3[,1], format = "%m/%d/%y")
  #duplicates <- which(colnames(current) %in% c(colnames(data_alfred), colnames(survey_series)))
  fredmd_4 <- current#[,-duplicates[-1]]
  fredmd_ts <- ts(fredmd_4[,-1], start = c(2004,6), end = c(2021, 6), frequency = 12)
 # fredmd_trans <- c(rep(0, ncol(fredmd_4)-2), survey_trans)
  data_fredmd <- nowcasting::Bpanel(base = fredmd_ts,
                        trans = rep(0, ncol(fredmd_4)-1),
                        aggregate = FALSE, na.prop = .75, NA.replace = T)
  data_big <- ts.intersect(data_fredmd, data_survey)

#54,97,
  mosumvar_ts <- ts(data_big[,-c(146)], start = time(data_big)[1], end = time(data_big)[nrow(data_big)], frequency = 12)
  #mosumvar_gdp <- diff(alfred_ts[,25],3)
  data_mosumvar <- nowcasting::Bpanel(base = mosumvar_ts,
                          trans = rep(0,ncol(mosumvar_ts)),
                          aggregate = FALSE, na.prop = 1, h = 12)
  
  return(data_mosumvar)
}
