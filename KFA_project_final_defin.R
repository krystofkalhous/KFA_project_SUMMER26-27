# Use of AI
# Tool: Claude (Anthropic), model Claude Opus 5.5, used through claude.ai, October 2026. -> Used to compile the code and for discussion of the results

### 1.) Setup: packages, paths, constants ---------------------------------------

# loading packages
library(ncdf4)
library(mgcv)

# We are in: 
getwd() # -> "/Users/krystofkalhous/Desktop/Projekt KFA"

# We are working with the datasets analysed in the thesis: Nonlinearity of teleconnections in the climate system (Petr Skala; 2025)
## which is available at: https://dspace.cuni.cz/handle/20.500.11956/203052

# We obtained the data from:
# NOTE: GISTEMP is updated every month (past values change slightly too), we downloaded the file on 8.9.2026
# 1.) GISTEMP (gistemp1200_GHCNv4_ERSSTv5.nc): https://data.giss.nasa.gov/gistemp/data_v4.html
# 2.) NAO_FILE (nao_jones.txt): https://climexp.knmi.nl/data/inao.dat
# 3.) NINO_FILE (nino34.txt): https://climexp.knmi.nl/data/iersst_nino3.4a.dat
# 4.) NINO_REL (nino34_rel.txt): https://climexp.knmi.nl/data/iersst_nino3.4a_rel.dat
# 5.) SKALA_NAO, SKALA_NINO: the nonlinearity classes of Skala (2025), figs. 5.4 and 5.2 (bottom
#     panels), digitised from the thesis PDF (0 = 0-1 sigma ... 4 = above 4 sigma)

BASE      <- getwd()                                                            # On my machine: "/Users/krystofkalhous/Desktop/Projekt KFA"
DATA      <- file.path(BASE, "Data")                                            # On my machine: "/Users/krystofkalhous/Desktop/Projekt KFA/Data"

GISTEMP   <- file.path(DATA, "gistemp1200_GHCNv4_ERSSTv5.nc")                   
NAO_FILE  <- file.path(DATA, "nao_jones.txt")                                   # CRU/Jones, Gibraltar-Stykkisholmur
NINO_FILE <- file.path(DATA, "nino34.txt")                                      # ERSSTv5 NINO3.4
NINO_REL  <- file.path(DATA, "nino34_rel.txt")                                  # same, minus the 20S-20N mean


SKALA_NAO  <- file.path(DATA, "skala_nao_nonlin_class.csv")                     # (!!!) optional: comparison with Skala (section 7). I haven't yet asked Skala for his original grid values: these classes were digitised from the colours of his published figures (see 5.) above), with the help of Claude Opus 5.5, so single cells may be misread. Once his data are obtained, the comparison can be rerun unchanged.
SKALA_NINO <- file.path(DATA, "skala_nino34_nonlin_class.csv")

# stop right here with a clear message if the working directory is wrong
stopifnot(all(file.exists(GISTEMP, NAO_FILE, NINO_FILE, NINO_REL)))

# Time window of our analysis is consistent with Skala (2025)
WIN_START <- as.Date("1880-01-01")
WIN_END   <- as.Date("2023-12-01")
N_MONTHS  <- 1728
ALPHA     <- 0.05 # chosen statistical significance threshold

# What to run?

RUN_PILOT  <- FALSE  # section 5: checks before the global run (1-2 minutes)
RUN_GLOBAL <- FALSE # section 6: every grid cell -> OUT_CSV (MacBook Air 2021 (M1): about 4 hours with RHO_METHOD = "reml", 1.5 hours with "acf")
RUN_RESULTS <- TRUE # section 7: BH correction, summary, comparison with Skala, fields for the maps (seconds)

# How the AR(1) parameter rho is chosen in every cell (section 3, ar1_fit):
#   "reml": the rho that minimises the REML score (Wood 2017, sec. 7.7.2) - the main analysis
#   "acf" : the lag-1 autocorrelation of the residuals of a fit with rho = 0 - our first run, kept
#           for comparison (it underestimates rho over the oceans; see ar1_fit)
RHO_METHOD <- "reml"
stopifnot(RHO_METHOD %in% c("reml", "acf"))

# Output: a results file and a map folder per method, so the two runs never mix
OUT_CSV   <- file.path(BASE, paste0("gam_grid_results_", RHO_METHOD, ".csv"))  # gam_grid_results_reml.csv or _acf.csv
MAPDIR    <- file.path(BASE, if (RHO_METHOD == "reml") "maps" else "maps_acf") # gam_maps_final.py draws maps/ (the main analysis)


### 2.) Reading the data -------------------------------------------------------

#' floor_month - converts any date to the first day of the given month. This ensures 
#' comparability between the elements of different time series with elements measured
#' on a monthly basis

floor_month <- function(d) {
  return(as.Date(format(d, "%Y-%m-01")))
  }

#' read_climexp - loads the Climate Explorer data into format preferred for the analysis
#' i.e. converts the data into "appropriately" ordered data frame

read_climexp <- function(path, missing = -999.9) {
  raw <- read.table(path, comment.char = "#", header = FALSE)
  
  if (ncol(raw) >= 13) {
    vals <- as.matrix(raw[, 2:13])                 # months only
    vals[vals <= missing + 1e-6] <- NA
    out <- data.frame(year  = rep(raw[[1]], each = 12),
                      month = rep(1:12, times = nrow(raw)),
                      value = as.vector(t(vals)))
  } else if (ncol(raw) == 2) {
    yr <- floor(raw[[1]])
    m  <- round((raw[[1]] - yr) * 12) + 1 # converts fractional part of the year into the corresponding month number
    over <- m == 13 # converts "13th month" of year n to 1st month of year n+1
    yr[over] <- yr[over] + 1
    m[over] <- 1
    v <- raw[[2]]
    v[v <= missing + 1e-6] <- NA
    out <- data.frame(year = yr, month = m, value = v)
  } else {
    stop("unrecognised Climate Explorer format: ", ncol(raw), " columns")
  }
  
  out$date <- as.Date(sprintf("%04d-%02d-01", out$year, out$month))
  return(out)
}


#' open_gistemp - opens the data from data.giss.nasa.gov, ie the file at GISTEMP
#' and converts it to a list

open_gistemp <- function(path = GISTEMP) {
  nc     <- nc_open(path)
  tunits <- ncatt_get(nc, "time", "units")$value
  origin <- as.Date(sub(".*since\\s+", "", tunits))
  list(nc    = nc,
       lon   = ncvar_get(nc, "lon"),
       lat   = ncvar_get(nc, "lat"),
       dates = floor_month(as.Date(ncvar_get(nc, "time"), origin = origin)))
}

#' to_file_lon - helper function to convert between the two common
#' longitude conventions (i.e. -180° to 180° and 0° to 360°)
#' The GISTEMP file uses −180° to 180° longitude convention. 

to_file_lon <- function(x, lon) {
  if (max(lon) > 180) {
    (x + 360) %% 360
  } else {
    ((x + 180) %% 360) - 180
  }
}

#' build_template - builds the predictor table once for the whole window: one row per month with
#' date, nao, nino34 and t (the same at every grid cell). Returns the table and keep = the positions
#' of the window months on the GISTEMP time axis.
#' *) t represents elapsed time in years since the beginning of the analysis period. It
#' allows us to account for a smooth long-term temperature trend (in the GAM model we include "s(t)").

build_template <- function(gdates) {
  idx  <- read_climexp(NAO_FILE)
  nino <- read_climexp(NINO_FILE)
  
  keep  <- which(gdates >= WIN_START & gdates <= WIN_END)
  wdate <- gdates[keep]
  
  tmpl <- data.frame(
    date   = wdate,
    nao    = idx$value[match(wdate, idx$date)], # these two lines provide more robust matching
    nino34 = nino$value[match(wdate, nino$date)]
  )
  
  # tmpl$t <- as.numeric(tmpl$date - min(tmpl$date)) / 365.25 # "+ 0.25" is accounting for leap years on average
  # We have decided to use: 
  tmpl$t <- (seq_len(nrow(tmpl)) - 1) / 12 # later s(t) in the model accounts for long term temperature changes
  
  stopifnot(sum(complete.cases(tmpl)) == N_MONTHS)
  list(tmpl = tmpl, keep = keep)
}


### 3.) Model ------------------------------------------------------------------

# Now we proceed with the construction of the GAM model, we find it sensible to approach
# the task at hand as in the chapter: "6.12.3 Testing a parametric term against a smooth alternative"
# of Generalized Additive Models: An Introduction with R (2nd edition; 2017) by Simon N. Wood.
#
# The method answers a question: "do we need a cubic spline term, or would a straight line be good enough?", i.e.
# it allows us to detect nonlinearity in the data, with the benefit of being able to assess its significance and
# shape.
#
# From Wood (2017, p. 313): "summary(gam(y ~ x + s(x,bs="tp",m=c(2,0)))) will automatically return the p-value
# for testing whether the non-linear components of the cubic thin place [sic] spline are required. m=c(2,0)
# specifies that a spline with a second order penalty and no null space basis is required."
#
# We include here the MODEL for a direct comparability with Skala (2025), which pools the months. Further we fit the data with
# MODEL_SEASON that takes into account the seasonal variability in teleconnections.


MODEL <- temp ~ nao + s(nao, bs = "tp", m = c(2, 0)) +
  nino34 + s(nino34, bs = "tp", m = c(2, 0)) +
  s(t)

# Seasonal model: a separate straight line for each season (season:nao = one NAO slope per
# season), so the smooths test only curvature beyond season-specific straight lines. When the
# linear response differs between seasons, pooling all months bends the all-month average
# even if every season responds linearly (index extremes are mostly winter months); this model
# removes that effect. Needs a column "season" (see season_of below).
MODEL_SEASON <- temp ~ season + season:nao + s(nao, bs = "tp", m = c(2, 0)) +
  season:nino34 + s(nino34, bs = "tp", m = c(2, 0)) +
  s(t)

# MODEL_SEASON with a larger basis for the index smooths (k = 20): checks that the basis size does
# not limit the fitted curves of the main model (s(t) stays at k = 10: k = 20 for the trend worsened
# the AR(1) fit over the ocean).
# -> We conclude that k = 10 is enough: doubling the basis changes the edf by at most 0.34 at the
#    contrast points (pilot, section 5). In the global run (rho by REML) an index edf exceeded 6 in
#    only 8 of the 10 674 cells; there, k = 20 changed it by at most 1.8.
MODEL_SEASON_K20 <- temp ~ season + season:nao + s(nao, bs = "tp", m = c(2, 0), k = 20) +
  season:nino34 + s(nino34, bs = "tp", m = c(2, 0), k = 20) +
  s(t)

# Further we include straight line fits (ie. models without the smooth term) to directly quantify the strength of the teleconnection.
# (Where the curve bends, the straight-line coefficient of MODEL is only the linear part of the split, not the
# slope of the best straight line.)
# Straight lines only
MODEL_LIN <- temp ~ nao + nino34 + s(t)

# One straight line per season (the seasonal counterpart of MODEL_LIN).
MODEL_SEASON_LIN <- temp ~ season + season:nao + season:nino34 + s(t)


# for convenience:
SM_NAO  <- "s(nao)"
SM_NINO <- "s(nino34)"
SM_T    <- "s(t)"


## ------------------------------------------------------------ Helper functions

#' rho1 - lag-1 autocorrelation of a single series (helper function)

rho1 <- function(x) {
  return(as.numeric(acf(x, lag.max = 1, plot = FALSE)$acf[2]))
}

#' ar_whiten - Serves us as a check. Removes the AR(1) carry-over from a residual series:
#' u_t = (e_t - rho * e_(t-1)) / sqrt(1 - rho^2), the first value stays as it is.
#' This is the transformation bam() applies internally when rho is given; if rho was
#' estimated well, the result has lag-1 autocorrelation close to 0.

ar_whiten <- function(e, rho) {
  n <- length(e)
  return(c(e[1], (e[-1] - rho * e[-n]) / sqrt(1 - rho^2)))
}

#' season_of - meteorological season (DJF, MAM, JJA, SON) of each date
#' naming convention: DJF == December, January, February, MAM == March, April, May, etc.

season_of <- function(date) {
  s <- c("DJF", "DJF", "MAM", "MAM", "MAM", "JJA",
         "JJA", "JJA", "SON", "SON", "SON", "DJF")[as.integer(format(date, "%m"))]
  return(factor(s, levels = c("DJF", "MAM", "JJA", "SON")))
}


#' ar1_fit - fits a model with AR(1) residuals; how rho is chosen is set by RHO_METHOD (section 1):
#'   "reml": the rho that minimises the REML score of bam(), as in Wood (2017, sec. 7.7.2: "we can
#'           easily search for a REML optimal rho"). bam() adds the AR(1) term to the score, so the
#'           scores at different rho can be compared. Brent's method (optimize() on [0, 0.95]; near
#'           1 the AR(1) fit becomes numerically unstable) finds the minimum in about eight fits,
#'           and the fit with the lowest score is kept.
#'   "acf" : (our first run) the lag-1 autocorrelation of the residuals of a fit with rho = 0,
#'           clipped to [0, 0.95], then one fit with that rho. Only two fits, but it underestimates
#'           rho where the persistence is strong: the fit with rho = 0 lets the smooths soak up slow
#'           variations, so its residuals look less autocorrelated than they are. Over the oceans
#'           (rho above about 0.6) the whitened residuals kept a positive lag-1 autocorrelation
#'           (median +0.08), and too many cells came out significant (NINO3.4 nonlinearity in the
#'           main model, BH: 1514 cells with "acf", 1112 with "reml").
#' weighted = TRUE: the column w of d gives the weights. Returns list(model, rho), or NULL if the
#' fit fails.

RHO_TOL <- 0.005   # precision of the REML search for rho

ar1_fit <- function(formula, d, weighted = FALSE) {
  fit <- function(r) tryCatch(if (weighted) bam(formula, data = d, method = "fREML", weights = w, rho = r)
                              else          bam(formula, data = d, method = "fREML", rho = r),
                              error = function(e) NULL)
  if (RHO_METHOD == "acf") {
    m0 <- fit(0)
    if (is.null(m0)) return(NULL)
    r <- max(min(rho1(residuals(m0)), 0.95), 0)   # clip to [0, 0.95]: near 1 the AR(1) fit becomes numerically unstable; negative values are physically implausible. Deviance residuals: they include the weights, if any
    m <- fit(r)
    if (is.null(m)) return(NULL)
    return(list(model = m, rho = r))
  }
  best  <- list(score = Inf, model = NULL, rho = NA_real_)
  score <- function(r) {
    m <- fit(r)
    if (is.null(m)) return(Inf)
    if (m$gcv.ubre < best$score) best <<- list(score = m$gcv.ubre, model = m, rho = r)
    return(m$gcv.ubre)
  }
  optimize(score, c(0, 0.95), tol = RHO_TOL)
  if (is.null(best$model)) return(NULL)
  return(list(model = best$model, rho = best$rho))
}

#' fit_ar - fits a model (MODEL by default) to one grid cell with AR(1) residuals, rho chosen as
#' set by RHO_METHOD (ar1_fit). bam() whitens the residuals with rho, so that the smoothing
#' parameters, edf and p-values account for the autocorrelation.
#' Returns list(model, rho, rho_left), or NULL if the fit fails. rho_left is the lag-1
#' autocorrelation left in the whitened residuals: within about +-0.05 (= 2/sqrt(1728)) if the
#' AR(1) correction worked.
#'
#' Why: monthly temperature residuals are autocorrelated (a warm month tends to follow a warm
#' month), so 1728 months carry less independent information than 1728 independent values.
#' Ignoring this makes the p-values too small and lets REML choose too wiggly curves (too high
#' edf). bam(rho = ) builds AR(1) errors into the fit and so corrects the edf as well as the
#' p-values. gam() has no such option; gamm() has, but is far too slow for ~10^4 cells. bam()
#' needs rho to be given, hence ar1_fit.

fit_ar <- function(d, formula = MODEL) {
  R <- ar1_fit(formula, d)
  if (is.null(R)) return(NULL)
  list(model = R$model, rho = R$rho,
       rho_left = rho1(ar_whiten(residuals(R$model, type = "response"), R$rho)))
}


#' fit_season_w - MODEL_SEASON, each month weighted by the inverse of its noise variance,
#' plus AR(1) residuals. The noise variance is modelled as a seasonal level times a slow drift
#' over the decades: the noise differs between seasons and early decades are noisier where observations were sparse. 
#' Noisier months get less weight. Only the noise is modelled; the temperature itself is not rescaled, so all results
#' stay in kelvins.
#' 
#' Two steps:
#' (1) unweighted fit, rho = 0; a gamma GAM of its squared residuals on season + s(t) gives
#'     the noise variance of every month (if that GAM fails: one variance per season);
#' (2) weighted fit(s) with rho chosen as set by RHO_METHOD (ar1_fit).
#' 
#' bam() applies the weights before the AR(1) whitening, so it is the weighted (standardised)
#' residuals that are assumed to follow AR(1). Needs columns "season" and "t".
#' Returns list(model, rho, rho_left, sd_season = residual SD of each season from step (1),
#' drift = noise variance around 1900 relative to 2014, w = the weights (mean 1),
#' noise_sd = modelled noise SD of every month (K)), or NULL if a fit fails.

fit_season_w <- function(d) {
  m0 <- tryCatch(bam(MODEL_SEASON, data = d, method = "fREML", rho = 0), error = function(e) NULL)
  if (is.null(m0)) return(NULL)
  e0   <- residuals(m0, type = "response")
  v    <- tapply(e0, d$season, var)                                  # noise variance of each season
  d$e2 <- e0^2                                                       # squared residuals
  vf   <- tryCatch(gam(e2 ~ season + s(t, k = 5), family = Gamma(link = "log"), data = d, # Gamma distribution is suitable for modelling positive, continuous values.; k=5 binds over ocean -> bit of a hack but we don't want it to fit well the trend we are looking for (only slow changes in observation density, so the weights don’t follow ENSO-driven ups and downs)
                       method = "REML"), error = function(e) NULL)    # noise variance: season x slow drift
  if (is.null(vf)) {
    nv    <- as.numeric(v[as.character(d$season)])                   # fallback: one variance per season
    drift <- NA
  } else {
    nv    <- as.numeric(fitted(vf))                                  # noise variance of every month
    lv    <- as.numeric(predict(vf, data.frame(season = d$season[1], t = c(1900, 2014) - 1880)))
    drift <- exp(lv[1] - lv[2])                                      # t = years since 1880
  }
  d$w <- 1 / nv
  d$w <- d$w / mean(d$w)                                             # normalize weights; mean weight 1
  R <- ar1_fit(MODEL_SEASON, d, weighted = TRUE)                     # rho as set by RHO_METHOD, with the weights
  if (is.null(R)) return(NULL)
  list(model = R$model, rho = R$rho,
       rho_left  = rho1(ar_whiten(sqrt(d$w) * residuals(R$model, type = "response"), R$rho)),
       sd_season = sqrt(v), drift = drift, w = d$w, noise_sd = sqrt(nv))
}

#' season_lines - strength and significance of the dependence. Fits MODEL_SEASON_LIN (one
#' straight line per season, no smooths) with the same weights and AR(1) as the weighted
#' seasonal model (uses the output of fit_season_w i.e. “The seasonal straight lines are fitted with the same weights and
#' AR(1) parameter as the main model, so both describe the noise identically.”).
#' 
#' Returns the slope of each season (DJF, MAM, JJA, SON), Kelvin per index unit, and, for each index, a joint
#' test of "no dependence at all": all four seasonal slopes = 0 (Wald chi-squared test, 4 df,
#' from the fitted covariance matrix).
#' 
#' The nonlinear part is tested separately, by the smooth in MODEL_SEASON.
#' 
#' Returns list(model, slopes = 2 x 4 matrix (rows nao, nino34; columns DJF, MAM, JJA, SON),
#' p_dep = c(nao, nino34)), or NULL if the fit fails.

season_lines <- function(d, W) {
  d$w <- W$w
  m <- tryCatch(bam(MODEL_SEASON_LIN, data = d, method = "fREML", weights = w, rho = W$rho),
                error = function(e) NULL)
  if (is.null(m)) return(NULL)
  b <- coef(m)
  V <- vcov(m) # mgcv's default: the Bayesian covariance matrix (freq = FALSE). For unpenalized coefficients like these slopes it is "rather similar" to the frequentist one, vcov(m, freq = TRUE) (mgcv help, ?summary.gam); at 4 contrast points the standard errors differ by less than 0.1 %
  slopes <- matrix(NA_real_, 2, nlevels(d$season),
                   dimnames = list(c("nao", "nino34"), levels(d$season)))
  p_dep  <- c(nao = NA_real_, nino34 = NA_real_)
  for (ix in rownames(slopes)) {
    j <- paste0("season", levels(d$season), ":", ix)                            # the four seasonal slopes
    slopes[ix, ] <- b[j]
    wald <- drop(b[j] %*% solve(V[j, j], b[j]))                                 # Wald statistic
    p_dep[[ix]] <- pchisq(wald, df = length(j), lower.tail = FALSE)
  }
  list(model = m, slopes = slopes, p_dep = p_dep)
}

### 4.) Load the data ----------------------------------------------------------

G <- open_gistemp()
B <- build_template(G$dates)
tmpl <- B$tmpl
keep <- B$keep

# Check that the data loaded properly
cat(sprintf("grid %d x %d, template %d months, %s to %s\n",
            length(G$lon), length(G$lat), nrow(tmpl),
            format(min(tmpl$date)), format(max(tmpl$date))))

#' cell_series - gets the data for one grid point with given (lon., lat.)

cell_series <- function(target_lat, target_lon) {
  i <- which.min(abs(G$lon - to_file_lon(target_lon, G$lon)))
  j <- which.min(abs(G$lat - target_lat))                                       # latitude has only one convention
  v <- ncvar_get(G$nc, "tempanomaly", start = c(i, j, 1), count = c(1, 1, -1))
  d <- tmpl
  d$temp <- v[keep]
  list(d = d, lon = G$lon[i], lat = G$lat[j], complete = !anyNA(d$temp))        # complete = TRUE if no month is missing
}


### 5.) PILOT: checks before the global run ------------------------------------
# Runs fit_cell() - the same function as the global loop - on 8 contrast points and checks the
# alignment, that every cell is fitted, the basis size, the AR(1) correction and the run time.

## The analysis of one cell (this code is shared by the pilot and the global run)

SEAS     <- levels(season_of(tmpl$date))                                        # "DJF" "MAM" "JJA" "SON"
IDX_MEAN <- c(nao = mean(tmpl$nao), nino34 = mean(tmpl$nino34))                 # all months, 1880-2023;
IDX_SD   <- c(nao = sd(tmpl$nao),   nino34 = sd(tmpl$nino34))                   # the same for every cell
IX_TERM  <- c(nao = "nao", nino = "nino34")                                     # column suffix in the results -> model term
IX_NAME  <- c(nao = "NAO", nino = "NINO3.4")                                    # for printed headings
K20_EDF  <- 6                                                                   # global run: refit with k = 20 only where an index edf of the model exceeds this (by default k = 10)


#' index_rows - the part of a fitted model that belongs to one index (nao, nino) (its straight-line
#' coefficient plus its smooth), evaluated at the index values x with everything else (season,
#' trend, the other index) held fixed. 
#' m: fitted model; term: "nao" or "nino34"; season: for
#' the seasonal models the season whose straight line is used ("DJF", ...), NULL for the pooled
#' ones. Returns list(X = one row per value of x, j = which coefficients).

index_rows <- function(m, term, x, season = NULL) {
  nd  <- m$model[rep(1, length(x)), ]                                           # any data row as a template;
  nd[[term]] <- x                                                               # only the chosen index varies
  lin <- term                                                                   # name of the straight-line coefficient
  if (!is.null(season)) {
    nd$season <- factor(season, levels = levels(m$model$season))
    lin <- paste0("season", season, ":", term)                                  # e.g. "seasonDJF:nao"
  }
  X  <- predict(m, newdata = nd, type = "lpmatrix")                             # type = "lpmatrix" asks mgcv to return the linear-predictor design matrix
  cf <- names(coef(m))
  j  <- cf == lin | startsWith(cf, paste0("s(", term, ")."))                    # this index's coefficients -> i.e. identifies the linear X coefficient and all the X smooth coefficients, while excluding the intercept and the Y coefficients
  return(list(X = X[, j, drop = FALSE], j = j))
}

#' resp4 - the response to a change of the index by 4 SD, read off a fitted model at three
#' index values: mean - 2 SD, mean and mean + 2 SD (IDX_MEAN, IDX_SD):
#'   full : f(mean + 2 SD) - f(mean - 2 SD): strongly positive minus strongly negative phase;
#'          for a straight line = slope x 4 SD
#'   neg  : f(mean) - f(mean - 2 SD): the half from the negative phase to the mean
#'   pos  : f(mean + 2 SD) - f(mean): the half from the mean to the positive phase (neg + pos = full)
#'   bend : pos - neg = f(mean + 2 SD) - 2 f(mean) + f(mean - 2 SD): 0 for a straight line, so the
#'          size of the nonlinearity in K. 
#'          
#' The full difference compares only the two ends, so it misses a curve that bends between them
#' (e.g. both phases colder than the mean); the halves and the bend show it.
#' 
#' Standard errors from the covariance matrix corrected for the uncertainty of the smoothing
#' parameters (by unconditional = TRUE) (Wood 2017, sec. 6.11.1).
#' Where the smooth is shrunk to a straight line (edf ~ 0) the bend
#' is 0 with SE ~ 0; the test of the smooth is then the relevant information.
#' m, term, season: as in index_rows(). Returns a matrix: rows full, neg, pos, bend; columns est
#' and se (K).

resp4 <- function(m, term, season = NULL) {
  R   <- index_rows(m, term, IDX_MEAN[[term]] + c(-2, 0, 2) * IDX_SD[[term]], season) # IDX_MEAN and IDX_SD are defined at the start of this section
  K   <- rbind(full = c(-1,  0, 1),                                             # weights of f(mean - 2 SD), f(mean), f(mean + 2 SD)
               neg  = c(-1,  1, 0),
               pos  = c( 0, -1, 1),
               bend = c( 1, -2, 1))
  cmb <- K %*% R$X                                                              # each row: one combination of the coefficients
  V   <- vcov(m, unconditional = TRUE)[R$j, R$j]
  return(cbind(est = drop(cmb %*% coef(m)[R$j]),
               se  = sqrt(rowSums((cmb %*% V) * cmb))))
}

#' CELL_COLS - the results of one cell, in a fixed order (one row of the global results file):
#'   status            : 0 = fitted, 1 = a fit failed, 2 = cell not complete (not fitted),
#'                       3 = unexpected error (set by the global loop)
#'   pooled MODEL (all months, the set-up for comparison with Skala (2025)):
#'     edf_*, p_*      : edf of the nonlinear part (0 = straight line) and its test
#'     rho, rho_left   : AR(1) parameter (chosen as set by RHO_METHOD) and the lag-1 autocorrelation left after whitening
#'                       (within about +-0.05 -> AR(1) adequate)
#'   straight line through all months (MODEL_LIN):
#'     lin4_*, p_lin_* : slope x 4 SD (K) and its test
#'   main model (MODEL_SEASON, weighted, AR(1)):
#'     edf_*_sw, p_*_sw : nonlinear part and its test - what the nonlinearity maps show
#'     edf_*_k20       : the same edf with k = 20 (NA unless refitted): close to edf_*_sw -> k = 10 suffices
#'     bend_*, bend_se_* : size of the nonlinearity (in K) and its standard error (resp4())
#'     lin4_<S>_*      : straight line of season S x 4 SD (season_lines; K)
#'     curve4_<S>_*    : the 4-SD change read off the curve in season S (K)
#'     p_dep_*         : joint test of "no dependence at all" (Wald) (all four seasonal slopes = 0)
#'     rho_sw, rho_left_sw : AR(1) as above
#'     drift           : noise variance around 1900 relative to 2014 (if NA: then variance per season used)
#'     sd_DJF, sd_JJA  : residual noise level in winter and summer (K)

CELL_COLS <- c("status",
               "edf_nao", "p_nao", "edf_nino", "p_nino", "rho", "rho_left",
               "lin4_nao", "p_lin_nao", "lin4_nino", "p_lin_nino",
               "edf_nao_sw", "p_nao_sw", "edf_nao_k20", "bend_nao", "bend_se_nao",
               "edf_nino_sw", "p_nino_sw", "edf_nino_k20", "bend_nino", "bend_se_nino",
               paste0("lin4_", SEAS, "_nao"), paste0("curve4_", SEAS, "_nao"), "p_dep_nao",
               paste0("lin4_", SEAS, "_nino"), paste0("curve4_", SEAS, "_nino"), "p_dep_nino",
               "rho_sw", "rho_left_sw", "drift", "sd_DJF", "sd_JJA")

#' fit_cell - the whole analysis of one grid cell: the pooled MODEL with AR(1) (fit_ar), the
#' straight lines through all months (MODEL_LIN, same rho), the main model (fit_season_w), its
#' seasonal straight lines (season_lines) and, if needed, the main model with k = 20.
#' d: the cell's data (output of cell_series()$d). k20_always: refit with k = 20 in every cell (the
#' pilot); otherwise only where an index edf of the main model exceeds K20_EDF (the global run).
#' keep_models: also return the fitted models (for plots; never in the global run - memory).
#' Returns list(row = one-row data frame with the columns CELL_COLS, models = list(A, W, L, m_lin)
#' or NULL).

fit_cell <- function(d, k20_always = FALSE, keep_models = FALSE) {
  row <- as.data.frame(as.list(setNames(rep(NA_real_, length(CELL_COLS)), CELL_COLS)))
  row$status <- 2                                           # coverage rule: complete cells only
  if (anyNA(d$temp)) return(list(row = row, models = NULL))
  row$status <- 1
  d$season <- season_of(d$date)
  A  <- fit_ar(d)                                           # pooled MODEL with AR(1)
  W  <- fit_season_w(d)                                     # main model
  L  <- if (is.null(W)) NULL else season_lines(d, W)        # seasonal straight lines, p_dep
  ml <- if (is.null(A)) NULL else tryCatch(bam(MODEL_LIN, data = d, method = "fREML", rho = A$rho),
                                           error = function(e) NULL)    # straight lines, all months pooled
  if (is.null(A) || is.null(W) || is.null(L) || is.null(ml)) return(list(row = row, models = NULL))

  st <- summary(A$model)$s.table                            # pooled
  sw <- summary(W$model)$s.table                            # main model
  pl <- summary(ml)$p.table                                 # straight lines, all months pooled
  sk <- NULL                                                # basis check: main model with k = 20
  if (k20_always || max(sw[c(SM_NAO, SM_NINO), "edf"]) > K20_EDF) {
    d$w <- W$w                                              # same weights and rho as the main model
    sk  <- tryCatch(summary(bam(MODEL_SEASON_K20, data = d, method = "fREML", weights = w,
                                rho = W$rho))$s.table, error = function(e) NULL)
  }

  for (ix in names(IX_TERM)) {
    term <- IX_TERM[[ix]]
    sm   <- if (ix == "nao") SM_NAO else SM_NINO
    sc4  <- 4 * IDX_SD[[term]]                                                  # K per index unit -> K per 4 SD
    row[[paste0("edf_", ix)]]        <- st[sm, "edf"]
    row[[paste0("p_", ix)]]          <- st[sm, "p-value"]
    row[[paste0("lin4_", ix)]]       <- pl[term, "Estimate"] * sc4
    row[[paste0("p_lin_", ix)]]      <- pl[term, "Pr(>|t|)"]
    row[[paste0("edf_", ix, "_sw")]] <- sw[sm, "edf"]
    row[[paste0("p_", ix, "_sw")]]   <- sw[sm, "p-value"]
    if (!is.null(sk)) row[[paste0("edf_", ix, "_k20")]] <- sk[sm, "edf"]
    b <- resp4(W$model, term, season = SEAS[1])["bend", ]   # the same in every season
    row[[paste0("bend_", ix)]]       <- b[["est"]]
    row[[paste0("bend_se_", ix)]]    <- b[["se"]]
    for (sn in SEAS) {
      row[[paste0("lin4_", sn, "_", ix)]]   <- L$slopes[term, sn] * sc4
      row[[paste0("curve4_", sn, "_", ix)]] <- resp4(W$model, term, season = sn)["full", "est"]
    }
    row[[paste0("p_dep_", ix)]]      <- L$p_dep[[term]]
  }
  row$rho         <- A$rho
  row$rho_left    <- A$rho_left
  row$rho_sw      <- W$rho
  row$rho_left_sw <- W$rho_left
  row$drift       <- W$drift
  row$sd_DJF      <- W$sd_season[["DJF"]]
  row$sd_JJA      <- W$sd_season[["JJA"]]
  row$status      <- 0
  return(list(row = row, models = if (keep_models) list(A = A, W = W, L = L, m_lin = ml) else NULL))
}

#' print_table - prints a results table with rounded numbers: p-values (columns whose names
#' start with "p_") to 2 significant digits, all other numbers to k decimals.

print_table <- function(df, k = 2) {
  for (v in names(df)) if (is.numeric(df[[v]]))
    df[[v]] <- if (startsWith(v, "p_")) signif(df[[v]], 2) else round(df[[v]], k)
  print(df, row.names = FALSE)
}

## -------------------------------------------------------- Helper for the pilot

#' ccf_check - cross-correlation between an index and a temperature series. The peak must be
#' at lag 0, otherwise the two series are offset by whole months and everything downstream is
#' wrong (same check as in Skala (2025), fig. 4.7). Returns the lag of the peak (invisibly).

ccf_check <- function(x, y, label) {
  cc   <- ccf(x, y, lag.max = 24, plot = TRUE, na.action = na.pass,
              main = paste("CCF", label, "- peak must be at lag 0"))
  peak <- cc$lag[which.max(abs(cc$acf))]
  cat(sprintf("CCF peak, %s: lag %d months %s\n", label, peak,
              if (peak != 0) "<-- ALIGNMENT PROBLEM" else "OK"))
  invisible(peak)
}


if (RUN_PILOT) {
  # 1. Alignment: one check per index file (the Czechia one is Skala's fig. 4.7)
  CZ  <- cell_series(49, 15)
  BOX <- cell_series(-1, -145)
  ccf_check(CZ$d$nao,     CZ$d$temp,  "NAO vs Czechia")
  ccf_check(BOX$d$nino34, BOX$d$temp, "NINO3.4 vs box")
  
  # 2. The contrast points through the global code (with the k = 20 refit): land and ocean,
  #    strong and weak links, curved and straight responses (chosen by agreement with Skala 2025 and after the conversation with LLM (Claude Opus 5.5) about the spatio-temporal characteristics of studied teleconnections)
  pts  <- list(c(49, 15), c(55, -61), c(65, -21), c(49, -105), c(-1, -145), c(-49, -173),
               c(53, -85), c(7, -175))
  labs <- c("Czechia", "Labrador", "Iceland", "Canadian Prairies", "NINO3.4 box", "SW Pacific",
            "E Canada", "C tropical Pacific")
  t0  <- Sys.time()
  res <- do.call(rbind, Map(function(p, lab)
    data.frame(label = lab, fit_cell(cell_series(p[1], p[2])$d, k20_always = TRUE)$row),
    pts, labs))
  sec <- as.numeric(difftime(Sys.time(), t0, units = "secs")) / length(pts)
  
  # 3. Three compact tables: fits and AR(1); nonlinearity with the k check; strength by season
  print_table(res[, c("label", "status", "rho_sw", "rho_left_sw", "drift")], k = 3)
  print_table(res[, c("label", "edf_nao_sw", "edf_nao_k20", "p_nao_sw", "bend_nao",
                      "edf_nino_sw", "edf_nino_k20", "p_nino_sw", "bend_nino")])
  print_table(res[, c("label", paste0("lin4_", SEAS, "_nao"), "p_dep_nao",
                      paste0("lin4_", SEAS, "_nino"), "p_dep_nino")])
  
  # 4. Go / no-go (the time per cell includes the k = 20 refit, so the estimate is an upper bound)
  n_cells <- 10700                                                              # estimate for rough timing of global run
  n_cores <- max(1, parallel::detectCores() - 1)                                # as N_CORES in section 6 -> We are working on MacBook Air 2021 (M1 chip); leaves one core free, so the Mac stays usable
  dk      <- max(abs(c(res$edf_nao_k20 - res$edf_nao_sw, res$edf_nino_k20 - res$edf_nino_sw)), na.rm = TRUE)
  cat(sprintf("fitted %d of %d | max |edf(k=20) - edf| = %.2f | %.1f s per cell -> about %.1f h for %d cells on %d cores\n",
              sum(res$status == 0), nrow(res), dk, sec, sec * n_cells / 3600 / n_cores, n_cells, n_cores))
}

### 6.) Global loop ------------------------------------------------------------
#
# fit_cell() for every grid cell (90 latitude rows x 180 cells): one row of results per cell in
# OUT_CSV - lat, lon and the columns CELL_COLS (see section 5). Incomplete cells (any month
# missing) are not fitted (as in Skala 2025, sec. 4.1) and get status 2.
#   - The GISTEMP file is read only here, in the main R process, one latitude row at a time; the
#     cells of the row are then fitted in parallel. (mclapply copies the running R process into
#     workers; reading the netCDF file inside the workers can silently return wrong data.)
#   - Each latitude row is appended to OUT_CSV as soon as it is done. With RESUME = TRUE an
#     interrupted run continues where it stopped (rows already in OUT_CSV are skipped).
#   - A cell whose analysis stops with an unexpected error gets status 3; the run goes on.
#   - RUN_INFO records the data version, software versions, settings and timing of the run.
#   - mclapply works on macOS and Linux; on Windows the loop runs on one core.
#
#   On Mac, run from terminal using: caffeinate -i Rscript KFA_project_final_checked.R; (set RUN_PILOT <- FALSE, RUN_GLOBAL <- TRUE)

RESUME      <- TRUE     # FALSE: delete OUT_CSV and start again
GLOBAL_ROWS <- NULL     # NULL = all 90 latitude rows; e.g. 70 (49 N, Czechia's row) to test the loop first
N_CORES     <- if (.Platform$OS.type == "windows") 1 else max(1, parallel::detectCores() - 1)
RUN_INFO    <- sub("\\.csv$", "_info.txt", OUT_CSV)                 # e.g. gam_grid_results_reml_info.txt
OUT_COLS    <- c("lat", "lon", CELL_COLS)

#' read_row - temperatures of one latitude row in the analysis window: matrix of 180 cells
#' (longitudes) x 1728 months.

read_row <- function(j) {
  return(ncvar_get(G$nc, "tempanomaly", start = c(1, j, min(keep)), count = c(-1, 1, length(keep))))
}

#' na_row - an empty result row (all NA) with the given status.

na_row <- function(status) {
  r <- as.data.frame(as.list(setNames(rep(NA_real_, length(CELL_COLS)), CELL_COLS)))
  r$status <- status
  return(r)
}

if (RUN_GLOBAL) {
  
  stopifnot(all(diff(keep) == 1))                 # the window months are consecutive in the file
  
  # -- Which cells are complete: one quick pass over the whole window --------------
  TA       <- ncvar_get(G$nc, "tempanomaly", start = c(1, 1, min(keep)), count = c(-1, -1, length(keep)))
  complete <- apply(TA, c(1, 2), function(z) !anyNA(z))           # longitude x latitude, TRUE = no month missing
  rm(TA); invisible(gc())
  rows <- if (is.null(GLOBAL_ROWS)) seq_along(G$lat) else GLOBAL_ROWS
  
  # -- Resume or start again ------------------------------------------------------
  if (file.exists(OUT_CSV) && RESUME) {
    old <- read.csv(OUT_CSV)
    if (!identical(names(old), OUT_COLS))
      stop("OUT_CSV has other columns (an older version of the analysis?): delete it or set RESUME <- FALSE")
    done <- as.numeric(names(which(table(old$lat) == length(G$lon))))    # latitude rows written in full
    if (nrow(old) != length(done) * length(G$lon))                       # a row cut off mid-write: drop it
      write.csv(old[old$lat %in% done, ], OUT_CSV, row.names = FALSE)
    cat(sprintf("\nresuming: %d of the %d requested latitude rows are already in OUT_CSV\n",
                sum(G$lat[rows] %in% done), length(rows)))
    rows <- rows[!(G$lat[rows] %in% done)]
  } else if (file.exists(OUT_CSV)) {
    file.remove(OUT_CSV)
  }
  
  # -- Record of the run ----------------------------------------------------------
  cat(file = RUN_INFO, append = RESUME && file.exists(RUN_INFO), sep = "",
      sprintf("\n=== run started %s ===\n", format(Sys.time())),
      "GISTEMP: ", ncatt_get(G$nc, 0, "history")$value, "\n",
      "window: ", format(WIN_START), " to ", format(WIN_END), " (", N_MONTHS, " months)\n",
      "indices: ", basename(NAO_FILE), ", ", basename(NINO_FILE), "\n",
      R.version.string, ", mgcv ", as.character(packageVersion("mgcv")), "\n",
      sprintf("complete cells: %d of %d | rows to do now: %d | cores: %d | k = 20 refit where an index edf > %g\n",
              sum(complete), length(complete), length(rows), N_CORES, K20_EDF),
      if (RHO_METHOD == "reml") sprintf("AR(1): rho chosen by REML on [0, 0.95], tolerance %g\n", RHO_TOL)
      else "AR(1): rho = lag-1 autocorrelation of the residuals of a fit with rho = 0\n")

  # -- The loop over latitude rows ------------------------------------------------
  cat(sprintf("\nglobal run: %d latitude rows, %d cells to fit, %d cores\n",
              length(rows), sum(complete[, rows]), N_CORES))
  t0     <- Sys.time()
  n_todo <- sum(complete[, rows])
  n_done <- 0
  for (j in rows) {
    slab  <- read_row(j)                          # read here, in the main process
    cells <- which(complete[, j])                 # the complete cells of this row
    fits  <- parallel::mclapply(cells, function(i) {
      d <- tmpl
      d$temp <- slab[i, ]
      tryCatch(fit_cell(d)$row, error = function(e) na_row(3))
    }, mc.cores = N_CORES)
    
    out <- na_row(2)[rep(1, length(G$lon)), ]     # incomplete cells: status 2
    for (k in seq_along(cells)) {                 # a worker that crashed returns no proper row
      ok <- is.data.frame(fits[[k]]) && identical(names(fits[[k]]), CELL_COLS)
      out[cells[k], ] <- if (ok) fits[[k]] else na_row(3)
    }
    out <- data.frame(lat = G$lat[j], lon = G$lon, out)
    write.table(out, OUT_CSV, sep = ",", row.names = FALSE,
                col.names = !file.exists(OUT_CSV), append = file.exists(OUT_CSV))
    
    n_done <- n_done + length(cells)
    el     <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
    cat(sprintf("lat %6.1f: %3d fitted, %d failed | %6.1f min elapsed, about %.0f min left\n",
                G$lat[j], sum(out$status == 0), sum(out$status %in% c(1, 3)), el,
                if (n_done > 0) el / n_done * (n_todo - n_done) else NA))
  }
  
  # -- Summary --------------------------------------------------------------------
  res_all <- read.csv(OUT_CSV)
  st  <- table(factor(res_all$status, levels = 0:3,
                      labels = c("fitted", "fit failed", "incomplete", "error")))
  fit <- res_all$status == 0
  k20 <- fit & !is.na(res_all$edf_nao_k20)
  msg <- sprintf(paste0("run finished %s, %.2f h | OUT_CSV: %d rows | %s\n",
                        "k = 20 refits: %d, max |edf(k = 20) - edf| = %.2f | ",
                        "rho_left_sw within +-0.05 in %.0f %% of fitted cells\n"),
                 format(Sys.time()), as.numeric(difftime(Sys.time(), t0, units = "hours")),
                 nrow(res_all), paste(names(st), st, collapse = ", "), sum(k20),
                 if (any(k20)) with(res_all[k20, ], max(abs(c(edf_nao_k20 - edf_nao_sw,
                                                              edf_nino_k20 - edf_nino_sw)))) else NA,
                 100 * mean(abs(res_all$rho_left_sw[fit]) <= 0.05))
  cat("\n", msg, sep = "")
  cat(msg, file = RUN_INFO, append = TRUE)
}


### 7.) Results and maps -----------------------------------------------------------
#
# Reads OUT_CSV (one row per cell) and
#   1. corrects the p-values for testing ~10 700 cells at once: Benjamini-Hochberg (BH; Benjamini
#      and Hochberg 1995), separately for each test and index, over the fitted cells. Calling the
#      cells with q < ALPHA significant keeps the expected share of false positives among them (the
#      false discovery rate) at or below ALPHA. Neighbouring cells are positively correlated
#      (GISTEMP's 1200 km smoothing): BH is proven to hold under a form of positive dependence (PRDS;
#      Benjamini and Yekutieli 2001), and Wilks (2016) recommends it for gridded climate fields.
#   2. puts our p-values on Skala's significance scale and compares them with his digitised classes,
#      cell by cell: for the pooled model (his set-up: all months pooled) and for the main model
#      (after the seasonality correction). Skala (2025) did not correct for multiple testing
#      formally (he judges as physically relevant only regions where several neighbouring points
#      are similarly significant, p. 34), so the like-for-like comparison uses uncorrected per-cell
#      values; the BH version is printed next to it.
#   3. prints the summary numbers for the report;
#   4. writes RES_Q_CSV (OUT_CSV plus the new columns) and one 90 x 180 field per mapped quantity to
#      MAPDIR, in Skala's format (rows = GISTEMP latitudes, columns = longitudes, NaN = no result),
#      which gam_maps_final.py draws. OUT_CSV itself is not changed (the resume check of section 6 needs it).

RES_Q_CSV <- sub("\\.csv$", "_q.csv", OUT_CSV)      # e.g. gam_grid_results_reml_q.csv
SKALA_SIG <- 2   # Skala's classes >= 2 sigma counted as significant: his 5 % level (S = 1.645) lies inside class 1, so 2 sigma is the nearest class boundary above it (p <= 0.023)
MIN_NB    <- 4   # coherence: a significant cell with at least 4 of its 8 neighbours significant too; our criterion: Skala (2025, p. 34) asks for "multiple neighbouring points" without giving a number

#' sigma_class - a p-value on Skala's significance scale: the one-sided normal score S = qnorm(1 - p).
#' His score is S = (statistic - surrogate mean) / surrogate SD (Skala 2025, eq. 2.16), with the
#' surrogate distribution assumed Gaussian (p. 17), and his MI test is one-sided: p = 0.05 <-> S = 1.645
#' (p. 33). (His eq. 2.15, as printed, is the two-sided formula p = 2(1 - Phi(S)), which would give
#' S = 1.96; we follow p. 33.) In his classes: 0 = 0-1 sigma, 1 = 1-2, 2 = 2-3, 3 = 3-4,
#' 4 = above 4 sigma (S below 0 -> class 0).
#'
#' Puts our p-values on his scale, so the maps can be compared class by class. The two tests still
#' ask different questions: his MI test against IAAFT surrogates detects any departure from a linear,
#' stationary, Gaussian system; ours only whether the mean response to the index curves.

sigma_class <- function(p) {
  return(pmin(pmax(floor(qnorm(p, lower.tail = FALSE)), 0), 4))
}

#' agree_class - our significance (TRUE/FALSE) next to Skala's, cell by cell:
#' 0 = neither, 1 = GAM only, 2 = Skala only, 3 = both; NA where either has no result. (this "encoding": "ours + 2 * his" allows us to have 4 possible values 0(FF),1(TF),2(FT),3(TT))

agree_class <- function(ours, his) {
  return(ifelse(is.na(ours) | is.na(his), NA, ours + 2 * his))
}

#' agree_stats - the same comparison summed over all cells where both have a result: counts, the
#' share of cells where the two agree, and Cohen's kappa (agreement beyond chance: 0 = no more than
#' chance, 1 = perfect).

agree_stats <- function(ours, his) {
  k  <- !is.na(ours) & !is.na(his)
  o  <- ours[k]; h <- his[k]
  po <- mean(o == h)                                                   # observed agreement
  pe <- mean(o) * mean(h) + (1 - mean(o)) * (1 - mean(h))              # agreement expected by chance
  return(c(cells = length(o), neither = sum(!o & !h), GAM_only = sum(o & !h),   
           Skala_only = sum(!o & h), both = sum(o & h), agree_pct = 100 * po,
           kappa = if (pe < 1) (po - pe) / (1 - pe) else NA))   # NA if neither side flags any cell
}

#' read_skala / skala_cells - Skala's digitised class map (NULL if the file is not there), and its
#' value at each of the cells (lat, lon) (NA where he had no result).

read_skala <- function(path) {
  if (!file.exists(path)) {
    message("Skala's class map not found: ", path, " -> no comparison for this index")
    return(NULL)
  }
  return(as.matrix(read.csv(path, header = FALSE, na.strings = c("nan", "NaN", "NA"))))
}

skala_cells <- function(F, lat, lon) {
  if (is.null(F)) return(rep(NA_real_, length(lat)))
  return(F[cbind(match(lat, G$lat), match(lon, G$lon))])
}

#' as_field - one results column -> 90 x 180 matrix in GISTEMP order (rows = latitudes, columns =
#' longitudes), NA where there is no result.

as_field <- function(res, v) {
  m <- matrix(NA_real_, length(G$lat), length(G$lon))
  m[cbind(match(res$lat, G$lat), match(res$lon, G$lon))] <- v
  return(m)
}

#' coherent - how many of the significant cells (logical 90 x 180 field S) have at least MIN_NB of
#' their 8 neighbours significant too, i.e. lie inside a patch rather than alone. Longitudes wrap
#' around; there are no neighbours beyond the poles.
#' Criterion loosely adapted from Skala (2025, p. 34) (He counts only patches of neighbouring significant cells as physically relevant.)
#' 
#' Caveat: our greater spatial coherence (compared with Skala 2025) need not mean a stronger physical
#' signal. Each cell is fitted on its own (the splines smooth along the index and time, not in space),
#' but the GAM p-value carries no Monte Carlo noise, while his score comes from a noisy MI estimate
#' and only 100 IAAFT surrogates, which scatters his significant cells. GISTEMP's 1200 km smoothing
#' makes neighbouring cells similar for both methods. Our digitising of his figures may also add
#' isolated misread cells to his map.

coherent <- function(S) {
  S[is.na(S)] <- FALSE
  nb <- matrix(0, nrow(S), ncol(S))
  for (dj in -1:1) for (di in -1:1) {
    if (dj == 0 && di == 0) next
    r  <- seq_len(nrow(S)) + dj
    ok <- r >= 1 & r <= nrow(S)
    sh <- matrix(FALSE, nrow(S), ncol(S))
    sh[ok, ] <- S[r[ok], ((seq_len(ncol(S)) + di - 1) %% ncol(S)) + 1]
    nb <- nb + sh
  }
  return(c(cells = sum(S), in_patches = sum(S & nb >= MIN_NB)))
}

if (RUN_RESULTS) {

  if (!file.exists(OUT_CSV)) stop("OUT_CSV not found: run the global loop (section 6) first")
  res <- read.csv(OUT_CSV)
  stopifnot(identical(names(res), OUT_COLS))
  fit <- res$status == 0
  st  <- table(factor(res$status, levels = 0:3, labels = c("fitted", "fit failed", "incomplete", "error")))
  cat(sprintf("\n=== results: %d cells in OUT_CSV | %s ===\n", nrow(res), paste(names(st), st, collapse = ", ")))
  if (nrow(res) < length(G$lat) * length(G$lon))
    cat("NOTE: OUT_CSV does not cover the whole grid (partial run) - the numbers are preliminary\n")

  # -- 1. Multiple testing: BH, separately for each test and index ------------------
  P_COLS <- c("p_nao_sw", "p_nino_sw",       # nonlinearity, main model (the main maps)
              "p_dep_nao", "p_dep_nino",     # any dependence, one straight line per season
              "p_nao", "p_nino",             # nonlinearity, pooled (comparison with Skala)
              "p_lin_nao", "p_lin_nino")     # straight line through all months
  for (v in P_COLS) {
    q <- rep(NA_real_, nrow(res))
    q[fit] <- p.adjust(res[[v]][fit], method = "BH")
    res[[sub("^p_", "q_", v)]] <- q              # e.g. p_nao_sw -> q_nao_sw
  }

  # -- 2. Skala (2025) scale and the comparison with his classes -----------------------------
  for (ix in names(IX_TERM)) {
    his <- skala_cells(read_skala(if (ix == "nao") SKALA_NAO else SKALA_NINO), res$lat, res$lon)
    res[[paste0("skala_", ix)]]        <- his
    res[[paste0("sigma_", ix)]]        <- sigma_class(res[[paste0("p_", ix)]])          # pooled
    res[[paste0("sigma_", ix, "_sw")]] <- sigma_class(res[[paste0("p_", ix, "_sw")]])   # main model
    res[[paste0("agree_", ix)]]        <- agree_class(res[[paste0("sigma_", ix)]] >= SKALA_SIG,
                                                      his >= SKALA_SIG)
    res[[paste0("agree_", ix, "_sw")]] <- agree_class(res[[paste0("sigma_", ix, "_sw")]] >= SKALA_SIG,
                                                      his >= SKALA_SIG)
  }
  write.csv(res, RES_Q_CSV, row.names = FALSE)

  # -- 3. Summary numbers -------------------------------------------------------------
  r     <- res[fit, ]
  share <- function(x) round(100 * mean(x, na.rm = TRUE), 1)
  rowname <- function(m) { rownames(m) <- IX_NAME[rownames(m)]; m }

  cat(sprintf("\n--- nonlinearity: %% of the %d fitted cells where a straight line is rejected ---\n", sum(fit)))
  print(rowname(t(sapply(names(IX_TERM), function(ix) c(
    pooled_raw_p = share(r[[paste0("p_", ix)]] < ALPHA),
    pooled_BH_q  = share(r[[paste0("q_", ix)]] < ALPHA),
    main_raw_p   = share(r[[paste0("p_", ix, "_sw")]] < ALPHA),
    main_BH_q    = share(r[[paste0("q_", ix, "_sw")]] < ALPHA))))))

  cat("\n--- what the seasonality correction changes (cells with BH q < ALPHA) ---\n")
  print(rowname(t(sapply(names(IX_TERM), function(ix) {
    a <- r[[paste0("q_", ix)]] < ALPHA; b <- r[[paste0("q_", ix, "_sw")]] < ALPHA
    c(both = sum(a & b), pooled_only = sum(a & !b), main_only = sum(!a & b))
  }))))

  cat("\n--- size of the nonlinearity where significant (main model, BH q < ALPHA) ---\n")
  print(rowname(t(sapply(names(IX_TERM), function(ix) {
    s <- r[[paste0("q_", ix, "_sw")]] < ALPHA
    b <- abs(r[[paste0("bend_", ix)]][s])
    c(cells = sum(s), median_edf = round(median(r[[paste0("edf_", ix, "_sw")]][s]), 2),
      median_abs_bend_K = round(median(b), 2), p90_abs_bend_K = round(unname(quantile(b, 0.9)), 2))
  }))))

  cat("\n--- dependence (any season, BH q < ALPHA) and strength (median |K per 4 SD| in those cells) ---\n")
  print(rowname(t(sapply(names(IX_TERM), function(ix) {
    s <- r[[paste0("q_dep_", ix)]] < ALPHA
    c(dependent_pct = share(s),
      all_months = round(median(abs(r[[paste0("lin4_", ix)]][s])), 2),
      sapply(setNames(SEAS, SEAS), function(sn) round(median(abs(r[[paste0("lin4_", sn, "_", ix)]][s])), 2)))
  }))))

  cat(sprintf("\n--- comparison with Skala (2025): his classes >= %d sigma vs ours; cells where both have a result ---\n",
              SKALA_SIG))
  cmp <- do.call(rbind, lapply(names(IX_TERM), function(ix) {
    his  <- r[[paste0("skala_", ix)]] >= SKALA_SIG
    ours <- list("pooled, >= 2 sigma" = r[[paste0("sigma_", ix)]] >= SKALA_SIG,
                 "pooled, BH q"       = r[[paste0("q_", ix)]] < ALPHA,
                 "main, >= 2 sigma"   = r[[paste0("sigma_", ix, "_sw")]] >= SKALA_SIG,
                 "main, BH q"         = r[[paste0("q_", ix, "_sw")]] < ALPHA)
    data.frame(index = IX_NAME[[ix]], ours = names(ours),
               t(sapply(ours, agree_stats, his = his)), row.names = NULL)
  }))
  print_table(cmp)

  # Coherence: do the significant cells form patches (Skala 2025, p. 34)? Our maps both uncorrected
  # (>= 2 sigma, like-for-like with Skala's map) and BH-corrected (our main results).
  cat(sprintf("\n--- coherence: significant cells inside patches (>= %d of 8 neighbours significant) ---\n", MIN_NB))
  print(do.call(rbind, lapply(names(IX_TERM), function(ix) {
    S <- list("pooled, >= 2 sigma" = as_field(res, res[[paste0("sigma_", ix)]] >= SKALA_SIG),
              "pooled, BH q"       = as_field(res, res[[paste0("q_", ix)]] < ALPHA),
              "main, >= 2 sigma"   = as_field(res, res[[paste0("sigma_", ix, "_sw")]] >= SKALA_SIG),
              "main, BH q"         = as_field(res, res[[paste0("q_", ix, "_sw")]] < ALPHA),
              "Skala, >= 2 sigma"  = as_field(res, res[[paste0("skala_", ix)]] >= SKALA_SIG))
    st <- t(sapply(S, coherent))                                   # cells, in_patches per map
    data.frame(index = IX_NAME[[ix]], map = names(S), st,
               pct_in_patches = round(100 * st[, "in_patches"] / pmax(st[, "cells"], 1)),
               row.names = NULL)
  })), row.names = FALSE)

  # AR(1) check: share of fitted cells whose whitened residuals keep at most +-0.05 lag-1 autocorrelation
  cat(sprintf("\nAR(1): rho_left within +-0.05 in %.1f %% of fitted cells for the main model (median rho %.2f), %.1f %% for the pooled model (median rho %.2f)\n",
              share(abs(r$rho_left_sw) <= 0.05), median(r$rho_sw),
              share(abs(r$rho_left) <= 0.05), median(r$rho)))

  # Where is the nonlinearity: in cells where the AR(1) correction works, or where it fails?
  # (This check is why we choose rho by REML. In our first run, rho from the residual ACF, AR(1) was
  # adequate in only 47 % of the cells, the leftover autocorrelation was mostly positive (median
  # +0.055, i.e. rho too low, so p-values too small), and NINO3.4 nonlinearity was found in 21.8 % of
  # the cells where AR(1) failed against 5.5 % where it worked. The NAO was not concentrated in the
  # problem cells (9.3 % against 12.5 %): it is close to random from month to month, so leftover
  # autocorrelation barely affects its p-values. With rho by REML: 64 %, median +0.020, 18.8 % against 5.7 %.)
  ok <- factor(abs(r$rho_left_sw) <= 0.05, levels = c(FALSE, TRUE), labels = c("AR(1) fails", "AR(1) adequate"))
  cat("\n--- nonlinearity where the AR(1) check fails or passes (main model) ---\n")
  print(rbind("NINO3.4: % of cells with BH q < ALPHA" = round(100 * tapply(r$q_nino_sw < ALPHA, ok, mean), 1),
              "NAO: % of cells with raw p < ALPHA"    = round(100 * tapply(r$p_nao_sw  < ALPHA, ok, mean), 1)))
  cat("leftover lag-1 autocorrelation (main model):\n"); print(summary(r$rho_left_sw))

  # -- 4. The fields for the maps (gam_maps_final.py) -----------------------------------
  # Headerless 90 x 180 CSVs in Skala's format, NaN = no result. Per index (nao, nino):
  #   sigma_*, sigma_*_sw   : nonlinearity on Skala's scale, pooled and main model (classes 0-4)
  #   agree_*, agree_*_sw   : agreement with Skala (0 neither, 1 GAM only, 2 Skala only, 3 both)
  #   edf_*, q_*            : nonlinear edf and its BH q, pooled model
  #   edf_*_sw, q_*_sw      : the same, main model
  #   bend_*                : size of the nonlinearity in K, main model
  #   lin4_*, q_lin_*       : straight line through all months (K per 4 SD) and its BH q
  #   lin4_<S>_*, q_dep_*   : straight line of each season (K per 4 SD) and the BH q of the joint test
  # and for the whole grid: status, rho_left_sw (the AR(1) check).
  dir.create(MAPDIR, showWarnings = FALSE)
  FIELDS <- c("status", "rho_left_sw")
  for (ix in names(IX_TERM))
    FIELDS <- c(FIELDS, paste0(c("sigma_", "sigma_", "agree_", "agree_", "edf_", "q_", "edf_", "q_",
                                 "bend_", "lin4_", "q_lin_", "q_dep_"), ix,
                               c("", "_sw", "", "_sw", "", "", "_sw", "_sw", "", "", "", "")),
                paste0("lin4_", SEAS, "_", ix))
  for (f in FIELDS)
    write.table(as_field(res, res[[f]]), file.path(MAPDIR, paste0("gam_", f, ".csv")),
                sep = ",", row.names = FALSE, col.names = FALSE, na = "NaN")
  cat(sprintf("\n%d fields written to %s; all results with the new columns in %s\n",
              length(FIELDS), MAPDIR, RES_Q_CSV))
}

### End ------------------------------------------------------------------------
# Keep this as the very last line of the script (after the results section):
nc_close(G$nc)
