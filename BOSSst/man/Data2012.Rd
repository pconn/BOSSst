\docType{data}
\name{Data2012}
\alias{Data2012}
\title{Data to perform spatio-temporal abundance estimation for 2012 surveys in the eastern Bering Sea}
\description{
  The list \code{Data2012} includes 15 elements:
  \itemize{
    \item{\code{Area_hab}} A vector giving the relative area of effective habitat in each modeled grid cell.  The maximum value is 1.0, which indicates the entire cell is composed of sea water/ice.  Cells that include land have values less than 1.0.
    \item{\code{Area_trans}} A vector giving the proportion of area sampled in each surveyed grid cell. Note that this dataset includes 100 'pseudo-zeroes' - randomly selected grid cells with an ice concentration of zero, where counts were set to zero.  \code{Area_trans} was set to 0.9 for 'pseudo-zeroes'.
    \item{\code{Beta_psi}} Species misclassification estimates on the multinomial logit scale (from McClintock multiple observer study)
    \item{\code{Coords}} Matrix giving projected xy-coordinates of each grid cell. Not needed for this analysis, but might be useful for other analyses.  The projection used is "+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0"
    \item{\code{Dat}} Data frame consisting of 3324 rows (one for each `on effort' thermal detection), and 4 columns, corresponding to  
    ``Transect", ``Photo", ``Obs", and ``Group".  The ``Transect" column indexes
    the day and cell that was surveyed (see the `Mapping' element for a link to `Transect' to grid cell and day ).  The ``Photo" column is an artifact from a previous analysis and is just a vector of 1s (all animals in the present dataset were photographed).  ``Obs" gives the observation type (1 to 13, indexing species and certainty type - `spotted - guess', `spotted - likely', `spotted-certain', `ribbon-guess', `ribbon-likely', `ribbon-certain', `bearded-guess', `bearded-likely', `bearded-certain', `ringed-guess',`ringed-likely',`ringed-certain', `unknown seal species'), and ``Group" gives the number of seals that made up the thermal hotspot (usually 1, but sometimes 2-3, most often mother-pup pairs)
    \item{\code{DayHour}} A 2 column matrix holding the day (first column; day = 1 corresponds to 10 Apr 2012, and hour (2nd column; solar hour) for each surveyed grid cell.  These pseudo-zeroes have hour set to 1.0.
    \item{\code{Det_priors}} A list giving an estimate of total detection probability (factoring in availability<1, incomplete detection of sensors, and escape behavior) and an associated variance-covariance matrix for each species (bd - bearded, sd - spotted, rn - ribbon, rd - ringed) and in each surveyed cell.  A VC matrix isn't provided for ringed seals since we're missing data on many of the components and they aren't the focus of this paper.
    \item{\code{Hab_cov}} A data frame giving covariate values (columns) for each
    cell and each day surveyed.  For instance, the first 1331 rows give habitat values for the first day of the survey, rows 1332-2662 give values for the second day, etc.
    \item{\code{K}} This is a spatial weights matrix for spatial modeling, but is not used in the present paper
    \item{\code{Knot_loc}}  This gives the location of spatial knots potentially used in spatial models.  Not used in the present paper.
    \item{\code{Mapping}}  A two column matrix giving the spatial cell surveyed (first column) and the day surveyed for each surveyed cell.  Note the last 100 records are pseudo-zeroes.
    \item{\code{Prop_photo}} Proportion of surveyed cells that were photographed.  This is a relic of a previous analysis approach and is set to 1.0 here since we calculate effort as the area of the photographed footprint (no unphotographed hot spots are modeled)
    \item{\code{Psi}} A 3 dimensional array holding MCMC draws from McClintock et al.'s species misidentification analysis.  The first dimension is true species, the second dimension is observation type, and the third indexes posterior sample (for 10000 MCMC samples)
    \item{\code{Surveyed}} A vector giving the indexes of \code{Area_trans} that were actually surveyed (i.e. omitting pseudo-zeroes)
    \item{\code{VC_psi}} An estimated variance covariance matrix for species misclassification parameters in multinomial logit space. 
  }
}
\author{
  Paul Conn \email{paul.conn@noaa.gov}
}
\keyword{data}

