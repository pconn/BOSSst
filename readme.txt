A number of steps are needed to get data together for a spatio-temporal analysis of BOSS 
count data.

These include:
-Assemble analysis grid and spatial covariates using the "bass" package.  The most recent 
version of the R code, "create_Bering_BOSS_grids_2012_2013.R" produces a grid for 
2012-2013 for which cells with any ice in 2012-2013 are included.  For the eastern 
Bering, running this code results in "AlaskaBeringData2012_2013_14Dec2015.Rdat"

-Effort and hot spot data were previously processed via the R script 
STabundance/inst/pep_geo_boss_hotspots.R, which resulted in the file 
"BOSS_2012Effort_22Apr14.Rdata."  To produce 2012 and 2013 effort and hotspots datasets 
from finalized BOSS data, a combination of database queries and R programming was used.  
From the database, a combination of files are produced (see ./data_from_JML/readme_boss_data_products).
To formulate effort, we use a SpatialPointsDataFrame (in "boss_geo_sp.rda") that has a record
for each "on effort" photograph.  We then 1) use a state space model to improve estimates of the
aircraft's location at each photo, 2) employ a filter to remove effort during turns (where roll is problematic),
and 3) calculate a table of area surveyed by georeferencing remaining photographs and calculating the 
cumulative area covered (thereby accounting for overlapping photographs).  This is done via the script
/BOSSst/inst/Calculate_area_surveyed_crawl.R.  

-Next, we format remaining data products and convert them into a form suitable for application of 
MCMC estimation functions.  For 2012, this is done with /BOSSst/inst/format_Bering2012_all.R.  Requisite
data for hierarchical modeling is stored in "BOSS_data_2012.Rda"

-Finally, we run models, using e.g. run_Bering2012_all.R
