####################
# Standardize Philadelphia PD data
####################

####
# SETUP PACKAGES
####

require(optparse, quietly=T)
require(lubridate, quietly=T)
require(stringr, quietly=T)
require(sp, quietly=T)
require(rgdal, quietly=T)


####
# READ COMMAND LINE PARAMS
####

option_list = list(
  make_option(c("--verbose"), type="integer", default=2,
              help="Adjust level of status messages [default %default]",
              metavar="number"),
  make_option(c("--toWebMercator"), action="store_true", default=FALSE,
              help="Reproject points to web mercator [default %default]"),
  make_option(c("--reprojectTo"), type="character", default="none", 
              help="Reproject points to match a shapefile [default %default]")
)

# parse the options
parser = OptionParser(usage = "%prog [options] csvfile outputdir", option_list=option_list)
arguments = parse_args(parser, positional_arguments=TRUE)

# extract the input file name or quit
if(length(arguments$args) == 2) {
  kInputCSV = arguments$args[1]
  kOutputDir = arguments$args[2]
} else {
  stop("Run with --help for assistance.")
}
# extract the remaining options
opt = arguments$options


####
# SETUP PARAMETERS
####

kVerbose = opt$verbose
kToWebMercator = opt$toWebMercator
kReprojectTo = opt$reprojectTo

####
# SETUP/CLEANUP ENVIRONMENT
####

# clean up
remove(arguments)
remove(opt)

####
# FUNCTIONS
####

PrintStatus = function (visible.at, ...) {
  # Outputs status messages for logging purposes
  # 
  # Args:
  #   visible.at: works with kVerbose constant to determine what to print
  #   ...: any number of other variables
  if(visible.at <= kVerbose) print(paste0("Status --  ", ...))
}



####
# MAIN
####

# read the data file into a data frame
PrintStatus(1, "Reading csv file ", kInputCSV)
crimes.raw = read.csv(file.path(kInputCSV), stringsAsFactors=FALSE, comment.char="")

if(kVerbose > 1) {
  PrintStatus(2, "Raw data information:")
  head(crimes.raw)
  str(crimes.raw)
  summary(crimes.raw)
}

# standardize the data

# setup a new variable to store cleaned data
crimes.cleaned = data.frame(id=crimes.raw$DC_KEY, stringsAsFactors=FALSE)

# parse datetimes
crimes.cleaned$datetime = parse_date_time(crimes.raw$DISPATCH_DATE_TIME, orders=c("ymdhms"))

# set event classes
crimes.cleaned$class = crimes.raw$TEXT_GENERAL_CODE

# set xy coordinates
crimes.cleaned$pointx = crimes.raw$POINT_X
crimes.cleaned$pointy = crimes.raw$POINT_Y

if(kVerbose > 1) {
  PrintStatus(2, "Cleaned data information:")
  head(crimes.cleaned)
  str(crimes.cleaned)
  summary(crimes.cleaned)
}

# omit rows with NAs
crimes.cleaned.complete = na.omit(crimes.cleaned)

# calculate percent dropped
crimes.percentdropped = 100 * (nrow(crimes.cleaned)-nrow(crimes.cleaned.complete)) / nrow(crimes.cleaned)

PrintStatus(1, "Dropped ", format(crimes.percentdropped, digits=1), "% of events with incomplete data.")

if(kVerbose > 1) {
  PrintStatus(2, "Cleaned data, complete cases information:")
  summary(crimes.cleaned.complete)
}

####
# REPROJECT
####

if(kToWebMercator == TRUE | kReprojectTo != "none") {
  if(kReprojectTo != "none") {
    PrintStatus(1, "Projecting to match shapefile...")
    target.shapefile = readOGR(dirname(kReprojectTo),gsub("\\.shp", "", basename(kReprojectTo)))
    reprojection.crs = CRS(proj4string(target.shapefile))
  } else {
    PrintStatus(1, "Projecting to web mercator...")
    # project to spherical mercator EPSG:3857 http://spatialreference.org/ref/sr-org/6864/
    reprojection.crs = CRS("+proj=merc +lon_0=0 +lat_ts=0.0 +k=1 +x_0=0 +y_0=0 +a=6378137 +b=6378137 +units=m +wktext +nadgrids=@null +no_defs ")
  }
  
  # fetch just the locations
  points.xy = crimes.cleaned.complete[,c("pointx","pointy")]
  
  # create a point set
  points.sp.lonlat = SpatialPoints(points.xy, proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs"))
  if(kVerbose > 1) {
    summary(points.sp.lonlat)
  }
  
  # reproject
  points.sp.webmercator = spTransform(points.sp.lonlat, CRS=reprojection.crs)
  if(kVerbose > 1) {
    summary(points.sp.webmercator)
  }
  
  # overwrite x and y
  crimes.cleaned.complete$pointx = as.data.frame(points.sp.webmercator)$pointx
  crimes.cleaned.complete$pointy = as.data.frame(points.sp.webmercator)$pointy
  
  
  if(kVerbose > 1) {
    PrintStatus(2, "Sample projected rows:")
    crimes.cleaned.complete[sample(nrow(crimes.cleaned.complete),10),]
  }
}
####
# SAVE FILES
####

# save one csv of all incidents
PrintStatus(1, "Writing all crimes to CSV...")
write.csv(crimes.cleaned.complete, file.path(kOutputDir,"alltime-allcrimes.csv"))

# extract unique class values
classes = unique(str_trim(crimes.cleaned.complete$class))
if(kVerbose > 1) {
  PrintStatus(2, "Classes found:")
  print(classes)
}

# save one csv for each class value
for(c in classes) {
  # filename
  c.filename = paste0("alltime-",gsub(" |-","",tolower(c)),".csv")
  PrintStatus(1, "Writing ", c.filename)
  write.csv(crimes.cleaned.complete[str_trim(crimes.cleaned.complete$class)==c,], file.path(kOutputDir,c.filename),
            row.names=FALSE)
}

####
# EXIT
####
PrintStatus(1, "Complete")