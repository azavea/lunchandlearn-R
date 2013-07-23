####################
# Visualize Philadelphia PD data
####################

####
# SETUP PACKAGES
####

require(optparse, quietly=T)
require(lubridate, quietly=T)
require(stringr, quietly=T)
require(sp, quietly=T)
require(spatstat, quietly=T)
require(rgdal, quietly=T)
require(raster, quietly=T)
require(RColorBrewer, quietly=T)
require(maptools)


####
# READ COMMAND LINE PARAMS
####

option_list = list(
  make_option(c("--verbose"), type="integer", default=2,
              help="Adjust level of status messages [default %default]",
              metavar="number"),
  make_option(c("--distance"), action="store_true", default=FALSE,
              help="Generate nearest crime raster [default %default]"),
  make_option(c("--density"), action="store_true", default=FALSE,
              help="Generate density raster [default %default]"),
  make_option(c("--aggregate"), type="character", default="none",
              help="Aggregate to polygons [default %default]"),
  make_option(c("--cellsize"), type="integer", default=250,
              help="Set raster cell size [default %default]",
              metavar="number")
)

# parse the options
parser = OptionParser(usage = "%prog [options] boundaryshp csvfile outputdir", option_list=option_list)
arguments = parse_args(parser, positional_arguments=TRUE)

# extract the input file name or quit
if(length(arguments$args) == 3) {
  kBoundaryFile = arguments$args[1]
  kInputCSV = arguments$args[2]
  kOutputDir = arguments$args[3]
} else {
  stop("Run with --help for assistance.")
}
# extract the remaining options
opt = arguments$options


####
# SETUP PARAMETERS
####

kVerbose = opt$verbose

kDistance = opt$distance
kDensity = opt$density
kAggregate = opt$aggregate

kCellSize = opt$cellsize

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

SaveRaster = function (raster, label, inputfile, outputdir) {
  inputfile.simple = gsub("\\.csv", "", basename(inputfile))
  output.file = file.path(outputdir, paste0(inputfile.simple, "-", label))
  
  # write geotiff
  writeRaster(raster, output.file, format="GTiff", overwrite=TRUE)
  
  # write KML
  # first project back to longlat
  raster.longlat =projectRaster(raster, crs="+proj=longlat +datum=WGS84", method='ngb')
  KML(raster.longlat, output.file, overwrite=TRUE, col=adjustcolor(brewer.pal(9, "YlOrRd"), alpha.f=0.5), blur=4)
}
  

SaveSHP = function (polygons, label, inputfile, outputdir) {
  inputfile.simple = gsub("\\.csv", "", basename(inputfile))
  output.file = file.path(outputdir, paste0(inputfile.simple, "-", label))
  
  # write shapefile
  writeOGR(polygons, dirname(output.file), basename(output.file), driver="ESRI Shapefile", overwrite_layer=TRUE)
  
}
  

####
# MAIN
####

####
# Define raster
####

PrintStatus(1, "Reading boundary shapefile...")
boundary = readOGR(dirname(kBoundaryFile),gsub("\\.shp", "", basename(kBoundaryFile)))
boundary.bbox = bbox(boundary)
boundary.width = boundary.bbox[1,2]-boundary.bbox[1,1]
boundary.height = boundary.bbox[2,2]-boundary.bbox[2,1]

PrintStatus(1, "Summary of boundary:")
summary(boundary)

PrintStatus(1, "Generate raster covering...")
# calculate raster width
raster.cols = ceiling(boundary.width / kCellSize)
raster.width = raster.cols * kCellSize
raster.rows = ceiling(boundary.height / kCellSize)
raster.height = raster.rows * kCellSize
raster.extent = extent(c(xmin=boundary.bbox[1,1],xmax=boundary.bbox[1,1]+raster.width,
                         ymin=boundary.bbox[2,1],ymax=boundary.bbox[2,1]+raster.height))
# make raster
raster.overall = raster(extent(raster.extent), 
                        nrow=raster.rows, ncol=raster.cols,
                        crs=CRS(proj4string(boundary)))

if(kVerbose > 1) {
  PrintStatus("Summary of raster covering:")
  print(raster.overall)
}

####
# Read in points
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

PrintStatus(1, "Building spatial point data...")

# fetch just the locations
points.xy = crimes.raw[,c("pointx","pointy")]

# create a point set
points.sp = SpatialPoints(points.xy, proj4string=CRS(proj4string(boundary)))
if(kVerbose > 1) {
  summary(points.sp)
}

# associate each point with the polygon it falls within
points.sp.polygonassignments = over(points.sp, as(boundary, "SpatialPolygons"))

# remove points that don't intersect a polygon
points.sp = points.sp[!is.na(points.sp.polygonassignments)]

# rasterize the points
points.raster = rasterize(points.sp, raster.overall, fun="count")

if(kVerbose > 1) {
  PrintStatus(2, "Raster count summary:")
  summary(getValues(points.raster))
}

####
# DISTANCE
####

if(kDistance == TRUE) {
  PrintStatus(1, "Calculating distances...")
  
  if(sum(is.na(values(points.raster))) > 0) {
    # raster contains at least 1 NA cell
    points.raster.distance = distance(points.raster)
    
    SaveRaster(points.raster.distance, "distance", kInputCSV, kOutputDir)

  } else {
    PrintStatus(1, "No NA cells present, skipping distance calculation.")
  }
}  

####
# DENSITY
####

if(kDensity == TRUE) {
  PrintStatus(1, "Calculating density...")
  
  points.ppp = ppp(coordinates(points.sp)[,1], coordinates(points.sp)[,2],
                   window=owin(xrange=bbox(raster.overall)[1,], yrange=bbox(raster.overall)[2,]))
  
  points.raster.density = raster(density(points.ppp, dimyx=c(nrow(raster.overall), ncol(raster.overall)), sigma=1000))
  projection(points.raster.density) = projection(raster.overall)
  
  SaveRaster(points.raster.density, "density", kInputCSV, kOutputDir)
}

####
# AGGREGATE
####

if(kAggregate != "none") {
  PrintStatus(1, "Calculating aggregation...")
  
  PrintStatus(1, "Reading aggregation shapefile...")
  aggregate.polygons = readOGR(dirname(kAggregate),gsub("\\.shp", "", basename(kAggregate)))
  
  PrintStatus(1, "Projecting to our analysis projection...")
  aggregate.polygons = spTransform(aggregate.polygons, CRS(proj4string(raster.overall)))
  
  PrintStatus(1, "Count points in each polygon...")
  polygon.counts = over(as(aggregate.polygons, "SpatialPolygons"), points.sp, fn=count)
  
  aggregate.polygons$counts = polygon.counts  
  
  PrintStatus(1, "Save as shapefile...")
  SaveSHP(aggregate.polygons, "counts", kInputCSV, kOutputDir)
  
}
  
####
# EXIT
####
PrintStatus(1, "Complete")