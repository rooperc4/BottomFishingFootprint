#' Function to produce summarized fishing effort data by gear type and year group while accounting for confidentiality concerns
#'
#' This function takes a raw fishing effort data frame and summarizes it for mapping, creating a footprint_data_object. The raw fishing effort data should
#' include rows for each individual fishing event (set), such as a trawl haul or longline set. For each set the year, latitude, longitude, gear type and a
#' vessel identifier must be included. A minimum number of vessels can be specified so that the output does not include grid cells where less than this
#' number of vessels fished a given gear type in a given year group. The default value is 3 (consistent with Canadian requirements for confidentiality),
#' The starting year of interest (usually the earliest year in the data) and how many years should be combined should be specified. Additionally,
#' the spatial resolution that is required and the units of that spatial resolution (defaults are 0.5 decimal degrees of latitude and longitude) need to
#' be specified. A map projection (in proj4 format) should also be given.
#'
#' The function outputs a list of eight data objects that include a data frame of raw fishing effort data (fishing_effort_data) summarized to the spatial
#' resolution #' provided for each fishing event with an additional column for the year group to which the event belongs. There should be one record in this
#' table for each of the fishing events that were input. This table may not be appropriate to share with the NPFC, as it will contain records for grid cells
#' where fewer than the minimum number of vessels fished. The second data object output by the function (fishing_footprint_table) is a table that summarizes
#' the fishing events by year group, latitude, longitude and gear type. This table produces the number of sets that occurred in each year by each gear type
#' at a grid cell with the spatial resolution specified and the center position given by the latitude and longitude where it was grouped. In addition, a column
#' specifies the number of vessels that are included at each grid cell. The table will not include grid cells where there were fewer than the minimum number of
#' vessels fishing. This table may be appropriate for sharing with the wider group, as it will meet confidentiality rules. The third item is a table of the number
#' of grid cells by gear type and year group where the minimum number of vessels was not met (not_enough_vessels). This table indicates how many grid cells for which data exists, but
#' cannot be reported. The fourth item is the year groups identified from the inputs (years) and the fifth item is the gear types represented in the data (gear_
#' types). The last three items are the spatial resolution and units of the data summary and the map projection of the data (map_projection). These last items
#' are for tracking inputs and outputs and are needed for creating the figures using the footprint_map function.
#'
#' @param Year The year in which the fishing event occurred
#' @param Longitude The longitude of the fishing event in decimal degrees
#' @param Latitude The latitude of the fishing event in decimal degrees
#' @param Gear_type The gear type of the fishing event (e.g. trawl, trap, longline, gillnet)
#' @param Vessel_ID The vessel identifier for the set
#' @param Minimum_vessels The minimum number of vessels to fish in in a grid cell in order to satisfy confidentiality rules (default = 3)
#' @param Start_year The year to start summmarizing data
#' @param Years_to_combine The number of years to combine in the summary of data
#' @param Spatial_resolution The spatial resolution of the output grid (default = 0.5). If using "ISEA_12" or "ISEA_13" this value should be NA.
#' @param Units The units of the spatial resolution, can be "dd" for decimal degrees, "m" for projected data, "ISEA_12" or "ISEA_13"
#' @param map_projection The map projection should be specified in proj4 format (e.g. "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs" for decimal degrees or "+proj=utm +zone=59 +ellps=WGS84 +datum=WGS84 +units=m +no_defs" for UTM
#'
#'
#' @keywords NPFC, fishing effort, fishing footprint
#' @export
#' @examples
#' Spatial_resolution_dd<-data.frame(Spatial_resolution=.5,Units="dd")
#' Start_year<-1991
#' Years_to_combine<-5
#' Minimum_vessels<-3
#' map_projection_dd<-"+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
#' example_dd<-footprint_data(effort_data$Year,effort_data$Longitude,effort_data$Latitude,effort_data$Gear_type,effort_data$Vessel,Minimum_vessels,Start_year,
#' Years_to_combine,Spatial_resolution_dd$Spatial_resolution,Spatial_resolution_dd$Units, map_projection=map_projection_dd)


footprint_data<-function(Year,Longitude,Latitude,Gear_type,Vessel_ID,Minimum_vessels=3,Start_year,Years_to_combine,Spatial_resolution=.5,Units="dd", map_projection=NULL){


  if(Units=="dd"){
    LONGITUDE<-round_any(Longitude,Spatial_resolution)
    LATITUDE<-round_any(Latitude,Spatial_resolution)}

  if(Units=="m"){
    sed1<-cbind(Longitude,Latitude) #Create the position data (lat-long pairs)
    seddata.project <- project(sed1, map_projection) #Project the position data to AEA
    LONGITUDE<-round_any(seddata.project[,1],Spatial_resolution)
    LATITUDE<-round_any(seddata.project[,2],Spatial_resolution)}

  if(Units=="ISEA_12"){
    dggs<- dgconstruct(projection="ISEA",aperture=3,res=12)
    dggs_effort_cell<-dgGEO_to_SEQNUM(dggs,Longitude,Latitude)$seqnum
    cellcenters<- dgSEQNUM_to_GEO(dggs,dggs_effort_cell)
    LONGITUDE<-cellcenters$lon_deg
    LATITUDE<-cellcenters$lat_deg
    #Spatial_resolution<-min_dist(Longitude,Latitude)
    Spatial_resolution<-.08}

  if(Units=="ISEA_13"){
    dggs<- dgconstruct(projection="ISEA",aperture=3,res=13)
    dggs_effort_cell<-dgGEO_to_SEQNUM(dggs,Longitude,Latitude)$seqnum
    cellcenters<- dgSEQNUM_to_GEO(dggs,dggs_effort_cell)
    LONGITUDE<-cellcenters$lon_deg
    LATITUDE<-cellcenters$lat_deg
   # Spatial_resolution<-min_dist(Longitude,Latitude)
    Spatial_resolution<-0.05}


  End_year<-round_any(max(Year,na.rm=TRUE),Years_to_combine,ceiling)
  label_years<-paste(seq(Start_year,End_year,Years_to_combine),seq(Start_year+(Years_to_combine-1),End_year,Years_to_combine),sep="-")

  fishing_effort<-data.frame(cbind(YEAR=Year,GEAR_TYPE=Gear_type,VESSEL=Vessel_ID,LONGITUDE=LONGITUDE,LATITUDE=LATITUDE))

  fishing_effort$YEAR_GROUP<-cut(as.numeric(as.character(fishing_effort$YEAR)),labels=label_years,breaks=seq(Start_year,End_year+Years_to_combine,Years_to_combine),right=FALSE)
  d1<-aggregate(as.numeric(fishing_effort$VESSEL),list(fishing_effort$YEAR_GROUP,fishing_effort$LONGITUDE,fishing_effort$LATITUDE,fishing_effort$GEAR_TYPE),FUN="unique")
  d1$NumberOfSets<-aggregate(as.numeric(fishing_effort$VESSEL),list(fishing_effort$YEAR_GROUP,fishing_effort$LONGITUDE,fishing_effort$LATITUDE,fishing_effort$GEAR_TYPE),FUN="length")$x
  d1$NumberOfVessels<-lengths(d1$x)
  fishing_footprint_table<-data.frame(YEARS=d1$Group.1,LONGITUDE=d1$Group.2,LATITUDE=d1$Group.3,GEAR_TYPE=d1$Group.4,NumberOfSets=d1$NumberOfSets,NumberOfVessels=d1$NumberOfVessels)
  fishing_footprint_table$YEARS<-factor(fishing_footprint_table$YEARS)
  did_not_meet_minimum_vessels<-table(fishing_footprint_table$YEARS[fishing_footprint_table$NumberOfVessels<Minimum_vessels],fishing_footprint_table$GEAR_TYPE[fishing_footprint_table$NumberOfVessels<Minimum_vessels],exclude=NULL)
  gears<-sort(unique(as.character(fishing_effort$GEAR_TYPE)))
  years<-sort(unique(as.character(fishing_effort$YEAR_GROUP)))
 # if(Units=="ISEA_13"|Units=="ISEA_12"){Spatial_resolution<-min(min(diff(sort(unique(LATITUDE)))),min(diff(sort(unique(LONGITUDE)))))}

  fishing_footprint_table<-subset(fishing_footprint_table,fishing_footprint_table$NumberOfVessels>=Minimum_vessels)
  return(list(fishing_effort_data=fishing_effort,fishing_footprint_table=fishing_footprint_table,not_enough_vessels=did_not_meet_minimum_vessels,years=years,gear_types=gears,spatial_resolution=Spatial_resolution,units=Units,map_projection=map_projection))
}


#' Function to produce raster layers of effort by gear type and year group
#'
#' This function takes a fishing effort data object (created by the footprint_data function) and makes a raster map of the effort data. The outputs are a raster
#' stack with layers for each year group and gear type combination and a .png figure for each year group and gear type combination. The .png object
#' is overlaid on a static google map of the eastern or western Pacific, based on the input longitude values of the data. The footprint_data_object must
#' be a list including a data frame of raw fishing effort data (Fishing_effort_data), spatially, temporally effort data compiled by gear type (fishing_
#' footprint_table), the number of grid cells where their are not enough vessels represented to plot the data due to confidentiality rules (not_enough_vessels),
#' gear type-year-year groups over which the data is summarized (years), gear types over which the data is summarized (gear_types), the spatial resolution
#' of the data summary, #' the spatial resolution units (units) and the map projection of the data (map_projection). These 8 items are ouput from the
#' footprint_data function.
#' #'
#' @param footprint_data_object A list of data objects output from a footprint_data function
#' @keywords NPFC, fishing effort, fishing footprint
#' @export
#' @examples
#' effort_rasters<-footprint_data_object<-example_dd


footprint_map<-function(footprint_data_object){
 # footprint_data_object<-example_dd
Spatial_resolution<-footprint_data_object$spatial_resolution
Units<-footprint_data_object$units
map_projection<-footprint_data_object$map_projection
fishing_effort<-footprint_data_object$fishing_effort_data
extent1<-extent(c(min(as.numeric(as.character(fishing_effort$LONGITUDE))),max(as.numeric(as.character(fishing_effort$LONGITUDE))),min(as.numeric(as.character(fishing_effort$LATITUDE))),max(as.numeric(as.character(fishing_effort$LATITUDE)))))
raster1<-raster(extent1,resolution=Spatial_resolution,crs=map_projection)
raster.stack<-stack(raster1)

fishing_footprint<-footprint_data_object$fishing_footprint_table
years<-footprint_data_object$years
gear_types<-footprint_data_object$gear_types

for(i in 1:length(years)){
  for(j in 1:length(gear_types)){
    fishing_footprint_sub<-subset(fishing_footprint,as.character(fishing_footprint$YEARS)==years[i])
    fishing_footprint_sub<-subset(fishing_footprint_sub,fishing_footprint_sub$GEAR_TYPE==gear_types[j])
    if(dim(fishing_footprint_sub)[1]>0){
      coords1<-data.frame(as.numeric(as.character(fishing_footprint_sub$LONGITUDE)),as.numeric(as.character(fishing_footprint_sub$LATITUDE)))
      currents.project<-SpatialPointsDataFrame(coords=coords1,data=fishing_footprint_sub["NumberOfSets"], proj4string=CRS(map_projection))
      raster2<-rasterize(currents.project,raster1)
      names(raster2)[2]<-paste("YEARS",years[i],gear_types[j],sep="_")
      raster.stack<-stack(raster.stack,raster2[[2]])

    }
    if(dim(fishing_footprint_sub)[1]==0){ print(paste("Not enough vessels to plot",years[i],gear_types[j],"data",sep=" "))}
  }}

data(map4)
data(Emporer)
res1<-paste(round(Spatial_resolution,2),Units,sep="_")
if(Units=="m"){raster.stack<-projectRaster(raster.stack,crs="+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")}
rtp <- rasterToPolygons(raster.stack,na.rm=FALSE)
minsets<-min(minValue(raster.stack))
maxsets<-max(maxValue(raster.stack))
for(i in 1:length(names(raster.stack))){
  year1<-unlist(strsplit(x=names(raster.stack)[i],"_"))[2]
  gear1<-unlist(strsplit(x=names(raster.stack)[i],"_"))[3]
  year1<-gsub(".","-",year1,fixed=TRUE)
  gear1<-gsub("."," ",gear1,fixed=TRUE)
  notenough<-footprint_data_object$not_enough_vessels[row.names(footprint_data_object$not_enough_vessels)==year1,colnames(footprint_data_object$not_enough_vessels)==gear1]
  if(extent1[1,]<0){map4<-map4}
  if(extent1[1,]>0){map4<-Emporer}
  r<-rtp[i]
  m1<-ggmap(map4)+geom_polygon(data = r,
                               aes(x = long, y = lat,group=group,fill = rep(r@data[[1]],each=5)),
                               size = 0,
                               alpha = .75) +
    scale_fill_gradientn("Number Of Sets", colors = inferno(30),na.value=NA,limits=c(minsets,maxsets))+xlab("Longitude")+ylab("Latitude")
  if(extent1[1,]<0){ m1<-m1+annotate("text",-140,55,label=year1,color="white",size=4)+annotate("text",-140,53,label=gear1,color="white",size=4)+annotate("text",-138,51,label=paste("Grid cells < minimum vessels = ",notenough,sep=""),color="white",size=3)+theme(panel.background=element_rect(fill="transparent",colour=NA),legend.position="bottom")}
  if(extent1[1,]>0){ m1<-m1+annotate("text",160,50,label=year1,color="white",size=4)+annotate("text",160,49,label=gear1,color="white",size=4)+annotate("text",175,49,label=paste("Grid cells < minimum vessels = ",notenough,sep=""),color="white",size=3)+theme(panel.background=element_rect(fill="transparent",colour=NA),legend.position="bottom")}
    print(m1)

    png(paste0(year1,"_",gear1,"_",res1,".png"),height=7,width=7,unit="in",res=300)
    plot(m1)
    dev.off()


}
return(raster.stack)}

#' This function estimates distance of one set of points to another
#'
#' This function takes two lat-long pairs and calculates the distance between
#' them in m or km
#' @param lat1 Latitude of point 1
#' @param lat2 Latitude of point 2
#' @param long1 Longitude of point 1
#' @param long2 Longitude of point 2
#' @param unit either "km" or "m"
#' @keywords distance
#' @export
#' @examples
#' dist_xy()

min_dist<-function(long1,lat1){
require(sp)
 positions<-cbind(long1,lat1) #
 positions<-unique(positions)
dist<-array(length(lat1))  #unit="m"
for(i in 1:length(lat1)){
  dist2<-spDists(positions,cbind(long1[i],lat1[i]),longlat=FALSE)
  dist[i]<-min(dist2[dist2>0])
}
  return(min(dist))
}

