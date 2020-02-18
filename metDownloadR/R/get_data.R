############################################################################################
#' @title  Pull data for a site found with metScanR

#' @author Robert Lee \email{rlee@battelleecology.org}\cr
#' Josh Roberti\cr

#' @description This function takes an input of the site metadata from a metScanR search,
#' as well as start and end dates to download data for, and downloads
#'
#' @param site_meta a metScanR list element.
#' @param start_date A YYYY-MM-DD hh:mm:ss formated start time
#' @param end_date A YYYY-MM-DD hh:mm:ss formated end time
#' @param token Optional, but required for Mesonet stations. See: \url{https://developers.synopticdata.com/mesonet/} for more information.
#'
#' @return Data frame data for the input site found with metScanR, if available.

#' @keywords USCRN, data, process quality, data quality, gaps, commissioning
#' @export
#' @examples
#' \dontrun{
#' # Get Data for one SCAN site near CPER
#'cper_sites=metScanR::getNearby(siteID="NEON:CPER", radius = 30)

#' uscrn_out=metDownloadR::getData(site_meta = cper_sites$USW00094074,
#'                               start_date = "2018-10-01",
#'                              end_date = "2018-10-31", 
#'                              temp_agg="monthly")

#' # Get all october 2018 data from sites within a 5 km radius of CPER
#' out=lapply(
#' metScanR::getNearby(siteID="NEON:CPER", radius = 30),
#' metDownloadR::getData,
#' start_date = "2018-10-01",
#' end_date = "2018-10-31", 
#' temp_agg="daily")
#' }
#' @seealso Currently none

# changelog and author contributions / copyrights
#   Robert Lee (2017-07-18)
#     original creation
#   Josh Roberti (2018-12-17)
#     logic changes, reordering of sites and option to only grab data for N sites
##############################################################################################


getData<-function(site_meta, start_date, end_date, temp_agg, token=NA, save_file=NULL){
  #browser()
  #set out to no data right off the bat:
  out<-"NO DATA"
  #check to make sure temp_agg is appropriate
  
  if(!(tolower(temp_agg) %in% c("monthly","daily","hourly","subhourly"))){
    stop("temp_agg must equal 'monthly', 'daily', 'hourly', or 'subhourly' ")
  }
  temp_agg=tolower(temp_agg)
  
  platform<-site_meta$platform
  identifiers<-site_meta$identifiers
  #make dataframe with site's platform and identifiers:
  stationMeta<-data.frame(platform,identifiers)
  
  #filter the resTracking to rows where the platform and temp_agg have a match (map be multiple)
  if(any(stationMeta$platform=="Mesonet")){
    resTrackFilter<-resTracking[which(resTracking$platform==site_meta$platform & resTracking$value=="nominal"),]
  }
  else{
    resTrackFilter<-resTracking[which(resTracking$platform==site_meta$platform & resTracking$value==temp_agg),]
  }
  #browser()
  #open station resolution mapping:
  #save(resTracking,file = "C:/Users/jroberti/Git/external-sites/metDownloadR/data/stationResolutionMap.rda")
  #load("../metDownloadR/data/stationResolutionMap.rda")
  #load resTracking from /data folder of R package
  #metDownloadR:::stationResolutionMap
  #get station's platform and IDs
  
  #merge resTrackFilter with station identifiers:
  
  if(nrow(resTrackFilter)>0){
    mergedMeta<-merge(stationMeta,resTrackFilter,by=intersect(names(stationMeta),names(resTrackFilter)))
    #browser()
    #update on 2019-02-04.  Just want to use mesonet repository because they'll have data we nned. ACIS only has few variables
    mergedMeta<-mergedMeta[grep("Mesonet",mergedMeta$Rpackage),]
    #fget unique R packages in case it's repeating:
    useFuncR<-unique(mergedMeta$Rpackage)
    #logic to pull data
    if(length(useFuncR)!=0){ 
      if(length(useFuncR)==1){
        useTheseMeta<-mergedMeta[1,]
      }else if(length(useFuncR)>=1){#example site: USW00094074 (has IDs that could be used with ACIS or USCRN)
        #find the one with the more returns:
        useFuncR<-useFuncR[which.max(table(mergedMeta$Rpackage))]
        useTheseMeta<-mergedMeta[grep(useFuncR,mergedMeta$Rpackage)[1],]
      }
      #older logic using platform:
      # if(platform=="COOP"){
      #   sid<-useTheseMeta$id
      #   out<-getSingleACIS(sid = sid, start_date = start_date, end_date = end_date)
      # }
      #grab the first row of the mergedData (just want to run data pull function once per station, not per ID)
      if(useTheseMeta$Rpackage=="ACIS"){
        #browser()
        sid<-useTheseMeta$id
        out<-getSingleACIS(sid = sid, start_date = start_date, end_date = end_date)
        #browser()
      }
      else if(useTheseMeta$Rpackage=="RNRCS"){
        if(useTheseMeta$idType=="SCAN"){
          #browser()
          idType<-site_meta$identifiers$idType[!grepl(x=site_meta$identifiers$idType, pattern = "Mesonet")]
          sid<-site_meta$identifiers[site_meta$identifiers$idType==idType, "id"]
          out<-RNRCS::grabNRCS.data(site_id = sid, network = idType, timescale = temp_agg, DayBgn = start_date, DayEnd =  end_date)
        }else if(useTheseMeta$idType=="BOR"){
          sid<-gsub(pattern = "BOR:", replacement = "", x = site_meta$identifiers[site_meta$identifiers$idType=="BOR", "id"])
          out<-RNRCS::grabBOR.data(site_id = sid, timescale = temp_agg, DayBgn = start_date, DayEnd =  end_date)
        }
      }
      else if(useTheseMeta$Rpackage=="USCRN"){
        #browser()
        #useTheseMeta<-mergedMeta[4,]
        sid<-useTheseMeta$id
        out<-metDownloadR::getUSCRNData(stationID = sid, timeScale = temp_agg, timeBgn = start_date, timeEnd = end_date)
      }
      else if(useTheseMeta$Rpackage=="Mesonet"){
        out="No data returned- please input a Mesonet API token.\nSee: https://developers.synopticdata.com/mesonet/ "
        if(!is.na(token)){
          #add test dates for grabbing preliminary data:
          #pingTimeBgn<-mesoTime(start_date)
          #pingTimeEnd<-mesoTime(as.character(as.POSIXct(start_date)+36*60*60))
          #assign mesonet ID:
          mesoId=useTheseMeta$id
          #collapse the time stamps into strings without separators; need this for mesonet API
          timeBgn<-mesoTime(start_date)
          timeEnd<-mesoTime(end_date)
          #print mesonetID to end-user:
          print(mesoId)
          # out<-metDownloadR::getMesonetData(token="16088448e1b149509e45e401196106f0",mesoId=mesoId,
          #                timeBgn=timeBgn,timeEnd=timeEnd)
          out<-metDownloadR::getMesonetData(token=token,mesoId=mesoId,timeBgn=timeBgn,timeEnd=timeEnd)
          #pingTimeBgn=pingTimeBgn,
          #pingTimeEnd=pingTimeEnd
          
        }
      }
      ### input helper function to clean missing data rows out.
      if(class(out)=="data.frame"){
        out=.remove.nas(out)
      }
    }
  }
  
  if(is.null(save_file)){
    return(out)
  }else{
    message("Saving data")
    saveRDS(object = out, file = save_file)    
  }
}

#function to convert time to mesonet format:
mesoTime<-function(x){
  x<-gsub("-|:| ","",x)
  return(x)
}
