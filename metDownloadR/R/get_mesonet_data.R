############################################################################################
#' @title  Pull data from meterologic stations in the mesonet network

#' @author Josh Roberti \email{jroberti@battelleecology.org}\cr

#' @description Provide a valid mesonet ID for a site in the mesonetwork; your personal
#' mesonetwork data token, and requested begin and end times in POSIX format,
#' and a data frame will be returned to the global environment.
#'
#' @param token Personalized Mesonetwork token. This is needed to download data; 
#' see \url{https://developers.synopticdata.com/mesonet/v2/auth/} for more information.
#' @param mesoId Mesonet ID for site of interest
#' @param timeBgn A YYYY-MM-DD hh:mm:ss formated start time, in UTC.
#' @param timeEnd A YYYY-MM-DD hh:mm:ss formated end time, in UTC.

#' @return Data frame of mesonet data

#' @keywords mesonet, weather, meteorology, climatology, ecology, data, metScanR

#' @examples
#' \dontrun{
# Example inputs:
#' token="16088448e1b149509e45e401196106f0"
#' mesoId="TR951"
#' timeBgn="2019-01-01 00:00:00"
#' timeEnd="2019-02-01 00:00:00"
#'
#' getMesonetData(timeBgn = timeBgn,
#' timeEnd = timeEnd, 
#' mesoId=mesoId, 
#' token=token)
#' }
#' 
#' @export


# changelog and author contributions / copyrights
#   Josh Roberti(2019-02-06)
#     original creation
#   Robert Lee (2019-02-12)
#     updated example and param descriptions
#   Josh Roberti (2019-02-18)
#     logic for pulling 1-day's worth of data to determine reported variables
##############################################################################################

getMesonetData<-function(token,mesoId,timeBgn,timeEnd){
  #pingTimeBgn,pingTimeEnd
  library(httr)
  library(jsonlite)
  #grab data via "ping" dates to check what types of data are available:
  #callTextPing<-paste0("https://api.synopticdata.com/v2/stations/timeseries?&token=",token,"&start=",pingTimeBgn,"&end=",pingTimeEnd,"&obtimezone=utc&output=json&stid=",mesoId)
  #data.call<-httr::GET(callTextPing)
  #daily:
  # data.raw<-rawToChar(data.call$content)
  # data.list <- jsonlite::fromJSON(data.raw)
  # if(length(data.list$STATION)!=0){
  #   browser()
  #   if(class(data.list$STATION)=="list"){
  #     data.list$STATION<-unlist(data.list$STATION,recursive = F)
  #   }
  # }
  
  #replace the last two digits of timeEnd with 0s:
  timeBgn<-substr(timeBgn,1,nchar(timeBgn)-2)
  timeEnd<-substr(timeEnd,1,nchar(timeEnd)-2)
  #timeEnd<-paste0(timeEnd,"00")
  
  callText<-paste0("https://api.synopticdata.com/v2/stations/timeseries?&token=",token,"&start=",timeBgn,"&end=",timeEnd,"&obtimezone=utc&output=json&stid=",mesoId)
  #get station data for a single station,; WBB:
  #hourly
  data.call<-httr::GET(callText)
  #browser()
  #daily:
  data.raw<-rawToChar(data.call$content)
  data.list <- jsonlite::fromJSON(data.raw)
  #check to see that data were actually pulled
  if(length(data.list$STATION)!=0){
    if(class(data.list$STATION)=="list"){
      data.list$STATION<-unlist(data.list$STATION,recursive = F)
    }
    #browser()
    #check to see which observations are lists, these will all be NULL/missing:
    for(i in 1:length(data.list$STATION$OBSERVATIONS)){
      #find the null entries:
      nullEntries<-which(unlist(lapply(data.list$STATION$OBSERVATIONS[[i]], function(x) is.null(x)))==T)
      #replace null entires with NA (will retain entry)
      if(length(nullEntries)!=0){
        data.list$STATION$OBSERVATIONS[[i]][nullEntries]<-NA
      }
      #unlist and put back into nested list.
      #data.list$STATION$OBSERVATIONS[[i]]<-unlist(data.list$STATION$OBSERVATIONS[[i]])
      #data.list$STATION$OBSERVATIONS[[i]]<-data.list$STATION$OBSERVATIONS[[i]]
      #browser()
    }
    #data.list$STATION$OBSERVATIONS[i]
    data.df <- data.frame(do.call(cbind,lapply(data.list$STATION$OBSERVATIONS, as.data.frame)))
    #assign names:
    names(data.df)<-names(data.list$STATION$OBSERVATIONS)
    #convert date_time from factor to POSIX:
    data.df$date_time<-as.POSIXct(data.df$date_time,format='%Y-%m-%dT%H:%M:%SZ',tz = "UTC")
    print(paste0("Data downloaded for ", mesoId))
  }
  else{
    data.df<-"NO DATA"
    print(paste0("NO DATA for ", mesoId))
  }
  #return the data:
  return(data.df)
}

