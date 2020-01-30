############################################################################################
#' @title  Pull data from meterologic stations in the USCRN network

#' @author Robert Lee \email{rlee@battelleecology.org}\cr

#' @description Provide a valid GHCND or WBAN ID for a site in the USCRN, the desired
#' timescale of measurements, and requested begin and end times in POSIX format,
#' and a data frame will be returned to the global environment.
#'
#' @param stationID The WBAN ID of a valid USCRN site.
#' @param TimeBgn A YYYY-MM-DD hh:mm:ss formated start time, in UTC.
#' @param TimeEnd A YYYY-MM-DD hh:mm:ss formated end time, in UTC.
#' @param timeScale The time scale of measurements to return.
#' Enter 'monthly', 'daily', 'hourly', or 'subhourly'.
#'
#' @return Data frame of USCRN data.

#' @keywords USCRN, data, process quality, data quality, gaps, commissioning
#' @export
#' @examples
#' \dontrun{
# Example inputs:
#' timeScale <- "subhourly"
#' stationID <- "USW00003047"
#' timeBgn <- as.POSIXct("2014-04-01 00:00:01", format="%Y-%m-%d %H:%M:%S", tz="UTC")
#' timeEnd <- as.POSIXct("2015-02-01 00:00:00", format = "%Y-%m-%d %H:%M:%S", tz="UTC")
#'
#' grabUSCRN(timeScale = timeScale, TimeBgn = TimeBgn, TimeEnd = TimeEnd, stationID = stationID)
#' }
#' @seealso Currently none

# changelog and author contributions / copyrights
#   Robert Lee (2017-07-18)
#     original creation
#
##############################################################################################


getUSCRNData = function (timeScale, stationID, timeBgn, timeEnd){
  ### FN START OPTIONS
  print(stationID)
  functionStart = Sys.time()
  options(stringsAsFactors = FALSE)
  library(data.table)
  library(stringr)

  ### Handle timescales
  if (!timeScale %in% c("monthly", "daily", "hourly", "subhourly")) {
    stop("Invalid 'timeScale'! Please enter one of the following: 'monthly', 'daily', 'hourly', 'subhourly'.")
  }
  timeBgn = as.POSIXct(timeBgn,  tz="UTC")
  timeEnd = as.POSIXct(timeEnd,  tz="UTC")
  years=seq(substr(timeBgn, 0, 4), substr(timeEnd, 0, 4))
  ref.seq=.make.time.seq(timeBgn, timeEnd, timeScale)

  #### Get table header info
  if (timeScale == "hourly") {
    header <- readLines("ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/hourly02/HEADERS.txt")
  }else {
    header <- readLines(paste0("ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/",
                               timeScale, "01/HEADERS.txt"))
  }
  header=unlist(strsplit(x = header[2], split = " "), recursive=F)

  # generate a master.df representing the data we want to get.
  # data returned thru the FTP will merge into this DF
  master.df=data.frame(matrix(nrow = length(ref.seq), ncol = length(header), data = NA))
  colnames(master.df)=header
  master.df$WBANNO=substr(stationID, (nchar(stationID) - 4), nchar(stationID))

  # add date time info, including unique date time rownames
  if(timeScale %in% c("hourly", "subhourly")){
    master.df$UTC_DATE=gsub(pattern = "-", replacement = "", x = substr(x=ref.seq, start=1, stop=10))
    master.df$UTC_TIME=gsub(pattern = ":", replacement = "", x = substr(x=ref.seq, start=12, stop=16))
    rownames(master.df)=paste0(master.df$UTC_DATE, "T", master.df$UTC_TIME)
  }else if(timeScale=="daily"){
    master.df$LST_DATE=gsub(pattern = "-", replacement = "", x = substr(x=ref.seq, start=1, stop=10))
    rownames(master.df)=master.df$LST_DATE
  }else if(timeScale=="monthly"){
    master.df$LST_YRMO=gsub(pattern = "-", replacement = "", x = substr(x=ref.seq, start=1, stop=7))
    rownames(master.df)=master.df$LST_YRMO
  }

  #### Generate the base FTP link
  baseLink=.make.base.link(timeScale)

  ##### Extract the station name given the ID
  station.name=.get.uscrn.name(substr(stationID, (nchar(stationID) - 4), nchar(stationID)))

  # do the downloading. Note that each timeScale is different in what it needs
  if(timeScale %in% c("subhourly", "hourly", "daily")){
    raw.data=lapply(years, function(x) .get.data(baseLink = baseLink, station.name=station.name, year = x, timeScale = timeScale, header=header))
    is.df.data=unlist(lapply(raw.data, class))=="data.frame"
    if(any(is.df.data)){
      raw.data=raw.data[is.df.data]
      ok.data=unlist(lapply(raw.data, ncol)==max(unlist(lapply(raw.data, ncol))))

      data.df=data.frame(do.call(rbind, raw.data[ok.data]))
      if(timeScale %in% c("subhourly", "hourly")){
        rownames(data.df)=paste0(data.df$UTC_DATE, "T", data.df$UTC_TIME) #Rownames are datetime
      }else{
        rownames(data.df)=data.df$LST_DATE # Row names are dates
      }
    }else{data.df=NULL} #if none of the data returned are vaild, write noting out
  }else{
    data.df=.get.monthly.data(station.name = station.name, header=header)
    rownames(data.df)=data.df$LST_YRMO #Rownames are year mon
  }

  # Put data in rows common to both master df and data
  common.rows=intersect(x = rownames(master.df), y = rownames(data.df))
  master.df[common.rows,]=data.df[common.rows,]

  return(master.df)
}