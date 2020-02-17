############################################################################################
#' @title  Download Data Served through ACIS?

#' @author Cody Flagg \email{cflagg@battelleecology.org}\cr

#' @description GETUSOMEDATA
#'
#' @param sid Parameter of class character. The (GHCND) site ID of the site data should be downloaded for.
#' @param start_date Parameter of class character. The year-month-day (e.g. "2017-01-02") of the last
#'  date to get data for.
#' @param end_date Parameter of class character. The year-month-day (e.g. "2017-01-02") of the last
#'  date to get data for.
#' @param vars Optional. A string of element codes, separated by commas only. Defaults to "mint,maxt,avgt,pcpn,snow,snwd,obst"
#' @return Writes results file to the specified directory.

#' @keywords process quality, data quality, gaps, commissioning

#' @examples Currently none
#' @seealso Currently none
#' @export
##############################################################################################

getSingleACIS <- function(sid, start_date,end_date,vars="mint,maxt,avgt,pcpn,snow,snwd,obst"){
  library(jsonlite)
  library(httr)
  library(plyr)
  # pull station meta-data!
  #vars <- "mint,maxt,avgt,pcpn,snow,snwd,obst"
  #browser()
  #params <- list(sid = sid, elems = vars, sdate = start_date, edate = end_date, output = "json")
  params <- list(sid = sid, elems = vars, sdate = start_date, edate = end_date, output = "json")
  #req <- httr::POST(url = "http://data.rcc-acis.org/MultiStnData", body = params, encode='json')
  req <- suppressMessages(httr::POST(url = "http://data.rcc-acis.org/StnData", body = params, encode='json'))
  dat <- suppressMessages(jsonlite::fromJSON(content(req, "text"))) ## much better than parsing a csv
  out <- as.data.frame(dat$data)
  proceed <- ncol(out)
  ## if for some reason there are no data -- keep going
  if (proceed == 0){
    return("NO DATA")
  } else {
    names(out) <- c("date", strsplit(vars, ",")[[1]])
    print(sid)
    out$sid <- sid
    out$state #dat$meta$statedat$meta$state
    return(out)
  }
}

