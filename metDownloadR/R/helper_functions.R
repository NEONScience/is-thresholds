## Helper functions

## The find the download name given WBAN
.get.uscrn.name=function(wban){
  USCRN=metScanR::getNetwork(network = "USCRN")
  names=unlist(lapply(USCRN, "[[", "namez"))
  ids=lapply(USCRN, "[[", "identifiers")
  temp=list(wban.list=unlist(lapply(ids, function(i) i$id[i$idType=="WBAN"])),
            download.names=tolower(gsub(pattern = " ", replacement = "_", x = names))
  )
  out=temp$download.names[which(temp$wban.list==wban)]
  return(out)
}

# Reference time seq to merge into
.make.time.seq=function(timeBgn, timeEnd, timeScale){
  time.info=list(
    timeScale=c("monthly", "daily", "hourly", "subhourly"),
    steps=c("1 month", "1 day", "1 hour", "5 min")
  )
  ref.seq=lubridate::round_date(
    x=seq.POSIXt(
      from = timeBgn,
      to = timeEnd,
      by = time.info$steps[time.info$timeScale==timeScale]),
    unit = time.info$steps[time.info$timeScale==timeScale]
  )
  return(ref.seq)
}

# Query for and return base links
.make.base.link=function(timeScale){
  uscrn <- readLines("https://www.ncdc.noaa.gov/crn/qcdatasets.html",
                     warn = F)
  dataLines <- unlist(strsplit(uscrn[grep(pattern = paste0("/",
                                                           timeScale), uscrn)], "\""))
  dataLinks <- unlist(dataLines[which(grepl("/data/", dataLines) &
                                        grepl("http", dataLines))])
  baseLink <- gsub(pattern = "http", replacement = "https",
                   x = dataLinks[which(min(nchar(dataLinks)) == nchar(dataLinks))]
  )
  return(baseLink)
}

# Return link to the file
.find.file=function(base, station.name){
  files=readLines(base)
  file=files[grep(pattern = station.name, x = files, ignore.case = T)]
  temp=unlist(stringr::str_split(string=file, pattern = " "), recursive = F)
  file.name=temp[length(temp)]
  out=paste0(base, file.name)
  return(out)
}

# Get data (for non-monthly data) ####
.get.data=function(baseLink, station.name, year, timeScale, header){
  raw.data=NA
  #browser()
  year<-2003
  if(any(grepl(pattern = paste0(">", year), x = readLines(baseLink)))){
    tails="01/"
    if(timeScale=="hourly"){tails="02/"}
    base=paste0("ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/", timeScale, tails, year, "/")
    file.link=.find.file(base, station.name)
    raw.data=read.delim(file =file.link, sep="", header = F, colClasses = "character")
    if(ncol(raw.data)==length(header)){
      colnames(raw.data)=header
    }
  }
  return(raw.data)
}

# Get monthly data ####
.get.monthly.data=function(station.name, header){
  base="ftp://ftp.ncdc.noaa.gov/pub/data/uscrn/products/monthly01/"
  file.link=.find.file(base, station.name)
  raw.data=read.delim(file =file.link, sep="", header = F, colClasses = "character")
  if(ncol(raw.data)==length(header)){
    colnames(raw.data)=header
  }
  return(raw.data)
}

# clean out NAs and Ms
.remove.nas=function(df){
  data.indx=!grepl(colnames(df), pattern = "date|*id|utc*|LST*|wbanno", ignore.case = T)
  missing.rows=unlist(lapply(rownames(df), function(r) all(df[r,data.indx]=="M"|is.na(df[r,data.indx]))))
  out=df
  if(any(missing.rows)){out=df[-which(missing.rows),]}
  if(nrow(out)==0){out="NO DATA"}
  return(out)
}
