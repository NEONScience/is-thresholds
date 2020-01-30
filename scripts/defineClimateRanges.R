##############################################################################################
#' @title Function to define climatological thresholds for NEON sites

#' @author
#' Josh Roberti \email{jroberti@battelleecology.org} \cr
#'
#' @description This function munges climatological data from sites located near a point of
#' interest over a specified time period.  It then quantifies climatological ranges for said
#' environmental variable (e.g., air temperature)

#' @param stack Name of database stack to pull from; int, cert, or prod [character]
#' @param meas Name of measurement stream to get data for, in the form NEON.DOM.SITE.DPL.PRNUM.001.TERMS.HOR.VER.TMI, e.g. NEON.D01.HOPB.DP0.20016.001.01379.101.100.000 [character]
#' @param startDate Earliest end date of data to pull, in the form YYYY-MM-DDTHH:MM:SS.SSSZ [character]
#' @param endDate Latest end date of data to pull, in the form YYYY-MM-DDTHH:MM:SS.SSSZ [character]

#' @return This function returns a data frame

#' @references
#' License:

#' @export

##############################################################################################

#source the functions needed for script:
 #devtools::install(path = '~/GitHub/metget/')
#
# # FUNCTION 1:
# ## quickly extract data from the metScanR DB outputs
# sf_getter <- function(input, first, second=NULL){
#   if (is.null(second)){
#     lapply(input, "[[", first)
#   } else {
#     lapply(lapply(input, "[[", first), "[[", second)
#   }
# }

#MAIN FUNCTION:
#this will be built off of metScanR logic

defineClimateRanges<-function(siteID="NEON:ORNL",radius=40,startDate="2010-01-01",
                              endDate=Sys.Date(),recrunchThresholds=F,overwriteData=F,numStations=3, save.dir, token){
  #browser()
  startSysTime<-Sys.time()
  options(stringsAsFactors = F)
  #first check to see if thresholds have already been gathered:
  strt=Sys.time()
  sitename<-gsub(".*:","",siteID)
  if(exists("vars")){
    varsFile<-gsub("\\s+","", vars)
  }else{
    varsFile="all"
  }
  filename<-paste0(save.dir, "/", sitename,"_",varsFile,"_thresholds.rds")
  #filename<-paste0("~/GitHub/NEON-FIU-document-IPT/dataQAQC/thresholds/",sitename,"_",varsFile,"_thresholds.rds")
  if(file.exists(filename) & overwriteData==F & recrunchThresholds==F){
      stop("File already exists.  If you want to re-scrape data AND re-generated thresholds, plesae set 'overwriteData' = T. \n
         If you'd like to just re-generated Thresholds, please set 'recrunchThresholds' = T ")
  }else if(file.exists(filename) & overwriteData==F & recrunchThresholds==T){
    datafile<-readRDS(filename)
    #browser()
    seasThresh.list<-createThresholds(data=datafile$DataRangesThresholds$data)
    out<-list(siteID,seasThresh.list,siteMeta.df)
    names(out)<-c("siteNEON","DataRangesThresholds","externalMetadata")
  }else{
    # library(metDownloadR)
    # library(rjson)
    #grep NEON site name for use in website name
    siteNEON<-gsub(".*\\:","",siteID)
    #set dynamic website url:
    website<-paste0("http://data.neonscience.org/api/v0/locations/",siteNEON)
    #get the data from the website:
    json_data <- rjson::fromJSON(file=website)
    #grab elevations
    elevNEON<-json_data$data$locationElevation
    latNEON<-json_data$data$locationDecimalLatitude
    lonNEON<-json_data$data$locationDecimalLongitude
    #set elevation bounds for filtering nearby stations
    #elevBounds<-c(elevNEON-50,elevNEON+50)
    #nearbySites<-siteFinder(siteID=siteID,radius=radius,vars=vars) #still need to add in the
    #MAKE SURE TO ADD IN ELEVATION!!
   # nearbySites2<-metScanR::siteFinder(siteID=siteID,radius=radius,elevMin = elevBounds[1],elevMax = elevBounds[2])
    nearbySites2<-metScanR::siteFinder(lat=latNEON,lon=lonNEON,radius=radius,network = c("mesonet"))
    #check to see how many stations were returned:
    lengthResult<-length(nearbySites2)
    if(lengthResult==0){
      stop("No nearby sites returned.  Please expand the search radius")
    }

    ####### [START] Reorganizing nearby sites based on distance from POI #######
    #get metadata for each site:
    siteMetaLocations.list<-lapply(nearbySites2, function(x) x$location[1:2])
    #make into dataframe
    siteMetaLocations.df<-do.call(rbind, siteMetaLocations.list)
    #add station names:
    siteMetaLocations.df$station<-names(siteMetaLocations.list)
    #create distance from POI
    siteMetaLocations.df$distFromSite_m<-apply(siteMetaLocations.df[,c("longitude_dec","latitude_dec")], 1,
                                      function(x) geosphere::distCosine(p1=data.frame(lonNEON,latNEON),p2=x))
    #add sequence so we can retain where the initial placement of the sites are in the nearbSites list:
    siteMetaLocations.df$index<-1:length(siteMetaLocations.df$distFromSite_m)
    #sort the metadata by distance:
    siteMetaLocations.df<-siteMetaLocations.df[with(siteMetaLocations.df,order(distFromSite_m)),]
    #reorder nearby sites so it looks for closer stations first:
    nearbySites2<-nearbySites2[siteMetaLocations.df$index]
    ####### [END] Reorganizing nearby sites based on distance from POI #######


    ###### [START] Grab all reported variables at a site ######################
    #elementsNearby<-lapply(lapply(nearbySites2, "[[", "elements"),"[[", "element")
    elementsNearby<-lapply(nearbySites2, "[[", "elements")
    #grep over the elements to find those that we need for NEON thresholds
    mesoElements<-metScanR:::metScanR_terms$master[which(metScanR:::metScanR_terms$master$Agency.Network=="Mesonet"),]
    #NEON threshold elements:
    elementsNEON<-c("air temperature","precipitation","pressure","wind speed","radiation","net radiation",
                    "shortwave radiation","longwave radiation","particulate","surface temperature",
                    "relative humidity","dew point","heat_flux","soil temperature","soil moisture")
    #go thru the NEON elements vector and grep these terms within the mesoElements DF:
    elementsMeso<-lapply(elementsNEON, function(x) mesoElements$ElementCode[grep(x,mesoElements$ElementName)])
    #name the list using elementsNEON:
    names(elementsMeso)<-elementsNEON
    #create empty list to store matched terms:
    elementsMatched<-list()
    for(i in 1:length(elementsNearby)){
      #grep elementsMeso at each site:
      elementsMatched[[i]]<-unlist(lapply(unlist(elementsMeso), function(x) elementsNearby[[i]]$element[grep(x,elementsNearby[[i]]$element)]))
      #grep the returned element codes in the mesoElements DF:
      elementsMatched[[i]]<-unique(names(elementsMeso)[unlist(lapply(elementsMatched[[i]], function(x) grep(x,elementsMeso)))])
    }
    sortedElements<-sort(table(unlist(elementsMatched)),decreasing = T)
    #create rankings of the elements:
    ranked.df<-data.frame(elements=names(sortedElements),rank=rev(1:length(sortedElements)))
    #sort by rank
    ranked.df<-ranked.df[with(ranked.df, order(rank)), ]
    #go thru the nearbySites list and start by grabbing sites in order of rank:
    elementsSearch<-lapply(ranked.df$elements, function(x) unlist(elementsMeso[grep(x,names(elementsMeso))]))
    #go thru and grep each element code of the
    pullData<-list()
    for(i in 1:length(elementsSearch)){
      sites<-names(elementsNearby)[unlist(lapply(elementsSearch[[i]], function(x) grep(x,lapply(elementsNearby,"[[","element"))))]
      #index of [[1]] will reflect the specific variable; don't change
      # if(i==10){
      #   browser()
      # }
      #elementDF<-lapply(elementsSearch[[i]], function(x) elementsNearby[grep(x,lapply(elementsNearby,"[[","element"))])[[1]]
      elementDF<-unlist(lapply(elementsSearch[[i]], function(x) elementsNearby[grep(x,lapply(elementsNearby,"[[","element"))]),recursive = F)
      #for each site that has variable of interest, go and grab the dates of the variable:
      dates<-list()
      for(j in 1:length(elementDF)){
        #  print(j)
        # if(j==1){
        #   browser()
        # }

        varGrep<-unlist(sapply(elementsSearch[[i]], function(x) grep(x,elementDF[[j]]$element)))
        #grep(elementsSearch[[i]],elementDF[[j]]$element)
        temp.df<-elementDF[[j]][varGrep,]
        dates[[j]]<-data.frame(site=sites[j],
                               dateBegin=temp.df$date.begin,
                               dateEnd=temp.df$date.end)
        #get record Length:
        presentDate<-grep("present",dates[[j]]$dateEnd)

        if(length(presentDate)!=0){
          dates[[j]]$dateEnd[presentDate]<-as.character(Sys.Date())
        }
        # if(dates[[j]]$dateEnd=="present"){
        #   #set dateEnd to System Date:
        #   dates[[j]]$dateEnd<-as.character(Sys.Date())
        #   if(dates[[j]]$dateEnd==17946){
        #     browser()
        #   }
        #}
        #subtract End Date from Start Date for record length
        dates[[j]]$recLenAbs<-as.numeric(as.Date(dates[[j]]$dateEnd)-as.Date(dates[[j]]$dateBegin))
        #get record length relative to the time requested by the user:
        dates[[j]]$recLenRel<-round(dates[[j]]$recLenAbs/as.numeric(as.Date(endDate)-as.Date(startDate)),3)

      }
      dates.out<-do.call(rbind,dates)
      #remove any rows where the length of record is 0:
      zeroRecordLen<-which(dates.out$recLenAbs==0)
      if(length(zeroRecordLen)!=0){
        dates.out<-dates.out[-zeroRecordLen,]
      }
      #elementDates<-lapply(elementsSearch[[i]], function(x) grep(x,lapply(unlist(elementDF,recursive = F),"[[","element")))
      pullData[[i]]<-dates.out
    }
    names(pullData)<-ranked.df$elements

    #browser()
    #Go thru each nested list and identify sites with multiple entries for the same variable.
    #Navigate to these entries and then keep the entry with the longest record:
    for(i in 1:length(pullData)){
      #find duplicate site entries per variable:
      dupSiteVars<-which(duplicated(pullData[[i]]$site)==T)
      if(length(dupSiteVars)!=0){
        dupSiteName<-pullData[[i]]$site[dupSiteVars]
        #find the longest record for the duplicated site-variable combo:
        for(j in 1:length(dupSiteName)){
          recordRemove<-which.min(pullData[[i]]$recLenRel[grep(dupSiteName[j],pullData[[i]]$site)])
          #remove Site Index:
          removeSiteInd<-grep(dupSiteName[j],pullData[[i]]$site)[recordRemove]
          #remove the removeSiteInd row from pullData[[i]]
          if(length(removeSiteInd)!=0){
            pullData[[i]]<-pullData[[i]][-removeSiteInd,]
          }
        }
      }

      #Go thru and rank the sites based on the length of the recLenRel for each variable (pullData[[i]]):
      pullData[[i]]$scoreRec<-rank(pullData[[i]]$recLenRel)
      pullData[[i]]<-pullData[[i]][with(pullData[[i]], order(-scoreRec)), ]
    }

    #check to make sure all nested DFs have data:
    haveData<-which(sapply(pullData, function(x) nrow(x))>0)
    #keep only those that have data
    pullData<-pullData[haveData]


    numVars<-list()
    keepThese.list<-list()
    for(i in 1:length(pullData)){
      for(j in 1:nrow(pullData[[i]])){
        #go thru each variable, and find how many variables in our list of interest that each site reports:
        numVars[[j]]<-unlist(sapply(pullData, function(x) grep(pullData[[i]]$site[j],x$site)))
        numVars[[j]]["rank"]<-mean(numVars[[j]])
        numVars[[j]]<-data.frame(t(cbind(numVars[[j]])))
        #numVars[[j]]$site<-pullData[[i]]$site[j]
        numVars[[j]]<-data.frame(site=pullData[[i]]$site[j],numVars[[j]])

      }
      numVars.df<-do.call(plyr::rbind.fill,numVars)
      keepThese.list[[i]]<-numVars.df[which.min(numVars.df$rank),]
    }
    keepThese.df<-do.call(plyr::rbind.fill,keepThese.list)
    #keep unique DF:
    keepThese.df<-unique(keepThese.df)
    #remove all NA cols:
    dontUse1<-which(apply(keepThese.df, 2, function(x) all(is.na(x))))
    if(length(dontUse1)!=0){
      keepThese.df<-keepThese.df[,-dontUse1]
    }
    #remove site and rank from search
    dontUse2<-grep("site|rank",names(keepThese.df))
    varsOnly<-keepThese.df[,-dontUse2]
    rownames(varsOnly)<-keepThese.df$site

    #find first instances of each variable:
    endSite<-max(apply(varsOnly, 2, function(x) min(which(!is.na(x)))))
    pullData<-rownames(varsOnly)[1:endSite]

    #filter to sites to be pulled
    nearbySites2<-nearbySites2[which(names(nearbySites2) %in% pullData)]

    #browser()

    #
    #
    # browser()
    #
    # #Come up with variable ranks based on availability:
    # browser()
    # varLengths<-sapply(pullData, function(x) nrow(x))
    # varScore<-sum(varLengths)/varLengths
    # for(i in 1:length(varScore)){
    #   pullData[[i]]$scoreVar<-varScore[i]
    # }
    # #combine all pullData into dataframe:
    # pullData.df<-data.frame(do.call(rbind,pullData),row.names = NULL)
    # #come up with overall score as a function of record length score and variable score
    # pullData.df$scoreTot<-pullData.df$scoreRec*pullData.df$scoreVar
    # #go thru and come up with a collective score for each site:
    # siteScores<-sort(tapply(pullData.df$scoreTot, pullData.df$site, sum),decreasing = T)
    # names(siteScores)
    # browser()


    ###### [END] Grab all reported variables at a site ######################

    ####### [START] Grab data via Robert's getData function #############
    #only get data from N number of sites so code doesn't run forever.  User can specify at function level:
    siteData.list<-list()
    keepThese<-list() #for keeping names of stations that have data:
    keepNames<-list()
    j<-1
    #browser()
    while(length(siteData.list)<numStations & j<=length(nearbySites2)){
      #browser()
      out<-metDownloadR::getData(site_meta = nearbySites2[[j]],start_date = startDate, end_date = endDate, temp_agg="daily", token=token)
      #only add it if it's a dataframe:
      if(class(out)=="data.frame"){
        k<-length(siteData.list)+1
        siteData.list[[k]]<-out
        keepThese[[k]]<-j
        keepNames[[k]]<-siteMetaLocations.df$station[j]
        print(k)
      }
      #remove NHLL ind
      print(length(siteData.list))
      j<-j+1
    }

    #browser()

    saveRDS(siteData.list,filename)
    #browser()
    end=Sys.time()
    message(difftime(end, strt, units = "secs"))
    print(paste0("Mesonet data for ", siteID, " saved to ", filename))
}

  #   #siteData.list<-lapply(nearbySites2, getData, start_date = startDate, end_date = endDate,temp_agg="daily")
  #   #Data will be cleaned (no NAs or Ms since Robert added cleaning helper function)
  #   ### [START] Data Filtering Round 1 ###
  #   #keep data frames, remove other ish:
  #   #keepThese<-which(unlist(lapply(siteData.list, function(x) class(x)))=="data.frame")
  #   #keepNames<-names(nearbySites2)[keepThese]
  #   #find the appropriate sites in the metScanR DB:
  #   #siteData.list<-siteData.list[keepThese]
  #   #name the list with station IDs of metScanR DB:
  #
  #
  #
  #   #INFORMATION BELOW IS COMMENTED OUT AS OF 2019-02-07 IN ORDER TO GENERATE AND THEN KEEP DATA FOR FUTURE ANALYSIS
  #   names(siteData.list)<-unlist(keepNames)
  #   #throw error if no sites are returned
  #   if(length(siteData.list)==0){
  #     stop("Nearby Sites found, but no data were returned for temporal period of interest.  Please try expanding your search.")
  #   }
  #   ####### [END] Grab data via Robert's getData function #############
  #
  #
  #   ####### STOPPING HERE GIVES YOU ALL THE DATA IN A NESTED DATA FRAME FORMAT #######
  #
  #
  #
  #   ### [START] Data Filtering Round 2 ###
  #   #find the number of columns from each dataframe
  #   ncols<-lapply(siteData.list, function(x) ncol(x))
  #   #group by same number of columns:
  #   colGroups<-names(table(unlist(ncols)))
  #   #go thru list of returned dataframes and rbind those with similar columns:
  #   climateData.list<-list()
  #   #browser()
  #   for(i in 1:length(colGroups)){
  #     climateData.list[[i]]<-do.call(rbind,siteData.list[grep(paste0("\\b",colGroups[i],"\\b"),unlist(ncols))])
  #   }
  #   ### [END] Data Filtering Round 2 ###
  #
  #
  #   ### [START] Grab Data Type(s) requested ###
  #   fuzzyTerms<-metScanR:::metScanR_terms$traceability
  #   #manually add ACIS names to the fuzzy terms database (they are not in metScanR fuzzy database yet)
  #   termsACIS<-tibble::tibble(mainTerms=c("Temperature","Temperature","Temperature","Precipitation","Snow","Snow","Temperature"),
  #                         subTerms=c("Minimum Air Temperature","Maximum Air Temperature","Average Air Temperature",
  #                                    "Precipitation","Snowfall","Snow Depth","Observed Air Temperature"),
  #                         subCodes=c("mint","maxt","avgt","pcpn","snow","snwd","obst"),
  #                         agency.network=rep("ACIS",7))
  #   fuzzyTerms<-rbind(fuzzyTerms,termsACIS)
  #
  #
  #
  #   #run below to filter to specific variables:
  #   # filter dataframes using user input (vars) and term traceability via fuzzy sheet:
  #   keepCols<-fuzzyTerms[grep(tolower(vars), tolower(fuzzyTerms$subTerms)),]
  #   #browser()
  #   #Add date column to the grep call:
  #   grepTerms<-c("Date","date","DATE",unique(unlist(c(keepCols[,c("subCodes","subTerms")]),use.names = F)))
  #   #remove periods from names, convert all to lowercase, and trim white space from start and end:
  #   dataNames<-lapply(climateData.list, function(x) trimws(gsub("\\.+"," ",names(x)),"both"))
  #   #browser()
  #   ##get names of dataframes within climateData
  #   #natch the grepCodes with the column names within the
  #   climateData.list2<-list()
  #   for(i in 1:length(dataNames)){
  #     keepDataCols<-unique(unlist(lapply(grepTerms, function(x) grep(x,tolower(dataNames[[i]])))))
  #     climateData.list2[[i]]<-climateData.list[[i]][,keepDataCols]
  #     #bin the data by months and create seasons:
  #     #partition the data into meteorological seasons
  #     date_col<-grep("date",tolower(names(climateData.list2[[i]])))
  #     print(tolower(names(climateData.list2[[i]])))
  #     #print(date_col)
  #     browser()
  #     climateData.list2[[i]]$month<-as.character(substr(climateData.list2[[i]][,date_col],6,7))
  #     #Initiate a season vector
  #     climateData.list2[[i]]$season<-NA
  #     climateData.list2[[i]]$season[which(climateData.list2[[i]]$month %in% c("03","04","05"))]<-"spring"
  #     climateData.list2[[i]]$season[which(climateData.list2[[i]]$month %in% c("06","07","08"))]<-"summer"
  #     climateData.list2[[i]]$season[which(climateData.list2[[i]]$month %in% c("09","10","11"))]<-"fall"
  #     climateData.list2[[i]]$season[which(climateData.list2[[i]]$month %in% c("12","01","02"))]<-"winter"
  #     #just keep min and maximums, date, and season
  #     #if(vars="air temperature"){
  #     grepKeep<-paste0(c(dataNames[[i]][keepDataCols],"season"),collapse="|")
  #     keepCols<-grep(grepKeep, tolower(names(climateData.list2[[i]])))
  #     #}
  #     #else{
  #
  #     #}
  #     #keepCols<-grep("min|max|date|season", tolower(names(climateData.list2[[i]])))
  #     climateData.list2[[i]]<-climateData.list2[[i]][,keepCols]
  #     #change min column to just 'min' and do the same for 'max' column:
  #     minCol<-grep("min",tolower(names(climateData.list2[[i]])))
  #     maxCol<-grep("max",tolower(names(climateData.list2[[i]])))
  #     if(length(minCol)!=0 & length(maxCol)!=0){
  #       colnames(climateData.list2[[i]])[minCol] <- "min"
  #       colnames(climateData.list2[[i]])[maxCol] <- "max"
  #     }
  #     #convert all names to lowercase:
  #     names(climateData.list2[[i]])<-tolower(names(climateData.list2[[i]]))
  #   }
  #   #now we should be able to combine everything into one massive dataframe since col names are identical among nested dfs:
  #   climData<-do.call(rbind,climateData.list2)
  #   #browser()
  #   #run the data thru the getThresholds function:
  #   seasThresh.list<-createThresholds(data=climData,vars=vars)
  #   #browser()
  #   #get metadata for each site:
  #   siteMetaDates.list<-lapply(siteData.list, function(x) data.frame(startDate=x[1,1],endDate=x[nrow(x),1]))
  #   siteMeta.df<-do.call(rbind, siteMetaDates.list)
  #   siteMeta.df$station<-names(siteMetaDates.list)
  #   siteMeta.df$var<-vars
  #   #browser()
  #   keepThese<-unlist(keepThese)
  #   siteMetadata1.df<-do.call(rbind,lapply(nearbySites2[keepThese], function(x) x$location[1:3]))
  #   siteMeta.df<-cbind(siteMeta.df,siteMetadata1.df)
  #   #distance from POI (NEON site)
  #   siteMeta.df$distFromSite_m<-apply(siteMeta.df[,c("longitude_dec","latitude_dec")], 1,
  #           function(x) geosphere::distCosine(p1=data.frame(lonNEON,latNEON),p2=x))
  #   # #elevation deviation from site:
  #   siteMeta.df$elevDiff_m<-elevNEON-siteMeta.df$elev
  #   #browser()
  #   #put everything into a list
  #   out<-list(siteID,seasThresh.list,siteMeta.df)
  #   names(out)<-c("siteNEON","DataRangesThresholds","externalMetadata")
  # }
  # #write the file out:
  # #sitename<-gsub(".*:","",siteID)
  # #varsFile<-gsub("\\s+","", vars)
  # #filename<-paste0("C:/Users/jroberti/Git/NEON-FIU-document-IPT/dataQAQC/thresholds/",sitename,"_",varsFile,"_thresholds.rds")
  # saveRDS(out,filename)
  # #browser()
  # end=Sys.time()
  # message(difftime(end, strt, units = "secs"))
}

#set thresholds function:
createThresholds<-function(data,seasons=c("winter","spring","summer","fall"),vars){
  #convert data (if applicable)
  if(vars=="air temperature"){
    data<-suppressWarnings(degF2degC(df = data,x = c("min","max")))
  }
  else if(vars=="precipitation"){
    data<-suppressWarnings(inch2mm(df = data,x = "pcpn"))
    #browser()
  }
  #browser()
  #check if min and max columns exist:
  grepColMin<-grep("min",names(data))
  grepColMax<-grep("max",names(data))
  if(length(grepColMin)!=0 & length(grepColMax)!=0){
    seasRange.list.min<-lapply(seasons, function(x) min(data$min[which(data$season==x)],na.rm = T))
    seasRange.list.max<-lapply(seasons, function(x) max(data$max[which(data$season==x)],na.rm = T))
  }
  else{
    #find column that isn't date or season:
    dataCol<-which(!grepl("date|season",names(data)))
    #convert dataCol to as.numeric:
    if(class(data[,dataCol])=="factor"){
      data[,dataCol]<-as.numeric(levels(data[,dataCol]))[data[,dataCol]]
    }
    else{
      data[,dataCol]<-as.numeric(data[,dataCol])
    }
    seasRange.list.min<-lapply(seasons, function(x) min(data[which(data$season==x),dataCol],na.rm = T))
    seasRange.list.max<-lapply(seasons, function(x) max(data[which(data$season==x),dataCol],na.rm = T))
  }
  seasRange.list<-Map(cbind, seasRange.list.min, seasRange.list.max)

  #put into dataframe:
  seasRange.df<-data.frame(do.call(rbind,seasRange.list),row.names = NULL)
  #assign names:
  names(seasRange.df)<-c("min","max")
  #name the seasonal thresholds
  seasRange.df$season<-seasons
  #move seasons column to the first column:
  seasRange.df<-seasRange.df[,c("season","min","max")]
  #create seasonal thresholds from the ranges of the neighboring seasons:
  #Min thresholds:
  #browser()
  minRangeThresh.df<-seasRange.df[,c("season","min")]
  #sort seasons from min to max on min values:
  minRangeThresh.df<-minRangeThresh.df[with(minRangeThresh.df, order(min)), ]
  #make new vector with min values adjusted down 1 row:
  minThreshVals<-c(minRangeThresh.df$min[1]-(0.25*abs(minRangeThresh.df$min[1])),minRangeThresh.df$min[1:3])
  #minRangeThresh.df$min2<-minThreshVals
  #subtract seasonal min vs projected threshold:
  seasThreshDiff<-minRangeThresh.df$min-minThreshVals
  #check to make sure the variable isn't something like precip, where the min for all seasons should be 0
  if(all(seasThreshDiff!=0)){
    swapThese<-which(seasThreshDiff<=4)
    if(any(swapThese==1)){
      minThreshVals[1]<-minRangeThresh.df$min[1]-5
      swapWith<-minThreshVals[swapThese[-1]-1]
      minThreshVals[swapThese[-1]]<-swapWith
      #browser()
    }
    else{
    swapWith<-minThreshVals[swapThese-1]
    minThreshVals[swapThese]<-swapWith
    }
  }
  minRangeThresh.df$min<-minThreshVals
  #find threshold diffs that are <=4 degrees, swap these with threshold that gives the next largest threshold difference:
  # swapThese<-which(seasThreshDiff<=4)
  # swapWith<-minRangeThresh.df$min2[swapThese-1]
  # minRangeThresh.df$min2[swapThese]<-swapWith
  #browser()
  #check to make sure that the thresholds for a given season are reasonable:
  #
  # for(i in 1:nrow(minRangeThresh.df)){
  #   if(minRangeThresh.df$min-minRangeThresh.df$min2)
  # }
  #
  # minSeasonMin<-seasRange.df$season[which.min(seasRange.df$min)]
  # maxSeasonMin<-seasRange.df$season[which.max(seasRange.df$min)]
  # minRangeThresh.df$min[grep("winter",minRangeThresh.df$season)]<-seasRange.df$min[grep("winter",seasRange.df$season)]-(0.25*diff(range(minRangeThresh.df$min)))
  # minRangeThresh.df$min[grep("spring",minRangeThresh.df$season)]<-seasRange.df$min[grep("winter",seasRange.df$season)]
  # minRangeThresh.df$min[grep("summer",minRangeThresh.df$season)]<-seasRange.df$min[grep("spring",seasRange.df$season)]
  # minRangeThresh.df$min[grep("fall",minRangeThresh.df$season)]<-seasRange.df$min[grep("winter",seasRange.df$season)]
  #Max thresholds:
  #browser()
  maxRangeThresh.df<-seasRange.df[,c("season","max")]
  #sort seasons from max to min on max values:
  maxRangeThresh.df<-maxRangeThresh.df[with(maxRangeThresh.df, order(-max)), ]
  #make new vector with max values adjusted down 1 row:
  maxThreshVals<-c(maxRangeThresh.df$max[1]+(0.25*abs(maxRangeThresh.df$max[1])),maxRangeThresh.df$max[1:3])
  #maxRangeThresh.df$max2<-maxThreshVals
  #subtract seasonal min vs projected threshold:
  seasThreshDiff<-maxThreshVals-maxRangeThresh.df$max
  #find threshold diffs that are <=4 degrees, swap these with threshold that gives the next largest threshold difference:
  swapThese<-which(seasThreshDiff<=4)
  if(any(swapThese==1)){
    browser()
  }
  swapWith<-maxThreshVals[swapThese-1]
  maxThreshVals[swapThese]<-swapWith
  maxRangeThresh.df$max<-maxThreshVals

  # maxRangeThresh.df$max[grep("winter",maxRangeThresh.df$season)]<-seasRange.df$max[grep("fall",seasRange.df$season)]
  # maxRangeThresh.df$max[grep("spring",maxRangeThresh.df$season)]<-seasRange.df$max[grep("summer",seasRange.df$season)]
  # maxRangeThresh.df$max[grep("fall",maxRangeThresh.df$season)]<-seasRange.df$max[grep("summer",seasRange.df$season)]
  # maxRangeThresh.df$max[grep("summer",maxRangeThresh.df$season)]<-seasRange.df$max[grep("summer",seasRange.df$season)]+(0.25*diff(range(maxRangeThresh.df$max)))
  #merge the data:
  rangeThresholds.df<-merge(minRangeThresh.df,maxRangeThresh.df)
  returnThese<-list(data=data, seasonalRanges=seasRange.df,seasonalThresholds=rangeThresholds.df)
  return(returnThese)
}


### function to convert degF to Celsius (temperature)
degF2degC<-function(df,x){
  x<-paste0(x,collapse = "|")
  colInd<-grep(x,names(df))
  df[,colInd]<-apply(df[,colInd], 2, function(x) as.numeric(x))
  #browser()
  df[,colInd]<-round((df[,colInd] -32) * (5/9),2)
  return(df)
}

### function to convert inch to mm (precip)
inch2mm<-function(df,x){
  #browser()
  x<-paste0(x,collapse = "|")
  colInd<-grep(x,names(df))
  if(length(colInd)==1){
    #check class:
    if(class(df[,colInd])=="factor"){
      df[,colInd]<-as.numeric(levels(df[,colInd]))[df[,colInd]]
    }
    else{
      df[,colInd]<-as.numeric(df[,colInd])
    }
  }
  else{
    if(class(df[,colInd])=="factor"){
      df[,colInd]<-apply(df[,colInd], 2, function(x) as.numeric(levels(x))[x])
    }
    else{
      df[,colInd]<-apply(df[,colInd], 2, function(x) as.numeric(x))
    }
  }
  #browser()
  df[,colInd]<-round(df[,colInd]*25.4,2)
  return(df)
}

#run the function:
#defineClimateRanges()
