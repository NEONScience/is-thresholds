if(dir.exists("~/Github")){
  baseDir="~/Github/"
}else{
  baseDir="~/GitHub/"
}

#MAKE SURE TO BUILD THE metDownloadR package first before running!
setwd(paste0(baseDir, "is-thresholds/metDownloadR"))
library(devtools)
build()
install()

#move back to is-thresholds directory:
setwd(paste0(baseDir, "is-thresholds/"))

#source the define climate ranges function:
base::source(paste0(baseDir, "is-thresholds/scripts/defineClimateRanges.R"))


#source the define climate ranges function:
base::source(paste0(baseDir, "/NEON-FIU-document-IPT/dataQAQC/defineClimateRanges.R"))
base::source(paste0(baseDir,"/is-thresholds/scripts/getClimateData.R"))
load(paste0(baseDir,"is-thresholds/metDownloadR/data/resTracking.rda"))

#sort by ID:
#saveRDS(NEONsites,"C:/Users/jroberti/Git/NEON-FIU-document-IPT/dataQAQC/neonSiteList.rds")
NEONsites<-readRDS(paste0(baseDir,"is-thresholds/data/neonSiteList.rds"))
NEONsites<-NEONsites[with(NEONsites,order(Site.ID)),]
sites<-NEONsites$Site.ID

#startSite<-grep("YELL",sites)
#useSites<-sites[startSite:length(sites)]
saveDir=paste0(baseDir,"is-thresholds/colocated_data/")
#open directory to see what I've already run, and then run for sites that haven't been run yet:
temp<-list.files(saveDir)

#filter temp by variable:
temp<-temp[grep("all",temp)]
#sites that have already run:
alreadyRun<-gsub("\\_.*","",temp)

#runThese<-sites[which(sites %in% alreadyRun ==F)]
#useSites<-runThese[-c(1,2,3)]
#runThese=runThese[c(1:10)]

soloSites<-read.csv(file = paste0(baseDir, "is-thresholds/data/solo sites.txt"), stringsAsFactors = F, header = F)
#runThese=soloSites[24:37,1] # monolith's split of sites
#runThese=soloSites[1:18, 1] # rhlee12 CyVerse split of sites

clustSites<-read.csv(file = paste0(baseDir, "is-thresholds/data/sites_clusters.txt"), stringsAsFactors = F, header = F)
#runThese<-clustSites[18:nrow(clustSites),1]
#runThese<-clustSites[4:nrow(clustSites),1]

#identify all sites that need to be run between 'soloSites' and the first column of clustSites:
allSites<-c(soloSites$V1,clustSites$V1)
#identify sites that still need to be run:
runThese<-sort(allSites[!allSites %in% alreadyRun])

#Synoptic tokens:
#token="a2d2b292a8b74cb496060f09501204c6"   #ROBERT'S TOKEN
token="16088448e1b149509e45e401196106f0"  #JOSH'S TOKEN


for(i in 1:length(runThese[1:6])){
  siteCode<-paste0("NEON:",runThese[i])
  #siteCode<-paste0("NEON:","CPER")
  result=try(getClimateData(siteID=siteCode,
                                 radius=100,
                                 startDate="2010-01-01 00:00:00",
                                 endDate=as.character(Sys.time()),
                                 recrunchThresholds=T,
                                 overwriteData=T,
                                 numStations=3,
                                 save.dir = saveDir,
                                 token = token))
  if(!class(result)=="try-error"){
    message(paste0("Thresholds Generated for ", runThese[i]))
    # git2r::pull(repo = "~/GitHub/is-thresholds/", name= "NEONScience/is-thresholds")
    # git2r::add(repo = "~/GitHub/is-thresholds/", path = "./data/")
    # git2r::commit(repo = "~/GitHub/is-thresholds/", message = paste0(siteCocde, " threshold data"))
    # git2r::push()
    # 
  }else{message(paste0("Thresholds Generation failed for ", runThese[i], "!"))}
}

