
# 00b_Temporal Data Cleaning.R

library(stringr)
library(igraph)
library(lubridate)
library(hms)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

library(magrittr); library(ggregplot)

#Load functions

group = c("F","V","KK","V","F","F","KK","V","V","KK","S","V","F","V")
years = c(2015,2015,2015,
          2016,2016,2017,2017,2017,
          2018, 2018,2019, 2019,2021,2021)

groupyears = paste0(group, years)

gy <- 1 # 2 # 16

edgelist.all = data.frame()

savePath = 'Greg Data/R.Data/'

library(fs)

dir_create("Greg Data/R.Data")

InputRoot <- "Data/Data All Cleaned/BehavioralDataFiles/"

EdgeListList <- list()

for (gy in 1:length(groupyears)){ #for all group & years
  
  print(paste("%%%%%%%%%%%%%%%%%% ",groupyears[gy], "%%%%%%%%%%%%%%%%%%"))
  
  if (years[gy] == 2018) {
    
    scans2018 <- read.csv(paste(InputRoot, 
                                "Group",groupyears[gy],"_scansamples_FULL_CLEANED.txt", sep = ""))
    
    meta_data <- read.csv(paste(InputRoot, 
                                "Group",groupyears[gy],"_GroupByYear.txt", sep = "")) #load meta data
    
    #Initialize SocialCaiptalData
    
    SocialCapitalData <- 
      meta_data[,c("id","sex","age","ordinal.rank","percofsex.dominanted","hrs.focalfollowed","focalcutoff_met")] %>% 
      mutate(group = group[gy],
             year = years[gy])
    
    ##Format data##
    
    scans2018 %<>% 
      mutate(date = lubridate::mdy(as.character(date))) %>% 
      mutate(year = lubridate::year(date),
             Q = lubridate::quarter(date),
             date = as.character(date)) %>%  #re-format to character after finding year and quarter
      mutate(#Add unique scan identifier
        observation.name = as.factor(paste(date, scan.num,sep="."))) %>% 
      mutate(#Add hurricane info
        isPost = 1)
    
    #Format time and create timeBlock column
    #IMPORTANT: MAKE SURE TIME COLUMNS ARE FORMATTED IN EXCEL IN FORMAT "13:55:00"
    scans2018$start.time = as_hms(as.character(scans2018$start.time))
    scans2018$stop.time = as_hms(as.character(scans2018$stop.time))
    
    scans2018$timeBlock = NA
    scans2018$timeBlock[which(scans2018$start.time <= as_hms("11:00:00"))] = "AM";
    scans2018$timeBlock[which(scans2018$start.time > as_hms("11:00:00"))] = "PM";
    
    #Format XEX names
    unique(scans2018$subject.ID) #check spelling of subject id
    scans2018$subject.ID=sub("'","",as.character(scans2018$subject.ID)) #Replace 'XEX byXEX names if needed
    scans2018$subject.ID=str_trim(scans2018$subject.ID,side="both") #Remove blanks
    unique(scans2018$prox.adult.IDs)
    scans2018$prox.adult.IDs=sub("'","",as.character(scans2018$prox.adult.IDs))
    
    #Clean up: rename and delete unused columns
    scans2018[,c("stop.time","observer.initials","cayo.map.code","nearest.adult.neighbour.ID","distance.nearest.neighbour","start.time")]=NULL
    names(scans2018)[3]="scan.number"; names(scans2018)[4]="focalID"; names(scans2018)[5]="focal.activity"; names(scans2018)[7]="in.proximity"
    
    #create column with equivalent activity to pre-hurricane
    scans2018[] = lapply(scans2018,str_trim)
    scans2018$focal.activity.isPost = as.character(scans2018$focal.activity) #preserve post-hurricane activity code
    #re-code activity to pre-hurricane for comparison
    scans2018$focal.activity = as.character(scans2018$focal.activity)
    scans2018$focal.activity[unique(c(which(scans2018$focal.activity=='G'),which(scans2018$focal.activity=='E'),which(scans2018$focal.activity=='E,P'),which(scans2018$focal.activity=='G,E')))]="social"
    scans2018$focal.activity[unique(c(which(scans2018$focal.activity=='R'),which(scans2018$focal.activity=='P')))]="rest"
    scans2018$focal.activity[unique(c(which(scans2018$focal.activity=='AG'),which(scans2018$focal.activity=='AR')))]="aggression"
    scans2018$focal.activity[unique(c(which(scans2018$focal.activity=='SR'),which(scans2018$focal.activity=='SG')))]="submit"
    scans2018$focal.activity[grep('T',scans2018$focal.activity)]="travel"
    scans2018$focal.activity[grep('F',scans2018$focal.activity)]="feed"
    scans2018$focal.activity[grep('D',scans2018$focal.activity)]="drink"
    scans2018$focal.activity[grep('SD',scans2018$focal.activity)]="sdb"
    scans2018$focal.activity[grep('N/A',scans2018$focal.activity)]="UNK"
    #unique(scans2018$focal.activity) #Check correct activity categories
    
    #Format in.proximity, count number of prox partners
    scans2018$partner.ID = as.character(scans2018$partner.ID); scans2018$partner.ID[which(scans2018$partner.ID=="N/A")]=NA
    scans2018$in.proximity = as.character(scans2018$in.proximity); scans2018$in.proximity[which(scans2018$in.proximity=="N/A")]=NA
    scans2018$num.prox = str_count(as.character(scans2018$in.proximity),",")+1
    scans2018$num.prox[is.na(scans2018$num.prox)]=0
    
    #Add social information
    scans2018$isProx=1; scans2018$isProx[which(scans2018$num.prox==0)]=0
    scans2018$isSocial=0; scans2018$isSocial[which(scans2018$focal.activity=="social")]=1
    scans2018$isSocialGive = 0; scans2018$isSocialGive[which(scans2018$focal.activity.isPost=="G")]=1
    scans2018$isSocialGet = 0; scans2018$isSocialGet[which(scans2018$focal.activity.isPost=="E")]=1
    scans2018$isAgg=0; scans2018$isAgg[which(scans2018$focal.activity=="aggression" | scans2018$focal.activity=="submit")]=1
    
    #Order columns
    col_order <- c("date","observation.name","focalID","group","year","scan.number","focal.activity","focal.activity.isPost","partner.ID","in.proximity","num.prox","isProx","isSocial","isSocialGive", "isSocialGet","isAgg", "Q","isPost","timeBlock")
    scans2018 <- scans2018[, col_order]
    prox_data =  scans2018#[,c("focalID","in.proximity")]
    
    prox_data %<>% rename(focal.monkey = focalID) %>% mutate(time = NA)
    
    # prox_data$in.proximity[prox_idx] <- paste(prox_data$focal.monkey[prox_idx], 
    #                                           prox_data$in.proximity[prox_idx], sep = ",")
    
  }else{ #if not 2018 (i.e. regular focal data)
    ############################################################################ 
    
    #Load data
    
    prox_data <- read.csv(paste0(InputRoot, "Group", groupyears[gy], "_ProximityGroups.txt"))
    meta_data <- read.csv(paste0(InputRoot, "Group", groupyears[gy], "_GroupByYear.txt")) #load meta data
    
    #cleaned_data = read.csv(paste("Group",groupyears[gy],"_CleanedData.txt", sep = ""))
    
    prox_data$time %<>% lubridate::as_datetime()
    
    #Set NA focal activity in scans to "Rest"
    
    prox_data %<>% mutate_at("focal.activity", ~replace_na("rest"))
    
    if (groupyears[gy] == "HH2016"){
      
      #Add HH dominance for subadults. 
      hh.dominance <- read.csv("HH_Dominance.csv");names(hh.dominance)[1]="id"
      hh.dominance$id = as.character(hh.dominance$id)
      hh.dominance$id[hh.dominance$id=="2.00E+09"]="2E9"
      hh.dominance$id[hh.dominance$id=="2.00E+08"]="2E8"
      hh.dominance$id[hh.dominance$id=="7.00E+00"]="7E0"
      hh.dominance$id[hh.dominance$id=="7.00E+03"]="7E3"
      hh.dominance$id[hh.dominance$id=="8.00E+02"]="8E2"
      
      meta_data[,c("ordinal.rank","percofsex.dominanted")]=hh.dominance[match(meta_data$id, hh.dominance$id),c("ordinal.rank","percofsex.domianted")]
      
    }
    
    if (groupyears[gy] == "KK2017"){
      
      kk.dominance <- # Add HH dominance for juveniles
        read.csv(paste0("Data/Data All Cleaned/Dominance/", "KK_dominance_withSubadults.csv")) %>% 
        rename(id = 1)
      
      meta_data[,c("ordinal.rank","percofsex.dominanted")] = 
        kk.dominance[match(meta_data$id, kk.dominance$id),
                     c("ordinal.rank","percofsex.domianted")]
      
    }
    
    if (groupyears[gy] == "V2019"){ #quick fix for now
      
      prox_idx = !is.na(prox_data$partners.activity..sequential.) #find indices where there are individuals in proximity
      
      #add focal monkey in proximity column
      
      prox_data$in.proximity[prox_idx] <- paste(prox_data$focal.monkey[prox_idx], 
                                                prox_data$in.proximity[prox_idx], sep = ",")
      
    }
    
  } #end of year clause (2018 vs. other)
  
  # #Format data with aggregate format
  # # Output the Master Edgelist of all possible pairs given the unique IDs.
  # unqIDs = meta_data$id#[meta_data$focalcutoff_met=="Y"]
  # edgelist = calcMasterEL(unqIDs);
  # df_obs_agg  = calcEdgeList(prox_data, edgelist); ##IMPORTANT NOTE: the structure of the proximity data is different for 
  # #2018 relative to other years (the focal ID is not counted in proximity for 2018 but it is for other years)
  # #This is not a problem because the focal ID column name is not the same so the function works properly as is.
  # names(df_obs_agg)=c("ID1", "ID2", "dyad_id","count")
  # 
  # library(magrittr)
  # 
  # df_obs_agg %<>% arrange(ego, alter)
  # 
  # prox_data
  
  prox_data
  
  library(purrr)
  
  # BlankEdgeList <- 
  #   expand.grid(focal.monkey = meta_data$id, in.proximity = meta_data$id)
  
  Edges <- 
    prox_data %>% 
    dplyr::select(From = focal.monkey, To = in.proximity, DateTime = time) %>% 
    mutate(Year = str_split(DateTime, "-") %>% map_chr(first)) %>% 
    mutate_at("To", ~str_split(.x %>% str_remove_all(" "), ",")) %>% 
    unnest(To) %>% 
    filter(From != To) %>% 
    arrange(DateTime)
  
  EdgeListList[[gy]] <- Edges
  
}

names(EdgeListList) <- groupyears

EdgeListList %<>% bind_rows(.id = "Rep")

EdgeListList %>% saveRDS("Greg Data/MetaEdges.rds")

