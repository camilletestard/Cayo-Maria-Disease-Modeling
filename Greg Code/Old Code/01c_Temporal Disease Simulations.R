
# 01c_Temporal Disease Simulations ####

library(dplyr)
library(igraph)
library(foreach)
library(doParallel)
library(stringr)

library(tidyverse); library(magrittr); library(fs)

EdgeListList <- readRDS("Greg Data/TimeEdges.rds")

EdgeListList %<>% 
  separate(DateTime, sep = " ", into = c("Date", "Time"))

EdgeListList %<>% mutate_at("Date", ~lubridate::ymd(.x))

EdgeListList %<>% 
  filter(!(str_detect(Rep, "K2018") & Date == "2017-01-18"))

pinf = 1
pinf = 0.5

Window <- 30
reps <- 1000

Reps <- EdgeListList$Rep %>% unique %>% sort

FocalRep <- Reps[1]

dir_create("Greg Data/Outputs/Temporal")

for(FocalRep in Reps[str_detect(Reps, "2018")]){
  
  SubEdgeList <- EdgeListList %>% filter(Rep == FocalRep)
  
  SubEdgeList %<>% mutate(NDate = Date - min(Date) + 1)
  
  FullDays <- max(SubEdgeList$Date) - min(SubEdgeList$Date) - Window
  
  library(tidygraph)
  
  BlankIndivs <- 
    SubEdgeList[,c("From", "To")] %>% unlist %>% unique %>% 
    data.frame(ID = .)
  
  GraphList <- 
    0:FullDays %>% 
    map(function(Day){
      
      SubSubEdgeList <- 
        SubEdgeList %>% 
        filter(NDate %in% (1:Window + Day)) %>% 
        dplyr::select(From, To) %>% 
        # graph_from_data_frame %>% #as_tbl_graph()
        tbl_graph(edges = ., 
                  nodes = BlankIndivs) %>% 
        get.adjacency(sparse = F)
      
    })
  
  r <- 1
  
  IndivList <- list()
  
  for (r in r:reps){
    
    t1 = Sys.time()
    
    print(r)
    
    N = nrow(BlankIndivs)
    
    Health = rep(0, N)
    
    Health[sample(1:length(Health), 1)] = 1
    
    Indivs <- data.frame(ID = BlankIndivs$ID, 
                         Infected = Health,
                         Time = Health - 1)
    
    s <- 1
    
    while(s < FullDays){
      
      Network <- GraphList[[s]]
      
      diag(Network) <- 0
      
      TransmissionMatrix <- array(0, dim = dim(Network))
      
      I2 <- which(Indivs$Infected > 0) # Identifying infected
      
      NI2 <- setdiff(1:nrow(Indivs), 
                     which(Indivs$Infected > 0)) # Identifying uninfected
      
      if(0){ # Removing resistant individuals?
        
        NI2 <- setdiff(NI2, Indivs[Indivs$Time > (-1), "ID"])
        
      }
      
      if(0){
        
        TransmissionMatrix[I2,] <- 
          
          # rbinom(length(Network[I2,]), 1, Network[I2,])* # Identifying if they interact
          
          as.numeric(runif(length(Network[I2,]), 0, 1) < pinf) # Identifying if they infect
        
        Infected <- which(colSums(TransmissionMatrix) > 0)# %>% as.numeric()
        
        NewlyInfected <- setdiff(Infected, 
                                 which(Indivs$Infected == 1))
        
      }
      
      if(1){
        
        NPairs <- which(Network[I2,]>0)
        
        TransmissionMatrix[I2,][NPairs] <- 
          
          # rbinom(length(NPairs), 1, Network[I2,][NPairs])* # Identifying if they interact
          
          as.numeric(runif(length(NPairs), 0, 1) < pinf) # Identifying if they infect
        
        Infected <- which(colSums(TransmissionMatrix) > 0)# %>% as.numeric()
        
        NewlyInfected <- setdiff(Infected, 
                                 which(Indivs$Infected == 1))
        
      }
      
      if(length(NewlyInfected) > 0){
        
        Indivs[NewlyInfected, "Time"] <- s
        Indivs[NewlyInfected, "Infected"] <- 1
        
      }
      
      Indivs[Indivs$Time < (s - Window), "Infected"] <- 0
      
      s <- s + 1
      
    }
    
    IndivList[[r]] <- Indivs
    
    print(Sys.time() - t1)
    
  }
  
  saveRDS(IndivList, file = paste("Greg Data/Outputs/Temporal/PI", pinf, "Window", Window,  FocalRep, ".rds", sep = "_"))
  
}
