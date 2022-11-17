
# 01c_Temporal Disease Simulations ####

library(dplyr)
library(igraph)
library(foreach)
library(doParallel)
library(stringr)

library(magrittr); library(fs)
EdgeListList <- readRDS("Greg Data/TimeEdges.rds")

EdgeListList %<>% 
  separate(DateTime, sep = " ", into = c("Date", "Time"))

EdgeListList %<>% mutate_at("Date", ~lubridate::ymd(.x))

Window <- 30

Reps <- EdgeListList$Rep %>% unique %>% sort

FocalRep <- Reps[1]

SubEdgeList <- EdgeListList %>% filter(Rep == FocalRep)

SubEdgeList %<>% mutate(NDate = Date - min(Date) + 1)

FullDays <- max(SubEdgeList$Date) - min(SubEdgeList$Date) - Window

library(tidygraph)

GraphList <- 
  0:FullDays %>% 
  map(function(Day){
    
    SubSubEdgeList <- 
      SubEdgeList %>% 
      filter(NDate %in% (1:Window + Day)) %>% 
      dplyr::select(From, To) %>% 
      graph_from_data_frame %>% #as_tbl_graph()
      get.adjacency(sparse = F)
    
  })
