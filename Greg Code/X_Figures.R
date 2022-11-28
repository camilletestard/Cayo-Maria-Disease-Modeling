
# X_Figures ####

library(tidyverse); library(cowplot)

theme_set(theme_cowplot())

dir_create("Figures")

EdgeListList <- readRDS("Greg Data/TimeEdges.rds")

EdgeListList %<>%
  separate(DateTime, sep = " ", into = c("Date", "Time"))

EdgeListList %<>% mutate_at("Date", ~lubridate::ymd(.x))

EdgeListList %<>% 
  mutate(Year = str_split(Date, "-") %>% map_chr(1))

EdgeListList %>% 
  count(Rep, Date, Year) %>% 
  ggplot(aes(Rep, Date)) + 
  geom_tile(aes(fill = n)) + 
  coord_flip() +
  facet_wrap(~Year, scales = "free")

DiffDF <- 
  EdgeListList %>% 
  group_by(Rep, Date, Year) %>% 
  count %>% group_by(Rep, Year) %>% 
  mutate(DateGap = c(0, diff(as.numeric(Date))))

DiffDF %>% 
  filter(DateGap > 0) %>% 
  # mutate_at("Gap", ~(.x + 1)) %>% 
  SinaGraph("Rep", "DateGap", Just = T) + 
  scale_y_log10()

# Plotting temporal windows ####

Reps <- EdgeListList$Rep %>% unique %>% sort

FocalRep <- Reps[1]

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

library(ggregplot); library(ggraph); library(patchwork)

a <- GraphList[[1]]

Samples <- 
  sample(1:length(GraphList), 10) %>% 
  sort

GraphList[Samples] %>% 
  map(function(a){
    
    a %>% 
      ggraph("stress") + 
      geom_edge_link0() +
      geom_node_point() +
      coord_fixed()
    
  }) %>% ArrangeCowplot()
