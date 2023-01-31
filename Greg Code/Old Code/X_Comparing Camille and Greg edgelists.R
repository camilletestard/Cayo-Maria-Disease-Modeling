
# X_Comparing Camille and Greg edgelists ####

library(tidyverse); library(ggregplot); library(magrittr)

MetaEdges <- readRDS("Greg Data/MetaEdges.rds")

MetaEdges %<>% # Cleaning the columns and names
  mutate_at("in.proximity", ~str_remove_all(.x, " ")) %>% 
  rename(From = focal.monkey, To = in.proximity)

AggregatedEdges <- # Getting the total observations for each dyad:rep
  MetaEdges %>% 
  group_by(From, To, Rep) %>% 
  summarise(Count = n())

Observations <-  # Getting the total observations for each ID:Rep
  AggregatedEdges %>% 
  pivot_longer(c("To", "From")) %>% 
  group_by(value, Rep) %>% 
  summarise_at("Count", ~sum(.x)) %>% 
  rename(Obs = Count)

AggregatedEdges %<>% # Joining the edge list with the total observations
  left_join(Observations, by = c("From" = "value", "Rep")) %>% 
  left_join(Observations, by = c("To" = "value", "Rep"), suffix = c(".From", ".To"))

AggregatedEdges %<>% mutate(Weight = Count/(Obs.From + Obs.To - Count))

AggregatedEdges$Weight %>% qplot

AggregatedEdges %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate(Group = substr(Rep, 1, 1), 
         Year = substr(Rep, str_count(Rep) - 3, str_count(Rep)))

AggregatedEdges %<>% mutate(Post = Year >= 2018)

# Importing Camille version ####

load("Data/R.Data/proximity_data.RData")

edgelist.all %>% nrow

edgelist.all %>% dplyr::select(ID1, ID2, Year = year, Group = group) %>% unique

AggregatedEdges %>% nrow
