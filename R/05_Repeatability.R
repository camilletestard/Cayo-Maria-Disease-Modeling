
# X_Repeatability_Observed.R ####

rm(list = ls())

library(dplyr); library(igraph); library(foreach); library(doParallel); library(stringr); library(tidyverse)
library(magrittr); library(fs);
library(MCMCglmm)

# load("Data/R.Data/BisonFittedNetworks.RData")
# load("Data/BisonFittedNetworks (2).RData")

# Importing data ####

load("Data/Intermediate/proximity_data.RData")

AggregatedEdges <- 
  edgelist.all %>% #
  # bind_rows(.id = "Rep")
  mutate_at("group", ~substr(.x, 1, 1)) %>% 
  rename(Pop = group, Year = year) %>% 
  mutate(Rep = paste0(Pop, Year))

AggregatedEdges %<>% rename(From = 1, To = 2)

AggregatedEdges %<>% mutate(Weight = count/total_samples)

AggregatedEdges %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate(Group = substr(Rep, 1, 1), 
         Year = substr(Rep, str_count(Rep) - 3, str_count(Rep)))

Reps <- AggregatedEdges$Rep %>% unique %>% sort

dir_create("Data/Outputs/Repeatability_Observed")

# Calculating repeatability metrics ####

reps = 1000 # number of times the simulation should be repeated

sims = 10000 # number of time steps/times each dyad should be allowed to potentially interact

MeanInf <- 0.15

InfSD <- 0.04

FocalRep <- Reps[1]

StrengthList <- list()

for(FocalRep in Reps){
  
  print(FocalRep)
  
  RepData <- AggregatedEdges %>% filter(Rep == FocalRep)
  
  Network <- 
    AdjMatrix <-
    RepData[,c("From", "To", "Weight")] %>% 
    arrange(From, To) %>% 
    graph.data.frame(directed = F) %>% 
    get.adjacency(sparse = FALSE, attr = 'Weight')
  
  StrengthDF <-
    Network %>% rowSums %>% data.frame %>% 
    rownames_to_column %>%
    rename(ID = 1, Strength = 2) %>% 
    bind_cols(
      
      Network %>% 
        DegreeGet %>% 
        rowSums %>% data.frame %>% 
        rownames_to_column %>%
        rename(ID = 1, Degree = 2) %>% 
        dplyr::select(Degree)
      
    )
  
  # }) %>% bind_cols()
  
  # colnames(StrengthDF) <- paste0("draw.", 1:1000)
  
  StrengthList[[FocalRep]] <- StrengthDF
  
}

IDDF <- 
  StrengthList %>% bind_rows(.id = "Rep") %>% 
  mutate(Population = substr(Rep, 1, 1), 
         Year = substr(Rep, 2, 5))

IDDF %<>%
  mutate(PostMaria = as.factor(as.numeric(Year %in% c(2018:2021))))

IDDF %<>% mutate_at("Year", as.numeric)

IDDF %>% saveRDS("Data/Outputs/Repeatability_ObservedDF.rds")

# Repeatability Models ####

IDDF %<>% filter(Year > 2013)

MC1 <- MCMCglmm(Strength ~ Year + Population + PostMaria, data = IDDF, random =~ ID)

MC2 <- MCMCglmm(Strength ~ Year + Population + PostMaria, data = IDDF, random =~ ID + ID:PostMaria)

MC3a <- MCMCglmm(Strength ~ Year + Population, # + PostMaria, 
                 data = IDDF %>% filter(PostMaria == 0), 
                 random =~ ID)

MC3b <- MCMCglmm(Strength ~ Year + Population, # + PostMaria, 
                 data = IDDF %>% filter(PostMaria == 1), 
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% saveRDS("Data/Outputs/RepeatabilityModels.rds")

ModelList %>% map(MCMCRep) %>% bind_rows(.id = "Model") %>% 
  mutate_at(3:5, ~round(as.numeric(.x), 3)) %>% 
  mutate_at("Model", ~c("ID Overall", "ID:Hurricane", "ID Before", "ID After")[as.numeric(.x)]) %>% 
  rename(Lower = lHPD, Upper = uHPD) %>% 
  write.csv("Data/Outputs/RepeatabilityValues.csv", row.names = F)

ModelList %>% 
  map(MCMCRep) %>% 
  bind_rows(.id = "Model") %>% filter(Component != "units") %>% 
  # mutate_at(3:5, ~round(as.numeric(.x), 2))
  mutate_at(3:5, ~as.numeric(.x)) %>% 
  ggplot(aes(factor(Model), Mode, 
             fill = Component)) + 
  geom_col(colour = "black", position = "stack") + 
  labs(x = "Model") + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  lims(y = c(NA, 1)) +
  geom_hline(yintercept = 1, lty = 2, alpha = 0.4) +
  scale_x_discrete(labels = c("Overall", "ID:Hurricane", "ID Before", "ID After"))

MC1 <- MCMCglmm(Degree ~ Year + Population + PostMaria, data = IDDF, random =~ ID)

MC2 <- MCMCglmm(Degree ~ Year + Population + PostMaria, data = IDDF, random =~ ID + ID:PostMaria)

MC3a <- MCMCglmm(Degree ~ Year + Population, # + PostMaria, 
                 data = IDDF %>% filter(PostMaria == 0), 
                 random =~ ID)

MC3b <- MCMCglmm(Degree ~ Year + Population, # + PostMaria, 
                 data = IDDF %>% filter(PostMaria == 1), 
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% 
  map(MCMCRep) %>% 
  bind_rows(.id = "Model") %>% filter(Component != "units") %>% 
  # mutate_at(3:5, ~round(as.numeric(.x), 2))
  mutate_at(3:5, ~as.numeric(.x)) %>% 
  ggplot(aes(factor(Model), Mode, 
             fill = Component)) + 
  geom_col(colour = "black", position = "stack") + 
  labs(x = "Model") + 
  theme(axis.text.x = element_text(angle = 45, 
                                   hjust = 1)) +
  lims(y = c(NA, 1)) +
  geom_hline(yintercept = 1, lty = 2, alpha = 0.4) +
  scale_x_discrete(labels = c("Overall", "ID:Hurricane", "ID Before", "ID After"))

# Running raw social models ####

load("Data/Intermediate/proximity_data.RData")

IndividualTraits <-
  edgelist.all %>%
  mutate(IsPost = ifelse(year < 2018, "Pre", "Post")) %>%
  mutate(Pop = paste0(group, year)) %>%
  dplyr::select(ID = ID1, Sex = ID1_sex, Rank = ID1_rank, Age = ID1_age, Pop, IsPost) %>%
  unique

IndividualTraits %<>% mutate_at("Sex", ~str_replace_all(.x, "FALSE", "F"))

TestDF <-
  IndividualTraits %>% 
  left_join(IDDF, by = c("ID", "Pop" = "Rep")) %>% 
  na.omit

TestDF %<>%
  filter(Age > 5) %>% 
  mutate_at("Age", ~.x/10)

TestDF %<>% mutate_at(c("Strength", "Degree"), ~log(.x + 1))

mdl.strength <- inla(Strength ~ IsPost * (Age + Sex + Rank) + 
                       f(ID, model = "iid") + f(Pop, model = "iid"),
                     # family = "poisson",
                     control.compute = list(dic = TRUE),
                     data = TestDF)

mdl.strength <- lmer(Strength ~ IsPost * (Age + Sex + Rank) + 
                       (1|ID) + (1|Pop),
                     # family = "poisson",
                     # control.compute = list(dic = TRUE),
                     data = TestDF)

mdl.strength %>% resid %>% qplot

qplot(fitted(mdl.strength), resid(mdl.strength))

mdl.degree <- inla(Degree ~ IsPost * (Age + Sex + Rank) + 
                     f(ID, model = "iid") + f(Pop, model = "iid"),
                   # family = "poisson",
                   control.compute = list(dic = TRUE),
                   data = TestDF)

mdl.degree <- lmer(Degree ~ IsPost * (Age + Sex + Rank) + 
                     (1|ID) + (1|Pop),
                   # family = "poisson",
                   # control.compute = list(dic = TRUE),
                   data = TestDF)

mdl.degree %>% resid %>% qplot

qplot(fitted(mdl.degree), resid(mdl.degree))

mdl.degree %>% plot



