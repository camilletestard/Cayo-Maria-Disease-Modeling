
# X_Repeatability_Observed.R ####

rm(list = ls())

library(dplyr); library(igraph); library(foreach); library(doParallel); library(stringr); library(tidyverse)
library(magrittr); library(fs)

# load("Data/R.Data/BisonFittedNetworks.RData")
# load("Data/BisonFittedNetworks (2).RData")
load("Greg Data/proximity_data.RData")

AggregatedEdges <- edgelist.all %>% #
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

dir_create("Greg Data/Outputs/Repeatability_Observed")

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
    rename(ID = 1, Strength = 2)
  
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

IDDF %>% saveRDS("Greg Data/Outputs/Repeatability_ObservedDF.rds")

library(lme4); library(lmerTest)

LM <- lmer(Strength ~ (1|ID) + Year + Population, data = IDDF)

library(rptR)

Rpt <- rpt(Strength ~ (1|ID) + Year + Population, data = IDDF, 
           grname = "ID",
           datatype = "Gaussian")

Rpt %>% summary

Rpt

Rpt2 <- rpt(Strength ~ (1|ID) + Year + Population + , data = IDDF, 
            grname = "ID",
            datatype = "Gaussian")

Rpt %>% summary

Rpt

# MCMC

library(MCMCglmm)

IDDF %<>% filter(Year > 2013)

MC1 <- MCMCglmm(Strength ~ Year + Population + PostMaria, data = IDDF, random =~ ID)

summary(MC1)

MC2 <- MCMCglmm(Strength ~ Year + Population + PostMaria, data = IDDF, random =~ ID + ID:PostMaria)

summary(MC2)

list(MC1, MC2) %>% map("DIC")

library(ggregplot)

MC1 %>% MCMCRep()

MC2 %>% MCMCRep()


MC3a <- MCMCglmm(Strength ~ Year + Population, # + PostMaria, 
                 data = IDDF %>% filter(PostMaria == 0), 
                 random =~ ID)

MC3b <- MCMCglmm(Strength ~ Year + Population, # + PostMaria, 
                 data = IDDF %>% filter(PostMaria == 1), 
                 random =~ ID)

list(MC1, MC2, MC3a, MC3b) %>% map(MCMCRep) %>% bind_rows(.id = "Model") %>% 
  mutate_at(3:5, ~round(as.numeric(.x), 3)) %>% 
  mutate_at("Model", ~c("ID Overall", "ID:Hurricane", "ID Before", "ID After")[as.numeric(.x)]) %>% 
  rename(Lower = lHPD, Upper = uHPD) %>% 
  write.csv("RepeatabilityValues.csv", row.names = F)

library(cowplot)

theme_set(theme_cowplot())

list(MC1, MC2, MC3a, MC3b) %>% map(MCMCRep) %>% 
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

ggsave("Figures/Repeatability.jpeg", units = "mm", height = 120, width = 150)

(FocalID <- IDDF$ID %>% unique %>% extract2(2))

(FocalID2 <- IDDF$ID %>% unique %>% extract(c(210, 219)))

OverlapDF <- 
  IDDF %>% group_by(ID) %>% summarise(MinYear = min(Year), MaxYear = max(Year)) %>% 
  filter(MinYear < 2020, MaxYear > 2020)

(FocalID2 <- c("8K1", "00V"))
(FocalID2 <- c("8K1", OverlapDF[10, ]$ID))

(FocalID2 <- c("8K1", OverlapDF[14, ]$ID))

(FocalID2 <- c("8K1", OverlapDF[15, ]$ID))

(FocalID2 <- c("8K1", OverlapDF[16, ]$ID))

MaintainedLines <- 
  IDDF %>% mutate_at("Year", as.numeric) %>% 
  ggplot(aes(Year, Strength+0.001)) +
  geom_line(aes(colour = as.factor(!ID %in% FocalID2),
                alpha = as.factor(ID %in% FocalID2),
                group = ID)) +
  geom_point() +
  scale_y_log10() + 
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.01, 1), name = paste(FocalID2, collapse = ", ")) +
  scale_colour_manual(values = c("red", "black"), name = paste(FocalID2, collapse = ", "))

(BrokenFocalID2 <- c("8K1", OverlapDF[17, ]$ID))

BrokenLines <-
  IDDF %>% mutate_at("Year", as.numeric) %>% 
  ggplot(aes(Year, Strength+0.001)) +
  geom_line(aes(colour = as.factor(!ID %in% BrokenFocalID2),
                alpha = as.factor(ID %in% BrokenFocalID2),
                group = ID)) +
  geom_point() +
  scale_y_log10() + 
  theme(legend.position = "none") +
  scale_alpha_manual(values = c(0.01, 1), name = paste(FocalID2, collapse = ", ")) +
  scale_colour_manual(values = c("red", "black"), name = paste(FocalID2, collapse = ", "))

library(patchwork)

MaintainedLines/BrokenLines

ggsave("Maintained_vs_Broken.jpeg", units = "mm", height = 180, width = 180)

