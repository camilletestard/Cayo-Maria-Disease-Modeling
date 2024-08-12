
# X_Results.R ####

{
  
  library(tidyverse); library(cowplot); library(fs); library(magrittr)
  library(patchwork); library(ggforce); library(colorspace); library(ggregplot)
  library(MCMCglmm); library(INLA)
  
  theme_set(theme_cowplot())
  
  dir_create("Figures")
  
}

# Inequalities in disease risk were exacerbated after Hurricane Maria

IndivDF <- readRDS("Data/Intermediate/FullIndividualInfectionData.rds")

# IndivDF

IndivDF %<>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  mutate_at("Age", ~.x/10) %>% 
  # mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H"))) %>% 
  mutate_at("Rank", ~.x/100) %>% 
  mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
  na.omit %>% 
  group_by(ID, Pop, Sex, Rank, Age, Hurricane) %>% 
  summarise(MeanInf = mean(Infected), 
            MeanTime = mean(Time)) %>% 
  ungroup

Model1 <- readRDS("Data/Intermediate/IndividualInfectionModelAdd.rds")

Model2 <- readRDS("Data/Intermediate/IndividualInfectionModelAddStrength.rds")

sink("Results.txt")

# We found strong sex differences in infection risk, 
# which were strongly exacerbated following the hurricane (Figure 2A, 2D).
# Males were on average XXX steps slower to get infected than females
# in the years pre-hurricane (Mean, CI, P, Figure 2A, 2D).

Model1$FinalModel %>% GetEstimates("Sex")
Model1$FinalModel %>% INLAPValue("SexM")

# This sex difference increased following the hurricane, 
# driving a greater inequality between males and females (Mean, CI, P, Figure 2A, 2D). 

Model1$FinalModel %>% GetEstimates("HurricanePost:Sex")
Model1$FinalModel %>% INLAPValue("HurricanePost:SexM")

# In our initial models, age did not affect epidemic risk before or after the hurricane 
# (Mean, CI, P, Model Figure, Figure 2B), 

Model1$FinalModel %>% GetEstimates("Age")
Model1$FinalModel %>% INLAPValue("Age")

# and there was no interaction with hurricane status (Mean, CI, P, Model Figure); 

Model1$FinalModel %>% GetEstimates("HurricanePost:Age")
Model1$FinalModel %>% INLAPValue("HurricanePost:Age")

# however, when we included connection Strength in the model as an explanatory variable, 
# this interaction became significant, signifying that – although neither slope was significantly different from zero – 
# the slopes of epidemic risk against age were significantly different pre- and post-hurricane 
# (Mean, CI, P, Figure 2A, 2C). 

Model2$FinalModel %>% GetEstimates("HurricanePost:Age")
Model2$FinalModel %>% INLAPValue("HurricanePost:Age")

# In contrast, there were strong differences in epidemic risk among individuals of different ranks, which did not change following the hurricane: 
# lower-rank individuals were slower to get infected than 
# medium- (XXX, CI, P) and high-rank individuals (XXX, CI, P, Figure 2A-B).

# Model1$FinalModel %>% GetEstimates("RankM")
# Model1$FinalModel %>% INLAPValue("RankM")
# 
# Model1$FinalModel %>% GetEstimates("RankH")
# Model1$FinalModel %>% INLAPValue("RankH")

Model1$FinalModel %>% GetEstimates("Rank")
Model1$FinalModel %>% INLAPValue("Rank")

# Connection strength had a strong accelerating effect on infection timestep 
# when included as an explanatory variable (Figure 2A): 
# an extra full-strength connection was associated with XXX earlier infection (CI, P, Figure 2A). 

Model2$FinalModel %>% GetEstimates("Strength")
Model2$FinalModel %>% INLAPValue("Strength")

# However, its inclusion did not remove any effects (Figure 2A), 
# demonstrating that strength alone did not fully explain the hurricane’s effect (or any other effect). 
# Instead, as outlined above, it revealed variation in the slope of epidemic risk on age. 

# Hurricane Maria substantially redistributed epidemic risk across individuals

ModelList <- readRDS("Data/Intermediate/InfectionRepeatabilityModels.rds")

RepDF <- 
  ModelList %>% 
  map(MCMCRep) %>% 
  bind_rows(.id = "Model") %>% filter(Component != "units") %>% 
  mutate_at(c("Mode", "lHPD", "uHPD"), ~round(as.numeric(.x), 3))

# To further explore how the hurricane altered the distribution of epidemic risk
# across the population, we quantified between-individual repeatability of infection timestep; 
# see Figure S3. We found that individuals’ epidemic risks were repeatable across the study period (R=XX; 95% Credibility Intervals Lower-Upper), 

RepDF %>% 
  filter(Model == 1, Component == "ID") %>% 
  ReportRep

# but that fitting an ID:hurricane interaction effectively removed this repeatability (R=XX, Lower-Upper). 

ReportRep <- function(a) paste0(a$Mode, " (", a$lHPD, ", ", a$uHPD, ")")

RepDF %>% 
  filter(Model == 2, Component == "ID") %>% 
  ReportRep

# The variance accounted for by the ID:hurricane interaction itself was large (R=XX; Lower-Upper), 

RepDF %>% 
  filter(Model == 2, Component == "ID:Hurricane") %>% 
  ReportRep

# indicating that individuals’ rank order of network strength greatly changed after the hurricane compared to before. 
# That is, the same individuals were not necessarily high-risk before versus after the hurricane relative to the rest of the population. 
# Supporting this observation, the repeatability for the period after the hurricane was much higher than the overall period (R=XX) 

RepDF %>% 
  filter(Model == 4, Component == "ID") %>% 
  ReportRep

# while there was no detectable repeatability before (R=XX), 

RepDF %>% 
  filter(Model == 3, Component == "ID") %>% 
  ReportRep

# demonstrating that the hurricane’s shuffling of connection strengths had a profound impact on the distribution of epidemic risk across the population (Figure S3).

sink()
