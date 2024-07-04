
# 02b_BISoN Individual Summarising ####

{
  
  library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
  library(magrittr); library(colorspace); library(lme4); library(lmerTest); library(patchwork)
  library(INLA); library(MCMCglmm)
  
  theme_set(theme_cowplot())
  
  dir_create("Intermediate/IndividualINLA")
  
  # Import the population timecourses ####
  
  FileList <- 
    "Data/Outputs/BISoN" %>% 
    dir_ls()
  
  names(FileList) <- "Data/Outputs/BISoN" %>%
    list.files
  
  IndivDFList <- 
    FileList %>% str_split("/") %>% map_chr(last) %>% str_split("_") %>% 
    map_chr(first) %>% unique %>% sort %>% 
    map(function(Pop){
      
      print(Pop)
      
      FileList[str_detect(FileList, Pop)] %>% 
        
        map(function(a){
          
          print(a)
          
          b <- a %>% readRDS
          
          # b %>% arrange(ID)
          
          b %>% return
          
        }) %>% bind_rows(.id = "File") %>% mutate(Pop = Pop)
      
    })
  
  names(IndivDFList) <-
    FileList %>% str_split("/") %>% map_chr(last) %>% str_split("_") %>% 
    map_chr(first) %>% unique %>% sort
  
  IndivDF <- IndivDFList %>% bind_rows(.id = "Rep")
  
  # Attaching individual traits ####
  
  load("Data/Intermediate/proximity_data.RData")
  
  IndividualTraits <-
    edgelist.all %>%
    mutate(Hurricane = factor(ifelse(year < 2018, "Pre", "Post"), levels = c("Pre", "Post"))) %>%
    mutate(Pop = paste0(group, year)) %>%
    dplyr::select(ID = ID1, Sex = ID1_sex, Rank = ID1_rank, Age = ID1_age, Pop, Hurricane) %>%
    unique
  
  IndivDF %<>% left_join(IndividualTraits)
  
  IndivDF2 <- 
    IndivDF %>% 
    mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
    filter(Age > 5) %>% 
    filter(Time > 0) %>% 
    mutate_at("Age", ~.x/10) %>% 
    mutate(S_I = str_split(File, "_") %>% map_chr(last) %>% str_remove(".rds") %>% as.numeric %>% multiply_by(10)) %>%
    na.omit
  
}

# IndivDF %>% saveRDS("Data/Intermediate/FullIndividualInfectionData.rds")

IndivDF <- readRDS("Data/Intermediate/FullIndividualInfectionData.rds")

# Analysing averaged infection timestep ####

TestDF <- 
  IndivDF %>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  mutate_at("Age", ~.x/10) %>% 
  mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H"))) %>% 
  mutate_at("Hurricane", ~factor(.x, levels = c("Pre", "Post"))) %>% 
  na.omit %>% 
  group_by(ID, Pop, Sex, Rank, Age, Hurricane) %>% 
  summarise(MeanInf = mean(Infected), 
            MeanTime = mean(Time),
            TimeSD = sd(Time),
            TimeSE = TimeSD/(1000^0.5)) %>% 
  ungroup

TestDF %<>% mutate_at("MeanTime", round)

INLAGaussian <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                       f(ID, model = "iid") + f(Pop, model = "iid"),
                     # family = "poisson",
                     control.compute = list(dic = TRUE),
                     data = TestDF)

INLALogGaussian <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                          f(ID, model = "iid") + f(Pop, model = "iid"),
                        # family = "poisson",
                        family = "lognormal",
                        control.compute = list(dic = TRUE),
                        data = TestDF)

INLACount <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "poisson",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLAGamma <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "gamma",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLANB <- inla(MeanTime ~ Hurricane * (Age + Sex + Rank) + 
                 f(ID, model = "iid") + f(Pop, model = "iid"),
               family = "nbinomial",
               control.compute = list(dic = TRUE),
               data = TestDF)

list(#INLAGaussian, 
  INLALogGaussian, INLAGamma, INLACount, INLANB) %>% Efxplot(Intercept = F)

list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB) %>% MDIC

list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB) %>% 
  saveRDS("Data/Intermediate/IndividualInfectionModelList.rds")

IM1 <- INLAModelAdd(Data = TestDF, 
                    Family = "poisson",
                    Response = "MeanTime",
                    Explanatory = c("Hurricane", "Age", "Sex", "Rank"), 
                    Add = paste0("Hurricane:", c("Age", "Sex", "Rank")),
                    AllModels = T,
                    Random = c("ID", "Pop"), RandomModel = "iid"
)

IM1 %>% saveRDS("Data/Intermediate/IndividualInfectionModelAdd.rds")

IM1$FinalModel %>% Efxplot(Intercept = F)

IM1$AllModels[[2]] %>% Efxplot(Intercept = F)

IM1$dDIC

# Repeatability (alternative formulation) ####

TestDF %<>% data.frame

TestDF %<>% 
  mutate(Year = substr(Pop, str_count(Pop) - 3, str_count(Pop))) %>% 
  mutate(Pop = substr(Pop, 1, 1))

MC1 <- MCMCglmm(MeanTime ~ Year + Pop, # + Hurricane, 
                data = TestDF, 
                # family = "poisson",
                random =~ ID)

MC2 <- MCMCglmm(MeanTime ~ Year + Pop, # + Hurricane, 
                data = TestDF, 
                # family = "poisson",
                random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(MeanTime ~ Year + Pop, 
                 data = TestDF %>% filter(Hurricane == "Pre"), 
                 # family = "poisson",
                 random =~ ID)

MC3b <- MCMCglmm(MeanTime ~ Year + Pop, 
                 data = TestDF %>% filter(Hurricane == "Post"), 
                 # family = "poisson",
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% saveRDS("Data/Intermediate/InfectionRepeatabilityModels.rds")

ModelList %>% map(~MCMCRep(.x)) %>% bind_rows(.id = "Model") %>% 
  mutate_at(3:5, ~round(as.numeric(.x), 3)) %>% 
  mutate_at("Model", ~c("ID Overall", "ID:Hurricane", "ID Before", "ID After")[as.numeric(.x)]) %>% 
  rename(Lower = lHPD, Upper = uHPD)# %>% 
# write.csv("Data/Outputs/InfectionRepeatabilityValues.csv", row.names = F)

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

# MCMC Repeatability Model ####

TestDF %<>% data.frame

TestDF %<>% mutate_at("Pop", ~substr(.x, 1, 1))

MC1 <- MCMCglmm(MeanTime ~ Hurricane * (Age + Sex + Rank), 
                data = TestDF, 
                family = "poisson",
                random =~ ID)

MC2 <- MCMCglmm(MeanTime ~ Hurricane * (Age + Sex + Rank), 
                data = TestDF, 
                family = "poisson",
                random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(MeanTime ~ Age + Sex + Rank, 
                 data = TestDF %>% filter(Hurricane == "Pre"), 
                 family = "poisson",
                 random =~ ID)

MC3b <- MCMCglmm(MeanTime ~ Age + Sex + Rank, 
                 data = TestDF %>% filter(Hurricane == "Post"), 
                 family = "poisson",
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% saveRDS("Data/Intermediate/InfectionRepeatabilityModels.rds")

ModelList %>% map(~MCMCRep(.x)) %>% bind_rows(.id = "Model") %>% 
  mutate_at(3:5, ~round(as.numeric(.x), 3)) %>% 
  mutate_at("Model", ~c("ID Overall", "ID:Hurricane", "ID Before", "ID After")[as.numeric(.x)]) %>% 
  rename(Lower = lHPD, Upper = uHPD) %>% 
  write.csv("Data/Outputs/InfectionRepeatabilityValues.csv", row.names = F)

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

# One Big Model ####

LMER1 <- lmer(Time ~ Hurricane * (Age + Sex + Rank) + S_I + (1|ID) + (1|Pop) + (1|Rep),
              data = IndivDF2)

LMER1 %>% summary

LMER1 %>% 
  saveRDS(file = glue::glue("Data/Intermediate/FullIndividualLMER.rds"))

IM1 <- inla(Time ~ Hurricane * (Age + Sex + Rank) + S_I + f(ID, model = "iid") + f(Pop, model = "iid") + f(Rep, model = "iid"),
            data = IndivDF2)

IM1 %>% Efxplot


# Looping through reps ####

IndivDF$Rep <- IndivDF$File %>% str_split("_") %>% map_chr(2) %>% as.numeric

NIterations <- IndivDF$Rep %>% max

i <- 1

i <- "Data/Intermediate/IndividualLMER" %>% dir_ls %>% length %>% add(1)

for(i in i:NIterations){
  
  print(i)
  
  TestDF <- IndivDF2 %>% filter(Rep == i)
  
  # TestDF$Time %>% qplot
  
  # TestDF <- 
  #   IndivDF %>% 
  #   mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  #   filter(Age > 5) %>% 
  #   filter(Time > 0) %>% 
  #   mutate_at("Age", ~.x/10) %>% 
  #   na.omit %>% 
  #   group_by(ID, Pop, Sex, Rank, Age, Hurricane) %>% 
  #   summarise(MeanInf = mean(Infected), 
  #             MeanTime = mean(Time),
  #             TimeSD = sd(Time),
  #             TimeSE = TimeSD/(1000^0.5)) %>% 
  #   mutate_at("MeanTime", ~log(.x))
  
  # MCMC1 <- MCMCglmm(Time ~ Hurricane * (Age + Sex + Rank), 
  #                   random =~ ID + Pop,
  #                   data = TestDF)
  # 
  # MCMC1 %>% 
  #   saveRDS(file = glue::glue("Data/Intermediate/IndividualMCMC_{i}.rds"))
  
  # t1 <- Sys.time()
  # 
  # LMER1 <- glmer(Time ~ Hurricane * (Age + Sex + Rank) + (1|ID) + (1|Pop),
  #                family = "poisson",
  #                data = TestDF)
  # 
  # print(Sys.time() - t1)
  # 
  # LMER1 %>% 
  #   saveRDS(file = glue::glue("Data/Intermediate/IndividualLMER/{i}.rds"))
  # 
  # t1 <- Sys.time()
  
  INLACount <- inla(Time ~ Hurricane * (Age + Sex + Rank) + 
                      f(ID, model = "iid") + f(Pop, model = "iid"),
                    family = "poisson",
                    # control.compute = list(dic = TRUE),
                    data = TestDF)
  
  INLACount %>%
    saveRDS(file = glue::glue("Data/Intermediate/IndividualINLA/{i}.rds"))
  
  # print(Sys.time() - t1)
  
  # INLANB <- inla(Time ~ Hurricane * (Age + Sex + Rank) + 
  #                  f(ID, model = "iid") + f(Pop, model = "iid"),
  #                family = "nbinomial",
  #                control.compute = list(dic = TRUE),
  #                data = TestDF)
  
  # INLA1 <- inla(Time ~ Hurricane * (Age + Sex + Rank) + f(ID, model = "iid") + f(Pop, model = "iid"),
  #               data = TestDF %>% mutate_if(is.character, as.factor))
  
}

# Importing estimates ####

IMList <- 
  "Data/Intermediate/IndividualINLA" %>% dir_ls %>% 
  map(readRDS)

MCMCEstimates <- 
  IMList %>% map("summary.fixed") %>% map(c(as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "Rep") %>% 
  rename(Var = rowname)

MCMCEstimateSummaries <- 
  MCMCEstimates %>% 
  group_by(Var) %>% 
  summarise(Mean = mean(mean), 
            Lower = HPDinterval(as.mcmc(mean))[1], 
            Upper = HPDinterval(as.mcmc(mean))[2])

# LM below ####

LMList <-
  "Data/Intermediate/IndividualLMER" %>% dir_ls %>% extract(1:1000) %>% 
  map(readRDS)

MCMCEstimates <- 
  LMList %>% 
  map(c(summary, coef, as.data.frame, rownames_to_column)) %>% 
  bind_rows(.id = "rep") %>% 
  rename(Var = rowname)

MCMCEstimateSummaries <- 
  MCMCEstimates %>% 
  group_by(Var) %>% 
  summarise(Mean = mean(Estimate), 
            Lower = HPDinterval(as.mcmc(Estimate))[1], 
            Upper = HPDinterval(as.mcmc(Estimate))[2])

MCMCEstimates %>% 
  ggplot() + 
  geom_hline(yintercept = 0, lty = 2, colour = "grey") +
  geom_sina(aes(Var, Estimate, colour = Var)) +
  geom_errorbar(data = MCMCEstimateSummaries, 
                aes(x = Var, y = Mean, ymin = Lower, ymax = Upper),
                colour = "black") +
  geom_point(data = MCMCEstimateSummaries, 
             aes(x = Var, y = Mean),
             colour = "black") +
  coord_flip()
