
# 02b_BISoN Individual Summarising ####

library(tidyverse); library(ggregplot); library(ggforce); library(cowplot); library(fs)
library(magrittr); library(colorspace); library(lme4); library(lmerTest); library(patchwork)

theme_set(theme_cowplot())

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

# IndivDF$S_I <- IndivDF$File %<>% str_split("_") %>% map_chr(last) %>% str_remove(".rds")
# 
# IndivDF$Rep <- 
#   IndivDF$File %>% 
#   str_split("_") %>% 
#   map_chr(2)
# 
# IndivDF$Rep %<>% as.numeric

# IndivDF %>% saveRDS("Data/Intermediate/IndividualInfectionStatus.rds")

# Running model iterations ####

load("Data/Intermediate/proximity_data.RData")

IndividualTraits <-
  edgelist.all %>%
  mutate(IsPost = ifelse(year < 2018, "Pre", "Post")) %>%
  mutate(Pop = paste0(group, year)) %>%
  dplyr::select(ID = ID1, Sex = ID1_sex, Rank = ID1_rank, Age = ID1_age, Pop, IsPost) %>%
  # bind_rows(edgelist.all %>%
  #             dplyr::select(ID = ID2, Sex = ID2_sex, Rank = ID2_rank, Age = ID2_age, Group = group, Year = year)) %>%
  unique

# IDList <- 
#   "Data/Input" %>% dir_ls(regex = "GroupByYear") %>% map(read.delim)
# 
# IDList %>% 
#   bind_rows(.id = "File")

IndivDF %<>% left_join(IndividualTraits)

IndivDF

# IndivDF %<>% mutate(IsPost = ifelse(Year < 2018, "Pre", "Post"))

IndivDF$Rep <- IndivDF$File %>% str_split("_") %>% map_chr(2) %>% as.numeric

NIterations <- IndivDF$Rep %>% max

i <- 1

dir_create("Data/Intermediate/IndividualMCMC")

IndivDF2 <- 
  IndivDF %>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  # mutate_at("Time", ~log(.x)) %>% 
  mutate_at("Age", ~.x/10) %>% 
  mutate(S_I = str_split(File, "_") %>% map_chr(last) %>% str_remove(".rds") %>% as.numeric %>% multiply_by(10)) %>%
  na.omit

# Mean Model ####

TestDF <- 
  IndivDF %>% 
  mutate_at("Sex", ~substr(.x, 1, 1)) %>% 
  filter(Age > 5) %>% 
  filter(Time > 0) %>% 
  mutate_at("Age", ~.x/10) %>% 
  na.omit %>% 
  group_by(ID, Pop, Sex, Rank, Age, IsPost) %>% 
  summarise(MeanInf = mean(Infected), 
            MeanTime = mean(Time),
            TimeSD = sd(Time),
            TimeSE = TimeSD/(1000^0.5)) %>% 
  # mutate_at("MeanTime", ~log(.x))
  ungroup

LMER1 <- lmer(MeanTime ~ IsPost * (Age + Sex + Rank) + (1|ID) + (1|Pop),
              data = TestDF)

library(INLA)

TestDF %<>% mutate_at("MeanTime", round)

INLAGaussian <- inla(MeanTime ~ IsPost * (Age + Sex + Rank) + 
                       f(ID, model = "iid") + f(Pop, model = "iid"),
                     # family = "poisson",
                     control.compute = list(dic = TRUE),
                     data = TestDF)# %>% ungroup %>% mutate_at("MeanTime", ~c(scale(.x))))

INLALogGaussian <- inla(MeanTime ~ IsPost * (Age + Sex + Rank) + 
                          f(ID, model = "iid") + f(Pop, model = "iid"),
                        # family = "poisson",
                        family = "lognormal",
                        control.compute = list(dic = TRUE),
                        data = TestDF# %>% mutate_at("MeanTime", log)
)

INLACount <- inla(MeanTime ~ IsPost * (Age + Sex + Rank) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "poisson",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLAGamma <- inla(MeanTime ~ IsPost * (Age + Sex + Rank) + 
                    f(ID, model = "iid") + f(Pop, model = "iid"),
                  family = "gamma",
                  control.compute = list(dic = TRUE),
                  data = TestDF)

INLANB <- inla(MeanTime ~ IsPost * (Age + Sex + Rank) + 
                 f(ID, model = "iid") + f(Pop, model = "iid"),
               family = "nbinomial",
               control.compute = list(dic = TRUE),
               data = TestDF)

list(#INLAGaussian, 
  INLALogGaussian, INLAGamma, INLACount, INLANB) %>% Efxplot(Intercept = F)

list(INLAGaussian, INLALogGaussian, INLAGamma, INLACount, INLANB) %>% MDIC

LMER1 %>% 
  saveRDS(file = glue::glue("Data/Intermediate/MeanLMER.rds"))

IM1 <- INLAModelAdd(Data = TestDF, 
                    Family = "poisson",
                    Response = "MeanTime",
                    Explanatory = c("IsPost", "Age", "Sex", "Rank"), 
                    Add = paste0("IsPost:", c("Age", "Sex", "Rank")),
                    Random = c("ID", "Pop"), RandomModel = "iid"
)

IM1$FinalModel %>% Efxplot


# One Big Model ####

LMER1 <- lmer(Time ~ IsPost * (Age + Sex + Rank) + S_I + (1|ID) + (1|Pop) + (1|Rep),
              data = IndivDF2)

LMER1 %>% summary

LMER1 %>% 
  saveRDS(file = glue::glue("Data/Intermediate/FullIndividualLMER.rds"))

IM1 <- inla(Time ~ IsPost * (Age + Sex + Rank) + S_I + f(ID, model = "iid") + f(Pop, model = "iid") + f(Rep, model = "iid"),
            data = IndivDF2)

IM1 %>% Efxplot


# Looping through reps ####

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
  #   group_by(ID, Pop, Sex, Rank, Age, IsPost) %>% 
  #   summarise(MeanInf = mean(Infected), 
  #             MeanTime = mean(Time),
  #             TimeSD = sd(Time),
  #             TimeSE = TimeSD/(1000^0.5)) %>% 
  #   mutate_at("MeanTime", ~log(.x))
  
  # MCMC1 <- MCMCglmm(Time ~ IsPost * (Age + Sex + Rank), 
  #                   random =~ ID + Pop,
  #                   data = TestDF)
  # 
  # MCMC1 %>% 
  #   saveRDS(file = glue::glue("Data/Intermediate/IndividualMCMC_{i}.rds"))
  
  # t1 <- Sys.time()
  # 
  # LMER1 <- glmer(Time ~ IsPost * (Age + Sex + Rank) + (1|ID) + (1|Pop),
  #                family = "poisson",
  #                data = TestDF)
  # 
  # print(Sys.time() - t1)
  # 
  # LMER1 %>% 
  #   saveRDS(file = glue::glue("Data/Intermediate/IndividualLMER/{i}.rds"))
  # 
  # t1 <- Sys.time()
  
  INLACount <- inla(Time ~ IsPost * (Age + Sex + Rank) + 
                      f(ID, model = "iid") + f(Pop, model = "iid"),
                    family = "poisson",
                    # control.compute = list(dic = TRUE),
                    data = TestDF)
  
  INLACount %>%
    saveRDS(file = glue::glue("Data/Intermediate/IndividualINLA/{i}.rds"))

  # print(Sys.time() - t1)
  
  # INLANB <- inla(Time ~ IsPost * (Age + Sex + Rank) + 
  #                  f(ID, model = "iid") + f(Pop, model = "iid"),
  #                family = "nbinomial",
  #                control.compute = list(dic = TRUE),
  #                data = TestDF)
  
  # INLA1 <- inla(Time ~ IsPost * (Age + Sex + Rank) + f(ID, model = "iid") + f(Pop, model = "iid"),
  #               data = TestDF %>% mutate_if(is.character, as.factor))
  
}

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

# Averaging for Camille's models (DEFUNCT) ####

IndivDF %<>% 
  mutate(S_I_Category = cut(S_I, breaks = c(quantile(S_I, 0:3/3)), 
                            labels = c("Low", "Med", "High"), 
                            include.lowest = T))

IndivDF %<>% 
  group_by(ID, Pop, S_I_Category) %>% 
  summarise(MeanInf = mean(Infected), 
            MeanTime = mean(Time),
            TimeSD = sd(Time),
            TimeSE = TimeSD/(1000^0.5))

# IndivDF %<>% 
#   group_by(ID, Pop) %>% 
#   summarise(MeanInf = mean(Infected), 
#             MeanTime = mean(Time),
#             TimeSD = sd(Time),
#             TimeSE = TimeSD/(1000^0.5))

# IDLevels <- IndivDF %>% arrange(MeanTime) %>% mutate_at("ID", ~paste0(.x, "_", Pop)) %>% pull(ID)

# IndivDF %>% arrange(MeanTime) %>% mutate_at("ID", ~paste0(.x, "_", Pop)) %>% 
#   mutate_at("ID", ~factor(.x, levels = IDLevels)) %>% 
#   ggplot(aes(ID, MeanTime)) +
#   geom_errorbar(aes(ymin = MeanTime - TimeSE, ymax = MeanTime + TimeSE, colour = Pop), 
#                 position = position_dodge(w = 0.4)) +
#   geom_point(aes(colour = Pop), 
#              position = position_dodge(w = 0.4)) +
#   theme(legend.position = "none")

IndivMeanDF %>% saveRDS("Data/Outputs/IndividualTimesteps.rds")

