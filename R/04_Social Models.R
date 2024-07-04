
# 03b_Individual Proximity Differences.R ####

# Initially written by Camille, rewritten and pipelined by Greg

{
  
  #Load libraries
  #data wrangling
  library(dplyr); library(forcats)
  
  #Modelling
  library(lme4); library(lmerTest); library(broom.mixed); library(emmeans)
  
  #For plotting
  library(ggplot2); library(sjPlot); library(sjmisc); library(cowplot)
  
  theme_set(theme_cowplot())
  
  dir_create("Data/Intermediate/SocialLMER")
  
  load("Data/Intermediate/BisonMetrics.RData")
  
}

#Set group year list

group = c("F","KK","F","HH","F","V","R","KK","R","V","F","HH","F","KK","V","V","KK","S","V","F","V","TT","V","F")

years = c(2013, 2013,2014,2014,2015,2015,2015,2015,
          2016,2016,2016,2016,2017,2017,2017,
          2018,2018, 2019, 2019,2021,2021,2022,2022,2022)

groupyears = paste0(group,years)

# DataList <- list()

# for (gy in 1:length(groupyears)){
#   
#   # data_path<-'~/Documents/GitHub/Cayo-Maria-Survival/Data/Data All Cleaned/BehavioralDataFiles/'
#   data_path <- 'Data/Input/'
#   data <- read.csv(paste0(data_path,"Group", groupyears[gy], "_GroupByYear.txt")) #load meta data
#   
#   data = data[data$hrs.focalfollowed>0,]
#   
#   data$group = group[gy]
#   data$year = years[gy]; 
#   data$Hurricane = ifelse(years[gy]<2018, "pre","post")
#   data$Hurricane.year = ifelse(years[gy]<2018,"pre",paste(data$Hurricane, data$year,sep='.'))
#   data$id = as.factor(data$id)
#   data$sex[data$sex == "FALSE"] = "F"; data$sex = as.factor(data$sex);
#   data$id.year = paste(data$id, data$year,sep='.')
#   
#   DataList[[gy]] <- data
#   
# }

load("Data/Intermediate/proximity_data.RData")

IndividualTraits <-
  edgelist.all %>%
  mutate(Hurricane = factor(ifelse(year < 2018, "Pre", "Post"), levels = c("Pre", "Post"))) %>%
  mutate(Pop = paste0(group, year)) %>%
  dplyr::select(ID = ID1, Sex = ID1_sex, Rank = ID1_rank, Age = ID1_age, Pop, Hurricane) %>%
  unique

IndividualTraits %<>% mutate_at("Sex", ~substr(.x, 1, 1))

# STEP 1: CREATE IMPUTED DATASSET ####

# Extract change in proximity from each draw to create imputed data set:

imputed.data = list(); imputed.data.prepost = list(); i = 1; gy = 1;

imputed.density = list()

StrengthModelList <- DegreeModelList <- list()

num_iter <- 10

i <- "Data/Intermediate/SocialLMER" %>% dir_ls %>% length %>% add(1)

i <- 1

for (i in i:num_iter){
  
  # print(i)
  # 
  # data.all <- data.frame(); density.iter = data.frame(matrix(NA, nrow = num_iter)); names(density.iter)="dens"
  # 
  # density.all <- data.frame();
  # 
  # for (gy in 1:length(groupyears)){
  # 
  #   data <- DataList[[gy]]
  # 
  #   strength <- node_strength_all[[groupyears[gy]]][i,]
  #   degree <- node_degree_all[[groupyears[gy]]][i,]
  #   node.id <- node_ids[[groupyears[gy]]]
  #   data$prox.strength<-as.numeric(strength[match(data$id, node.id)])
  #   data$prox.degree<-as.numeric(degree[match(data$id, node.id)])
  #   data.final <- data[,c("id", "sex", "age", "ordinal.rank", "group", "year", "Hurricane", "Hurricane.year", "id.year",
  #                         "prox.strength", "prox.degree")]
  # 
  #   data.final %<>% na.omit
  # 
  #   data.all <- rbind(data.all, data.final)
  # 
  #   density.iter$dens = density_all[[groupyears[gy]]]
  #   density.iter$group = group[gy]
  #   density.iter$year = years[gy]
  #   density.iter$Hurricane = ifelse(years[gy]<2018, "pre","post")
  #   density.iter$Hurricane.year = ifelse(years[gy]<2018,"pre",paste(data$Hurricane, data$year,sep='.'))
  # 
  #   density.all <- rbind(density.all, density.iter)
  # 
  # }
  # 
  # data.all %<>% mutate_at("prox.strength", log10)
  # 
  # #Standardize
  # density.all$std.dens <-(density.all$dens - mean(density.all$dens))/sd(density.all$dens)
  # data.all$std.prox.strength = (data.all$prox.strength-mean(data.all$prox.strength))/sd(data.all$prox.strength)
  # data.all$std.prox.degree = (data.all$prox.degree-mean(data.all$prox.degree))/sd(data.all$prox.degree)
  # 
  # #Create year factor
  # data.all$year.factor = as.factor(data.all$year)
  # 
  # #Adjust levels
  # data.all$group<-as.factor(data.all$group)
  # data.all = data.all %>%
  #   mutate(Hurricane = fct_relevel(Hurricane,
  #                               "pre", "post"))# %>%
  #   # mutate(Hurricane.year = fct_relevel(Hurricane.year,
  #   #                                  "pre", "post.2018", "post.2019", "post.2021", "post.2022"))
  # #"post.2018","pre","post.2019","post.2021", "post.2022"))
  # 
  # density.all = density.all %>%
  #   mutate(Hurricane = fct_relevel(Hurricane,
  #                               "pre", "post")) #%>%
  #   # mutate(Hurricane.year = fct_relevel(Hurricane.year,
  #   #                                  "pre", "post.2018","post.2019","post.2021", "post.2022"))
  # 
  # #"post.2018","pre","post.2019","post.2021", "post.2022"))
  # 
  # 
  # imputed.data[[i]] = data.all
  # imputed.density[[i]] = density.all
  # 
  # 
  # #If only consider IDs present both pre and post
  # id.Hurricane = table(data.all$id, data.all$Hurricane)
  # id.year = table(data.all$id, data.all$year);
  # id_pre = row.names(as.data.frame(which(id.Hurricane[,"pre"]>0))); id_post = row.names(as.data.frame(which(id.Hurricane[,"post"]>0)))
  # id.PreAndPost = as.data.frame(row.names(as.data.frame(which(id.Hurricane[,"pre"]>0 & id.Hurricane[,"post"]>0)))); names(id.PreAndPost)="id"
  # id.PreAndPost$group = data.all$group[match(id.PreAndPost$id, data.all$id)]
  # table(id.PreAndPost$group)
  # data.prepost = subset(data.all, id %in% id.PreAndPost$id)
  # 
  # imputed.data.prepost[[i]] <- data.prepost
  # 
  # data.all %<>% rename_all(CamelConvert) %>% 
  #   mutate(Pop = Group, ID = Id) %>% 
  #   rename(Rank = Ordinal.rank)
  # 
  # # Running internal Bayesian models 
  # 
  # # mdl.strength <- lmer(std.prox.strength ~ 1 + Hurricane*(ordinal.rank + age + sex) + (1|id/group),
  # #                      data = data.all)
  # 
  # # mdl.strength <- MCMCglmm(std.prox.strength ~ 1 + Hurricane*(ordinal.rank + age + sex), 
  # #                          random = ~ id + group,
  # #                          data = data.all)
  
  Strengths <-
    node_strength_all %>%  map(~.x[i,]) %>% map(c(as.data.frame, reshape2::melt)) %>%
    bind_rows(.id = "Pop") %>%
    rename(ID = variable, Strength = value)
  
  Degrees <-
    node_degree_all %>%  map(~.x[i,]) %>% map(c(as.data.frame, reshape2::melt)) %>%
    bind_rows(.id = "Pop") %>%
    rename(ID = variable, Degree = value)
  
  TestDF <-
    list(IndividualTraits,
         Strengths,
         Degrees) %>%
    reduce(~left_join(.x, .y, by = c("ID", "Pop"))) %>% 
    na.omit %>% 
    mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H")))
  
  TestDF %<>%
    filter(Age > 5) %>% 
    mutate_at("Age", ~.x/10)
  
  TestDF %<>% mutate_at("Strength", log)
  
  # mdl.strength <- lmer(Strength ~ Hurricane * (Age + Sex + Rank) +
  #                      (1|ID) + (1|Pop),
  #                    data = TestDF)
  
  mdl.strength <- inla(Strength ~ Hurricane * (Age + Sex + Rank) + 
                         f(ID, model = "iid") + f(Pop, model = "iid"),
                       # family = "poisson",
                       control.compute = list(dic = TRUE),
                       data = TestDF)
  
  StrengthModelList[[i]] <- mdl.strength
  
  # mdl.strength %>% saveRDS(glue::glue("Data/Intermediate/SocialLMERStrength/{i}.rds"))
  
  # mdl.degree <- lmer(std.prox.degree ~ 1 + Hurricane*(ordinal.rank + age + sex) +
  #                      (1|id/group),
  #                    data = data.all)
  
  # mdl.degree <- MCMCglmm(std.prox.degree ~ 1 + Hurricane*(ordinal.rank + age + sex),
  #                        random = ~ id + group,  
  #                        data = data.all)
  
  mdl.degree <- inla(Degree ~ Hurricane * (Age + Sex + Rank) + 
                       f(ID, model = "iid") + f(Pop, model = "iid"),
                     # family = "poisson",
                     control.compute = list(dic = TRUE),
                     data = TestDF)
  
  DegreeModelList[[i]] <- mdl.degree
  
  # mdl.degree %>% saveRDS(glue::glue("Data/Intermediate/SocialLMERDegree/{i}.rds"))
  
  # pairwise comparisons
  
  # xx.strength <- as.data.frame(summary(emmeans(mdl.strength, pairwise ~ Hurricane*ordinal.rank))$contrasts)
  
}

StrengthModelList[1:10] %>% Efxplot +
  
  DegreeModelList[1:10] %>% Efxplot

# Analysing means ####

Strengths <-
  node_strength_all %>% map(colMeans) %>% map(c(reshape2::melt, rownames_to_column)) %>%
  bind_rows(.id = "Pop") %>%
  rename(ID = rowname, Strength = value)

Degrees <-
  node_degree_all %>% map(colMeans) %>% map(c(reshape2::melt, rownames_to_column)) %>%
  bind_rows(.id = "Pop") %>%
  rename(ID = rowname, Degree = value)

TestDF <-
  list(IndividualTraits,
       Strengths,
       Degrees) %>%
  reduce(~left_join(.x, .y, by = c("ID", "Pop"))) %>% 
  na.omit %>% 
  mutate_at("Rank", ~factor(.x, levels = c("L", "M", "H")))

TestDF %<>%
  filter(Age > 5) %>% 
  mutate_at("Age", ~.x/10)

TestDF %<>% mutate_at("Strength", log)# %>% 
  # mutate_at(c("Degree", "Strength"), ~-.x) 

# mdl.strength <- inla(Strength ~ Hurricane * (Age + Sex + Rank) + 
#                        f(ID, model = "iid") + f(Pop, model = "iid"),
#                      # family = "poisson",
#                      control.compute = list(dic = TRUE),
#                      data = TestDF)

IM2 <- INLAModelAdd(Data = TestDF, 
                    # Family = "poisson",
                    Response = "Strength",
                    Explanatory = c("Hurricane", "Age", "Sex", "Rank"), 
                    Add = paste0("Hurricane:", c("Age", "Sex", "Rank")),
                    AllModels = T,
                    Random = c("ID", "Pop"), RandomModel = "iid"
)

# mdl.degree <- inla(Degree ~ Hurricane * (Age + Sex + Rank) + 
#                      f(ID, model = "iid") + f(Pop, model = "iid"),
#                    # family = "poisson",
#                    control.compute = list(dic = TRUE),
#                    data = TestDF)

IM3 <- INLAModelAdd(Data = TestDF, 
                    # Family = "poisson",
                    Response = "Degree",
                    Explanatory = c("Hurricane", "Age", "Sex", "Rank"), 
                    Add = paste0("Hurricane:", c("Age", "Sex", "Rank")),
                    AllModels = T,
                    Random = c("ID", "Pop"), RandomModel = "iid"
)

list(IM1, IM2, IM3) %>% map("FinalModel") %>% 
  Efxplot(ModelNames = c("Infection", "Strength", "Degree"), 
          Intercept = F)

IM1$FinalModel %>% Efxplot(Intercept = F) + 
  IM2$FinalModel %>% Efxplot(Intercept = F) + 
  IM3$FinalModel %>% Efxplot(Intercept = F)

# Repeatability ####

TestDF %<>% 
  mutate(Year = substr(Pop, str_count(Pop) - 3, str_count(Pop))) %>% 
  mutate(Pop = substr(Pop, 1, 1))

# MC1 <- MCMCglmm(Strength ~ Hurricane * (Age + Sex + Rank) + Pop, data = TestDF, random =~ ID)
# 
# MC2 <- MCMCglmm(Strength ~ Hurricane * (Age + Sex + Rank) + Pop, data = TestDF, random =~ ID + ID:Hurricane)
# 
# MC3a <- MCMCglmm(Strength ~ (Age + Sex + Rank) + Pop,
#                  data = TestDF %>% filter(Hurricane == "Pre"), 
#                  random =~ ID)
# 
# MC3b <- MCMCglmm(Strength ~ (Age + Sex + Rank) + Pop,
#                  data = TestDF %>% filter(Hurricane == "Post"), 
#                  random =~ ID)

# TestDF %<>% mutate(Hurricane = as.numeric(Hurricane == "Pre")) %>% droplevels

MC1 <- MCMCglmm(Strength ~ Year + Pop, # + Hurricane,
                data = TestDF, random =~ ID)

MC2 <- MCMCglmm(Strength ~ Year + Pop, # + Hurricane,
                data = TestDF, random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(Strength ~ Year + Pop, # + Hurricane, # + PostMaria,
                 data = TestDF %>% filter(Hurricane == "Pre") %>% droplevels,
                 random =~ ID)

MC3b <- MCMCglmm(Strength ~ Year + Pop,
                 data = TestDF %>% filter(Hurricane == "Post") %>% droplevels,
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% map("DIC")

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


MC1 <- MCMCglmm(Degree ~ Year + Pop, # + Hurricane,
                data = TestDF, random =~ ID)

MC2 <- MCMCglmm(Degree ~ Year + Pop, # + Hurricane,
                data = TestDF, random =~ ID + ID:Hurricane)

MC3a <- MCMCglmm(Degree ~ Year + Pop, # + Hurricane, # + PostMaria,
                 data = TestDF %>% filter(Hurricane == "Pre") %>% droplevels,
                 random =~ ID)

MC3b <- MCMCglmm(Degree ~ Year + Pop,
                 data = TestDF %>% filter(Hurricane == "Post") %>% droplevels,
                 random =~ ID)

ModelList <- list(MC1, MC2, MC3a, MC3b)

ModelList %>% map("DIC")

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

