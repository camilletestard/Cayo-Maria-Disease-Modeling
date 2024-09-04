---
title: "sna_hurricane"
author: "Alba Motes Rodrigo"
date: "2023-11-22"
output: html_document
---
#load packages and data
#packages

library(igraph); library(tibble); library(dplyr);library(ggplot2);library(foreach); library(doParallel); library(stringr); library(tidyverse);library(magrittr); library(fs); library(lme4)

# #data
# load("/Users/alba/Desktop/Postdoc_UNIL/Hurricane_project/BisonFittedNetworks.RData")
# load("/Users/alba/Desktop/Postdoc_UNIL/Hurricane_project/proximity_data-2.RData")

#Camille Load data
load("~/Documents/GitHub/Cayo-Maria-Disease-Modeling/Data/R.Data/BisonFittedNetworks.RData")
load("~/Documents/GitHub/Cayo-Maria-Disease-Modeling/Data/R.Data/proximity_data.RData")

#Remember that all networks are fully connected so modularity and clustering cannot be calculated as they are based on edge densities.


#create storage object for densities
obs_dens<-matrix(NA,ncol=1, nrow = length(unique(edgelist.all$groupyear)))
row.names(obs_dens)<-unique(edgelist.all$groupyear)
colnames(obs_dens)<-c("density")

#create groupyear combination variable
edgelist.all$groupyear<-paste(edgelist.all$group, edgelist.all$year, sep="")

#calculate observed network densities
for (g in unique(edgelist.all$groupyear)){
  z<-subset(edgelist.all, groupyear==g, select = c("ID1", "ID2", "count"))
  z<-z[z$count!= 0,]
  w<-graph_from_data_frame(z[,c(1,2)], directed = F, vertices = NULL)
  E(w)$weight <- z[,3]
  obs_dens[g,1]<-edge_density(w)}

#Create thresholded networks
posterior.filtered<-vector(mode = "list", length=length(unique(edgelist.all$groupyear)))
names(posterior.filtered)<-unique(edgelist.all$groupyear)

# set to 0 all edges with a value below threshold
s=1
for(s in unique(edgelist.all$groupyear)){
  
  df = posterior.el[[s]][,3:1002]
  df.filtered <- df
  df.filtered[df<0.0001]=0
  posterior.filtered[[s]] = posterior.el[[s]]
  posterior.filtered[[s]][,3:1002] = df.filtered
  
}

#create storage object with one value per network
modul_networks<-matrix(, ncol = length(names(posterior.el)), nrow = 1000)
colnames(modul_networks)<-names(posterior.el)

dense_networks<-matrix(, ncol = length(names(posterior.el)), nrow = 1000)
colnames(dense_networks)<-names(posterior.el)

trans_networks<-matrix(, ncol = length(names(posterior.el)), nrow = 1000)
colnames(trans_networks)<-names(posterior.el)

#calculate metrics for each network
for (g in 1:length(names(posterior.filtered))){
  for (d in 3:1002){
    #create each graph
    x<-graph_from_data_frame(posterior.filtered[[g]][,c(1,2)], directed = F, vertices = NULL)
    E(x)$weight <- posterior.filtered[[g]][,d]
    #calculate modularity
    #Greedy optimization for community assignment
    cfg<-cluster_fast_greedy(x, modularity = T, membership = T, weights = posterior.filtered[[g]][,d])
    modul_networks[d-2,g] <-modularity(cfg)
    #substract 0 edges
    y<-delete_edges(x,E(x)[E(x)$weight == 0])
    dense_networks[d-2,g]<- edge_density(y)
    trans_networks[d-2,g]<-transitivity(y, type="global")
  }
}

####################################
#Prepare dataframes for modelling
####################################

preyears<-as.character(seq(2013,2017))
#prepare modularity data for model
modul_networkst<- modul_networks %>% t()%>% as.data.frame() %>% rownames_to_column()
colnames(modul_networkst)[2:1001]<-paste("draw", seq(1:1000), sep=".")
colnames(modul_networkst)[1]<-"groupyear"
modul_networkst <- modul_networkst %>% 
  pivot_longer(
    cols = 2:1001, 
    names_to = "draw",
    values_to = "modularity"
)
modul_networkst$year<-sub('.*(\\d{4}).*', '\\1', modul_networkst$groupyear)
modul_networkst$prepost<- ifelse(modul_networkst$year %in% preyears, "pre", "post")
modul_networkst$unqid <-paste0(modul_networkst$groupyear, modul_networkst$draw)

#prepare transitivity data for model
trans_networkst<- trans_networks %>% t()%>% as.data.frame() %>% rownames_to_column()
colnames(trans_networkst)[2:1001]<-paste("draw", seq(1:1000), sep=".")
colnames(trans_networkst)[1]<-"groupyear"
trans_networkst <- trans_networkst %>% 
  pivot_longer(
    cols = 2:1001, 
    names_to = "draw",
    values_to = "transitivity"
)
trans_networkst$year<-sub('.*(\\d{4}).*', '\\1', trans_networkst$groupyear)
trans_networkst$prepost<- ifelse(trans_networkst$year %in% preyears, "pre", "post")
trans_networkst$unqid <-paste0(trans_networkst$groupyear, trans_networkst$draw)

#prepare density data for model
dense_networkst<- dense_networks %>% t()%>% as.data.frame() %>% rownames_to_column()
colnames(dense_networkst)[2:1001]<-paste("draw", seq(1:1000), sep=".")
colnames(dense_networkst)[1]<-"groupyear"
dense_networkst <- dense_networkst %>% 
  pivot_longer(
    cols = 2:1001, 
    names_to = "draw",
    values_to = "density"
)
dense_networkst$year<-sub('.*(\\d{4}).*', '\\1', dense_networkst$groupyear)
dense_networkst$prepost<- ifelse(dense_networkst$year %in% preyears, "pre", "post")
dense_networkst$unqid <-paste0(dense_networkst$groupyear, dense_networkst$draw)

#Combine all network metrics
network_metrics = merge(merge(dense_networkst, modul_networkst[,c("unqid","modularity")], by = "unqid", all = TRUE), trans_networkst[,c("unqid","transitivity")], by = "unqid", all = TRUE)
network_metrics$group = str_sub(network_metrics$groupyear,1,1)
network_metrics$prepost=factor(network_metrics$prepost, levels=c("pre","post"))

#Plot network metrics before after
#modularity
ggplot(network_metrics, aes(x=modularity, color=prepost))+#quite normal
  geom_histogram()+
  facet_grid(~group)
#density
ggplot(network_metrics, aes(x=density, color=prepost))+#quite normal
  geom_histogram()+
  facet_grid(~group)
#transitivity
ggplot(network_metrics, aes(x=transitivity, color=prepost))+#quite normal
  geom_histogram()+
  facet_grid(~group)

cor.test(network_metrics$density, network_metrics$transitivity)
cor.test(network_metrics$density, network_metrics$modularity)

#Run models
library(lmerTest)
library(broom.mixed)

#test running model on one iteration only
test=network_metrics[network_metrics$draw=="draw.2",]
mdl<-lmer(modularity ~ prepost + density + (1|group), data = test)
summary(mdl)

mdl<-lmer(transitivity ~ prepost + density + (1|group), data = test)
summary(mdl)

#Run model across all iterations using the mice package. 
df_list <- split(network_metrics, network_metrics$draw)
imp<- miceadds::datalist2mids(df_list, progress=T)

mdl <- with(imp, lme4::lmer(modularity ~ prepost + density + (1|group)))
mdl.pool<-summary(mice::pool(mdl))

mdl <- with(imp, lme4::lmer(transitivity ~ prepost + density + (1|group)))
mdl.pool<-summary(mice::pool(mdl))
#Note: I get the same results either way.

# install.packages("pdp")
# install.packages("ggplot2")
library(pdp)
library(ggplot2)

# Assume you have a linear model with response variable 'y',
# and covariates 'x1' and 'x2'
mdl<-lm(modularity ~ density, data = test)
test$modul.residuals = mdl$residuals

mdl<-lm(transitivity ~ density, data = test)
test$clust.residuals = mdl$residuals

ggplot(test, aes(y=modul.residuals, x=prepost))+
  geom_violin()+
  geom_hline(yintercept=0, linetype=2, label="As expected given density")+
  geom_point(size = 3, alpha=0.5)+
  theme_light()

ggplot(test, aes(y=clust.residuals, x=prepost))+
  geom_violin()+
  geom_point()+
  theme_light()

  

#####################
# Network resilience

library(EpiModel)
library(intergraph)

#create each graph
d=1
x<-graph_from_data_frame(posterior.filtered[[g]][,c(1,2)], directed = F, vertices = NULL)
E(x)$weight <- posterior.filtered[[g]][,d]

# Simulate epidemic spread
resilience_result <- EpiModel::resilience(x)
print(resilience_result)

