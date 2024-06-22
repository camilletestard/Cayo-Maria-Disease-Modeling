
# 01_Disease Simulations.R

rm(list = ls())

library(dplyr); library(igraph); library(foreach); library(doParallel); library(stringr); library(tidyverse)
library(magrittr); library(fs)

# load("Data/R.Data/BisonFittedNetworks.RData")
# load("Data/BisonFittedNetworks (2).RData")

load("Data/R.Data/BisonFittedNetworks.RData")

AggregatedEdges <- posterior.el %>% bind_rows(.id = "Rep")

AggregatedEdges %<>% rename(From = 2, To = 3)

AggregatedEdges %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate(Group = substr(Rep, 1, 1), 
         Year = substr(Rep, str_count(Rep) - 3, str_count(Rep)))

Reps <- AggregatedEdges$Rep %>% unique %>% sort

dir_create("Greg Data/Outputs/BISoN/Random")

reps = 1000 # number of times the simulation should be repeated

sims = 10000 # number of time steps/times each dyad should be allowed to potentially interact

MeanInf <- 0.15
InfSD <- 0.04

FocalRep <- Reps[1]

for(FocalRep in Reps){
  
  print(FocalRep)
  
  RepData <- AggregatedEdges %>% filter(Rep == FocalRep)

  r <- dir_ls("Greg Data/Outputs/BISoN/Random",
              regex = "Greg Data/Outputs/BISoN/Random/" %>%
                paste0(FocalRep)) %>%
    str_split("_") %>%
    map_chr(2) %>%
    as.numeric %>% max %>% add(1)

  if(r < 1) r <- 1
  
  if(r < 1001){
    
    IndivList <- list()
    
    for (r in r:reps){
      
      t1 = Sys.time()
      
      print(r)
      
      Network <- AdjMatrix <-
        RepData[,c("From", "To", paste0("draw.", r))] %>% 
        rename(Weight = 3) %>% 
        graph.data.frame(directed = F) %>% 
        get.adjacency(sparse = FALSE, attr = 'Weight')
      
      N = length(colnames(AdjMatrix))
      
      Health = rep(0, N)
      
      Health[sample(1:length(Health), 1)] = 1
      
      Indivs <- data.frame(ID = colnames(AdjMatrix), 
                           Infected = Health,
                           Time = Health - 1)
      
      Network <- AdjMatrix
      
      P_I <- rnorm(1, MeanInf, InfSD)
      
      if(P_I < 0) P_I <- 0.000001      
      
      s <- 1
      
      while(s < sims & sum(Indivs$Infected) < nrow(Indivs)){
        
        TransmissionMatrix <- array(0, dim = dim(Network))
        
        I2 <- which(Indivs$Infected > 0) # Identifying infected
        
        NI2 <- setdiff(1:nrow(Indivs), 
                       which(Indivs$Infected > 0)) # Identifying uninfected
        
        if(1){
          
          TransmissionMatrix[I2,] <- 
            
            rbinom(length(Network[I2,]), 1, Network[I2,])* # Identifying if they interact
            rbinom(length(Network[I2,]), 1, P_I) # Identifying if they infect
          
          Infected <- which(colSums(TransmissionMatrix) > 0)# %>% as.numeric()
          
          NewlyInfected <- setdiff(Infected, 
                                   which(Indivs$Infected == 1))
          
        }
        
        if(length(NewlyInfected) > 0){
          
          Indivs[NewlyInfected, "Time"] <- s
          Indivs[NewlyInfected, "Infected"] <- 1
          
        }
        
        s <- s + 1
        
      }
      
      saveRDS(Indivs, file = paste0("Greg Data/Outputs/BISoN/Random/",FocalRep, "_", 
                                    
                                    str_pad(r, width = 4, side = "left", pad = "0"), 
                                    
                                    "_", P_I, ".rds"))
      
      print(Sys.time() - t1)
      
    }
    
  }
  
}


