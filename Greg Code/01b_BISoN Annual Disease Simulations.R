
# 01_Disease Simulations.R

rm(list = ls())

library(dplyr)
library(igraph)
library(foreach)
library(doParallel)
library(stringr)

library(magrittr); library(fs)

load("Data/R.Data/BisonFittedNetworks.RData")

AggregatedEdges <- posterior.el %>% bind_rows(.id = "Rep")

AggregatedEdges %<>% rename(From = 2, To = 3)

AggregatedEdges %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate(Group = substr(Rep, 1, 1), 
         Year = substr(Rep, str_count(Rep) - 3, str_count(Rep)))

# AggregatedEdges %<>% mutate_at(vars(contains("draw")), ~ifelse(.x < 0.0001, 0, .x))

Reps <- AggregatedEdges$Rep %>% unique %>% sort

dir_create("Greg Data/Outputs/BISoN")

reps = 1000 #number of times the simulation should be repeated

#at each simulation time step, two individuals will interact according to their probability of proxomity (edge weight). In each time step, all possible dyads are considered.

sims = 10000 #number of time steps/times each dyad should be allowed to potentially interact/scans

pinfs = c(0.01, 0.1, 0.2)

pinf = 0.2

FocalRep <- Reps[1]

Reps <- Reps[-1]

for(FocalRep in Reps){
  
  ## ---------------------------------------------------------------------------------------------------------------------------
  
  V.data <- AggregatedEdges %>% filter(Rep == FocalRep)
  
  r <- 1
  
  IndivList <- list()
  
  for (r in r:reps){
    
    t1 = Sys.time()
    
    print(r)
    
    df <- V.data[,c("From", "To", paste0("draw.", r))] %>% rename(Weight = 3)
    
    mygraph <- graph.data.frame(df, directed = F)
    
    AdjMatrix <- get.adjacency(mygraph, sparse = FALSE, attr = 'Weight')
    
    N = length(colnames(AdjMatrix))
    
    Health = rep(0, N)
    
    Health[sample(1:length(Health), 1)] = 1
    
    Indivs <- data.frame(ID = colnames(AdjMatrix), 
                         Infected = Health,
                         Time = Health - 1)
    
    # if(!Indivs[which(Indivs$Infected == 1), "Unconnected"]){
    
    Network <- AdjMatrix
    
    # NPairs <- Network[which(Health == 1),]>0
    
    s <- 1
    
    while(s < sims & sum(Indivs$Infected) < nrow(Indivs)){
      
      TransmissionMatrix <- array(0, dim = dim(Network))
      
      I2 <- which(Indivs$Infected > 0) # Identifying infected
      
      NI2 <- setdiff(1:nrow(Indivs), 
                     which(Indivs$Infected > 0)) # Identifying uninfected
      
      if(0){
        
        TransmissionMatrix[I2,] <- 
          
          rbinom(length(Network[I2,]), 1, Network[I2,])* # Identifying if they interact
          
          as.numeric(runif(length(Network[I2,]), 0, 1) < pinf) # Identifying if they infect
        
        Infected <- which(colSums(TransmissionMatrix) > 0)# %>% as.numeric()
        
        NewlyInfected <- setdiff(Infected, 
                                 which(Indivs$Infected == 1))
        
      }
      
      if(1){
        
        NPairs <- which(Network[I2,]>0)
        
        TransmissionMatrix[I2,][NPairs] <- 
          
          rbinom(length(NPairs), 1, Network[I2,][NPairs])* # Identifying if they interact
          
          as.numeric(runif(length(NPairs), 0, 1) < pinf) # Identifying if they infect
        
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
    
    Indivs
    
    IndivList[[r]] <- Indivs
    
    print(Sys.time() - t1)
    
  }
  
  saveRDS(IndivList, file = paste0("Greg Data/Outputs/BISoN/", FocalRep, ".rds"))
  
}


