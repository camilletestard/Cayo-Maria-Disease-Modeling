
# 01_Disease Simulations.R

# knitr::purl("Greg Code/Disease-simulations.Rmd")

rm(list = ls())

## ---- warning=FALSE, message=F----------------------------------------------------------------------------------------------

library(dplyr)
library(igraph)
library(foreach)
library(doParallel)
library(stringr)

library(magrittr); library(fs)

AggregatedEdges <- readRDS("Greg Data/MetaEdges.rds")

AggregatedEdges %<>% rename(From = focal.monkey, To = in.proximity)

Observations <- AggregatedEdges %>% filter(From == To) %>% rename(Obs = Count)

AggregatedEdges %<>% # Joining the edge list with the total observations
  left_join(Observations %>% dplyr::select(-To), by = c("From", "Rep")) %>% 
  left_join(Observations %>% dplyr::select(-From), by = c("To", "Rep"), suffix = c(".From", ".To"))

AggregatedEdges %<>% mutate(Weight = Count/(Obs.From + Obs.To - Count))

AggregatedEdges$Weight %>% qplot

AggregatedEdges %<>% 
  mutate_at("Rep", str_trim) %>% 
  mutate(Group = substr(Rep, 1, 1), 
         Year = substr(Rep, str_count(Rep) - 3, str_count(Rep)))

AggregatedEdges %<>% filter(!Year == 2018)

AggregatedEdges %<>% mutate(Post = Year >= 2018)

AggregatedEdges %<>% filter(From != To)

Reps <- AggregatedEdges$Rep %>% unique %>% sort

dir_create("Greg Data/Outputs")

FocalRep <- Reps[1]

# Reps <- Reps[!str_detect(Reps, 2018)]

for(FocalRep in Reps){
  
  ## ---------------------------------------------------------------------------------------------------------------------------
  
  V.data <- AggregatedEdges %>% filter(Rep == FocalRep)
  
  library(magrittr)
  
  V.data %<>% dplyr::select(c("From", "To", "Weight", "Year"))
  
  #create different data frames for each year and combine into one list V.data.list
  
  yearsV <- V.data$Year %>% unique
  
  library(purrr)
  
  V.data.list <- yearsV %>% 
    map(~V.data %>% filter(Year == .x) %>% 
          dplyr::select(-Year))
  
  names(V.data.list) <- yearsV
  
  ## ---------------------------------------------------------------------------------------------------------------------------
  
  (Mins <- V.data.list %>% map(~min(.x$Weight)))
  
  (Maxes <- V.data.list %>% map(~max(.x$Weight)))
  
  list_names = yearsV
  
  ## ---------------------------------------------------------------------------------------------------------------------------
  
  list_df = V.data.list
  
  reps = 1000 #number of times the simulation should be repeated
  
  #at each simulation time step, two individuals will interact according to their probability of proxomity (edge weight). In each time step, all possible dyads are considered.
  
  sims = 10000 #number of time steps/times each dyad should be allowed to potentially interact/scans
  simu_res = list()#storing object for simulation results
  
  props_res = data.frame(matrix(NA, nrow = sims, ncol = reps))
  
  props_res_list = list()
  
  pinfs = c(0.01, 0.1, 0.2)
  
  d = 1
  # pinf = 0.01
  pinf = 0.2
  
  r <- 1
  
  df <- list_df[[d]]
  
  mygraph <- graph.data.frame(df, directed = F)
  
  my_mat <- get.adjacency(mygraph, sparse = FALSE, attr = 'Weight')*3
  
  Zeroes <- as.numeric(colSums(my_mat) == 0)
  
  N = length(colnames(my_mat))
  
  IndivList <- list()
  
  r <- 1
  
  for (r in 1:reps){
    
    t1 = Sys.time()
    
    print(r)
    
    #create vector to store the health status of each individual at each simulation step (sim)
    
    health = rep(0, N)
    
    #infect a random individual, change its status to infected in health vector
    
    # health[floor(runif(1, min = 1, max = N + 1))] = 1
    
    health[sample(1:length(health), 1)] = 1
    
    Indivs <- data.frame(ID = colnames(my_mat), 
                         Infected = health,
                         Unconnected = Zeroes,
                         Time = health - 1)
    
    if(!Indivs[which(Indivs$Infected == 1), "Unconnected"]){
      
      Network <- my_mat
      
      NPairs <- Network[which(health == 1),]>0
      
      s <- 1
      
      while(s < sims & sum(Network[which(Indivs$Infected == 1), -which(Indivs$Infected == 1)]) > 0){
        
        TransmissionMatrix <- array(0, dim = dim(Network))
        
        I2 <- which(Indivs$Infected > 0) # Identifying infected
        
        NI2 <- setdiff(1:nrow(Indivs), 
                       which(Indivs$Infected > 0)) # Identifying uninfected
        
        NPairs <- which(Network[I2,]>0)
        
        TransmissionMatrix[I2,][NPairs] <- 
          rbinom(length(NPairs), 1, Network[I2,][NPairs])* # Identifying if they interact
          as.numeric(runif(length(NPairs), 0, 1) < pinf) # Identifying if they infect
        
        Infected <- which(colSums(TransmissionMatrix) > 0)# %>% as.numeric()
        
        NewlyInfected <- setdiff(Infected, 
                                 which(Indivs$Infected == 1))
        
        if(length(NewlyInfected) > 0){
          
          Indivs[NewlyInfected, "Time"] <- s
          Indivs[NewlyInfected, "Infected"] <- 1
          
        }
        
        # #change the status of the newly infected individuals in the health status vector
        # 
        # health_unique <- unique(health_update)
        # 
        # for (h in health_unique){
        #   health[h] = 1
        # }
        # 
        # #calculate proportion of infected individuals in each simulation step (sims)
        # props[s]<- sum(health)/length(health)  
        
        s <- s + 1
        
      }
      
      Indivs
      
      # }
      
    }else{
      
      print("Unconnected!")
      
    }
    
    IndivList[[r]] <- Indivs
    
    print(Sys.time() - t1)
    
  }
  
  saveRDS(IndivList, file = paste0("Greg Data/Outputs/", FocalRep, ".rds"))
  
}

#name data frames in list with year names
names(simu_res_pinf)<-list_names

for (i in 1:6){
  names(simu_res_pinf[[i]]) <- pinfs
}


## ---------------------------------------------------------------------------------------------------------------------------
list_df=V.data.list
list_names=yearsV
pinfs=c(0.01, 0.1, 0.2)
simu_res=list()
sims=10000
simu_res_ids<-foreach(d=1:length(list_df))%:% foreach(pinf=pinfs)%dopar%{
  #for each year...
  df <- list_df[[d]]
  #create adjancecy matrix
  mygraph <- igraph::graph.data.frame(df)
  my_mat<-igraph::get.adjacency(mygraph, sparse = FALSE, attr='weight')
  N=length(colnames(my_mat))
  #next line optional but ensures lower tri is NAs
  my_mat[lower.tri(my_mat)] <- NA
  #empty vectors for simulation results
  #create vector to store the health status of each individual at each simulation step (sim)
  inf_ids_df=data.frame(matrix(NA, nrow = N, ncol = sims))
  health=rep(0, N)
  #infect a random individual, change its status to infected in health vector
  health[floor(runif(1, min=1, max=N+1))]=1
  for(s in 1:sims){
    #create new vector where infected individuals ID are stored after each iteraction over the matrix
    health_update=c()
    #iterate over upper triangle matrix
    for(i in 1:(N-1)){
      for(j in (i+1):N){
        #if one of the two individuals considered is infected..
        if(sum(health[i],health[j])==1){
          #evaluate whether they interact according to edge weight
          if(runif(n=1)<my_mat[i,j]){
            #if they interact, evaluate whether it gets infected. If it does...
            if(runif(n=1)<pinfs){
              #update its status to infected in the updates' vector
              health_update <- append(health_update,c(i,j))
            }
          }
        }
      }
    }
    #change the status of the newly infected individuals in the health status vector
    health_unique<- unique(health_update)
    for (h in health_unique){
      health[h]=1
    }
    inf_ids_df[, s]<-health
  }
  #store simulation steps at which thresholds were reached in each repetition of the simulation into list
  #inf_ids_df$ids<- colnames(my_mat)
  simu_res[d]<- inf_ids_df
  
}
#name data frames in list with year names
names(simu_res_ids)<-list_names
for (i in 1:6){
  names(simu_res_ids[[i]]) <- pinfs
}


#extract vectors of ids for individual present each year
ids_year=list()
for (d in 1:length(list_df)){
  df <- list_df[[d]]
  #create adjancecy matrix
  mygraph <- igraph::graph.data.frame(df)
  my_mat<-igraph::get.adjacency(mygraph, sparse = FALSE, attr='weight')
  ids_year[[d]]=colnames(my_mat)
}

#paste in each health_status data frame (simu_res_ids) the ids of the individuals present in that year
for (i in 1:length(simu_res_ids)){
  for (j in 1:length(simu_res_ids[[i]])){
    x=simu_res_ids[[i]][[j]]
    x=cbind(x, ids=ids_year[[i]])
    simu_res_ids[[i]][[j]]=x
  }
}


## ---------------------------------------------------------------------------------------------------------------------------
network_char<- data.frame(densities=c(), modularity_fg=c(), betweenness=c()) 
list_df<-V.data.list.raw

for(d in 1:length(list_df)){
  df <- list_df[[d]]
  #create adjancency matrix from edgelist
  mygraph <- igraph::graph_from_data_frame(df, directed=F)#create graph from data frame
  #density
  density_res<-sum(df$weight>0)/nrow(df)
  network_char[d, "densities"]=density_res
  #clustering
  modularity_res_fg<-cluster_fast_greedy(mygraph, modularity = T, membership = T)
  modularity_res_fg<-max(modularity_res_fg$modularity)
  network_char[d, "modularity_fg"]=modularity_res_fg
  network_char[d, "year"]=names(list_df[d])
}


## ---------------------------------------------------------------------------------------------------------------------------
list_df<-V.data.list
betweenness_data=list()
seeds_betweenness=list()
for(d in 1:length(list_df)){
  df <- list_df[[d]]
  #create adjancency matrix from edgelist
  mygraph <- igraph::graph_from_data_frame(df, directed=F)#create graph from data frame
  betweenness_data[[d]]<-betweenness(mygraph, directed = F,v = V(mygraph), normalized = FALSE)
  seeds_betweenness[[d]]<-subset(betweenness_data[[d]],betweenness_data[[d]] < quantile(betweenness_data[[d]], 0.1))
}
names(betweenness_data)<-list_names
names(seeds_betweenness)<-list_names


## ---------------------------------------------------------------------------------------------------------------------------
list_df=V.data.list
reps=1000#number of times the simulation should be repeated
#at each simulation time step, two individuals will interact according to their probability of proxomity (edge weight). In each time step, all possible dyads are considered.
sims=10000#number of time steps/times each dyad should be allowed to potentially interact/scans
simu_res = list()#storing object for simulation results
props_res=data.frame(matrix(NA, nrow = sims, ncol = reps))
props_res_list= list()
list_names=yearsV
pinfs=c(0.01, 0.1, 0.2)


simu_seeded_res_pinf<-foreach(d=1:length(list_df)) %:% foreach(pinf=pinfs)%dopar%{
  #for each year...
  df <- list_df[[d]]
  #create adjancecy matrix
  mygraph <- graph.data.frame(df)
  my_mat<-get.adjacency(mygraph, sparse = FALSE, attr='weight')
  #calculate group size
  N=length(colnames(my_mat))
  #next line optional but ensures lower tri is NAs
  my_mat[lower.tri(my_mat)] <- NA
  #empty vectors for simulation results
  value100 <- rep(NA, reps)
  value90 <- rep(NA, reps)
  value50 <- rep(NA, reps)
  value60 <- rep(NA, reps)
  value70 <- rep(NA, reps)
  value80 <- rep(NA, reps)
  maxprop <- rep(NA, reps)
  props <- rep(NA, sims)
  for (r in 1:reps){
    #create vector to store the health status of each individual at each simulation step (sim)
    health=rep(0, N)
    names(health)<-names(betweenness_data[[d]])
    #infect a random individual, change its status to infected in health vector
    health[names(sample(seeds_betweenness[[d]], size=1))]=1
    for(s in 1:sims){
      #create new vector where infected individuals ID are stored after each iteraction over the matrix
      health_update=c()
      #iterate over upper triangle matrix
      for(i in 1:(N-1)){
        for(j in (i+1):N){
          #if one of the two individuals considered is infected..
          if(sum(health[i],health[j])==1){
            #evaluate whether they interact according to edge weight
            if(runif(n=1)<my_mat[i,j]){
              #if they interact, evaluate whether it gets infected. If it does...
              if(runif(n=1)<pinf){
                #update its status to infected in the updates' vector
                health_update <- append(health_update,c(i,j))
              }
            }
          }
        }
      }
      #change the status of the newly infected individuals in the health status vector
      health_unique<- unique(health_update)
      for (h in health_unique){
        health[h]=1
      }
      #calculate proportion of infected individuals in each simulation step (sims)
      props[s]<- sum(health)/length(health)  
    }
    #write down at which simulation step different thresholds of proportions of individuals infected were reached for each repetition of the simulation (reps)
    value100[r] <- which(props == 1)[1]
    value90[r] <- which(props >= 0.9)[1]
    value80[r] <- which(props >= 0.8)[1]
    value70[r] <- which(props >= 0.7)[1]
    value60[r] <- which(props >= 0.6)[1]
    value50[r] <- which(props >= 0.5)[1]
    maxprop[r] <- max(props)
    #store proportion of infected individuals at each simulation step for each repetition of the simulation (reps) 
    props_res[, r]<-props
  }
  #store simulation steps at which thresholds were reached in each repetition of the simulation into list
  simu_res[[d]]<- list(data.frame("time_to_100"=value100,
                                  "time_to_90"=value90, 
                                  "time_to_50"=value50, 
                                  "time_to_80"=value80, 
                                  "time_to_70"=value70, 
                                  "time_to_60"=value60,
                                  "max_infected"=maxprop), props_res) 
  
}
#name data frames in list with year names
names(simu_seeded_res_pinf)<-list_names

for (i in 1:6){
  names(simu_seeded_res_pinf[[i]]) <- pinfs
}


## ---------------------------------------------------------------------------------------------------------------------------
unregister_dopar <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
#unregister_dopar()


## ---------------------------------------------------------------------------------------------------------------------------
save.image(file="Data/R.Data/macaque_disease_simulation.RData")
load(file="Data/R.Data/macaque_disease_simulation.RData")

