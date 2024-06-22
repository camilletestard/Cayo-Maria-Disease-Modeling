#Binary Network generation based on observed probabilities 
```{r}
sims=100
set.seed(1)
sim_res=list()
network_sims<-foreach(d=1:length(list_df))%dopar%{
  df <- list_df[[d]]
  for (s in 1:sims){
    new.vec=c()
    for(i in 1:nrow(df)){
      new.vec[i]<-ifelse(runif(n=1)<df[i,"weight"], 1, 0)
    }
    df<-cbind(df, new.vec)
  }
  sim_res[[d]]<-df
} 
names(network_sims)<-list_names
```

#Plot simulations V group
Each line is a simulation. Bold lines are the means per years. Dashed lines are after the hurricane.
```{r}
par(mfrow=c(1,6))
for(i in 1:length(simu_res_pinf)){
  hist(simu_res_pinf[[i]]$time_to_50, col=c("chocolate3"), breaks=20, ylim=c(0,5), xlim=c(0, 3000), xaxt="n", xlab=names(simu_res_pinf)[i], ylab="", main="",cex.axis=1)
  axis(1, at = seq(0, 3000, by = 500), las=2, cex.axis=1)
  clip(0,10, 0, 5)
  abline(v=mean(simu_res_pinf[[i]]$time_to_50, na.rm=T), col="orange", lty=1, lwd=3)
  abline(v=median(simu_res_pinf[[i]]$time_to_50, na.rm=T), col="red", lty=3, lwd=3)
}
```

```{r}
pinf_data<-do.call(rbind, simu_res_pinf)#complete data sets from each simulation
pinf_data<- tibble::rownames_to_column(pinf_data)
pinf_data[c('year', 'rep')] <- str_split_fixed(pinf_data$rowname, pattern='.', n=2)
```


# Plots results of simulation per year
```{r}
ydata_sub90 <- subset(all_data, metric=="90")
ydata_sub50 <- subset(all_data, metric=="50")

cols90 <- c("mean", "sd", "skew")

par(mfrow=c(2,3), mar=c(4,4,2,2), cex=1.25)
for (i in cols90){
  stripchart(ydata_sub90[,i]~ydata_sub90$sp, pch=19, col=colours[as.numeric(ydata_sub90$col)], method="jitter", jitter=0.1, ylab="Species", xlab=paste(i, "simulation cycles\n\ until 90% infected", sep=" "))
}

for (i in cols90){
  stripchart(ydata_sub50[,i]~ydata_sub50$sp, pch=19, col=colours[as.numeric(ydata_sub90$col)], method="jitter", jitter=0.1, ylab="Species", xlab=paste(i, "simulation cycles\n\ until 50% infected", sep=" "))
}
```