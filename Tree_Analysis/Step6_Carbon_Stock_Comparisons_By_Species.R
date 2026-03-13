path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/floodplain-experiment-repo"
setwd(path_to_repo)

library(tidyr)
library(dplyr)
library(reshape2)
library(brms)
library(bayestestR)

################################################################################
# load metadata

# treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# ecosystem carbon pools
pools = c("total.cwd.carbon","large.snag.carbon","live.tree.carbon")
df.paths = c("Clean_Data_By_Species/CWD_Carbon_Stocks_By_Species.csv",
             "Clean_Data_By_Species/Snag_Carbon_Stocks_By_Species.csv",
             "Clean_Data_By_Species/Woody_Biomass_Carbon_Stocks_By_Species.csv")
n.pool = length(pools)

################################################################################
# statistical analyses on coarse woody debris C stocks

## read in data frames
df.list = list()
for (i in 1:n.pool) { df.list[[pools[i]]] = read.csv(paste("Tree_Analysis/",df.paths[i],sep="")) }

## clean up live tree data

# make full species name column
df.list[["live.tree.carbon"]] = unite(df.list[["live.tree.carbon"]], col='species', c('genus','species'), sep=' ')

# make wide df
df.list[["live.tree.carbon"]] = pivot_wider(df.list[["live.tree.carbon"]], 
                                            id_cols = c(treatment, strip, plot),
                                            names_from = species, 
                                            values_from = bmC.ha3)
df.list[["live.tree.carbon"]][is.na(df.list[["live.tree.carbon"]])] = 0

# reshape dataframe
df.list[["live.tree.carbon"]] = melt(df.list[["live.tree.carbon"]], 
                                     id.vars=c("treatment","strip","plot"), 
                                     variable.name="species",
                                     value.name="live.tree.carbon")


# make simplified data lists
mean.list = list()
data.list = list()
for (i in 1:n.pool) {
  pool.i = pools[i]
  df.i = df.list[[pool.i]]
  spp.i = levels(factor(df.i$species))
  n.spp.i = length(spp.i)
  mean.list[[pool.i]] = list()
  data.list[[pool.i]] = list()
  for (j in 1:n.spp.i) {
    spp.ij = spp.i[j]
    spp.ij.id = which(df.i$species == spp.ij)
    mean.list[[pool.i]][[spp.ij]] = mean(df.i[spp.ij.id,pool.i])
    data.list[[pool.i]][[spp.ij]] = list(treatment = factor(df.i[spp.ij.id,"treatment"]), 
                                         strip = factor(df.i[spp.ij.id,"strip"]),
                                         y = df.i[spp.ij.id,pool.i]/mean.list[[pool.i]][[spp.ij]]) 
  }
}

# fit simple model for C stocks by species and pool
seed = 3141; n.iter=10000; n.chain=10
model.list = list()
rhat.df = data.frame(matrix(nrow=0, ncol=6))
df.int = data.frame(matrix(nrow=0, ncol=10))
colnames(rhat.df) = c("pool","species","species.mean","treatment","full.treatment.name","Rhat")
colnames(df.int) = c("pool","species","species.mean","treatment","full.treatment.name",
                     "posterior.mean","X5","X95","X25","X75")
for (i in 1:n.pool) {
  pool.i = pools[i]
  print(pool.i)
  df.i = df.list[[pool.i]]
  spp.i = unique(df.i$species)
  n.spp.i = length(spp.i)
  model.list[[pool.i]] = list()
  for (j in 1:n.spp.i) {
    if (i == 1 & j > 1) {
      rhat.df = read.csv("Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Rhat_Statistic_10Chains.csv",header=T)
      df.int = read.csv("Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Means_Intervals_10Chains.csv",header=T)
    } else if (i > 1) {
      rhat.df = read.csv("Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Rhat_Statistic_10Chains.csv",header=T)
      df.int = read.csv("Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Means_Intervals_10Chains.csv",header=T)
    }
    spp.ij = spp.i[j]
    print(spp.ij)
    m.ij = brm(data = data.list[[pool.i]][[spp.ij]], 
               y ~ treatment + (1|treatment:strip),
               family = gaussian(link="log"),
               prior=prior(normal(0,1), class=b),
               chains=n.chain, seed=seed, iter=n.iter)
    #control = list(adapt_delta = 0.99,
    #               max_treedepth = 12)
    print("Chains finished")
    model.list[[pool.i]][[spp.ij]] = m.ij
    
    
    # evaluate convergence with r-hat statistic
    model.sum.ij = summary(m.ij)
    rhat.post.ij = data.frame(matrix(nrow=n.t, ncol=6))
    colnames(rhat.post.ij) = c("pool","species","species.mean","treatment","full.treatment.name","Rhat")
    rhat.post.ij$pool = pool.i
    rhat.post.ij$species = spp.ij
    rhat.post.ij$species.mean = mean.list[[pool.i]][[spp.ij]]
    rhat.post.ij$treatment = trt.letters
    rhat.post.ij$full.treatment.name = trt.names
    rhat.post.ij$Rhat = model.sum.ij$fixed$Rhat
    rhat.df = rbind(rhat.df, rhat.post.ij)
    #rhat.df[which(rhat.df$pool == pool.i & rhat.df$species == spp.ij),] = rhat.post.ij
    write.csv(rhat.df, "Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Rhat_Statistic_10Chains.csv", row.names=F)
    print("Rhat file written")
    
    # posterior intervals
    post.ij = data.frame(matrix(nrow=n.t, ncol=10))
    colnames(post.ij) = c("pool","species","species.mean","treatment","full.treatment.name",
                          "posterior.mean","X5","X95","X25","X75")
    hdi.90 = hdi(m.ij, ci=0.90)
    hdi.50 = hdi(m.ij, ci=0.50)
    post.ij$pool = pool.i
    post.ij$species = spp.ij
    post.ij$species.mean = mean.list[[pool.i]][[spp.ij]]
    post.ij$treatment = trt.letters
    post.ij$full.treatment.name = trt.names
    post.ij$posterior.mean = fixef(m.ij)[,"Estimate"]
    post.ij[1:n.t,"X5"] = hdi.90[,"CI_low"]
    post.ij[1:n.t,"X95"] = hdi.90[,"CI_high"]
    post.ij[1:n.t,"X25"] = hdi.50[,"CI_low"]
    post.ij[1:n.t,"X75"] = hdi.50[,"CI_high"]
    df.int = rbind(df.int, post.ij)
    #df.int[which(df.int$pool == pool.i & df.int$species == spp.ij),] = post.ij
    write.csv(df.int, "Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Means_Intervals_10Chains.csv", row.names=F)
    print("Intervals written")
  }
}

################################################################################
# Run PCA's on species-specific C stock data for each pool and 
# make one combined figure

