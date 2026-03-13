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

# variable names
vars = read.csv("Metadata/Plot_Level_Carbon_Richness_Variables.csv")[,"variable"]
var.labels = read.csv("Metadata/Plot_Level_Carbon_Richness_Variables.csv", stringsAsFactors=F)[,"label"]

# variables to leave out
leave.out = c("small.fwd.carbon","intermediate.fwd.carbon","dead.stem.carbon","large.snag.carbon","tc.stock",
              "abg.live.tree.carbon","bg.live.tree.carbon","total.live.tree.carbon",
              "abg.live.stem.carbon","bg.live.stem.carbon","total.live.stem.carbon")
omit.var.id = which(vars %in% leave.out)
vars = vars[-omit.var.id]
var.labels = var.labels[-omit.var.id]
n.v = length(vars)

# model types
models = c("simple","strip.random")
model.labels = c("Fixed effects only","Strip random effects")
n.m = length(models)

################################################################################
## statistical analyses on live and dead C stock components by plot

# read in biomass C stock estimates
stock.richness.df = read.csv("Tree_Analysis/Clean_Data_By_Plot/Clean_Veg_Soil_C_Stocks_Richness_by_Plot.csv")

# make simplified data lists
all.lists = list()
all.means = list()
for (i in 1:n.v) {
  v = vars[i]
  all.means[[v]] = mean(stock.richness.df[,v])
  all.lists[[v]] = list(treatment = factor(stock.richness.df$treatment), 
                        strip = factor(stock.richness.df$strip), 
                        plot = factor(stock.richness.df$plot),
                        y = stock.richness.df[,v]/all.means[[v]])
  hist(all.lists[[v]]$y, main=v, xlim=c(0,6))
}

# fit simple models
seed = 3141; n.iter=10000; n.chain=10
stock.model.list = list()
for (i in 1:n.v) {
  v.i = vars[i]
  stock.model.list[[v.i]] = list()
  for (j in 1:n.m) {
    m.j = models[j]
    if (m.j == "simple") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                                y ~ treatment,
                                family = gaussian(link="log"),
                                prior=prior(normal(0,1), class=b),
                                chains=n.chain, seed=seed, iter=n.iter,
                                control = list(adapt_delta = 0.99,
                                               max_treedepth = 12))
    } else if (m.j == "strip.random") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ treatment + (1|treatment:strip),
                         family=gaussian(link="log"),
                         prior=prior(normal(0,1), class=b),
                         chains=n.chain, seed=seed, iter=n.iter,
                         control = list(adapt_delta = 0.99,
                                        max_treedepth = 12))
    }
    model.fit.ij = add_criterion(model.fit.ij,"loo")
    model.fit.ij = add_criterion(model.fit.ij,"waic")
    stock.model.list[[v.i]][[m.j]] = model.fit.ij
  }
}

# save Rhat statistic for each variable and model
rhat.df = data.frame(matrix(nrow=0,ncol=7))
colnames(rhat.df) = c("variable","variable.label","model","model.label",
                      "treatment","full.treatment.name","Rhat")
for (i in 1:n.v) {
  v.i = vars[i]
  for (j in 1:n.m) {
    m.j = models[j]
    model.sum.ij = summary(stock.model.list[[v.i]][[m.j]])
    model.df.ij = data.frame(matrix(nrow=6,ncol=7))
    colnames(model.df.ij) = c("variable","variable.label","model","model.label",
                              "treatment","full.treatment.name","Rhat")
    model.df.ij$variable = v.i
    model.df.ij$variable.label = var.labels[i]
    model.df.ij$model = m.j
    model.df.ij$model.label = model.labels[j]
    model.df.ij$treatment = trt.letters
    model.df.ij$full.treatment.name = trt.names
    model.df.ij$Rhat = model.sum.ij$fixed$Rhat
    rhat.df = rbind(rhat.df, model.df.ij)
  }
}
rhat.df$model.label = rep(0, nrow(rhat.df))
for (i in 1:n.m) { rhat.df$model.label[which(rhat.df$model == models[i])] = model.labels[i] }
write.csv(rhat.df, "Tree_Analysis/Posteriors/Stock_Rhat_Statistic.csv", row.names=F)

# compare models for each variale with WAIC and LOO
criteria = c("waic","loo")
criteria.labels = c("WAIC","LOO")
n.c = length(criteria)
comp.df = data.frame(matrix(nrow=0, ncol=13))
colnames(comp.df) = c("elpd_diff","se_diff","elpd","se_elpd","p","se_p","ic","se_ic",
                      "variable","variable.label","criterion","model","best.model")
for (i in 1:n.v) {
  v.i = vars[i]
  for (j in 1:n.c) {
    m1 = stock.model.list[[v.i]][[models[1]]]
    m2 = stock.model.list[[v.i]][[models[2]]]
    comp.ij = data.frame(loo_compare(m1, m2, criterion = criteria[j], model_names=models))
    colnames(comp.ij) = c("elpd_diff","se_diff","elpd","se_elpd","p","se_p","ic","se_ic")
    comp.ij$variable = v.i
    comp.ij$variable.label = var.labels[i]
    comp.ij$criterion = criteria.labels[j]
    comp.ij$model = row.names(comp.ij)
    comp.ij$best.model = c("Yes","No")
    comp.df = rbind(comp.df, comp.ij)
  }
}

# write information criterion comparisons to file
comp.df$model.label = rep(0, nrow(comp.df))
for (i in 1:n.m) { comp.df$model.label[which(comp.df$model == models[i])] = model.labels[i] }
write.csv(comp.df, "Tree_Analysis/Posteriors/Stock_Model_Information_Criteria.csv", row.names=F)

# get posterior intervals and write to file
df.int = data.frame(matrix(nrow=0, ncol=12))
int.cols = c("model","model.label","variable","variable.label","variable.mean","treatment",
             "full.treatment.name","posterior.mean","5","95","25","75")
colnames(df.int) = int.cols
for (i in 1:n.v) {
  v.i = vars[i]
  for (j in 1:n.m) {
    m.j = models[j]
    fit.ij = stock.model.list[[v.i]][[m.j]]
    df.ij = data.frame(matrix(nrow=n.t, ncol=12))
    colnames(df.ij) = int.cols
    hdi.90 = hdi(fit.ij, ci=0.90)
    hdi.50 = hdi(fit.ij, ci=0.50)
    df.ij$model = models[j]
    df.ij$model.label = model.labels[j]
    df.ij$variable = v.i
    df.ij$variable.label = var.labels[i]
    df.ij$variable.mean = all.means[[v.i]]
    df.ij$treatment = trt.letters
    df.ij$full.treatment.name = trt.names
    df.ij$posterior.mean = exp(fixef(fit.ij)[,"Estimate"])*all.means[[v.i]]
    df.ij[1:n.t,"5"] = exp(hdi.90[,"CI_low"])*all.means[[v.i]]
    df.ij[1:n.t,"95"] = exp(hdi.90[,"CI_high"])*all.means[[v.i]]
    df.ij[1:n.t,"25"] = exp(hdi.50[,"CI_low"])*all.means[[v.i]]
    df.ij[1:n.t,"75"] = exp(hdi.50[,"CI_high"])*all.means[[v.i]]
    df.int = rbind(df.int, df.ij)
  }
}
write.csv(df.int,"Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_10Chains_NaturalScale.csv", row.names=F)

