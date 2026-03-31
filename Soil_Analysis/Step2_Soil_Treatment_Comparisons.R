path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/floodplain-experiment-repo"
setwd(path_to_repo)

library(brms)
library(bayestestR)

### script that reads quadrat-level soil data, runs statistical analyses on the
### variables of interest, and then writes posterior intervals to a file

################################################################################
# metadata

# set seed 
set.seed(2718)

# treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# number of plots per treatment
trt.plot.strip.df = read.csv("Metadata/Treatments_Strips_Plots.csv")
n.p = length(unique(trt.plot.strip.df$Plot))

# treamtent, strip, plot, and quadrat id variables
id.cols = c("full.treatment.name","treatment","strip","plot","quadrat")
n.id = length(id.cols)

# model names
models = c("simple","strip.random","strip.plot.random")
model.labels = c("Fixed effects only","Strip random effects","Strip and plot random effects")
n.m = length(models)

# information criteria
criteria = c("waic","loo")
criteria.labels = c("WAIC","LOO")
n.c = length(criteria)

# variable labels
soil.var.labs = read.csv("Metadata/Quadrat_Level_Soil_Variable_Metadata.csv")
soil.var.list = list()
for (i in 1:nrow(soil.var.labs)) {
  v1 = soil.var.labs[i,"variable"]
  v2 = soil.var.labs[i,"label"]
  soil.var.list[[v1]] = v2
}

################################################################################
# run HMC with BRMS

# read in soil data
soil.df = read.csv("Soil_Analysis/Clean_Data/Soil_Data_by_Quadrat_June2023.csv", header=T)

# quadrat-level variables to omit from statistical analysis
leave.out = c("texture.class","volumetric.moisture","h.meq","elevation.mean","vegetation.height.mean","fpom")
soil.df = soil.df[,-which(colnames(soil.df) %in% leave.out)]

# variable names
vars = colnames(soil.df)[(n.id+1):ncol(soil.df)]
n.v = length(vars)

# normalize all columns by respective variable means for Bayes analysis
all.means = list()
all.lists = list()
for (i in 1:n.v) {
  v = vars[i]
  all.means[[v]] = mean(soil.df[,v])
  all.lists[[v]] = list(treatment = factor(soil.df$treatment),
                        strip = factor(soil.df$strip),
                        plot = factor(soil.df$plot),
                        y = soil.df[,v]/all.means[[v]])
}

# run three models for each variable and store results in list
seed = 3141; n.iter=20000; n.chain=20
soil.model.list = list()
select.vars = c("mg.meq")
n.v = length(select.vars)
for (i in 1:n.v) { 
  v.i = select.vars[i]
  soil.model.list[[v.i]] = list()
  for (j in 1:n.m) {
    m.j = models[j]
    if (m.j == "simple") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ 0 + treatment,
                         family=gaussian(link="log"),
                         prior=prior(normal(0,1), class=b),
                         chains=n.chain, seed=seed, iter=n.iter)
    } else if (m.j == "strip.random") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ 0 + treatment + (1|treatment:strip),
                         family=gaussian(link="log"),
                         prior=prior(normal(0,1), class=b),
                         chains=n.chain, seed=seed, iter=n.iter,
                         control = list(adapt_delta = 0.99))
    } else if (m.j == "strip.plot.random") {
      model.fit.ij = brm(data = all.lists[[v.i]], 
                         y ~ 0 + treatment + (1|treatment:strip) + (1|treatment:strip:plot),
                         family=gaussian(link="log"),
                         prior=prior(normal(0,1), class=b),
                         chains=n.chain, seed=seed, iter=n.iter,
                         control = list(adapt_delta = 0.99))
    }
    model.fit.ij = add_criterion(model.fit.ij,"loo")
    model.fit.ij = add_criterion(model.fit.ij,"waic")
    soil.model.list[[v.i]][[m.j]] = model.fit.ij
  }
}


# save Rhat statistic for each variable and model
#rhat.df = read.csv("Soil_Analysis/Posteriors/Soil_Rhat_Statistic.csv")
rhat.df = data.frame(matrix(nrow=0,ncol=7))
colnames(rhat.df) = c("variable","variable.label","model","model.label",
                      "treatment","full.treatment.name","Rhat")
for (i in 1:n.v) {
  v.i = select.vars[i]
  for (j in 1:n.m) {
    m.j = models[j]
    model.sum.ij = summary(soil.model.list[[v.i]][[m.j]])
    rhat.df.ij = data.frame(matrix(nrow=6,ncol=7))
    colnames(rhat.df.ij) = c("variable","variable.label","model","model.label",
                             "treatment","full.treatment.name","Rhat")
    rhat.df.ij$variable = v.i
    rhat.df.ij$variable.label = soil.var.list[[v.i]]
    rhat.df.ij$model = m.j
    rhat.df.ij$model.label = model.labels[j]
    rhat.df.ij$treatment = trt.letters
    rhat.df.ij$full.treatment.name = trt.names
    rhat.df.ij$Rhat = model.sum.ij$fixed$Rhat
    #rhat.df[rhat.df$variable == v.i & rhat.df$model == m.j,] = rhat.df.ij
    rhat.df = rbind(rhat.df, rhat.df.ij)
  }
}
write.csv(rhat.df, "Soil_Analysis/Posteriors/Soil_Rhat_Statistic.csv", row.names=F)

# compare models for each variale with WAIC, LOO, and kfold
#comp.df = read.csv("Soil_Analysis/Posteriors/Soil_Model_Information_Criteria.csv")
comp.df = data.frame(matrix(nrow=0, ncol=13))
colnames(comp.df) = c("elpd_diff","se_diff","elpd","se_elpd","p","se_p","ic","se_ic",
                      "variable","variable.label","criterion","model","best.model")
for (i in 1:n.v) {
  v.i = select.vars[i]
  for (j in 1:n.c) {
    m1 = soil.model.list[[v.i]][[models[1]]]
    m2 = soil.model.list[[v.i]][[models[2]]]
    m3 = soil.model.list[[v.i]][[models[3]]]
    comp.ij = data.frame(loo_compare(m1, m2, m3, criterion = tolower(criteria)[j], model_names=models))
    colnames(comp.ij) = c("elpd_diff","se_diff","elpd","se_elpd","p","se_p","ic","se_ic")
    comp.ij$variable = v.i
    comp.ij$variable.label = soil.var.list[[v.i]]
    comp.ij$criterion = criteria.labels[j]
    comp.ij$model = row.names(comp.ij)
    comp.ij$best.model = c("Yes","No","No")
    #comp.df[comp.df$variable == v.i & comp.df$criterion == toupper(criteria[j]),1:13] = comp.ij
    comp.df = rbind(comp.df, comp.ij)
  }
}
comp.df$model.label = 0
for (i in 1:n.m) { comp.df$model.label[which(comp.df$model == models[i])] = model.labels[i] }
write.csv(comp.df, "Soil_Analysis/Posteriors/Soil_Model_Information_Criteria.csv", row.names=F)

# get posterior intervals and write to file
#df.int = read.csv("Soil_Analysis/Posteriors/Soil_Posterior_Intervals_10Chains_NaturalScale.csv")
df.int = data.frame(matrix(nrow=0, ncol=12))
colnames(df.int) = c("model","model.label","variable","variable.label","variable.mean",
                     "treatment","full.treatment.name","posterior.mean","X5","X95","X25","X75")
for (i in 1:n.v) {
  v.i = select.vars[i]
  for (j in 1:n.m) {
    m.j = models[j]
    fit.ij = soil.model.list[[v.i]][[m.j]]
    df.ij = data.frame(matrix(nrow=n.t, ncol=12))
    colnames(df.ij) = c("model","model.label","variable","variable.label","variable.mean",
                        "treatment","full.treatment.name","posterior.mean","X5","X95","X25","X75")
    hdi.90 = hdi(fit.ij, ci=0.90)
    hdi.50 = hdi(fit.ij, ci=0.50)
    df.ij$model = models[j]
    df.ij$model.label = model.labels[j]
    df.ij$variable = v.i
    df.ij$variable.label = soil.var.list[[v.i]]
    df.ij$variable.mean = all.means[[v.i]]
    df.ij$treatment = trt.letters
    df.ij$full.treatment.name = trt.names
    df.ij$posterior.mean = exp(fixef(fit.ij, parts="treatment")[,"Estimate"])*all.means[[v.i]]
    df.ij[1:n.t,"X5"] = exp(hdi.90[,"CI_low"])*all.means[[v.i]]
    df.ij[1:n.t,"X95"] = exp(hdi.90[,"CI_high"])*all.means[[v.i]]
    df.ij[1:n.t,"X25"] = exp(hdi.50[,"CI_low"])*all.means[[v.i]]
    df.ij[1:n.t,"X75"] = exp(hdi.50[,"CI_high"])*all.means[[v.i]]
    df.int[which(df.int$variable == v.i & df.int$model == m.j),] = df.ij
    #df.int = rbind(df.int, df.ij)
  }
}
df.int$model.label = rep(0, nrow(df.int))
for (i in 1:n.m) { df.int$model.label[which(df.int$model == models[i])] = model.labels[i] }
write.csv(df.int, "Soil_Analysis/Posteriors/Soil_Posterior_Intervals_10Chains_NaturalScale.csv", row.names=F)













