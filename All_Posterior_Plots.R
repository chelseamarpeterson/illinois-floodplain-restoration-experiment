path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/floodplain-experiment-repo"
setwd(path_to_repo)

library(ggplot2)
library(patchwork)
library(RColorBrewer)
library(dplyr)
library(tidyr)
library(reshape2)
library(ggtext)
library(scales)

################################################################################
# load metadata

# treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# model names
models = c("simple","strip.random","strip.plot.random")
model.labels = c("Fixed effects only","Strip random effects","Strip and plot random effects")
n.m = length(models)

# criteria names
criteria = c("WAIC","LOO")
n.c = length(criteria)

# variable labels
soil.rhat.df = read.csv("Soil_Analysis/Posteriors/Soil_Rhat_Statistic.csv")
soil.vars = unique(soil.rhat.df$variable)
soil.var.labs = unique(soil.rhat.df$variable.label)
n.s.v = length(soil.vars)

# stock and richness variables
stock.rhat.df = read.csv("Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Rhat_Statistic.csv")
stock.vars = unique(stock.rhat.df$variable)
stock.var.labs = unique(stock.rhat.df$variable.label)
n.c.v = length(stock.vars)

################################################################################
# evaluate model convergence (Figures B1-B3)

# soil carbon variables: plot rhat's for each model and variable
soil.carbon.vars = c("tc.percent","tic.percent","toc.percent","maoc.percent","poc.percent",
                     "tc.stock","tic.stock","toc.stock","maoc.stock","poc.stock")
soil.carbon.labs = 0
for (i in 1:10) { soil.carbon.labs[i] = soil.rhat.df[soil.rhat.df$variable == soil.carbon.vars[i],"variable.label"][1] }
soil.carbon.rhat.df = soil.rhat.df[soil.rhat.df$variable %in% soil.carbon.vars,]
p.soil.carbon.rhat.comp = ggplot(soil.carbon.rhat.df, 
                          aes(x=(Rhat-1)*100, 
                              y=factor(full.treatment.name, levels=trt.names),
                              color=factor(model.label, levels=model.labels),
                              shape=factor(model.label, levels=model.labels))) + 
                          geom_vline(xintercept=0) +     
                          geom_point() + 
                          scale_shape_discrete(solid = F) +
                          labs(y="",
                               x="(R-hat statistic - 1)*100",
                               color="Model", 
                               shape="Model") +
                          geom_vline(xintercept=1) +
                          scale_x_continuous(breaks=seq(0,1,0.2),
                                             limits=c(-0.01,1)) + 
                          facet_wrap(.~factor(variable.label, levels=soil.carbon.labs), 
                                     ncol=2, dir="v")
p.soil.carbon.rhat.comp
ggsave("Supp_Figures/FigureB1_Soil_Carbon_Rhat_Score_Comparison.jpeg", 
       plot=p.soil.carbon.rhat.comp, width=18, height=16, units="cm",dpi=600)

# non-carbon soil variables: plot rhat's for each model and variable
soil.noncarbon.vars = soil.vars[-which(soil.vars %in% soil.carbon.vars)]
soil.noncarbon.labs = 0
for (i in 1:length(soil.noncarbon.vars)) { soil.noncarbon.labs[i] = soil.rhat.df[soil.rhat.df$variable == soil.noncarbon.vars[i],"variable.label"][1] }
soil.noncarbon.labs[soil.noncarbon.labs == "Root and wood fragment density (g/cm3)"] = "Root & wood\nfragment density (g/cm3)"
soil.rhat.df[soil.rhat.df$variable.label == "Root and wood fragment density (g/cm3)","variable.label"] = "Root & wood\nfragment density (g/cm3)"
soil.noncarbon.rhat.df = soil.rhat.df[soil.rhat.df$variable %in% soil.noncarbon.vars,]
p.soil.noncarbon.rhat.comp = ggplot(soil.noncarbon.rhat.df, 
                                    aes(x=(Rhat-1)*100, 
                                        y=factor(full.treatment.name, levels=trt.names),
                                        color=factor(model.label, levels=model.labels),
                                        shape=factor(model.label, levels=model.labels))) + 
                                    geom_vline(xintercept=0) +     
                                    geom_point() + 
                                    scale_shape_discrete(solid = F) +
                                    labs(y="",
                                         x="(R-hat statistic - 1)*100",
                                         color="Model",
                                         shape="Model") +
                                    geom_vline(xintercept=1) +
                                    scale_x_continuous(breaks=seq(0,1,0.2),
                                                       limits=c(-0.01,1)) + 
                                    facet_wrap(.~factor(variable.label, levels=soil.noncarbon.labs), ncol=4)
p.soil.noncarbon.rhat.comp
ggsave("Supp_Figures/FigureB2_Non_Carbon_Soil_Rhat_Score_Comparison.jpeg", 
       plot=p.soil.noncarbon.rhat.comp, width=26, height=24, units="cm",dpi=600)

# carbon stocks: plot rhat's for each model and variable
stock.rhat.df$variable.label[stock.rhat.df$variable.label == "Both tree and herbaceous layer"] = "Both tree & herbaceous\nlayer richness"
stock.rhat.df$variable.label[stock.rhat.df$variable.label == "Live F. pennsylvanica (>= 2.5 cm)"] = "Live F. pennsylvanica\n(>= 2.5 cm)"
stock.rhat.df$variable.label[stock.rhat.df$variable.label == "Dead F. pennsylvanica (>= 2.5 cm)"] = "Dead F. pennsylvanica\n(>= 2.5 cm)"
stock.rhat.df$variable.label[stock.rhat.df$variable.label == "Fine woody debris (< 7.6 cm)"] = "Fine woody debris\n(< 7.6 cm)"
stock.rhat.df$variable.label[stock.rhat.df$variable.label == "Coarse woody debris (>= 7.6 cm)"] = "Coarse woody debris\n(>= 7.6 cm)"
stock.var.labs[stock.var.labs == "Both tree and herbaceous layer"] = "Both tree & herbaceous\nlayer richness"
stock.var.labs[stock.var.labs == "Live F. pennsylvanica (>= 2.5 cm)"] = "Live F. pennsylvanica\n(>= 2.5 cm)"
stock.var.labs[stock.var.labs == "Dead F. pennsylvanica (>= 2.5 cm)"] = "Dead F. pennsylvanica\n(>= 2.5 cm)"
stock.var.labs[stock.var.labs == "Fine woody debris (< 7.6 cm)"] = "Fine woody debris\n(< 7.6 cm)"
stock.var.labs[stock.var.labs == "Coarse woody debris (>= 7.6 cm)"] = "Coarse woody debris\n(>= 7.6 cm)"
p.stock.rhat.comp = ggplot(stock.rhat.df, 
                           aes(x=(Rhat-1)*100, 
                               y=factor(full.treatment.name, levels=trt.names),
                               color=factor(model.label, levels=model.labels),
                                            shape=factor(model.label, levels=model.labels))) + 
                           geom_vline(xintercept=0) +     
                           geom_point() + 
                           scale_shape_discrete(solid = F) +
                           labs(y="",
                                x="(R-hat statistic - 1)*100",
                                color="Model",
                                shape="Model") +
                           geom_vline(xintercept=1) +
                           scale_x_continuous(breaks=seq(0,1,0.2),
                                              limits=c(-0.01,1)) + 
                           facet_wrap(.~factor(variable.label, levels=stock.var.labs), ncol=4)
p.stock.rhat.comp
ggsave("Supp_Figures/FigureB3_Carbon_Richness_Rhat_Score_Comparison.jpeg", 
       plot=p.stock.rhat.comp, width=26, height=24, units="cm",dpi=600)

################################################################################
# count the number of variables for which each model type minimizes the loo v. waic (Tables B1 & B2)

# read in soil data model comparison dataframe
comp.df.soil = read.csv("Soil_Analysis/Posteriors/Soil_Model_Information_Criteria.csv")
comp.df.stocks = read.csv("Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Model_Information_Criteria.csv")

# soil
soil.min.ic.count.df = data.frame(matrix(nrow=n.c, ncol=n.m))
colnames(soil.min.ic.count.df) = model.labels
row.names(soil.min.ic.count.df) = criteria
for (i in 1:n.c) {
  criterion.i = criteria[i]
  soil.min.ic.count.df[i,] = 0
  for (j in 1:n.s.v) {
    var.j = soil.vars[j]
    row.ij.ind = which(comp.df.soil$variable == var.j & comp.df.soil$criterion == criterion.i)
    row.ij = comp.df.soil[row.ij.ind,]
    min.ic.model = row.ij[which.min(row.ij$ic), "model.label"]
    soil.min.ic.count.df[i,min.ic.model] = soil.min.ic.count.df[i,min.ic.model] + 1
  }
}

# carbon stocks
stock.min.ic.count.df = data.frame(matrix(nrow=n.c, ncol=2))
colnames(stock.min.ic.count.df) = model.labels[1:2]
row.names(stock.min.ic.count.df) = criteria
for (i in 1:n.c) {
  criterion.i = criteria[i]
  stock.min.ic.count.df[i,] = 0
  for (j in 1:n.c.v) {
    var.j = stock.vars[j]
    row.ij.ind = which(comp.df.stocks$variable == var.j & comp.df.stocks$criterion == criterion.i)
    row.ij = comp.df.stocks[row.ij.ind,]
    min.ic.model = row.ij[which.min(row.ij$ic), "model.label"]
    stock.min.ic.count.df[i,min.ic.model] = stock.min.ic.count.df[i,min.ic.model] + 1
  }
}

################################################################################
# compare model information criteria (Figures B4-B6)

# soil carbon properties: plot comparison of information criteria
soil.carbon.labs = 0
for (i in 1:10) { soil.carbon.labs[i] = comp.df.soil[comp.df.soil$variable == soil.carbon.vars[i],"variable.label"][1] }
comp.df.soil.carbon = comp.df.soil[comp.df.soil$variable %in% soil.carbon.vars,]
p.carbon.soil.ic.score.comp = ggplot(comp.df.soil.carbon, 
                              aes(x=ic, y=factor(model.label, levels=model.labels), 
                                  color=criterion,
                                  linetype=criterion,
                                  shape=factor(best.model))) + 
                              geom_point(position=position_dodge(0.75),size=1.75) + 
                              geom_errorbarh(aes(xmin=ic-se_ic, xmax=ic+se_ic,
                                                 y=factor(model.label, levels=model.labels), 
                                                 color=criterion,
                                                 linetype=criterion),
                                             position=position_dodge(0.75), height=0.3) +
                              facet_wrap(.~factor(variable.label, levels=soil.carbon.labs), 
                                         scales="free_x", ncol=2, dir="v") + 
                              theme(text = element_text(size=14)) + 
                              labs(y="", x="Score",
                                   color="Information criterion",
                                   linetype="Information criterion",
                                   shape="Best model")
p.carbon.soil.ic.score.comp
ggsave("Supp_Figures/FigureB4_Soil_Carbon_IC_Score_Comparison.jpeg", 
       plot=p.carbon.soil.ic.score.comp, width=20, height=18, units="cm",dpi=600)

# soil non-carbon properties: plot comparison of information criteria
soil.noncarbon.vars = soil.vars[-which(soil.vars %in% soil.carbon.vars)]
soil.noncarbon.labs = 0
for (i in 1:length(soil.noncarbon.vars)) { soil.noncarbon.labs[i] = comp.df.soil[comp.df.soil$variable == soil.noncarbon.vars[i],"variable.label"][1] }
soil.noncarbon.labs[soil.noncarbon.labs == "Root and wood fragment density (g/cm3)"] = "Root & wood\nfragment density (g/cm3)"
comp.df.soil[comp.df.soil$variable.label == "Root and wood fragment density (g/cm3)","variable.label"] = "Root & wood\nfragment density (g/cm3)"
comp.df.soil.noncarbon = comp.df.soil[-which(comp.df.soil$variable %in% soil.carbon.vars),]
p.noncarbon.soil.ic.score.comp = ggplot(comp.df.soil.noncarbon, 
                                        aes(x=ic, y=factor(model.label, levels=model.labels), 
                                            color=criterion,
                                            linetype=criterion,
                                            shape=factor(best.model))) + 
                                        geom_point(position=position_dodge(0.75),size=1.75) + 
                                        geom_errorbarh(aes(xmin=ic-se_ic, xmax=ic+se_ic,
                                                           y=factor(model.label, levels=model.labels), 
                                                           color=criterion,
                                                           linetype=criterion),
                                                       position=position_dodge(0.75), height=0.3) +
                                        facet_wrap(.~factor(variable.label, levels=soil.noncarbon.labs), 
                                                   scales="free_x", ncol=4) + 
                                        theme(text = element_text(size=14)) + 
                                        labs(y="", x="Score",color="Information criterion",
                                             linetype="Information criterion",shape="Best model")
p.noncarbon.soil.ic.score.comp
ggsave("Supp_Figures/FigureB5_NonCarbon_Soil_IC_Score_Comparison.jpeg", 
       plot=p.noncarbon.soil.ic.score.comp, width=32, height=30, units="cm",dpi=600)

# soil: plot comparison of information criteria
comp.df.stocks$variable.label[comp.df.stocks$variable.label == "Both tree and herbaceous layer"] = "Both tree & herbaceous\nlayer richness"
comp.df.stocks$variable.label[comp.df.stocks$variable.label == "Live F. pennsylvanica (>= 2.5 cm)"] = "Live F. pennsylvanica\n(>= 2.5 cm)"
comp.df.stocks$variable.label[comp.df.stocks$variable.label == "Dead F. pennsylvanica (>= 2.5 cm)"] = "Dead F. pennsylvanica\n(>= 2.5 cm)"
comp.df.stocks$variable.label[comp.df.stocks$variable.label == "Fine woody debris (< 7.6 cm)"] = "Fine woody debris\n(< 7.6 cm)"
comp.df.stocks$variable.label[comp.df.stocks$variable.label == "Coarse woody debris (>= 7.6 cm)"] = "Coarse woody debris\n(>= 7.6 cm)"
p.stock.ic.score.comp = ggplot(comp.df.stocks, 
                               aes(x=ic, y=factor(model.label, levels=model.labels), 
                                   color=criterion,
                                   linetype=criterion,
                                   shape=factor(best.model))) + 
                               geom_point(position=position_dodge(0.75),size=1.75) + 
                               geom_errorbarh(aes(xmin=ic-se_ic, xmax=ic+se_ic,
                                                  y=factor(model.label, levels=model.labels), 
                                                  color=criterion,
                                                  linetype=criterion),
                                              position=position_dodge(0.75), height=0.3) +
                               facet_wrap(.~factor(variable.label, levels=stock.var.labs), 
                                          ncol=4, scales="free_x") + 
                               theme(text = element_text(size=14)) + 
                               labs(y="", x="Score",color="Information criterion",
                                    linetype="Information criterion",shape="Best model")
p.stock.ic.score.comp
ggsave("Supp_Figures/FigureB6_Stock_Richness_IC_Score_Comparison.jpeg", 
       plot=p.stock.ic.score.comp, width=32, height=30, units="cm",dpi=600)

###############################################################################
# plot HDIs for each model and variable (Figures B3 & B4)

# read posterior distributions for soil properties, carbon stocks, and species richness
soil.hdi.df = read.csv("Soil_Analysis/Posteriors/Soil_Posterior_Intervals_10Chains_NaturalScale.csv")
stock.hdi.df = read.csv("Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_10Chains_NaturalScale.csv")

# soil carbon: plot posterior HDIs for all soil variables and models
soil.carbon.labs = 0
for (i in 1:10) { soil.carbon.labs[i] = soil.hdi.df[soil.hdi.df$variable == soil.carbon.vars[i],"variable.label"][1] }
hdi.df.soil.carbon = soil.hdi.df[soil.hdi.df$variable %in% soil.carbon.vars,]
p.soil.carbon.hdi.comp = ggplot(hdi.df.soil.carbon, 
                         aes(y=factor(full.treatment.name, levels=trt.names),
                             x=posterior.mean,
                             color=factor(model.label, levels=model.labels),
                             shape=factor(model.label, levels=model.labels))) +
                         geom_point(position=position_dodge(0.8),size=1.5) +
                         geom_errorbarh(aes(xmin=`X5`,xmax=`X95`,
                                            y=factor(full.treatment.name, levels=trt.names),
                                            color=factor(model.label, levels=model.labels)),
                                        position=position_dodge(0.8),height=0.4,linewidth=0.4) +
                         geom_errorbarh(aes(xmin=`X25`,xmax=`X75`,
                                            y=factor(full.treatment.name, levels=trt.names),
                                            color=factor(model.label, levels=model.labels)),
                                        position=position_dodge(0.8),height=0,linewidth=0.8) +
                         facet_wrap(.~factor(variable.label, levels=soil.carbon.labs),
                                    ncol=2, scales="free_x", dir = "v") +
                         theme(text = element_text(size=12),
                               panel.spacing.x = unit(1,"lines")) + 
                         labs(y="", x="Posterior estimate", color="Model", shape="Model")
p.soil.carbon.hdi.comp
ggsave("Supp_Figures/FigureB7_Soil_Carbon_HDI_Comparison.jpeg", 
       plot=p.soil.carbon.hdi.comp, width=20, height=18, units="cm",dpi=600)

# soil non-carbon: plot posterior HDIs for all soil variables and models
soil.noncarbon.vars = soil.vars[-which(soil.vars %in% soil.carbon.vars)]
soil.noncarbon.labs = 0
for (i in 1:length(soil.noncarbon.vars)) { soil.noncarbon.labs[i] = soil.hdi.df[soil.hdi.df$variable == soil.noncarbon.vars[i],"variable.label"][1] }
soil.noncarbon.labs[soil.noncarbon.labs == "Root and wood fragment density (g/cm3)"] = "Root & wood\nfragment density (g/cm3)"
soil.hdi.df[soil.hdi.df$variable.label == "Root and wood fragment density (g/cm3)","variable.label"] = "Root & wood\nfragment density (g/cm3)"
hdi.df.soil.noncarbon = soil.hdi.df[-which(soil.hdi.df$variable %in% soil.carbon.vars),]
p.soil.noncarbon.hdi.comp = ggplot(hdi.df.soil.noncarbon, 
                                   aes(y=factor(full.treatment.name, levels=trt.names),
                                       x=posterior.mean,
                                       color=factor(model.label, levels=model.labels),
                                       shape=factor(model.label, levels=model.labels))) +
                                   geom_point(position=position_dodge(0.8),size=1.5) +
                                   geom_errorbarh(aes(xmin=`X5`,xmax=`X95`,
                                                      y=factor(full.treatment.name, levels=trt.names),
                                                      color=factor(model.label, levels=model.labels)),
                                                  position=position_dodge(0.8),height=0.4,linewidth=0.4) +
                                   geom_errorbarh(aes(xmin=`X25`,xmax=`X75`,
                                                      y=factor(full.treatment.name, levels=trt.names),
                                                      color=factor(model.label, levels=model.labels)),
                                                  position=position_dodge(0.8),height=0,linewidth=0.8) +
                                   facet_wrap(.~factor(variable.label, levels=soil.noncarbon.labs),
                                              ncol=4, scales="free_x", dir = "h") + 
                                   theme(text = element_text(size=12),
                                         panel.spacing.x = unit(1,"lines")) + 
                                   labs(y="", x="Posterior estimate", color="Model", shape="Model")
p.soil.noncarbon.hdi.comp
ggsave("Supp_Figures/FigureB8_NonCarbon_Soil_HDI_Comparison.jpeg", 
       plot=p.soil.noncarbon.hdi.comp, width=32, height=30, units="cm",dpi=600)

# soil: plot posterior HDIs for all soil variables and models
stock.hdi.df$variable.label[stock.hdi.df$variable.label == "Both tree and herbaceous layer"] = "Both tree & herbaceous\nlayer richness"
stock.hdi.df$variable.label[stock.hdi.df$variable.label == "Live F. pennsylvanica (>= 2.5 cm)"] = "Live F. pennsylvanica\n(>= 2.5 cm)"
stock.hdi.df$variable.label[stock.hdi.df$variable.label == "Dead F. pennsylvanica (>= 2.5 cm)"] = "Dead F. pennsylvanica\n(>= 2.5 cm)"
stock.hdi.df$variable.label[stock.hdi.df$variable.label == "Fine woody debris (< 7.6 cm)"] = "Fine woody debris\n(< 7.6 cm)"
stock.hdi.df$variable.label[stock.hdi.df$variable.label == "Coarse woody debris (>= 7.6 cm)"] = "Coarse woody debris\n(>= 7.6 cm)"
p.stock.hdi.comp = ggplot(stock.hdi.df, 
                          aes(y=factor(full.treatment.name, levels=trt.names),
                              x=posterior.mean,
                              color=factor(model.label, levels=model.labels),
                              shape=factor(model.label, levels=model.labels))) +
                          geom_point(position=position_dodge(0.8),size=1.75) +
                          geom_errorbarh(aes(xmin=`X5`,xmax=`X95`,
                                             y=factor(full.treatment.name, levels=trt.names),
                                             color=factor(model.label, levels=model.labels)),
                                         position=position_dodge(0.8),height=0.4,linewidth=0.4) +
                          geom_errorbarh(aes(xmin=`X25`,xmax=`X75`,
                                             y=factor(full.treatment.name, levels=trt.names),
                                             color=factor(model.label, levels=model.labels)),
                                         position=position_dodge(0.8),height=0,linewidth=0.8) +
                          facet_wrap(.~factor(variable.label, levels=stock.var.labs),
                                     ncol=4, scales="free_x") + 
                          theme(text = element_text(size=12)) +
                          labs(y="", x="Posterior estimate", color="Model", shape="Model")
p.stock.hdi.comp
ggsave("Supp_Figures/FigureB9_Stock_Richness_HDI_Comparison.jpeg", 
       plot=p.stock.hdi.comp, width=32, height=30, units="cm",dpi=600)

################################################################################
# make combined plot for carbon concentrations and stocks (Figure 3)

# IPCC estimates
ipcc.df = read.csv("Carbon_Calculations/IPCC_Carbon_Estimates.csv")

# isolate soil data for best model
soil.hdi.df.best = soil.hdi.df[which(soil.hdi.df$model == "strip.plot.random"),]
stock.hdi.df.best = stock.hdi.df[which(stock.hdi.df$model == "strip.random"),]

# Georgiu MAOC estimate
#soc.df = data.frame(matrix(nrow=1,ncol=2))
#soc.df$soc.type = "Georgiou et al. 2022 (MAOC)"
#soc.df$soc.value = 60.961376

# stacked plot for carbon fractions
c.conc.labs = c("tic.percent","poc.percent","maoc.percent")
c.conc.labs.new = c("TIC","POC","MAOC")
for (i in 1:3) { soil.hdi.df.best$variable.label[soil.hdi.df.best$variable == c.conc.labs[i]] = c.conc.labs.new[i] }
s.palette <- c(brewer.pal(9,"Greys")[5],brewer.pal(11,"BrBG")[c(2,1)])
p.c.concentrations = ggplot(soil.hdi.df.best[soil.hdi.df.best$variable.label %in% c.conc.labs.new,], 
                            aes(y=factor(full.treatment.name, levels=trt.names), 
                                x=posterior.mean, 
                                fill=factor(variable.label, levels=c.conc.labs.new))) + 
                            geom_bar(stat="identity",
                                     position="stack") +
                            labs(x="Concentration (% [g C/g soil])",
                                 y="",fill="Soil carbon fraction") + 
                            scale_fill_manual(values=s.palette) +
                            theme(text = element_text(size=14),
                                  plot.margin=unit(c(1,1,1,1),"lines")) +
                            coord_cartesian(xlim = c(0,4.93), clip="off") +
                            scale_x_continuous(breaks=seq(0,5)) +
                            geom_label(x=4.8,y=6,label="a",
                                       color="black",fill=alpha("white",0.9),
                                       label.r=unit(0,"pt"),label.size=0,
                                       size=8,fontface="bold") + 
                            guides(fill="none")
p.c.concentrations
ipcc.soil.vars = c("Annual crops","Revegetated cropland","Natural wetland")
c.stock.labs = c("tic.stock","poc.stock","maoc.stock")
c.stock.labs.new = c("TIC","POC","MAOC")
for (i in 1:3) { soil.hdi.df.best$variable.label[soil.hdi.df.best$variable == c.stock.labs[i]] = c.stock.labs.new[i] }
p.c.stocks = ggplot(soil.hdi.df.best[soil.hdi.df.best$variable.label %in% c.stock.labs.new,], 
                    aes(y=factor(full.treatment.name, levels=trt.names), 
                        x=posterior.mean, 
                        fill=factor(variable.label, levels=c.stock.labs.new))) +
                    geom_bar(stat="identity",position="stack") + 
                    scale_fill_manual(values=s.palette) +
                    labs(fill="",y="",x="Stock (Mg C/ha)",title="") +
                    geom_vline(data=ipcc.df, 
                               aes(xintercept=soc.value, 
                                   color=factor(soc.type, levels=ipcc.soil.vars),
                                   linetype=factor(soc.type, levels=ipcc.soil.vars)), 
                               linewidth=1.5) +
                    scale_color_manual(values=c("red","yellow1","royalblue1")) + 
                    scale_linetype_manual(values=c("solid","dotdash","dashed")) +
                    guides(color=guide_legend(title="IPCC soil organic carbon"),
                           linetype=guide_legend(title="IPCC soil organic carbon")) +
                    theme(text=element_text(size=14), 
                          axis.text.y=element_blank(),
                          legend.key.size=unit(0.7,'cm'),
                          plot.margin=unit(c(1,1,1,1),"lines"),
                          legend.key=element_rect(fill="darkgrey")) +
                    coord_cartesian(xlim=c(0,143), clip="off") +
                    scale_x_continuous(breaks=seq(0,125,25)) +
                    geom_label(x=138.8,y=6,label="b",
                               color="black",fill=alpha("white",0.9),
                               label.r=unit(0,"pt"),label.size=0,
                               size=8,fontface="bold")
p.c.stocks
p.c.all = p.c.concentrations + p.c.stocks
p.c.all =  p.c.all + theme(plot.margin = margin(0,0,0,0, "cm"))
p.c.all
ggsave("Main_Figures/Figure3_Soil_Carbon_Concentrations_and_Stocks.jpeg", 
       plot=p.c.all, width=30, height=12, units="cm",dpi=600)

################################################################################
# plot non-carbon chemical and physical variables (Figure 4)

# plot select variables 
phys.chem.vars = c("pH","P (ppm)","TN (%)",
                   "NO3-N (ppm)","NH4-N (ppm)","C:N Ratio",
                   "Ca (ppm)","K (ppm)","Mg (ppm)",
                   "Sand (%)","Silt (%)","Clay (%)",
                   "Temperature (C)","Gravitational moisture (%)",
                   "Bulk density (g/cm3)",">= 4.75 mm","2-4.75 mm",
                   "0.250-2 mm","0.053-0.250 mm","<0.053 mm")
phys.chem.labs = c("pH","P (ppm)","TN (%)",
                   "NO3-N (ppm)","NH4-N (ppm)","C:N ratio",
                   "Ca (ppm)","K (ppm)","Mg (ppm)",
                   "Sand (%)","Silt (%)","Clay (%)",
                   "Temperature (C)","Gravitational moisture (%)",
                   "Bulk density (g/cm3)","\u2265 4.75 mm (%)","2-4.75 mm (%)",
                   "0.250-2 mm (%)","0.053-0.250 mm (%)","< 0.053 mm (%)")
df.plot = soil.hdi.df.best[which(soil.hdi.df.best$variable.label %in% phys.chem.vars),]
df.plot$label.new = 0
for (i in 1:length(phys.chem.vars)) { df.plot$label.new[which(df.plot$variable.label == phys.chem.vars[i])] = phys.chem.labs[i] }
p.chem.phys = ggplot(df.plot, 
                     aes(y=factor(full.treatment.name, levels=trt.names), 
                         x=posterior.mean)) + 
                     geom_errorbar(aes(xmin=`X5`,xmax=`X95`), 
                                   width=0.3, color="black", 
                                   position=position_dodge(width=0.5),
                                   linewidth=0.35) +
                     geom_errorbar(aes(xmin=`X25`,xmax=`X75`), 
                                   width=0, color="black", 
                                   position=position_dodge(width=0.5),
                                   linewidth=0.8) +
                     geom_point(size=1.2, position=position_dodge(width=0.5)) +
                     facet_wrap(.~factor(label.new, levels=phys.chem.labs), 
                                ncol=4, scales="free_x") +
                     theme(panel.spacing.x = unit(0.6, "cm"),
                           text = element_text(size=12)) +
                     labs(y="",x="Posterior estimate") + 
                     scale_x_continuous(breaks=breaks_extended(n=4))
p.chem.phys
ggsave("Main_Figures/Figure4_All_Except_CEC_MWD.jpeg", 
       plot=p.chem.phys,width=24,height=24,units="cm",dpi=600)

################################################################################
# make combined plot for C stocks in biomass, debris, and soil + vegetation,
# along with species richness (Figure 6)

# define color palletes
v.palette <- brewer.pal(11,"RdYlGn")[c(7,9,11)]
d.palette <- brewer.pal(9,"YlOrBr")[seq(3,8)]
r.palette <- tail(brewer.pal(8,"BuPu"),3)

# update ipcc.df
ipcc.df = ipcc.df[1:2,]

#  live vegetation plot
stack.vars.l = c("Herbaceous biomass","Belowground woody biomass","Aboveground woody biomass")
ipcc.vars.l = c("Restored temperate forest","Natural temperate forest")
p.l = ggplot(stock.hdi.df.best[stock.hdi.df.best$variable.label %in% stack.vars.l,], 
             aes(y=factor(full.treatment.name, levels=trt.names), 
                 x=posterior.mean, 
                 fill=factor(variable.label, levels=stack.vars.l))) +
              geom_bar(stat="identity",position="stack") + 
              labs(fill="",y="",
                   x="Posterior mean stock (Mg C/ha)",title="") + 
              geom_vline(data=ipcc.df,
                         aes(xintercept=woody.value, 
                             color=factor(woody.type,levels=ipcc.vars.l),
                             linetype=factor(woody.type,levels=ipcc.vars.l)), 
                         linewidth=1.5) +
              scale_fill_manual(values=c("olivedrab3","olivedrab4",v.palette[3])) +
              scale_color_manual(values=c("red","royalblue1")) +
              scale_linetype_manual(values=c("solid","dashed")) + 
              guides(fill=guide_legend(order=1),
                     color=guide_legend(title="IPCC total woody biomass",order=2),
                     linetype=guide_legend(title="IPCC total woody biomass",order=2)) +
              theme(text=element_text(size=14), 
                    legend.key.size=unit(0.7,'cm'),
                    legend.background=element_rect(fill="transparent", color=NA),
                    plot.margin=unit(c(0,0,1,0),"lines"),
                    plot.background=element_rect(color="black",linewidth=1),
                    legend.key=element_rect(fill="darkgrey")) +
              coord_cartesian(xlim=c(0,170),clip="off") +
              geom_label(x=319,y=6.6,label="a",
                         color="black",fill=alpha("white",0.9),
                         label.r=unit(0,"pt"),label.size=0,
                         size=10,fontface="bold")
p.l
abg.df = stock.hdi.df.best[which(stock.hdi.df.best$variable.label == "Aboveground woody biomass"),]
bg.df = stock.hdi.df.best[which(stock.hdi.df.best$variable.label == "Belowground woody biomass"),]
round(bg.df$posterior.mean/(bg.df$posterior.mean+abg.df$posterior.mean)*100,1)

# debris plot
stack.vars.d.old = c("Herbaceous litter","Fine woody debris\n(< 7.6 cm)",
                     "Coarse woody debris\n(>= 7.6 cm)","Standing dead trees")
stack.vars.d.new = c("Herbaceous litter","Fine woody debris (< 7.6 cm)",
                     "Coarse woody debris (\u2265 7.6 cm)","Standing dead trees")
ipcc.vars.d = c("Restored temperate forest","Natural temperate forest")
stock.df.debris = stock.hdi.df.best[stock.hdi.df.best$variable.label %in% stack.vars.d.old,]
stock.df.debris$label.new = rep(0, nrow(stock.df.debris))
for (i in 1:n.t) { stock.df.debris$label.new[stock.df.debris$variable.label == stack.vars.d.old[i]] = stack.vars.d.new [i] }
p.d = ggplot(stock.df.debris, 
             aes(y=factor(full.treatment.name, levels=trt.names), 
                 x=posterior.mean,  
                 fill=factor(label.new, levels=stack.vars.d.new))) +
             geom_bar(stat="identity",position="stack") + 
             labs(fill="",y="",x="Posterior mean stock (Mg C/ha)",title="") + 
             geom_vline(data=ipcc.df, 
                        aes(xintercept=debris.value,
                            color=factor(debris.type,levels=ipcc.vars.d),
                            linetype=factor(debris.type,levels=ipcc.vars.d)), 
                        linewidth=1.5) +
             scale_color_manual(values=c("red","royalblue1")) +
             scale_linetype_manual(values=c("solid","dashed")) + 
             scale_fill_manual(values=d.palette) + 
             guides(fill=guide_legend(order=1),
                    color=guide_legend(title="IPCC litter\nand woody debris",order=2),
                    linetype=guide_legend(title="IPCC litter\nand woody debris",order=2)) +
             theme(text=element_text(size=14), 
                   legend.key.size=unit(0.7,'cm'),
                   legend.background=element_rect(fill="transparent", color=NA),
                   plot.margin=unit(c(0,0,1,0),"lines"),
                   plot.background=element_rect(color="black",linewidth=1),
                   legend.key=element_rect(fill="darkgrey")) +
             coord_cartesian(xlim = c(0,63), clip="off") +
             geom_label(x=121.6,y=6.6,label="b",
                        color="black",fill=alpha("white",0.9),
                        label.r=unit(0,"pt"),label.size=0,
                        size=10,fontface="bold")
p.d

# total ecosystem
stack.vars.e = c("Soil inorganic carbon","Soil organic carbon","Total dead vegetation","Total live vegetation")
ipcc.vars.e = c("Restored forested wetland","Natural forested wetland")
stock.soil.df = rbind(soil.hdi.df.best[soil.hdi.df.best$variable %in% c("tic.stock","toc.stock"),],
                      stock.hdi.df.best[stock.hdi.df.best$variable %in% c("total.dead.carbon","total.live.carbon"),])
stock.soil.df$variable.label[stock.soil.df$variable.label == "TIC"] = "Soil inorganic carbon"
stock.soil.df$variable.label[stock.soil.df$variable.label == "SOC (Mg/ha)"] = "Soil organic carbon"
p.e = ggplot(stock.soil.df, 
             aes(y=factor(full.treatment.name, levels=trt.names), 
                 x=posterior.mean,   
                 fill=factor(variable.label, levels=stack.vars.e))) +
             geom_bar(stat="identity",position="stack") + 
             labs(fill="",y="",x="Posterior mean stock (Mg C/ha)",title="") +
             scale_fill_manual(values=c(brewer.pal(9,"Greys")[5],s.palette[3],
                                        d.palette[4],v.palette[3]),
                               labels=stack.vars.e) +
              geom_vline(data=ipcc.df, 
                         aes(xintercept=total.value, 
                             color=factor(total.type,levels=ipcc.vars.e),
                             linetype=factor(total.type,levels=ipcc.vars.e)), 
                         linewidth=1.5) +
             scale_color_manual(values=c("red","royalblue1")) +
             scale_linetype_manual(values=c("solid","dashed")) +
             guides(fill=guide_legend(order=1),
                    color=guide_legend(title="IPCC total organic C",order=2),
                    linetype=guide_legend(title="IPCC total organic C",order=2)) +
             theme(text=element_text(size=14), 
                   legend.key.size=unit(0.7,'cm'),
                   legend.background=element_rect(fill="transparent", color=NA),
                   plot.margin=unit(c(0,0,1,0),"lines"),
                   plot.background=element_rect(color="black",linewidth=1),
                   legend.key=element_rect(fill="darkgrey")) +
             coord_cartesian(xlim = c(0,300), clip="off") +
             geom_label(x=563,y=6.6,label="c",
                        color="black",fill=alpha("white",0.9),
                        label.r=unit(0,"pt"),label.size=0,
                        size=10,fontface="bold")
p.e

# richness plot
stack.vars.r.old = c("Herbaceous layer richness","Tree layer richness","Both tree & herbaceous\nlayer richness") 
stack.vars.r.new = c("Herbaceous layer only","Tree layer only","Tree and herbaceous layer")
richness.df = stock.hdi.df.best[which(stock.hdi.df.best$variable.label %in% stack.vars.r.old),]
richness.df$label.new = rep(0, nrow(richness.df))
for (i in 1:n.t) { richness.df$label.new[which(richness.df$variable.label == stack.vars.r.old[i])] = stack.vars.r.new[i] }
p.r = ggplot(richness.df,
             aes(y=factor(full.treatment.name, levels=trt.names), 
                 x=posterior.mean,
                 fill=factor(label.new, levels=rev(stack.vars.r.new)))) + 
             geom_bar(stat="identity",position="stack") +
             labs(fill="",y="",x="Posterior mean richness",title="") + 
             scale_fill_manual(values=r.palette) +
             theme(text = element_text(size=14),
                   legend.key.size = unit(0.7,'cm'),
                   legend.background = element_rect(fill="transparent", color=NA),
                   plot.margin=unit(c(0,0,1,0),"lines"),
                   plot.background=element_rect(color="black",linewidth=1),
                   legend.key=element_rect(fill="darkgrey")) +
             coord_cartesian(xlim=c(0,20),clip="off") +
             geom_label(x=38.6,y=6.6,label="d",
                        color="black",fill=alpha("white",0.9),
                        label.r=unit(0,"pt"),label.size=0,
                        size=10,fontface="bold")
p.r

# combine all plots
p.c.stocks = (p.l+theme(plot.margin=unit(c(0,0,0,0),"pt"))+p.d)/(p.e+theme(plot.margin=unit(c(0,16,0,0),"pt"))+p.r)
p.c.stocks
ggsave("Main_Figures/Figure6_Soil_Vegetation_Ecosystem_Cstocks_Richness.jpeg", 
       plot=p.c.stocks, width=40, height=23, units="cm", dpi=1000)

################################################################################
# summarize Frax penn results for Table B3

live.frax.stock.df = stock.hdi.df.best[which(stock.hdi.df.best$variable == "snag.frax.live.carbon"),]
dead.frax.stock.df = stock.hdi.df.best[which(stock.hdi.df.best$variable == "snag.frax.dead.carbon"),]

cbind(live.frax.stock.df[,"full.treatment.name"],
      signif(live.frax.stock.df[,c("posterior.mean","X5","X95")],digits=3))

cbind(dead.frax.stock.df[,"full.treatment.name"],
      signif(dead.frax.stock.df[,c("posterior.mean","X5","X95")],digits=3))

################################################################################
# print C stock results for comparing with IPCC values

soc.stock.df = stock.hdi.df.best[which(stock.hdi.df.best$variable == "SOC"),]
cbind(soc.stock.df[,"full.treatment.name"],
      signif(soc.stock.df[,c("posterior.mean","X5","X95")]-88,digits=3))


################################################################################
# plot C stocks by species for snags, coarse woody debris, live trees,
# and herbaceous layer species

# ecosystem carbon pools
pools = c("total.cwd.carbon","large.snag.carbon","live.tree.carbon")
df.paths = c("Clean_Data_By_Species/CWD_Carbon_Stocks_By_Species.csv",
             "Clean_Data_By_Species/Snag_Carbon_Stocks_By_Species.csv",
             "Clean_Data_By_Species/Woody_Biomass_Carbon_Stocks_By_Species.csv")
n.pool = length(pools)

## read in data frames
df.list = list()
for (i in 1:n.pool) { df.list[[pools[i]]] = read.csv(paste("Tree_Analysis/",df.paths[i],sep="")) }

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


# calculate means
mean.df.all = data.frame(matrix(nrow=0,ncol=4))
colnames(mean.df.all) = c("treatment","species","mean","pool")
pool.order = rev(pools)
pool.labels = c("Live trees (\u2265 2.5 cm)",
                "Standing dead trees (\u2265 2.5 cm)",
                "Coarse woody debris (\u2265 7.6 cm)")
for (i in 1:n.pool) {
  pool.i = pool.order[i]
  colnames(df.list[[pool.i]])[5] = "value"
  df.sum.i = df.list[[pool.i]] %>%
    group_by(treatment, species) %>%
    summarize(mean = mean(value))
  df.sum.i$pool = pool.labels[i]
  mean.df.all = rbind(mean.df.all, df.sum.i)
}
mean.df.all$full.treatment.name = 0
for (i in 1:n.t) { mean.df.all$full.treatment.name[mean.df.all$treatment == trt.letters[i]] = trt.names[i] }

# plot means for live trees
uni.spp = unique(mean.df.all$species)
p.species = ggplot(mean.df.all, 
            aes(y=factor(full.treatment.name,levels=trt.names),
                x=mean,
                fill=forcats::fct_rev(species))) + 
            geom_bar(stat="identity",
                     position="stack") +
            labs(x="Empirical plot-level mean C stock (Mg/ha)", 
                 y="", fill="Species") +
            facet_wrap(.~factor(pool, levels=pool.labels), 
                       ncol=1, scales="free_x") + 
            guides(fill = guide_legend(reverse = TRUE))
ggsave("Supp_Figures/FigureB10_Woody_C_Stocks_By_Species.jpeg", 
       plot=p.species, width=14, height=16, units="cm")

################################################################################
# plot stacked means for herbaceous species groups

sp.groups = c("mixed.biomass.c.stock","h.japonicus.c.stock","p.arundinacea.c.stock")
sp.labels = c("Mixed biomass","H. japonicus biomass","P. arundinacea biomass")
herbaceous.biomass.c.stock.hdi.df = stock.hdi.df.best[stock.hdi.df.best$variable %in% sp.groups,]
p.herbC = ggplot(herbaceous.biomass.c.stock.hdi.df, 
                 aes(x=posterior.mean,
                     y=factor(full.treatment.name, levels=trt.names), 
                     fill=factor(variable.label, levels=sp.labels))) + 
          geom_bar(stat="identity",position="stack") +
          labs(y="", x="Posterior mean C stock (Mg/ha)", 
               title = "",
               fill="Species group") + 
          theme(text = element_text(size=12))
p.herbC
ggsave("Supp_Figures/FigureB11_Herbaceous_C_Stocks_By_Species.jpeg", 
       plot=p.herbC, width=18, height=10, units="cm")