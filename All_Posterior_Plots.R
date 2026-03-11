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
#soil.var.df = read.csv("Metadata/Quadrat_Level_Soil_Variables.csv")
soil.rhat.df = read.csv("Soil_Analysis/Posteriors/Soil_Rhat_Statistic_Updated.csv")
soil.vars = unique(soil.rhat.df$variable)
soil.var.labs = unique(soil.rhat.df$variable.label)
n.s.v = length(soil.vars)

# stock and richness variables
#stock.var.df = read.csv("Metadata/Plot_Level_Carbon_Richness_Variables.csv")
#leave.out = c("small.fwd.carbon","int.fwd.carbon","dead.stem.carbon.min","snag.carbon.min",
#              "abg.live.tree.carbon","bg.live.tree.carbon","total.live.tree.carbon",
#              "abg.live.stem.carbon","bg.live.stem.carbon","total.live.stem.carbon")
#omit.var.id = which(stock.var.df$variable %in% leave.out)
#stock.vars = stock.var.df[-omit.var.id,"variable"]
#stock.var.labs = stock.var.df[-omit.var.id,"label"]
#n.c.v = length(stock.vars)

################################################################################
# evaluate model convergence

# read in rhat statistic dataframes
#stock.rhat.df = read.csv("Tree_Analysis/Posteriors/Stock_Rhat_Statistic.csv")

# soil carbon variables: plot rhat's for each model and variable
soil.carbon.vars = c("tc.percent","tc.stock","toc.percent","toc.stock","tic.percent",
                     "tic.stock","poc.percent","poc.stock","maoc.percent","maoc.stock")
soil.carbon.labs = rep(0,10)
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
                          labs(y="",x="(R-hat statistic - 1)*100",
                               color="Model",shape="Model") +
                          geom_vline(xintercept=1) +
                          scale_x_continuous(breaks=seq(0,1,0.2),
                                             limits=c(-0.01,1)) + 
                          facet_wrap(.~factor(variable.label, levels=soil.carbon.labs), ncol=2)
p.soil.carbon.rhat.comp
ggsave("Supp_Figures/FigureB1_Soil_Carbon_Rhat_Score_Comparison.jpeg", 
       plot=p.soil.carbon.rhat.comp, width=18, height=16, units="cm",dpi=600)


# non-carbon soil variables: plot rhat's for each model and variable
soil.noncarbon.vars = soil.vars[-which(soil.vars %in% soil.carbon.vars)]
soil.noncarbon.labs = rep(0,length(soil.noncarbon.vars))
for (i in 1:length(soil.noncarbon.vars)) { 
  soil.noncarbon.labs[i] = soil.rhat.df[soil.rhat.df$variable == soil.noncarbon.vars[i],"variable.label"][1] 
}
soil.noncarbon.labs[which(soil.noncarbon.labs == "Root and wood fragment density (g/cm3)")] = "Root & wood\nfragment density (g/cm3)"
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
                                    labs(y="",x="(R-hat statistic - 1)*100",
                                         color="Model",shape="Model") +
                                    geom_vline(xintercept=1) +
                                    scale_x_continuous(breaks=seq(0,1,0.2),
                                                       limits=c(-0.01,1)) + 
                                    facet_wrap(.~factor(variable.label, levels=soil.noncarbon.labs), ncol=4)
p.soil.noncarbon.rhat.comp
ggsave("Supp_Figures/FigureB2_Non_Carbon_Soil_Rhat_Score_Comparison.jpeg", 
       plot=p.soil.noncarbon.rhat.comp, width=26, height=22, units="cm",dpi=600)

# carbon stocks: plot rhat's for each model and variable
p.stock.rhat.comp = ggplot(stock.rhat.df, 
                           aes(x=(Rhat-1)*100, 
                               y=factor(full.treatment.name, levels=trt.names),
                               color=factor(model.label, levels=model.labels),
                                            shape=factor(model.label, levels=model.labels))) + 
                           geom_vline(xintercept=0) +     
                           geom_point() + 
                           scale_shape_discrete(solid = F) +
                           labs(y="",x="(R-hat statistic - 1)*100",
                                color="Model",shape="Model") +
                           geom_vline(xintercept=1) +
                           scale_x_continuous(breaks=seq(0,1,0.2),limits=c(-0.01,1)) + 
                           facet_wrap(.~factor(variable.label, levels=stock.var.labs), ncol=4)
p.stock.rhat.comp

################################################################################
# count the number of variables for which each model type minimizes the loo v. waic (Tables B1 & B2)

# read in soil data model comparison dataframe
comp.df.soil = read.csv("Soil_Analysis/Posteriors/Soil_Model_Information_Criteria_Updated.csv")
#comp.df.stocks = read.csv("Tree_Analysis/Posteriors/Stock_Model_Information_Criteria.csv")

# soil
soil.min.ic.count.df = data.frame(matrix(nrow=n.c, ncol=n.m))
colnames(soil.min.ic.count.df) = model.labels
row.names(soil.min.ic.count.df) = criteria
for (i in 1:n.c) {
  criterion.i = criteria[i]
  soil.min.ic.count.df[i,] = rep(0,3)
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
  stock.min.ic.count.df[i,] = rep(0,2)
  for (j in 1:n.c.v) {
    var.j = stock.vars[j]
    row.ij.ind = which(comp.df.stocks$variable == var.j & comp.df.stocks$criterion == criterion.i)
    row.ij = comp.df.stocks[row.ij.ind,]
    min.ic.model = row.ij[which.min(row.ij$ic), "model.label"]
    stock.min.ic.count.df[i,min.ic.model] = stock.min.ic.count.df[i,min.ic.model] + 1
  }
}

################################################################################
# compare model information criteria (Figures B1 & B2)

# soil carbon properties: plot comparison of information criteria
soil.carbon.vars = c("tc.percent","toc.percent","tic.percent","poc.percent","maoc.percent",
                     "tc.stock","toc.stock","tic.stock","poc.stock","maoc.stock")
soil.carbon.labs = rep(0,10)
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
                              facet_wrap(.~factor(variable.label, levels=soil.var.labs), 
                                         scales="free_x", ncol=2) + 
                              theme(text = element_text(size=14)) + 
                              labs(y="", x="Score",color="Information criterion",
                                   linetype="Information criterion",shape="Best model")
p.carbon.soil.ic.score.comp
ggsave("Supp_Figures/FigureB4_Soil_Carbon_IC_Score_Comparison.jpeg", 
       plot=p.carbon.soil.ic.score.comp, width=20, height=18, units="cm",dpi=600)

# soil non-carbon properties: plot comparison of information criteria
soil.noncarbon.vars = soil.vars[-which(soil.vars %in% soil.carbon.vars)]
soil.noncarbon.labs = rep(0,length(soil.noncarbon.vars))
for (i in 1:length(soil.noncarbon.vars)) { soil.noncarbon.labs[i] = comp.df.soil[comp.df.soil$variable == soil.noncarbon.vars[i],"variable.label"][1] }
soil.noncarbon.labs[which(soil.noncarbon.labs == "Root and wood fragment density (g/cm3)")] = "Root & wood\nfragment density (g/cm3)"
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
       plot=p.noncarbon.soil.ic.score.comp, width=30, height=28, units="cm",dpi=600)

# soil: plot comparison of information criteria
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
ggsave("Figures/FigureB2_Stock_IC_Score_Comparison.jpeg", 
       plot=p.stock.ic.score.comp, width=35, height=30, units="cm",dpi=600)

###############################################################################
# plot HDIs for each model and variable (Figures B3 & B4)

# read posterior distributions for soil properties, carbon stocks, and species richness
soil.hdi.df = read.csv("Soil_Analysis/Posteriors/Soil_Posterior_Intervals_10Chains_NaturalScale_Updated.csv")
#stock.hdi.df = read.csv("Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_10Chains_NaturalScale.csv")

# soil carbon: plot posterior HDIs for all soil variables and models
soil.carbon.vars = c("tc.percent","tc.stock","toc.percent","toc.stock","tic.percent",
                     "tic.stock","poc.percent","poc.stock","maoc.percent","maoc.stock")
soil.carbon.labs = rep(0,10)
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
                                    ncol=2, scales="free_x", dir = "h") + 
                         theme(text = element_text(size=12)) + 
                         labs(y="", x="Posterior estimate", color="Model", shape="Model")
p.soil.carbon.hdi.comp
ggsave("Supp_Figures/FigureB7_Soil_Carbon_HDI_Comparison.jpeg", 
       plot=p.soil.carbon.hdi.comp, width=20, height=20, units="cm",dpi=600)

# soil non-carbon: plot posterior HDIs for all soil variables and models
soil.noncarbon.vars = soil.vars[-which(soil.vars %in% soil.carbon.vars)]
soil.noncarbon.labs = rep(0,length(soil.noncarbon.vars))
for (i in 1:length(soil.noncarbon.vars)) { soil.noncarbon.labs[i] = soil.hdi.df[soil.hdi.df$variable == soil.noncarbon.vars[i],"variable.label"][1] }
soil.noncarbon.labs[which(soil.noncarbon.labs == "Root and wood fragment density (g/cm3)")] = "Root & wood\nfragment density (g/cm3)"
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
                                   theme(text = element_text(size=12)) + 
                                   labs(y="", x="Posterior estimate", color="Model", shape="Model")
p.soil.noncarbon.hdi.comp
ggsave("Supp_Figures/FigureB8_NonCarbon_Soil_HDI_Comparison.jpeg", 
       plot=p.soil.noncarbon.hdi.comp, width=30, height=30, units="cm",dpi=600)

# soil: plot posterior HDIs for all soil variables and models
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
ggsave("Figures/FigureB4_Stock_HDI_Comparison.jpeg", 
       plot=p.stock.hdi.comp, width=32, height=32, units="cm",dpi=600)

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
soil.hdi.df.best$lab = rep(0, nrow(soil.hdi.df.best))
for (i in 1:3) { soil.hdi.df.best$lab[which(soil.hdi.df.best$variable == c.conc.labs[i])] = c.conc.labs.new[i] }
s.palette <- c(brewer.pal(9,"Greys")[5],brewer.pal(11,"BrBG")[c(2,1)])
p.c.concentrations = ggplot(soil.hdi.df.best[which(soil.hdi.df.best$lab %in% c.conc.labs.new),], 
                            aes(y=factor(full.treatment.name, levels=trt.names), 
                                x=posterior.mean, 
                                fill=factor(lab, levels=c.conc.labs.new))) + 
                            geom_bar(stat="identity",
                                     position="stack") +
                            labs(x="Concentration (% [g C/g soil])",
                                 y="",fill="Soil carbon fraction") + 
                            scale_fill_manual(values=s.palette) +
                            theme(text = element_text(size=14),
                                  plot.margin=unit(c(1,1,1,1),"lines")) +
                            coord_cartesian(xlim = c(0,5.5), clip="off") +
                            scale_x_continuous(breaks=seq(0,5)) +
                            geom_label(x=5.25,y=1,label="a",
                                       color="black",fill=alpha("white",0.9),
                                       label.r=unit(0,"pt"),label.size=0,
                                       size=8,fontface="bold") + 
                            guides(fill="none")
p.c.concentrations
ipcc.soil.vars = c("Annual crops","Revegetated cropland","Natural wetland")
p.c.stocks = ggplot(stock.hdi.df.best[which(stock.hdi.df.best$variable %in% c.conc.labs.new),], 
                    aes(y=factor(full.treatment.name, levels=trt.names), 
                        x=posterior.mean, 
                        fill=factor(variable, levels=c.conc.labs.new))) +
                    geom_bar(stat="identity",position="stack") + 
                    scale_fill_manual(values=s.palette) +
                    labs(fill="",y="",x="Stock (Mg C/ha)",title="") +
                    geom_vline(data=ipcc.df, 
                               aes(xintercept=soc.value, 
                                   color=factor(soc.type,levels=ipcc.soil.vars),
                                   linetype=factor(soc.type,levels=ipcc.soil.vars)), 
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
                    scale_x_continuous(breaks=seq(0,125,25))+
                    geom_label(x=136,y=1,label="b",
                               color="black",fill=alpha("white",0.9),
                               label.r=unit(0,"pt"),label.size=0,
                               size=8,fontface="bold")
p.c.stocks
p.c.all = p.c.concentrations + p.c.stocks
p.c.all =  p.c.all + theme(plot.margin = margin(0, 0, 0, 0, "cm"))
p.c.all
ggsave("Main_Figures/Figure3_Soil_Carbon_Concentrations_and_Stocks.jpeg", 
       plot=p.c.all, width=30, height=12, units="cm",dpi=600)

################################################################################
# plot non-carbon chemical and physical variables (Figure 4)

# plot select variables 
phys.chem.vars = c("Temperature (C)",
                   "Gravitational moisture (%)",
                   "Bulk density (g/cm3)",
                   "Root & wood\nfragment density (g/cm3)",
                   "POC:MAOC ratio",
                   "TN (%)",
                   "C:N Ratio",
                   "Sand (%)","Silt (%)","Clay (%)",
                   "pH","P (ppm)","NO3-N (ppm)","NH4-N (ppm)",
                   "K (ppm)","Ca (ppm)","Mg (ppm)","CEC (meq/100 g)",
                   "Mean-weight diameter (mm)",">= 4.75 mm","2-4.75 mm",
                   "0.250-2 mm","0.053-0.250 mm","<0.053 mm")
phys.chem.labs = c("Temperature (C)",
                   "Gravitational moisture (%)",
                   "Bulk density (g/cm3)",
                   "Root & wood\nfragment density (g/cm3)",
                   "POC:MAOC ratio",
                   "TN (%)",
                   "C:N ratio",
                   "Sand (%)","Silt (%)","Clay (%)",
                   "pH","P (ppm)","NO3-N (ppm)","NH4-N (ppm)",
                   "K (ppm)","Ca (ppm)","Mg (ppm)","CEC (meq/100 g)",
                   "Mean-weight diameter (mm)","\u2265 4.75 mm (%)","2-4.75 mm (%)",
                   "0.250-2 mm (%)","0.053-0.250 mm (%)","< 0.053 mm (%)")
df.plot = soil.hdi.df.best[which(soil.hdi.df.best$variable.label %in% phys.chem.vars),]
df.plot$label.new = rep(0, nrow(df.plot))
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
       plot=p.chem.phys,width=24,height=26,units="cm",dpi=600)

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
stack.vars.l = c("Herbaceous biomass",
                 "Belowground woody biomass",
                 "Aboveground woody biomass")
ipcc.vars.l = c("Restored temperate forest","Natural temperate forest")
p.l = ggplot(stock.hdi.df.best[which(stock.hdi.df.best$variable.label %in% stack.vars.l),], 
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
              coord_cartesian(xlim=c(0,180),clip="off") +
              geom_label(x=338,y=6.6,label="a",
                         color="black",fill=alpha("white",0.9),
                         label.r=unit(0,"pt"),label.size=0,
                         size=10,fontface="bold")
p.l
#abg.df = stock.hdi.df.best[which(stock.hdi.df.best$variable.label == "Aboveground woody biomass"),]
#bg.df = stock.hdi.df.best[which(stock.hdi.df.best$variable.label == "Belowground woody biomass"),]
#round(bg.df$posterior.mean/(bg.df$posterior.mean+abg.df$posterior.mean)*100,1)

# debris plot
stack.vars.d.old = c("Herbaceous litter","Fine woody debris (< 7.6 cm)",
                     "Coarse woody debris (>= 7.6 cm)","Standing dead trees")
stack.vars.d.new = c("Herbaceous litter","Fine woody debris (< 7.6 cm)",
                     "Coarse woody debris (\u2265 7.6 cm)","Standing dead trees")
ipcc.vars.d = c("Restored temperate forest","Natural temperate forest")
stock.df.debris = stock.hdi.df.best[which(stock.hdi.df.best$variable.label %in% stack.vars.d.old),]
stock.df.debris$label.new = rep(0, nrow(stock.df.debris))
for (i in 1:n.t) { stock.df.debris$label.new[which(stock.df.debris$variable.label == stack.vars.d.old[i])] = stack.vars.d.new [i] }
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
             coord_cartesian(xlim = c(0,30.5), clip="off") +
             geom_label(x=58.9,y=6.6,label="b",
                        color="black",fill=alpha("white",0.9),
                        label.r=unit(0,"pt"),label.size=0,
                        size=10,fontface="bold")
p.d

# total ecosystem
stack.vars.e = c("Soil inorganic carbon","Soil organic carbon","Total dead vegetation","Total live vegetation")
ipcc.vars.e = c("Restored forested wetland","Natural forested wetland")
p.e = ggplot(stock.hdi.df.best[which(stock.hdi.df.best$variable.label %in% stack.vars.e),], 
             aes(y=factor(full.treatment.name, levels=trt.names), 
                 x=posterior.mean,   
                 fill=factor(variable.label, levels=stack.vars.e))) +
             geom_bar(stat="identity",position="stack") + 
             labs(fill="",y="",x="Posterior mean stock (Mg C/ha)",title="") +
             geom_vline(data=ipcc.df, 
                        aes(xintercept=total.value, 
                            color=factor(total.type,levels=ipcc.vars.e),
                            linetype=factor(total.type,levels=ipcc.vars.e)), 
                         linewidth=1.5) +
             scale_fill_manual(values=c(brewer.pal(9,"Greys")[5],
                                        s.palette[3],d.palette[4],
                                        v.palette[3]),
                               labels=c("TIC","SOC","Litter and woody debris","Living biomass")) +
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
             geom_label(x=562,y=6.6,label="c",
                        color="black",fill=alpha("white",0.9),
                        label.r=unit(0,"pt"),label.size=0,
                        size=10,fontface="bold")
p.e

# richness plot
stack.vars.r.old = c("Herbaceous layer richness","Tree layer richness","Both tree and herbaceous layer") 
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
ggsave("Figures/Figure6_Veg_Ecosystem_Cstocks_Richness.jpeg", 
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

# read in posterior intervals
post.df.spp = read.csv("Tree_Analysis/Posteriors/Vegetation_Carbon_Stocks_Species_Means_Intervals_5Chains.csv")

# plot stacked means for live woody, snag/CWD area, and biomass stocks by species
type.rows = c("live.trees","dead.trees","coarse.woody.debris")
type.labs = c("Live trees (\u2265 2.5 cm)",
              "Standing dead trees (\u2265 2.5 cm)",
              "Coarse woody debris (\u2265 7.6 cm)")
post.df.woodC = post.df.spp[which(post.df.spp$type %in% type.rows),]
post.df.woodC$type.lab = rep(0, nrow(post.df.woodC))
for (i in 1:3) { post.df.woodC$type.lab[which(post.df.woodC$type == type.rows[i])] = type.labs[i] }
post.df.woodC$species = trimws(post.df.woodC$species)
all.sp = sort(unique(post.df.woodC$species))
p.woodyC = ggplot(post.df.woodC, 
                  aes(x=posterior.mean, 
                      y=factor(treatment, levels=trt.names),
                      fill=factor(species, levels=all.sp))) + 
                  geom_bar(stat="identity") +
                  facet_wrap(.~factor(type.lab, levels=type.labs), 
                             scales="free_x", ncol=1) +
                  labs(x="Posterior mean C stock (Mg/ha)", y="", fill="Species") +
                  guides(fill=guide_legend(ncol=1)) + 
                  theme(text = element_text(size=14))
p.woodyC
ggsave("Figures/FigureB1_Woody_C_Stocks_By_Species.jpeg", 
       plot=p.woodyC, width=17, height=19, units="cm")

# plot stacked means for herbaceous species grounds
post.df.herbC = post.df.spp[which(post.df.spp$type == "herbaceous.biomass"),]
spp.herbC = levels(factor(post.df.herbC$species))
p.herbC = ggplot(post.df.herbC, 
                 aes(x=posterior.mean,
                     y=factor(treatment, levels=trt.names), 
                     fill=factor(species, levels=spp.herbC))) + 
                 geom_bar(stat = "identity") + 
                 labs(y="", x="Posterior mean C stock (Mg/ha)", 
                      title = "",
                      fill="Species group") + 
                 theme(text = element_text(size=14))
p.herbC
ggsave("Figures/FigureB2_Herbaceous_C_Stocks_By_Species.jpeg", 
       plot=p.herbC, width=18, height=10, units="cm")

