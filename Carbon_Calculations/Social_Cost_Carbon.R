setwd("C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment")

library(reshape2)
library(patchwork)
library(dplyr)
library(tidyverse)

### script that uses social cost of carbon to estimate absolute and relative 
### carbon benefits of each treatment

################################################################################
# load data

# constants and conversions
ch4.100y = 27
n2o.100y = 273
c.mol.mass = 12.011
n.mol.mass = 14.007
o.mol.mass = 15.999
h.mol.mass = 1.008
y.sr = 25
baseline.cstock = data.frame(mean=55.9, lower=39.1, upper=75.6)
ch4.mol.mass = c.mol.mass + 4*h.mol.mass # 16.043
n2o.mol.mass = 2*n.mol.mass + o.mol.mass # 44.013
co2.mol.mass = c.mol.mass + 2*o.mol.mass # 44.009

# treatments
trt.df = read.csv("floodplain-experiment-repo/Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# statistics
stats = c("mean","lower","upper")
n.s = length(stats)

# social cost of carbon estimates 
scc.df = read.csv("floodplain-experiment-repo/Carbon_Calculations/SCC_Estimates.csv")
colnames(scc.df) = c("stat","per.CO2","per.CO2.C")
scc.df$stat = stats
rownames(scc.df) = stats

# establishment costs
est.df = read.csv("floodplain-experiment-repo/Carbon_Calculations/Treatment_Establishment_Costs.csv")
colnames(est.df) = c("trt","cost.2019","cost.2023")
est.df$cost.2023 = as.numeric(est.df$cost.2023)

# ecosystem carbon estimates
stock.df = read.csv("floodplain-experiment-repo/Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_10Chains_NaturalScale.csv")
ecoC.df.total = subset(stock.df, model == "strip.random" & variable == "total.organic.carbon")[,c("full.treatment.name","posterior.mean","X5","X95")]
colnames(ecoC.df.total) = c("trt",stats)
cbind(ecoC.df.total[1:6,"trt"],signif(ecoC.df.total[1:6,stats],3))

for (i in 1:n.s) { ecoC.df.total[,stats[i]] = ecoC.df.total[,stats[i]] - baseline.cstock[,stats[i]] }
cbind(ecoC.df.total[1:6,"trt"],signif(ecoC.df.total[1:6,stats],3))

# annual accumulation rates
ecoC.df.annual = ecoC.df.total
ecoC.df.annual[,stats] = ecoC.df.annual[,stats]/y.sr
cbind(ecoC.df.annual[1:6,"trt"],signif(ecoC.df.annual[1:6,stats],3))

# read in meta-analysis GHG emission estimates
ghg.meta = read.csv("floodplain-experiment-repo/Carbon_Calculations/He_2024_Meta_Analysis_GHG_Estimates.csv")
ghg.meta = subset(ghg.meta, Unit == "kg/ha/y")
colnames(ghg.meta) = c("Ecosystem change","molecule","mean","se","unit")
ghg.meta$lower = ghg.meta$mean - 1.645*ghg.meta$se
ghg.meta$upper = ghg.meta$mean + 1.645*ghg.meta$se

################################################################################
## combine datasets to estimate carbon benefit

# melt ecosystem carbon estimates by statistic
ecoC.df.total.melt = melt(ecoC.df.total, 
                          id.vars=c("trt"), 
                          variable.name="stat",
                          value.name="stock")
ecoC.df.annual.melt = melt(ecoC.df.annual, 
                           id.vars=c("trt"), 
                           variable.name="stat",
                           value.name="rate")

# estimate carbon benefit of each treatment
ecoC.df.total.melt$carbon.benefit = 0
tot.crop.ch4 = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Methane")[1,stats]*y.sr/1000 # (kg/ha/y)*(25 y)*(1 Mg/1000 kg) = Mg CH4/ha
tot.crop.n2o = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Nitrous oxide")[1,stats]*y.sr/1000 # (kg/ha/y)*(25 y)*(1 Mg/1000 kg) = Mg N2O/ha
for (i in 1:n.s) {
  s.i = stats[i]
  s.id = which(ecoC.df.total.melt$stat == s.i)
  ecoC.df.total.melt[s.id,"carbon.benefit"] = (ecoC.df.total.melt[s.id,"stock"]/c.mol.mass*co2.mol.mass - ch4.100y*tot.crop.ch4[1,s.i] - n2o.100y*tot.crop.n2o[1,s.i]) * scc.df[s.i,"per.CO2"]
}
cbind(subset(ecoC.df.total.melt, stat=="mean")[1:6,"trt"], signif(subset(ecoC.df.total.melt, stat=="mean")[1:6,"carbon.benefit"]/1000,3))
cbind(subset(ecoC.df.total.melt, stat=="lower")[1:6,"trt"], signif(subset(ecoC.df.total.melt, stat=="lower")[1:6,"carbon.benefit"]/1000,3))
cbind(subset(ecoC.df.total.melt, stat=="upper")[1:6,"trt"], signif(subset(ecoC.df.total.melt, stat=="upper")[1:6,"carbon.benefit"]/1000,4))

ecoC.df.annual.melt$carbon.benefit = 0
ann.crop.ch4 = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Methane")[1,stats]/1000 # (kg/ha/y)*(1 Mg/1000 kg) = Mg CH4/ha/y
ann.crop.n2o = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Nitrous oxide")[1,stats]/1000 # (kg/ha/y)*(1 Mg/1000 kg) = Mg N2O/ha/y
for (i in 1:n.s) {
  s.i = stats[i]
  s.id = which(ecoC.df.annual.melt$stat == s.i)
  ecoC.df.annual.melt[s.id,"carbon.benefit"] = (ecoC.df.annual.melt[s.id,"rate"]/c.mol.mass*co2.mol.mass - ch4.100y*ann.crop.ch4[1,s.i] - n2o.100y*ann.crop.n2o[1,s.i]) * scc.df[s.i,"per.CO2"]
}
cbind(subset(ecoC.df.annual.melt, stat=="mean")[1:6,"trt"], signif(subset(ecoC.df.annual.melt, stat=="mean")[1:6,"carbon.benefit"],3))
cbind(subset(ecoC.df.annual.melt, stat=="lower")[1:6,"trt"],signif(subset(ecoC.df.annual.melt, stat=="lower")[1:6,"carbon.benefit"],3))
cbind(subset(ecoC.df.annual.melt, stat=="upper")[1:6,"trt"],signif(subset(ecoC.df.annual.melt, stat=="upper")[1:6,"carbon.benefit"],3))

# estimate net benefit of each treatment by subtracting the establishment cost
ecoC.df.total.melt$net.benefit = 0
ecoC.df.total.melt$breakeven.scc = 0
for (i in 1:n.t) {
  cost.i = est.df[est.df$trt == trt.names[i],"cost.2023"]
  trt.id = which(ecoC.df.total.melt$trt == trt.names[i])
  ecoC.df.total.melt[trt.id,"net.benefit"] = ecoC.df.total.melt[trt.id,"carbon.benefit"] - cost.i
  for (j in 1:n.s) {
    s.j = stats[j]
    scc.ij = cost.i/(ecoC.df.total.melt[trt.id[j],"stock"]/c.mol.mass*co2.mol.mass - ch4.100y*tot.crop.ch4[1,s.j] - n2o.100y*tot.crop.n2o[1,s.j])
    if (scc.ij >= 0) {
      ecoC.df.total.melt[trt.id[j],"breakeven.scc"] = scc.ij
    } else {
      ecoC.df.total.melt[trt.id[j],"breakeven.scc"] = Inf
    }
  }
}
cbind(subset(ecoC.df.total.melt, stat=="mean")[1:6,"trt"], signif(subset(ecoC.df.total.melt, stat=="mean")[1:6,"net.benefit"]/1000,3))
cbind(subset(ecoC.df.total.melt, stat=="lower")[1:6,"trt"],signif(subset(ecoC.df.total.melt, stat=="lower")[1:6,"net.benefit"]/1000,3))
cbind(subset(ecoC.df.total.melt, stat=="upper")[1:6,"trt"],signif(subset(ecoC.df.total.melt, stat=="upper")[1:6,"net.benefit"]/1000,3))

cbind(subset(ecoC.df.total.melt, stat=="mean")[1:6,"trt"], signif(subset(ecoC.df.total.melt, stat=="mean")[1:6,"breakeven.scc"],3))
cbind(subset(ecoC.df.total.melt, stat=="lower")[1:6,"trt"],signif(subset(ecoC.df.total.melt, stat=="lower")[1:6,"breakeven.scc"],3))
cbind(subset(ecoC.df.total.melt, stat=="upper")[1:6,"trt"],signif(subset(ecoC.df.total.melt, stat=="upper")[1:6,"breakeven.scc"],3))

# estimate the breakeven year of each treatment 
ecoC.df.annual.melt$breakeven.year = 0
for (i in 1:n.t) {
  cost.i = est.df[est.df$trt == trt.names[i],"cost.2023"]
  trt.id = which(ecoC.df.total.melt$trt == trt.names[i])
  for (j in 1:n.s) {
    stat.j = stats[j]
    years.ij = cost.i/((ecoC.df.annual.melt[trt.id[j],"rate"]/c.mol.mass*co2.mol.mass - ch4.100y*ann.crop.ch4[1,stat.j] - n2o.100y*ann.crop.n2o[1,stat.j])*scc.df[stat.j,"per.CO2"])
    ecoC.df.annual.melt[trt.id[j],"breakeven.year"] = ceiling(years.ij + 1998)
  }
}
cbind(subset(ecoC.df.annual.melt, stat=="mean")[1:6,"trt"], round(subset(ecoC.df.annual.melt, stat=="mean")[1:6,"breakeven.year"]))
cbind(subset(ecoC.df.annual.melt, stat=="lower")[1:6,"trt"], round(subset(ecoC.df.annual.melt, stat=="lower")[1:6,"breakeven.year"]))
cbind(subset(ecoC.df.annual.melt, stat=="upper")[1:6,"trt"], round(subset(ecoC.df.annual.melt, stat=="upper")[1:6,"breakeven.year"]))

################################################################################
# Figure C5: Plot species richness v. social benefit of carbon 

# get species richness estimates
n.df = stock.df[stock.df$model == "strip.random" & stock.df$variable == "n.total",
                c("full.treatment.name","posterior.mean","X5","X95")]
colnames(n.df) = c("trt","mean","lower","upper")
n.df.melt = melt(n.df, id.vars=c("trt"), variable.name="stat",value.name="richness")

# join richness estimates 
ecoC.n.df.join = left_join(ecoC.df.total.melt, n.df.melt, by=c("trt","stat"))

# plot net carbon benefit v. richness
stock.df = pivot_wider(ecoC.n.df.join[,c("trt","stat","stock")],
                       names_from = "stat", values_from = "stock")
richness.df = pivot_wider(ecoC.n.df.join[,c("trt","stat","richness")],
                          names_from = "stat", values_from = "richness")
colnames(richness.df)[2:4] = paste("richness", stats, sep=".")
colnames(stock.df)[2:4] = paste("stock", stats, sep=".")
stock.richness.df = left_join(stock.df, richness.df, by=c("trt"))
scc.df$stat = c("Mean","Lower (5%)","Upper (95%)")
colnames(scc.df)[1] = "Social cost of carbon estimate"
p1 = ggplot(data=stock.richness.df) +
            geom_point(aes(y=richness.mean, 
                           x=stock.mean, 
                           color=factor(trt, levels=trt.names)),
                       size=1.5) + 
            geom_errorbar(aes(y=richness.mean, 
                              xmin=stock.lower, 
                              xmax=stock.upper,
                              color=factor(trt, levels=trt.names)),
                          orientation="y") + 
            geom_errorbar(aes(x=stock.mean, 
                              ymin=richness.lower, 
                              ymax=richness.upper,
                              color=factor(trt, levels=trt.names)),
                          orientation="x") + 
            guides(color="none") + 
            scale_y_continuous(breaks=seq(0,30,by=10),limits=c(0,30)) +
            scale_x_continuous(breaks=seq(0,300,by=100),limits=c(-1,350)) +
            labs(x="Total organic carbon stock\nbeyond row-crop basline (Mg/ha)",
                 y="Total species richness") +
            theme(text=element_text(size=12))
carbon.benefit.df = pivot_wider(ecoC.n.df.join[,c("trt","stat","carbon.benefit")],
                                names_from = "stat", values_from = "carbon.benefit")
colnames(carbon.benefit.df)[2:4] = paste("carbon.benefit", stats, sep=".")
carbon.benefit.richness.df = left_join(carbon.benefit.df, richness.df, by=c("trt"))
p2 = ggplot(data=carbon.benefit.richness.df) +
            geom_point(aes(y=richness.mean, x=carbon.benefit.mean/1000, 
                           color=factor(trt, levels=trt.names)),
                       size=1.5) + 
            geom_errorbar(aes(y=richness.mean, 
                              xmin=carbon.benefit.lower/1000, 
                              xmax=carbon.benefit.upper/1000,
                              color=factor(trt, levels=trt.names)),
                          orientation="y") + 
            geom_errorbar(aes(x=carbon.benefit.mean/1000, 
                              ymin=richness.lower, 
                              ymax=richness.upper,
                              color=factor(trt, levels=trt.names)),
                          orientation="x") + 
            guides(color="none") + 
            scale_y_continuous(breaks=seq(0,30,by=10),limits=c(0,30)) +
            scale_x_continuous(breaks=seq(-0,400,by=100),limits=c(-10,410)) +
            labs(x="Total social carbon benefit ($1,000/ha)",
                 y="") +
            theme(axis.text.y=element_blank(),
                  text=element_text(size=12))
net.benefit.df = pivot_wider(ecoC.n.df.join[,c("trt","stat","net.benefit")],
                                names_from = "stat", values_from = "net.benefit")
colnames(net.benefit.df)[2:4] = paste("net.benefit", stats, sep=".")
net.benefit.richness.df = left_join(net.benefit.df, richness.df, by=c("trt"))
p3 = ggplot(data=net.benefit.richness.df) +
            geom_point(aes(y=richness.mean, 
                           x=net.benefit.mean/1000, 
                           color=factor(trt, levels=trt.names)),
                       size=1.5) + 
            geom_errorbar(aes(y=richness.mean, 
                              xmin=net.benefit.lower/1000, 
                              xmax=net.benefit.upper/1000,
                              color=factor(trt, levels=trt.names)),
                          orientation="y") + 
            geom_errorbar(aes(x=net.benefit.mean/1000, 
                              ymin=richness.lower, 
                              ymax=richness.upper,
                              color=factor(trt, levels=trt.names)),
                          orientation="x") + 
            guides(color="none") +  
            scale_y_continuous(breaks=seq(0,30,by=10),limits=c(0,30)) +
            scale_x_continuous(breaks=seq(-200,300,by=100),limits=c(-200,364)) +
            labs(x="Relative economic benefit ($1,000/ha)",
                 y="Total species richness") +
            theme(text=element_text(size=12))
breakeven.df = pivot_wider(ecoC.n.df.join[,c("trt","stat","breakeven.scc")],
                             names_from = "stat", values_from = "breakeven.scc")
colnames(breakeven.df)[2:4] = paste("breakeven", stats, sep=".")
breakeven.richness.df = left_join(breakeven.df, richness.df, by=c("trt"))
p4 = ggplot(data=breakeven.richness.df) +
            geom_vline(data=scc.df, 
                       aes(xintercept=per.CO2, 
                           linetype=`Social cost of carbon estimate`)) +
            geom_point(aes(y=richness.mean, 
                           x=breakeven.mean, 
                           color=factor(trt, levels=trt.names)),
                       size=1.5) + 
            geom_errorbar(aes(y=richness.mean, 
                              xmin=breakeven.lower, 
                              xmax=breakeven.upper,
                              color=factor(trt, levels=trt.names)),
                          orientation="y") + 
            geom_errorbar(aes(x=breakeven.mean, 
                              ymin=richness.lower, 
                              ymax=richness.upper,
                              color=factor(trt, levels=trt.names)),
                          orientation="x") + 
            scale_y_continuous(breaks=seq(0,30,by=10),limits=c(0,30)) +
            scale_x_continuous(labels = scales::comma) +
            labs(x="Breakeven social cost of carbon ($/Mg CO2)",
                 y="",color="Treatment") +
            theme(axis.text.y=element_blank(),
                  text=element_text(size=12))
p.all = (p1 + p2)/(p3 + p4) + plot_layout(guides = "collect")
p.all
#ggsave("Manuscript/Supp_Figures/FigureC5_Social_Cost_Carbon.jpeg", 
#       plot=p.all,width=24,height=18,units="cm",dpi=1200)

################################################################################
# plot social cost of carbon as a function of GHG emissions

# methane
methane.line.df = data.frame(matrix(nrow=0, ncol=3))
colnames(methane.line.df) = c("trt","del.methane","mean")
for (i in 1:n.t) {
  max.x = ceiling(ecoC.df.total[i,"mean"]/c.mol.mass*co2.mol.mass/ch4.100y)
  x = seq(0, max.x+50, 1)
  n.x = length(x)
  methane.df.i = data.frame(matrix(nrow=n.x, ncol=3))
  colnames(methane.df.i) = c("trt","del.methane","mean")
  methane.df.i$trt = trt.names[i]
  methane.df.i$del.methane = x
  methane.df.i[1:n.x,"mean"] = pmax((ecoC.df.total[i,"mean"]/c.mol.mass*co2.mol.mass - ch4.100y*x)*scc.df["mean","per.CO2"],0)
  methane.line.df = rbind(methane.line.df, methane.df.i)
}

p.scc.methane = ggplot(methane.line.df) +
                       geom_vline(data=subset(ghg.meta, molecule == "Methane"), 
                                  aes(xintercept=mean/1000*y.sr,
                                      linetype=`Ecosystem change`),
                                  linewidth=0.75) +
                       geom_line(aes(x=del.methane,
                                     y=mean/1000,
                                     color=factor(trt, levels=trt.names)),
                                 linewidth=0.75,linetype="dotdash") +
                       xlim(0,30) + 
                       ylim(0,160) + 
                       geom_vline(xintercept=0, linewidth=1.25, color="darkgray") +
                       geom_hline(yintercept=0, linewidth=1.25, color="darkgray") +
                       labs(y="Social benefit of carbon ($1,000/ha)",
                            x="Cumulative methane emissions (Mg/ha)",
                            color="Treatment")
p.scc.methane

# nitrous oxide
nitrous.line.df = data.frame(matrix(nrow=0, ncol=3))
colnames(nitrous.line.df) = c("trt","del.nitrous","mean")
for (i in 1:n.t) {
  max.x = ceiling(ecoC.df.total[i,"mean"]/c.mol.mass*co2.mol.mass/n2o.100y)
  x = seq(0, max.x+5, 0.01)
  n.x = length(x)
  nitrous.df.i = data.frame(matrix(nrow=n.x, ncol=3))
  colnames(nitrous.df.i) = c("trt","del.nitrous","mean")
  nitrous.df.i$trt = trt.names[i]
  nitrous.df.i$del.nitrous = x
  nitrous.df.i[1:n.x,"mean"] = pmax((ecoC.df.total[i,"mean"]/c.mol.mass*co2.mol.mass - n2o.100y*x)*scc.df["mean","per.CO2"],0)
  nitrous.line.df = rbind(nitrous.line.df, nitrous.df.i)
}

p.scc.nitrous = ggplot(nitrous.line.df) +
                        geom_vline(data=subset(ghg.meta, molecule == "Nitrous oxide"), 
                                   aes(xintercept=mean/1000*y.sr,
                                       linetype=`Ecosystem change`),
                                   linewidth=0.75) +
                        geom_line(aes(x=del.nitrous,
                                      y=mean/1000,
                                      color=factor(trt, levels=trt.names)),
                                  linewidth=0.75,linetype="dotdash") +
                        xlim(0,3) + 
                        ylim(0,160) + 
                        geom_vline(xintercept=0, linewidth=1.25, color="darkgray") +
                        geom_hline(yintercept=0, linewidth=1.25, color="darkgray") +
                        theme(axis.text.y=element_blank()) +
                        labs(y="",
                             x="Cumulative nitrous oxide emissions (Mg/ha)",
                             color="Treatment")
p.scc.nitrous
p.scc.ghgs = p.scc.methane + p.scc.nitrous + plot_layout(guides = "collect")
p.scc.ghgs
#ggsave("Manuscript/Supp_Figures/FigureC6_Social_Cost_Carbon_Versus_GHGs.jpeg", 
#       plot=p.scc.ghgs, width=24, height=9, units="cm", dpi=1200)

################################################################################
# Fugre C7: 3D plot

library(plotly)
library(pracma)

# calculate social benefit of carbon as function of N2O and CH4
x.methane = seq(0,30,0.1)
y.nitrous = seq(0,3,0.01)
grid.ghg = meshgrid(x.methane, y.nitrous)
surface_list = list()
for (t in 1:n.t) {
  trt = trt.names[t]
  scc = grid.ghg$X
  for (i in 1:length(y.nitrous)) {
    for (j in 1:length(x.methane)) {
      trt.cstock = ecoC.df.total[t,"mean"]/c.mol.mass*co2.mol.mass
      scc[i,j] = (trt.cstock - ch4.100y*grid.ghg$X[i,j] - n2o.100y*grid.ghg$Y[i,j])*scc.df["mean","per.CO2"]
    }
  }
  surface_list[[trt]] = scc
}
colorscale_red <- list(c(0, 1), c("#F8766D", "#F8766D")) # Solid Red\
colorscale_gold <- list(c(0, 1), c("#B79F00", "#B79F00")) # Solid Red
colorscale_green <- list(c(0, 1), c("#00BA38", "#00BA38")) # Solid Green
colorscale_cyan <- list(c(0, 1), c("#00BFC4", "#00BFC4")) # Solid Purple
colorscale_blue <- list(c(0, 1), c("#619CFF", "#619CFF")) # Solid Red
colorscale_pink <- list(c(0, 1), c("#F564E3", "#F564E3")) # Solid Red
colorscale_gray <- list(c(0, 1), c("#808080", "#808080")) # Solid Red
colorscale_purple <- list(c(0, 1), c("#800080", "#800080")) # Solid Red
colorscale_orange <- list(c(0, 1), c("#FFA500", "#FFA500")) # Solid Red

p1  = plot_ly(showscale = FALSE) %>%
         add_surface(x = x.methane,
                     y = y.nitrous,
                     z = ~surface_list[[trt.names[1]]]/1000,
                     colorscale = colorscale_red,
                     name = "Balled-and-burlapped") %>%
         add_surface(x = x.methane,
                     y = y.nitrous,
                     z = ~surface_list[[trt.names[2]]]/1000,
                     colorscale = colorscale_gold,
                     name = "Balled-and-burlapped") %>%
         add_surface(x = x.methane,
                     y = y.nitrous,
                     z = ~surface_list[[trt.names[3]]]/1000,
                     colorscale = colorscale_green) %>%
         add_surface(x = c(0, 30),
                     y = c(0, 3),
                     z = matrix(0, nrow=2, ncol=2),
                    colorscale = colorscale_gray,
                    opacity = 0.9) %>%
         add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[1,"mean"]/1000*y.sr),
                     y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[1,"mean"]/1000*y.sr),
                     z = matrix(0, nrow=2, ncol=2),
                     colorscale = colorscale_purple) %>%
         add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[2,"mean"]/1000*y.sr),
                     y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[2,"mean"]/1000*y.sr),
                     z = matrix(0, nrow=2, ncol=2),
                     colorscale = colorscale_orange) %>%
        layout(scene = list(xaxis = list(title = "Methane (kg/ha)"),
                            yaxis = list(title = "Nitrous oxide (kg/ha)"),
                            zaxis = list(title = "Carbon benefit ($1,000/ha")))
p1
p2  = plot_ly(showscale = FALSE) %>%
        add_surface(x = x.methane,
                    y = y.nitrous,
                    z = ~surface_list[[trt.names[4]]]/1000,
                    colorscale = colorscale_cyan) %>%
        add_surface(x = x.methane,
                    y = y.nitrous,
                    z = ~surface_list[[trt.names[5]]]/1000,
                    colorscale = colorscale_blue) %>%
        add_surface(x = x.methane,
                    y = y.nitrous,
                    z = ~surface_list[[trt.names[6]]]/1000,
                    colorscale = colorscale_pink) %>%
        add_surface(x = c(0, 30),
                    y = c(0, 3),
                    z = matrix(0, nrow=2, ncol=2),
                    colorscale = colorscale_gray,
                    opacity = 0.9) %>%
        add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[1,"mean"]/1000*y.sr),
                    y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[1,"mean"]/1000*y.sr),
                    z = matrix(0, nrow=2, ncol=2),
                    colorscale = colorscale_purple) %>%
        add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[2,"mean"]/1000*y.sr),
                    y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[2,"mean"]/1000*y.sr),
                    z = matrix(0, nrow=2, ncol=2),
                    colorscale = colorscale_orange)  %>%
        layout(scene = list(xaxis = list(title = "Methane (kg/ha)"),
                            yaxis = list(title = "Nitrous oxide (kg/ha)"),
                            zaxis = list(title = "Carbon benefit ($1,000/ha")))
p2

################################################################################
# run simulate for social benefit of carbon
set.seed(3141)

n = 10000
ghg.meta.crop = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland")
ghg.crop.ch4 = subset(ghg.meta.crop, molecule == "Methane")[,stats]*y.sr/1000
ghg.crop.n2o = subset(ghg.meta.crop, molecule == "Nitrous oxide")[,stats]*y.sr/1000
baseline.cstock.range = runif(n=n, 
                              min=baseline.cstock$lower,
                              max=baseline.cstock$upper)
n2o.emission.range = runif(n=n,
                           min=ghg.crop.n2o$lower,
                           max=ghg.crop.n2o$upper)
ch4.emission.range = runif(n=n,
                           min=ghg.crop.ch4$lower,
                           max=ghg.crop.ch4$upper)
scc.range = runif(n=n,
                  min=scc.df["lower","per.CO2"],
                  max=scc.df["upper","per.CO2"])
sim.carbon.accrual = data.frame(matrix(nrow=n, ncol=n.t))
sim.carbon.seq = data.frame(matrix(nrow=n, ncol=n.t))
sim.carbon.benefit = data.frame(matrix(nrow=n, ncol=n.t))
sim.econ.benefit = data.frame(matrix(nrow=n, ncol=n.t))
sim.richness = data.frame(matrix(nrow=n, ncol=n.t))
colnames(sim.carbon.benefit) = trt.names
colnames(sim.carbon.accrual) = trt.names
colnames(sim.carbon.seq) = trt.names
colnames(sim.econ.benefit) = trt.names
colnames(sim.richness) = trt.names
for (i in 1:n.t) {
  trt.id = which(ecoC.df.total$trt == trt.names[i])
  cost.i = est.df[est.df$trt == trt.names[i],"cost.2023"]
  sim.richness[1:n,trt.names[i]] = runif(n=n,
                                         min=richness.df[richness.df$trt == trt.names[i],]$richness.lower,
                                         max=richness.df[richness.df$trt == trt.names[i],]$richness.upper)
  trt.total.oc.df = ecoC.df.total[trt.id,]
  trt.total.oc = runif(n=n, 
                       min=trt.total.oc.df$lower,
                       max=trt.total.oc.df$upper)
  carbon.accrual.C = trt.total.oc - baseline.cstock.range
  carbon.accrual.CO2 = carbon.accrual.C/c.mol.mass*co2.mol.mass
  sim.carbon.accrual[1:n,trt.names[i]] = carbon.accrual.C
  carbon.seq = carbon.accrual.CO2 - ch4.100y*ch4.emission.range - n2o.100y*n2o.emission.range
  sim.carbon.seq[1:n,trt.names[i]] = carbon.seq
  social.carbon.benefit = carbon.seq * scc.range
  sim.carbon.benefit[1:n,trt.names[i]] = social.carbon.benefit
  rel.econ.benefit = social.carbon.benefit - cost.i
  sim.econ.benefit[1:n,trt.names[i]] = rel.econ.benefit
}

# OC accrual
sim.carbon.accrual.melt = melt(sim.carbon.accrual)
p.oc.accrual = ggplot(sim.carbon.accrual.melt,
                      aes(y=factor(variable, level=trt.names),
                          x=value,
                          fill=factor(variable, levels=trt.names))) +
                      geom_vline(xintercept=0, color="black") +
                      geom_violin(color=NA, alpha=0.5) +
                      geom_boxplot(width=0.2) + 
                      theme(legend.position="none") +
                      scale_x_continuous(limits=c(-100,300),
                                         breaks=seq(-100,300,100)) +
                      labs(y="",x="Organic carbon accrual (Mg C/ha)") +
                      annotate("text", x = -Inf, y = Inf, label = "a",
                               hjust = -1, vjust = 1.5, 
                               size = 8, family = "sans", fontface = "plain") 
p.oc.accrual

# carbon sequestration
sim.carbon.seq.melt = melt(sim.carbon.seq)
p.oc.seq = ggplot(sim.carbon.seq.melt,
                  aes(y=factor(variable, level=trt.names),
                      x=value,
                      fill=factor(variable, levels=trt.names))) +
                  geom_vline(xintercept=0, color="black") +
                  geom_violin(color=NA, alpha=0.5) +
                  geom_boxplot(width=0.2) + 
                  theme(axis.text.y=element_blank(),
                        legend.position="none") +
                  scale_x_continuous(limits=c(-490,850),
                                     breaks=seq(-400,800,200)) +
                  labs(y="",x="Carbon sequestration (Mg CO2/ha)") +
                  annotate("text", x = -Inf, y = Inf, label = "b",
                           hjust = -1, vjust = 1.5, 
                           size = 8, family = "sans", fontface = "plain") 
p.oc.seq

# social carbon benefit
sim.carbon.benefit.melt = melt(sim.carbon.benefit)
p.scc = ggplot(sim.carbon.benefit.melt,
                aes(y=factor(variable, level=trt.names),
                    x=value/1000,
                    fill=factor(variable, levels=trt.names))) +
                geom_vline(xintercept=0, color="black") +
                geom_violin(color=NA, alpha=0.5) +
                geom_boxplot(width=0.2) +  
                theme(legend.position="none") +
                scale_x_continuous(limits=c(-205,400),
                                   breaks=seq(-200,400,100)) + 
                labs(y="",fill="",
                     x="Total social carbon benefit ($1,000/ha)") +
                annotate("text", x = -Inf, y = Inf, label = "c",
                         hjust = -1, vjust = 1.5, 
                         size = 8, family = "sans", fontface = "plain") 
p.scc

# relative economic benefit
sim.econ.benefit.melt = melt(sim.econ.benefit)
p.econ = ggplot(sim.econ.benefit.melt,
                aes(y=factor(variable, level=trt.names),
                    x=value/1000,
                    fill=factor(variable, levels=trt.names))) +
                geom_vline(xintercept=0, color="black") +
                geom_violin(color=NA, alpha=0.5) +
                geom_boxplot(width=0.2) +  
                theme(axis.text.y=element_blank(),
                      legend.position="none") +
                scale_x_continuous(limits=c(-325,325),
                                   breaks=seq(-300,300,100)) + 
                labs(y="",x="Relative economic benefit ($1,000/ha)") +
                annotate("text", x = -Inf, y = Inf, label = "d",
                         hjust = -1, vjust = 1.5, 
                         size = 8, family = "sans", fontface = "plain") 
p.econ

p.sim = (p.oc.accrual+p.oc.seq)/(p.scc+p.econ)
p.sim

ggsave("Manuscript/Supp_Figures/FigureC5_Social_Cost_Carbon_Simulation_by_Treatment.jpeg", 
       plot=p.sim, width=20, height=16, units="cm", dpi=600)
       
# richness
sim.n.carbon.all = data.frame(matrix(nrow=0, ncol=6))
colnames(sim.n.carbon.all) = c("trt","n","accrual","seq","benefit","econ")
for (i in 1:n.t) {
  sim.n.carbon = data.frame(matrix(nrow=n, ncol=6))
  colnames(sim.n.carbon) = c("trt","n","accrual","seq","benefit","econ")
  sim.n.carbon$trt = trt.names[i]
  sim.n.carbon$n = sim.richness[,trt.names[i]]
  sim.n.carbon$accrual = sim.carbon.accrual[,trt.names[i]]
  sim.n.carbon$seq = sim.carbon.seq[,trt.names[i]]
  sim.n.carbon$benefit = sim.carbon.benefit[,trt.names[i]]/1000
  sim.n.carbon$econ = sim.econ.benefit[,trt.names[i]]/1000
  sim.n.carbon.all = rbind(sim.n.carbon.all, sim.n.carbon)
}
sim.nc.all.melt = subset(melt(sim.n.carbon.all, id.vars=c("trt","n")), 
                         trt != "Reference")
p.acc.n = ggplot(subset(sim.nc.all.melt, variable == "accrual"),
                 aes(y=n,
                     x=value,
                     fill=factor(trt, levels=trt.names),
                     color=factor(trt, levels=trt.names))) +
                 geom_vline(xintercept=0, color="black") +
                 geom_point(alpha=0.02) +
                 stat_ellipse(geom = "polygon", 
                              alpha = 0.2, level = 0.90) +
                 theme(legend.position = "none") + 
                 labs(y="Total species richness",
                      x="Organic carbon accrual (Mg C/ha)") +
                 scale_y_continuous(limits=c(0,30)) +
                 scale_x_continuous(limits=c(-100,300),
                                    breaks=seq(-100,300,100)) +
                 annotate("text", x = -Inf, y = Inf, label = "a",
                           hjust = -1, vjust = 1.5, 
                           size = 8, family = "sans", fontface = "plain")
p.seq.n = ggplot(subset(sim.nc.all.melt, variable == "seq"),
                    aes(y=n,
                        x=value,
                        fill=factor(trt, levels=trt.names),
                        color=factor(trt, levels=trt.names))) +
                 geom_vline(xintercept=0, color="black") +
                 geom_point(alpha=0.02) +
                 stat_ellipse(geom = "polygon", 
                              alpha = 0.2, level = 0.90) +
                 theme(legend.position = "none",
                       axis.text.y = element_blank()) + 
                 labs(y="",
                      x="Carbon sequestration (Mg CO2/ha)") +
                 scale_y_continuous(limits=c(0,30)) +
                 scale_x_continuous(limits=c(-490,850),
                                    breaks=seq(-400,800,200)) +
                 annotate("text", x = -Inf, y = Inf, label = "b",
                          hjust = -1, vjust = 1.5, 
                          size = 8, family = "sans", fontface = "plain")
p.ben.n = ggplot(subset(sim.nc.all.melt, variable == "benefit"),
                 aes(y=n,
                     x=value,
                     fill=factor(trt, levels=trt.names),
                     color=factor(trt, levels=trt.names))) +
                  geom_vline(xintercept=0, color="black") +
                  geom_point(alpha=0.02) +
                  stat_ellipse(geom = "polygon", 
                               alpha = 0.2, level = 0.90) +
                  theme(legend.position = "none") + 
                  labs(y="Total species richness",
                       x="Total social carbon benefit ($1,000/ha)") +
                  scale_y_continuous(limits=c(0,30)) +
                  scale_x_continuous(limits=c(-205,400),
                                     breaks=seq(-200,400,100)) +
                  annotate("text", x = -Inf, y = Inf, label = "c",
                           hjust = -1, vjust = 1.5, 
                           size = 8, family = "sans", fontface = "plain")
p.econ.n = ggplot(subset(sim.nc.all.melt, variable == "econ"),
                 aes(y=n,
                     x=value,
                     fill=factor(trt, levels=trt.names),
                     color=factor(trt, levels=trt.names))) +
                  geom_vline(xintercept=0, color="black") +
                  geom_point(alpha=0.02) +
                  stat_ellipse(geom = "polygon", 
                               alpha = 0.2, level = 0.90) +
                  theme(legend.position = "none",
                        axis.text.y = element_blank()) + 
                  labs(y="",x="Relative economic benefit ($1,000/ha)") +
                  scale_y_continuous(limits=c(0,30)) +
                  scale_x_continuous(limits=c(-325,325),
                                     breaks=seq(-300,300,100)) +
                  annotate("text", x = -Inf, y = Inf, label = "d",
                           hjust = -1, vjust = 1.5, 
                           size = 8, family = "sans", fontface = "plain")
p.sim.n = (p.acc.n+p.seq.n)/(p.ben.n+p.econ.n)
p.sim.n

ggsave("Manuscript/Supp_Figures/FigureC6_Social_Cost_Carbon_Simulation_Richness_Versus.jpeg", 
       plot=p.sim.n, width=18, height=16, units="cm", dpi=600)
