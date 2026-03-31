path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment"
setwd(path_to_repo)

library(reshape2)
library(patchwork)
library(dplyr)
library(tidyverse)

### script that uses social cost of carbon to estimate absolute and relative 
### carbon benefits of each treatment

################################################################################
# load data

# constants and conversions
ch4.gwp100yr = 27
n2o.gwp100yr = 273
co2.molecular.mass = 44.009
c.molecular.mass = 12.011
n.molecular.mass = 14.007
years.since.restoration = 25
baseline.cstock = 55.9
ch4.molecular.mass = 16.043
n2o.molecular.mass = 44.013

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
ecoC.df.total[,stats] = ecoC.df.total[,stats] - baseline.cstock

# annual accumulation rates
ecoC.df.annual = ecoC.df.total
ecoC.df.annual[,stats] = ecoC.df.annual[,stats]/years.since.restoration

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
tot.crop.ch4 = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Methane")[1,stats]*years.since.restoration/1000 # (kg/ha/y)*(25 y)*(1 Mg/1000 kg) = Mg CH4/ha
tot.crop.n2o = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Nitrous oxide")[1,stats]*years.since.restoration/1000 # (kg/ha/y)*(25 y)*(1 Mg/1000 kg) = Mg N2O/ha
for (i in 1:n.s) {
  stat.i = stats[i]
  stat.id = which(ecoC.df.total.melt$stat == stat.i)
  ecoC.df.total.melt[stat.id,"carbon.benefit"] = (ecoC.df.total.melt[stat.id,"stock"]/c.molecular.mass*co2.molecular.mass - ch4.gwp100yr*tot.crop.ch4[1,stat.i] - n2o.gwp100yr*tot.crop.n2o[1,stat.i]) * scc.df[stat.i,"per.CO2"]
}
signif(subset(ecoC.df.total.melt, stat=="mean")[,"carbon.benefit"]/1000,3)
signif(subset(ecoC.df.total.melt, stat=="lower")[,"carbon.benefit"]/1000,3)
signif(subset(ecoC.df.total.melt, stat=="upper")[,"carbon.benefit"]/1000,3)

ecoC.df.annual.melt$carbon.benefit = 0
ann.crop.ch4 = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Methane")[1,stats]/1000 # (kg/ha/y)*(1 Mg/1000 kg) = Mg CH4/ha/y
ann.crop.n2o = subset(ghg.meta, `Ecosystem change` == "Cropland to wetland" & molecule == "Nitrous oxide")[1,stats]/1000 # (kg/ha/y)*(1 Mg/1000 kg) = Mg N2O/ha/y
for (i in 1:n.s) {
  stat.i = stats[i]
  stat.id = which(ecoC.df.annual.melt$stat == stat.i)
  ecoC.df.annual.melt[stat.id,"carbon.benefit"] = (ecoC.df.annual.melt[stat.id,"rate"]/c.molecular.mass*co2.molecular.mass - ch4.gwp100yr*ann.crop.ch4[1,stat.i] - n2o.gwp100yr*ann.crop.n2o[1,stat.i]) * scc.df[stat.i,"per.CO2"]
}
signif(subset(ecoC.df.annual.melt, stat=="mean")[,"carbon.benefit"],3)
signif(subset(ecoC.df.annual.melt, stat=="lower")[,"carbon.benefit"],3)
signif(subset(ecoC.df.annual.melt, stat=="upper")[,"carbon.benefit"],3)

# estimate net benefit of each treatment by subtracting the establishment cost
ecoC.df.total.melt$net.benefit = 0
ecoC.df.total.melt$breakeven.scc = 0
for (i in 1:n.t) {
  cost.i = est.df[est.df$trt == trt.names[i],"cost.2023"]
  trt.id = which(ecoC.df.total.melt$trt == trt.names[i])
  ecoC.df.total.melt[trt.id,"net.benefit"] = ecoC.df.total.melt[trt.id,"carbon.benefit"] - cost.i
  for (j in 1:n.s) {
    stat.j = stats[j]
    scc.ij = cost.i/(ecoC.df.total.melt[trt.id[j],"stock"]/c.molecular.mass*co2.molecular.mass - ch4.gwp100yr*tot.crop.ch4[1,stat.j] - n2o.gwp100yr*tot.crop.n2o[1,stat.j])
    if (scc.ij >= 0) {
      ecoC.df.total.melt[trt.id[j],"breakeven.scc"] = scc.ij
    } else {
      ecoC.df.total.melt[trt.id[j],"breakeven.scc"] = Inf
    }
  }
}
signif(subset(ecoC.df.total.melt, stat=="mean")[,"breakeven.scc"],3)
signif(subset(ecoC.df.total.melt, stat=="lower")[,"breakeven.scc"],3)
signif(subset(ecoC.df.total.melt, stat=="upper")[,"breakeven.scc"],3)

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
ggsave("Manuscript/Supp_Figures/FigureC5_Social_Cost_Carbon.jpeg", 
       plot=p.all,width=24,height=18,units="cm",dpi=1200)

################################################################################
# plot social cost of carbon as a function of GHG emissions

# methane
methane.line.df = data.frame(matrix(nrow=0, ncol=3))
colnames(methane.line.df) = c("trt","del.methane","mean")
for (i in 1:n.t) {
  max.x = ceiling(ecoC.df.total[i,"mean"]/c.molecular.mass*co2.molecular.mass/ch4.gwp100yr)
  x = seq(0, max.x+50, 1)
  n.x = length(x)
  methane.df.i = data.frame(matrix(nrow=n.x, ncol=3))
  colnames(methane.df.i) = c("trt","del.methane","mean")
  methane.df.i$trt = trt.names[i]
  methane.df.i$del.methane = x
  methane.df.i[1:n.x,"mean"] = pmax((ecoC.df.total[i,"mean"]/c.molecular.mass*co2.molecular.mass - ch4.gwp100yr*x)*scc.df["mean","per.CO2"],0)
  methane.line.df = rbind(methane.line.df, methane.df.i)
}

p.scc.methane = ggplot(methane.line.df) +
                       geom_vline(data=subset(ghg.meta, molecule == "Methane"), 
                                  aes(xintercept=mean/1000*years.since.restoration,
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
  max.x = ceiling(ecoC.df.total[i,"mean"]/c.molecular.mass*co2.molecular.mass/n2o.gwp100yr)
  x = seq(0, max.x+5, 0.01)
  n.x = length(x)
  nitrous.df.i = data.frame(matrix(nrow=n.x, ncol=3))
  colnames(nitrous.df.i) = c("trt","del.nitrous","mean")
  nitrous.df.i$trt = trt.names[i]
  nitrous.df.i$del.nitrous = x
  nitrous.df.i[1:n.x,"mean"] = pmax((ecoC.df.total[i,"mean"]/c.molecular.mass*co2.molecular.mass - n2o.gwp100yr*x)*scc.df["mean","per.CO2"],0)
  nitrous.line.df = rbind(nitrous.line.df, nitrous.df.i)
}

p.scc.nitrous = ggplot(nitrous.line.df) +
                        geom_vline(data=subset(ghg.meta, molecule == "Nitrous oxide"), 
                                   aes(xintercept=mean/1000*years.since.restoration,
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
ggsave("Manuscript/Supp_Figures/FigureC6_Social_Cost_Carbon_Versus_GHGs.jpeg", 
       plot=p.scc.ghgs, width=24, height=9, units="cm", dpi=1200)

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
      trt.cstock = ecoC.df.total[t,"mean"]/c.molecular.mass*co2.molecular.mass
      scc[i,j] = (trt.cstock - ch4.gwp100yr*grid.ghg$X[i,j] - n2o.gwp100yr*grid.ghg$Y[i,j])*scc.df["mean","per.CO2"]
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
         add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[1,"mean"]/1000*years.since.restoration),
                     y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[1,"mean"]/1000*years.since.restoration),
                     z = matrix(0, nrow=2, ncol=2),
                     colorscale = colorscale_purple) %>%
         add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[2,"mean"]/1000*years.since.restoration),
                     y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[2,"mean"]/1000*years.since.restoration),
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
        add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[1,"mean"]/1000*years.since.restoration),
                    y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[1,"mean"]/1000*years.since.restoration),
                    z = matrix(0, nrow=2, ncol=2),
                    colorscale = colorscale_purple) %>%
        add_surface(x = c(0, subset(ghg.meta, molecule == "Methane")[2,"mean"]/1000*years.since.restoration),
                    y = c(0, subset(ghg.meta, molecule == "Nitrous oxide")[2,"mean"]/1000*years.since.restoration),
                    z = matrix(0, nrow=2, ncol=2),
                    colorscale = colorscale_orange)  %>%
        layout(scene = list(xaxis = list(title = "Methane (kg/ha)"),
                            yaxis = list(title = "Nitrous oxide (kg/ha)"),
                            zaxis = list(title = "Carbon benefit ($1,000/ha")))
p2


