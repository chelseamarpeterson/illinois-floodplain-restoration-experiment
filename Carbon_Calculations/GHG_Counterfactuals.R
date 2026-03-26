path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/floodplain-experiment-repo"
setwd(path_to_repo)

library(tidyverse)
library(reshape2)
library(patchwork)

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
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# read in ecosystem C stock estimates
stats = c("posterior.mean","X5","X95","X25","X75")
n.s = length(stats)
stock.df = read.csv("Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_10Chains_NaturalScale.csv")
ecoC.df = stock.df[stock.df$model == "strip.random" & stock.df$variable.label == "Total organic",
                   c("full.treatment.name",stats)]

################################################################################
# individual treatment estimates

# estimate c stock relative to baseline
ecoC.accrual.df = ecoC.df
ecoC.accrual.df[,stats] = (ecoC.accrual.df[,stats]-baseline.cstock)/years.since.restoration

# estimate methane emissions necessary to offset co2 benefit
methane.offset.df = ecoC.accrual.df
methane.offset.df[,stats] = ecoC.accrual.df[,stats]*(co2.molecular.mass/c.molecular.mass)/ch4.gwp100yr*1000
methane.offset.df$molecule = "Methane"

# estimate nitrous oxide emissions necessary to offset co2 benefit
nitrous.offset.df = ecoC.accrual.df
nitrous.offset.df[,stats] = ecoC.accrual.df[,stats]*(co2.molecular.mass/c.molecular.mass)/n2o.gwp100yr*1000
nitrous.offset.df$molecule = "Nitrous oxide"

# combine methane and nitrous oxide dataframes
ghg.offset.df = rbind(methane.offset.df,nitrous.offset.df)

# read in meta-analysis GHG emission estimates
ghg.meta = read.csv("Carbon_Calculations/He_2024_Meta_Analysis_GHG_Estimates.csv")
ghg.stats = c("Average","Standard.deviation")
colnames(ghg.meta)[1:2] = c("Ecosystem change","molecule")

# plot results
p.ghg.offset = ggplot(ghg.offset.df) + 
                      geom_vline(data=ghg.meta[ghg.meta$Unit == "kg/ha/y",], 
                                 aes(xintercept=Average, 
                                     color=`Ecosystem change`,
                                     linetype=`Ecosystem change`), linewidth=1) +
                      geom_point(aes(y=factor(full.treatment.name, levels=trt.names),
                                     x=posterior.mean), size=2.5) + 
                      geom_errorbar(aes(y=factor(full.treatment.name, levels=trt.names),
                                        xmin=X25, xmax=X75), orientation="y", height=0.01, linewidth=1) +
                      geom_errorbar(aes(y=factor(full.treatment.name, levels=trt.names),
                                        xmin=X5, xmax=X95), orientation="y", width=0.2) +
                      facet_wrap(.~molecule, scales="free_x") +
                      scale_x_continuous(labels = scales::comma) + #limits = function(X5, X95) c(min(X5), max(X95))) +
                      labs(y="",x="Annual emissions needed to offset carbon accrual (kg/ha/yr)")
p.ghg.offset
ggsave("Supp_Figures/FigureC1_Greenhouse_Gas_Offsets_By_GHG_and_Treament.jpeg", 
       plot=p.ghg.offset, width=24, height=8, units="cm", dpi=600)

# make lines for combined methane and nitrous oxide emissions
n = 1000
ghg.offset.intercept.df = ecoC.accrual.df
ghg.offset.intercept.df[,stats] = ecoC.accrual.df[,stats]/n2o.gwp100yr
ghg.offset.slope = -ch4.gwp100yr/n2o.gwp100yr
ghg.offset.line.plot.df = data.frame(matrix(nrow=0,ncol=7))
colnames(ghg.offset.line.plot.df) = c("treatment","del.methane",stats)
for (i in 1:n.t) {
  trt.i = trt.names[i]
  line.plot.i = data.frame(matrix(nrow=n, ncol=7))
  colnames(line.plot.i) = c("treatment","del.methane",stats)
  line.plot.i$treatment = trt.i
  line.plot.i$del.methane = seq(0,1,length.out=n)
  for (j in 1:n.s) {
    line.plot.i[1:n,stats[j]] = ghg.offset.intercept.df[i,stats[j]] + ghg.offset.slope*line.plot.i[1:n,"del.methane"]
  }
  ghg.offset.line.plot.df = rbind(ghg.offset.line.plot.df, line.plot.i)
}

cols = c("Ecosystem change","molecule","Average")
ghg.meta = ghg.meta[ghg.meta$Unit == "kg/ha/y",]
cropland.to.wetland = pivot_wider(ghg.meta[ghg.meta$`Ecosystem change` == "Cropland to wetland",cols],
                                  names_from = "molecule", values_from = "Average")
restored.floodplain = pivot_wider(ghg.meta[ghg.meta$`Ecosystem change` == "Degraded to restored floodplain",cols],
                                  names_from = "molecule", values_from = "Average")
box.df = rbind(cropland.to.wetland, restored.floodplain)
p.ghg.abs = ggplot(ghg.offset.line.plot.df) + 
                   geom_rect(data=box.df, 
                             aes(xmin=0, xmax=Methane, 
                                 ymin=0, ymax=`Nitrous oxide`, 
                                 fill=`Ecosystem change`), 
                             alpha = 0.3, inherit.aes = FALSE) +
                   geom_line(aes(x=del.methane*1000,
                                 y=posterior.mean*1000,
                                 color=factor(treatment,levels=trt.names)),
                             linewidth=0.5,linetype="dashed") + 
                   xlim(0,350) + ylim(0,35) +
                   labs(color="Treatment",
                        x="Methane offset (kg/ha/yr)",
                        y="Nitrous oxide offset (kg/ha/yr)")
p.ghg.abs
ggsave("Supp_Figures/FigureC2_Combined_GHG_Thresholds_By_Treatment.jpeg", 
       plot=p.ghg.abs, width=16, height=10, units="cm", dpi=600)

################################################################################
# pairwise treatment comparisons

## add methane offsets to matrix
ch4.offsets = data.frame(matrix(nrow=n.t, ncol=n.t))
colnames(ch4.offsets) = trt.names
row.names(ch4.offsets) = trt.names
for (i in 1:(n.t-1)) {
  trt.i = trt.names[i]
  for (j in (i+1):n.t) {
    trt.j = trt.names[j]
    eco.cstock.i = subset(ecoC.df, full.treatment.name==trt.i)[,"posterior.mean"]
    eco.cstock.j = subset(ecoC.df, full.treatment.name==trt.j)[,"posterior.mean"]
    ch4.offsets[i,j] = (eco.cstock.i-eco.cstock.j)*(co2.molecular.mass/c.molecular.mass)/ch4.gwp100yr/years.since.restoration
  }
}

# swap negatives
for (i in 1:(n.t-1)) {
  trt.i = trt.names[i]
  for (j in (i+1):n.t) {
    trt.j = trt.names[j]
    ch4.offset.ij = ch4.offsets[i,j]
    if (ch4.offset.ij < 0) {
      ch4.offsets[i,j] = NA
      ch4.offsets[j,i] = -ch4.offset.ij
    }
  }
}

# create methane tile plot
ch4.offsets = rownames_to_column(ch4.offsets, var = "row.trt")
ch4.melt = melt(ch4.offsets, 
                id.vars = "row.trt",
                value.name="offset",
                variable.name="col.trt")
p.ch4 = ggplot(ch4.melt, 
               aes(x=factor(col.trt,levels=trt.names), 
                   y=factor(row.trt,levels=rev(trt.names)), 
                   fill=offset*1000)) +
               geom_tile(color="white", lwd=0.5) + 
               geom_text(aes(label=signif(offset*1000, 4)),
                         color="black", size=8) +
               scale_fill_gradient(low = "white", high = "blue") +
               coord_cartesian() +
               labs(x="", y="", fill="Methane offset\n(kg/ha/yr)") +
               theme(text=element_text(size=12)) + 
               geom_label(x=1, y=6, label="a",
                          color="black", fill=alpha("white",0),
                          label.r=unit(0,"pt"), label.size=0,
                          size=16, fontface="bold")
p.ch4

## add nitrous oxide offsets to matrix
n2o.offsets = data.frame(matrix(nrow=n.t, ncol=n.t))
colnames(n2o.offsets) = trt.names
row.names(n2o.offsets) = trt.names
for (i in 1:(n.t-1)) {
  trt.i = trt.names[i]
  for (j in (i+1):n.t) {
    trt.j = trt.names[j]
    eco.cstock.i = subset(ecoC.df, full.treatment.name==trt.i)[,"posterior.mean"]
    eco.cstock.j = subset(ecoC.df, full.treatment.name==trt.j)[,"posterior.mean"]
    n2o.offsets[i,j] = (eco.cstock.i-eco.cstock.j)*(co2.molecular.mass/c.molecular.mass)/n2o.gwp100yr/years.since.restoration
  }
}

# swap negatives
for (i in 1:(n.t-1)) {
  trt.i = trt.names[i]
  for (j in (i+1):n.t) {
    trt.j = trt.names[j]
    n2o.offset.ij = n2o.offsets[i,j]
    if (n2o.offset.ij < 0) {
      n2o.offsets[i,j] = NA
      n2o.offsets[j,i] = -n2o.offset.ij
    }
  }
}

# create nitrous oxide tile plot
n2o.offsets = rownames_to_column(n2o.offsets, var = "row.trt")
n2o.melt = melt(n2o.offsets, 
                id.vars = "row.trt",
                value.name="offset",
                variable.name="col.trt")
p.n2o = ggplot(n2o.melt, 
               aes(x=factor(col.trt,levels=trt.names), 
                   y=factor(row.trt,levels=rev(trt.names)), 
                   fill=offset*1000)) +
               geom_tile(color="white", lwd=0.5) + 
               geom_text(aes(label=signif(offset*1000, 3)),
                         color="black",size=8) +
               coord_cartesian() +
               scale_fill_gradient(low="white", high="darkorange") +
               labs(x="",y="",fill="Nitrous oxide\noffset (kg/ha/yr)") +
               theme(text=element_text(size=12)) + 
               geom_label(x=1,y=6,label="b",
                          color="black",fill=alpha("white",0),
                          label.r=unit(0,"pt"),label.size=0,
                          size=16,fontface="bold")
p.n2o
p.gdg.tiles = p.ch4 + p.n2o
ggsave("Supp_Figures/FigureC3_Greenhouse_Gas_Offsets_Pairwise_Comparisons.jpeg", 
       plot=p.gdg.tiles,, width=48, height=18, units="cm", dpi=600)

## calculate slopes and intercepts for combined methane and nitrous oxide threshold needed to offset CO2 benefit
line.df = data.frame(matrix(nrow=15,ncol=4))
colnames(line.df) = c("trt1","trt2","intercept","slope")
line.df$slope = -ch4.gwp100yr/n2o.gwp100yr
n = 1
for (i in 1:(n.t-1)) {
  trt.i = trt.names[i]
  for (j in (i+1):n.t) {
    trt.j = trt.names[j]
    line.df[n,"trt1"] = trt.i
    line.df[n,"trt2"] = trt.j
    eco.cstock.i = subset(ecoC.df, full.treatment.name==trt.i)[,"posterior.mean"]
    eco.cstock.j = subset(ecoC.df, full.treatment.name==trt.j)[,"posterior.mean"]
    line.df[n,"intercept"] = (eco.cstock.i-eco.cstock.j)*(co2.molecular.mass/c.molecular.mass)/n2o.gwp100yr
    n = n + 1
  }
}

# swap negatives
n = 1
for (i in 1:(n.t-1)) {
  trt.i = trt.names[i]
  for (j in (i+1):n.t) {
    trt.j = trt.names[j]
    intercept.ij = line.df[n,"intercept"]
    if (intercept.ij < 0) {
      line.df[n,"intercept"] = -intercept.ij
      line.df[n,"trt1"] = trt.j
      line.df[n,"trt2"] = trt.i
    }
    n = n + 1
  }
}

# create dataframe for each line
line.plot.df = data.frame(matrix(nrow=0, ncol=4))
colnames(line.plot.df) = c("trt1","trt2","del.methane","del.nitrous")
for (n in 1:nrow(line.df)) {
  trt.i = line.df[n,"trt1"]
  trt.j = line.df[n,"trt2"]
  line.plot.ij = data.frame(matrix(nrow=100, ncol=4))
  colnames(line.plot.ij) = c("trt1","trt2","del.methane","del.nitrous")
  line.plot.ij$trt1 = trt.i
  line.plot.ij$trt2 = trt.j
  line.plot.ij$del.methane = seq(0,28,length.out=100)
  line.plot.ij$del.nitrous = line.df[n,"intercept"] + line.df[n,"slope"]*line.plot.ij$del.methane
  line.plot.df = rbind(line.plot.df, line.plot.ij)
}

p.line.n2o.ch4 = ggplot(line.plot.df, 
                        aes(x=del.methane*1000/years.since.restoration,
                            y=del.nitrous*1000/years.since.restoration,
                            color=factor(trt2,levels=trt.names))) + 
                       geom_line(linetype="dashed",linewidth=0.5) +
                       ylim(-0.01,100) + xlim(0,1000) +
                       facet_wrap(~factor(trt1,levels=trt.names), ncol=3) +
                       labs(color="Treatment",
                            x="Methane offset (kg/ha/yr)",
                            y="Nitrous oxide offset (kg/ha/yr)") + 
                       geom_vline(xintercept=0,linewidth=0.5) +
                       geom_hline(yintercept=0,linewidth=0.5) + 
                       theme(text = element_text(size=10)) + 
                       scale_x_continuous(breaks=seq(0,1000,250),
                                          labels = scales::comma)
p.line.n2o.ch4
ggsave("Supp_Figures/FigureC4_Methane_Nitrous_Oxide_Threshold.jpeg", 
       plot=p.line.n2o.ch4,width=20,height=11,units="cm",dpi=600)

