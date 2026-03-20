path_to_repo= "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/floodplain-experiment-repo"
setwd(path_to_repo)

library(ggplot2)
library(reshape2)
library(patchwork)
library(dplyr)

### script that uses social cost of carbon to estimate absolute and relative 
### carbon benefits of each treatment

################################################################################
# load data

# treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# social cost of carbon estimates 
scc.df = read.csv("Carbon_Calculations/SCC_Estimates.csv")
colnames(scc.df) = c("stat","per.CO2","per.CO2.C")
scc.df$stat = c("mean","lower","upper")
rownames(scc.df) = c("mean","lower","upper")

# establishment costs
est.df = read.csv("Carbon_Calculations/Treatment_Establishment_Costs.csv")
colnames(est.df) = c("trt","cost.2019","cost.2023")
est.df$cost.2023 = as.numeric(est.df$cost.2023)

# ecosystem carbon estimates
stock.df = read.csv("Tree_Analysis/Posteriors/Carbon_Stocks_Richness_Means_Intervals_10Chains_NaturalScale.csv")
ecoC.df = stock.df[stock.df$model == "strip.random" & stock.df$variable.label == "Total ecosystem",
                   c("full.treatment.name","posterior.mean","X5","X95")]
colnames(ecoC.df) = c("trt","mean","lower","upper")

# constants and conversions
co2.molecular.mass = 44.01
c.molecular.mass = 12.01

################################################################################
## combine datasets to estimate carbon benefit

# melt ecosystem carbon estimates by statistic
ecoC.df.melt = melt(ecoC.df, 
                    id.vars=c("trt"), 
                    variable.name="stat",
                    value.name="stock")

# estimate carbon benefit of each treatment
ecoC.df.melt$carbon.benefit = 0
stats = c("mean","lower","upper")
n.s = length(stats)
for (i in 1:n.s) {
  stat.i = stats[i]
  stat.id = which(ecoC.df.melt$stat == stat.i)
  ecoC.df.melt[stat.id,"carbon.benefit"] = ecoC.df.melt[stat.id,"stock"]*scc.df[stat.i,"per.CO2.C"]
}

# estimate net benefit of each treatment by subtracting the establishment cost
ecoC.df.melt$net.benefit = 0
ecoC.df.melt$breakeven.scc = 0
for (i in 1:n.t) {
  cost.i = est.df[est.df$trt == trt.names[i],"cost.2023"]
  trt.id = which(ecoC.df.melt$trt == trt.names[i])
  ecoC.df.melt[trt.id,"net.benefit"] = ecoC.df.melt[trt.id,"carbon.benefit"] - cost.i
  ecoC.df.melt[trt.id,"breakeven.scc"] = cost.i/ecoC.df.melt[trt.id,"stock"] * c.molecular.mass/co2.molecular.mass
}

# get species richness estimates
n.df = stock.df[stock.df$model == "strip.random" & stock.df$variable == "n.total",
                c("full.treatment.name","posterior.mean","X5","X95")]
colnames(n.df) = c("trt","mean","lower","upper")
n.df.melt = melt(n.df, id.vars=c("trt"), variable.name="stat",value.name="richness")

# join richness estimates 
ecoC.n.df.join = left_join(ecoC.df.melt, n.df.melt, by=c("trt","stat"))

# plot net carbon benefit v. richness
mean.df = ecoC.n.df.join[ecoC.n.df.join$stat == "mean",]
scc.df$stat = c("Mean","Lower (5%)","Upper (95%)")
colnames(scc.df)[1] = "Social cost of carbon estimate"
p1 = ggplot(data=mean.df, 
            aes(y=richness, x=stock, 
                color=factor(trt, levels=trt.names))) +
            geom_point(size=4) + 
            guides(color="none") + 
            scale_y_continuous(breaks=seq(6,22,by=2),limits=c(6,22)) +
            scale_x_continuous(breaks=seq(100,300,by=50),limits=c(100,305)) +
            labs(x="Ecosystem carbon stock (Mg/ha)",
                 y="Total species richness") +
            theme(text=element_text(size=12))
p2 = ggplot(data=mean.df, 
            aes(y=richness, x=carbon.benefit/1000, 
                color=factor(trt, levels=trt.names))) +
            geom_point(size=4) + 
            guides(color="none") + 
            scale_y_continuous(breaks=seq(6,22,by=2),limits=c(6,22)) +
            scale_x_continuous(breaks=seq(50,250,by=50),limits=c(50,250)) +
            labs(x="Total carbon benefit ($1,000/ha)",
                 y="") +
            theme(axis.text.y=element_blank(),
                  text=element_text(size=12))
p3 = ggplot(data=mean.df, 
            aes(y=richness, x=net.benefit/1000, 
                color=factor(trt, levels=trt.names))) +
            geom_point(size=4) +
            guides(color="none") + 
            scale_y_continuous(breaks=seq(6,22,by=2),limits=c(6,22)) +
            scale_x_continuous(breaks=seq(-50,200,by=50),limits=c(-50,200)) +
            labs(x="Relative economic benefit ($1,000/ha)",
                 y="Total species richness") +
            theme(text=element_text(size=12))
p4 = ggplot() +
     geom_point(data=mean.df, 
                aes(y=richness, x=breakeven.scc, 
                    color=factor(trt, levels=trt.names)),
                size=4) +
     guides(color=guide_legend(title="Treatment")) + 
     scale_y_continuous(breaks=seq(6,22,by=2),limits=c(6,22)) +
     labs(x="Breakeven social cost of carbon ($/Mg CO2)",
          y="") +
     theme(axis.text.y=element_blank(),
           text=element_text(size=12)) +
     geom_vline(data=scc.df, 
                aes(xintercept=per.CO2, 
                linetype=`Social cost of carbon estimate`))
p.all = (p1 + p2)/(p3 + p4) + plot_layout(guides = "collect")
ggsave("Supp_Figures/FigureC5_Social_Cost_Carbon.jpeg", 
       plot=p.all,width=24,height=16,units="cm",dpi=600)

# print results for table 
cols = colnames(ecoC.n.df.join[3:7])
for (i in 1:length(cols)) {
  mean.df = ecoC.n.df.join[which(ecoC.n.df.join$stat == "mean"),c("trt",cols[i])]
  l.df = ecoC.n.df.join[which(ecoC.n.df.join$stat == "lower"),c("trt",cols[i])]
  u.df = ecoC.n.df.join[which(ecoC.n.df.join$stat == "upper"),c("trt",cols[i])]
  print(cols[i])
  for (j in 1:6) {
    trt = mean.df$trt[j]
    if (cols[i] == "stock") {
      print(paste(trt, ": ", 
                  signif(mean.df[j,cols[i]],3), " (", 
                  signif(l.df[j,cols[i]],3), "-", 
                  signif(u.df[j,cols[i]],3), ")", sep="")) 
    } else if (cols[i] %in% c("carbon.benefit","net.benefit")) {
      print(paste(trt, ": ", 
                  signif(mean.df[j,cols[i]]/1000,3), " (", 
                  signif(l.df[j,cols[i]]/1000,3), "-", 
                  signif(u.df[j,cols[i]]/1000,3), ")", sep=""))
    } else {
      print(paste(trt, ": ", 
                  signif(mean.df[j,cols[i]],3), " (", 
                  signif(l.df[j,cols[i]],3), "-", 
                  signif(u.df[j,cols[i]],3), ")", sep=""))
    }
  }
}
