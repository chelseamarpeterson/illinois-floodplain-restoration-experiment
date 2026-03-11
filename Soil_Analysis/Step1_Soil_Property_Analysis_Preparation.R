path_to_soil_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/floodplain-experiment-repo"
setwd(path_to_soil_folder)

library(tidyr)
library(dplyr)
library(reshape2)
library(tidyverse)

### script that brings together all files with quadrat-level soil data into
### one clean dataframe

################################################################################
### data upload and preparation

## treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

## read in all soil data

# bulk density
bd.data = read.csv("Soil_Analysis/Raw_Data/Bulk_Density_June2023.csv")[c(seq(2,4),8,
                                                                         seq(13,14),
                                                                         seq(19,23))]
bd.data = bd.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep=1, remove=T)
colnames(bd.data) = c("treatment","plot","quadrat","depth",
                      "mass.water","mass.coarse.roots.fragments","mass.dry.soil",
                      "corrected.volume","gravimetric.moisture","volumetric.moisture",
                      "bulk.density","root.fragment.density")

# temperature
temp.data = read.csv("Soil_Analysis/Raw_Data/Soil_Temperature_June2023.csv")[,seq(2,5)]
colnames(temp.data) = c("treatment","plot","quadrat","temperature")
temp.data$plot = as.character(temp.data$plot)

# total carbon
tc.data = read.csv("Soil_Analysis/Raw_Data/TC_Data_Updated_2026.csv")[,c(1,12,13)]
tc.data = tc.data %>% separate(Treatment_Plot_Quadrat, c("Treatment_Plot","Quadrat"), sep=" ", remove=T)
tc.data = tc.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep=1, remove=T)
colnames(tc.data) = c("treatment","plot","quadrat","tn.percent","tc.percent")

# total organic carbon
toc.data = read.csv("Soil_Analysis/Raw_Data/TOC_Data_Updated_2026.csv")[,c(1,4,5)]
toc.data = toc.data %>% separate(Treatment_Plot_Quadrat, c("Treatment_Plot","Quadrat"), sep=" ", remove=T)
toc.data = toc.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep=1, remove=T)
colnames(toc.data) = c("treatment","plot","quadrat","toc.n.percent","toc.percent")

# particulate organic carbon
poc.data = read.csv("Soil_Analysis/Raw_Data/POC_Data_Updated_2026.csv")[,c(2,8,9)]
poc.data = poc.data %>% separate(Treatment_Plot_Quadrat, c("Treatment_Plot","Quadrat"), sep = " ", remove = T)
poc.data = poc.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep = 1, remove = T)
colnames(poc.data) = c("treatment","plot","quadrat","poc.n.percent","poc.percent")

# soil chemistry data 
chem.data = read.csv("Soil_Analysis/Raw_Data/Waypoint_Results_June2023.csv")[,seq(2,25)]
chem.data = chem.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep=1, remove=T)
colnames(chem.data) = c("treatment","plot","quadrat","sand","silt","clay","texture.class",
                        "ph","p","k","ca","mg","som.percent","n.release","no3.n","nh4.n","cec",
                        "k.sat","ca.sat","mg.sat","h.sat","k.meq","ca.meq","mg.meq","h.meq")

# quadrat data
quad.data = read.csv("Soil_Analysis/Raw_Data/Quadrat_GIS_Data_PolygonMean_WGS1984Aux_2020Lidar.csv")[,c(seq(1,3),seq(5,6))]
colnames(quad.data) = c("treatment","plot","quadrat","elevation.mean","vegetation.height.mean")
quad.data$plot = as.character(quad.data$plot)

# assume second R1 Q1 is equal to R1 Q3
bd.ind1 = which((bd.data$treatment == "R" & bd.data$plot == "1") & (bd.data$quadrat == "Q1" & bd.data$depth == "0-15"))[2]
bd.ind2 = which((bd.data$treatment == "R" & bd.data$plot == "1") & (bd.data$quadrat == "Q3" & bd.data$depth == "0-15"))
bd.data[bd.ind2, colnames(bd.data)[seq(5,8)]] = bd.data[bd.ind1, colnames(bd.data)[seq(5,8)]]
bd.data = bd.data[-bd.ind1,]

# compute moisture and bulk density for 0-30 cm
bd.sum = bd.data %>% 
         group_by(treatment, plot, quadrat) %>%
         summarise(mass.water = sum(mass.water),
                   sum.soil = sum(mass.dry.soil),
                   sum.root.fragment = sum(mass.coarse.roots.fragments),
                   sum.vol = sum(corrected.volume))
bd.sum$gravimetric.moisture = bd.sum$mass.water/bd.sum$sum.soil
bd.sum$volumetric.moisture = bd.sum$mass.water/bd.sum$sum.vol
bd.sum$bulk.density = bd.sum$sum.soil/bd.sum$sum.vol
bd.sum$root.fragment.density = bd.sum$sum.root.fragment/bd.sum$sum.vol
bd.sum = bd.sum[,-which(colnames(bd.sum) %in% c("mass.water","sum.soil","sum.vol","sum.root.fragment"))]

# average the poc reps
poc.ave = poc.data[,c("treatment","plot","quadrat","poc.percent")] %>%
          group_by(treatment, plot, quadrat) %>%
          summarize(poc.percent = mean(poc.percent))

# add full treatment name
temp.data$full.treatment.name = rep(0, nrow(temp.data))
for (i in 1:n.t) { temp.data$full.treatment.name[which(temp.data$treatment == trt.letters[i])] = trt.names[i] }

# add column for treatment strip
trt.strip.plt.df = read.csv("Metadata/Treatments_Strips_Plots.csv")
colnames(trt.strip.plt.df) = c("treatment","strip","plot")
trt.strip.plt.df$plot = as.character(trt.strip.plt.df$plot)
trt.strip.plt.df$strip = as.character(trt.strip.plt.df$strip)
soil.data = left_join(temp.data[,c("full.treatment.name","treatment","plot",
                                   "quadrat","temperature")],
                      trt.strip.plt.df, by=c("treatment","plot"))
soil.data = soil.data[,c("full.treatment.name","treatment","strip","plot",
                         "quadrat","temperature")]

## combine all soil data frames

# bulk density
soil.data = left_join(soil.data, bd.sum, by=c("treatment","plot","quadrat"))

# chemistry
soil.data = left_join(soil.data, chem.data[,-which(colnames(chem.data) %in%
                                                     c("n.release","k.sat","ca.sat",
                                                       "mg.sat","h.sat"))],
                      by=c("treatment","plot","quadrat"))

# total carbon
soil.data = left_join(soil.data, tc.data, by=c("treatment","plot","quadrat"))

# total organic carbon
soil.data = left_join(soil.data, toc.data[,which(colnames(toc.data) != "toc.n.percent")], 
                      by=c("treatment","plot","quadrat"))

# particulate organic carbon
soil.data = left_join(soil.data, poc.ave[,which(colnames(poc.ave) != "poc.n.percent")], 
                      by=c("treatment","plot","quadrat"))

# ground elevation and vegetation height
soil.data = left_join(soil.data, quad.data, 
                      by=c("treatment","plot","quadrat"))

################################################################################
## clean up aggregate data 

# update columns
ag.data = read.csv("Soil_Analysis/Raw_Data/Aggregate_Masses_2023.csv")
ag.data = ag.data %>% separate(Treatment_Plot_Quadrat, c("Treatment_Plot","Quadrat"), sep = " ", remove = T)
ag.data = ag.data %>% separate(Treatment_Plot, c("Treatment","Plot"), sep = 1, remove = F)
key.cols = c("Size.class","Treatment_Plot","Treatment","Plot","Quadrat","Mass.of.sample.w.o.fragments..g.")
ag.data = ag.data[,colnames(ag.data) %in% key.cols]
colnames(ag.data) = c("size","treatment_plot","treatment","plot","quadrat","soil.mass")

# size classes, plots, and quadrats
sizes = c("fPOM",">4.75 mm",">2 mm",">250 um",">53 um","<=53 um")        
trt_plots = sort(unique(ag.data$treatment_plot))
quads = sort(unique(ag.data$quadrat))
n.s = length(sizes)
n.tp = length(trt_plots)

# make new dataframe to estimate mass of fraction <53 um
ag.df = data.frame(matrix(nrow=0, ncol=6))
ag.cols = c("size","treatment","plot","quadrat","soil.mass","rel.mass")
colnames(ag.df) = ag.cols
ag.data$rel.mass = rep(0, nrow(ag.data))
for (p in trt_plots) {
  for (q in quads) {
    # plot and quadrat id
    pq.id = which(ag.data$treatment_plot == p & ag.data$quadrat == q)
    
    # get total mass and estimate <53 mass
    pq.full.mass = ag.data[which(ag.data[pq.id,"size"] == "Full"),"soil.mass"]
    pq.53.mass = pq.full.mass - sum(ag.data[pq.id[2:length(pq.id)],"soil.mass"])
    
    # get subset of dataframe for plot p and quad q
    pq.df = ag.data[pq.id[2:length(pq.id)], ag.cols]
    
    # add row to p and q dataframe
    pq.row = data.frame(matrix(nrow=1, ncol=6))
    colnames(pq.row) = ag.cols
    pq.row[1, "size"] = "<=53 um"
    pq.row[1, ag.cols[2:4]] = pq.df[1, ag.cols[2:4]]
    pq.row[1, c("soil.mass","rel.mass")] = c(pq.53.mass, 0)
    pq.df = rbind(pq.df, pq.row)
    
    # estimate relative masses
    pq.df[,"rel.mass"] = pq.df[,"soil.mass"]/as.numeric(pq.full.mass)*100
    
    # append new dataframe
    ag.df = rbind(ag.df, pq.df)   
  }
}

# make all categorical variables into factors
ag.df$quadrat = as.factor(ag.df$quadrat)
ag.df$size = factor(ag.df$size, levels=sizes)

# convert to long form
ag.df$size.lab = rep("", nrow(ag.df))
size.labs = c("fPOM","g475mm","g2mm","g250um","g53um","l53um")
for (i in 1:n.s) {
  size.id = which(ag.df$size == sizes[i])
  ag.df[size.id,"size.lab"] = size.labs[i]
}
ag.new = ag.df[,c("treatment","plot","quadrat","rel.mass","size.lab")]
ag.wide = pivot_wider(ag.new,id_cols = c(treatment, plot, quadrat),
                      names_from = size.lab, values_from = rel.mass)

# compute mean weight diameter from aggregate size distribution
sizes = c(0, 0.053, 0.250, 2, 4.75, 8)
mean.diams = (sizes[2:6]-sizes[1:5])/2 + sizes[1:5]
ag.wide$mwd = (mean.diams[1]*ag.wide$l53um + mean.diams[2]*ag.wide$g53um + mean.diams[3]*ag.wide$g250um + mean.diams[4]*ag.wide$g2mm + mean.diams[5]*ag.wide$g475mm)/100

# add aggregate data to soil dataframe
colnames(ag.wide) = tolower(colnames(ag.wide))
soil.data = right_join(soil.data, ag.wide, by=c("treatment","plot","quadrat"))

################################################################################
## final quadrat-level data cleaning and writing to file

# calculate inorganic and mineral-associated organic c percents
soil.data$tic.percent = soil.data$tc.percent - soil.data$toc.percent
soil.data$maoc.percent = soil.data$toc.percent - soil.data$poc.percent

# visualize data
plot(seq(1,90),soil.data$bulk.c.percent,type="n",ylim=c(0,6))
points(seq(1,90),soil.data$tc.percent,col="black")
points(seq(1,90),soil.data$toc.percent,col="blue")
points(seq(1,90),soil.data$poc.percent,col="red")
points(seq(1,90),soil.data$maoc.percent,col="green")
points(seq(1,90),soil.data$tic.percent,col="purple")

# check for C data consistency
sum(soil.data$tc.percent >= soil.data$toc.percent)
sum(soil.data$toc.percent >= soil.data$maoc.percent)
sum(soil.data$toc.percent >= soil.data$poc.percent)
sum(soil.data$toc.percent >= soil.data$tic.percent)
sum(soil.data$tic.percent >= 0)
sum(soil.data$maoc.percent >= 0)

# look at histograms
hist(soil.data$tc.percent)
hist(soil.data$toc.percent)
hist(soil.data$maoc.percent)
hist(soil.data$poc.percent)
hist(soil.data$tic.percent)

# calculate POC:MAOC ratio and C:N ratio
soil.data$poc.maoc.ratio = soil.data$poc.percent/soil.data$maoc.percent
soil.data$cn.ratio = soil.data$tc.percent/soil.data$tn.percent

# scale moisture to percentage
soil.data$gravimetric.moisture = soil.data$gravimetric.moisture * 100

# estimate all C stocks
soil.data$tc.stock = soil.data$tc.percent * soil.data$bulk.density * 30
soil.data$toc.stock = soil.data$toc.percent * soil.data$bulk.density * 30
soil.data$maoc.stock = soil.data$maoc.percent * soil.data$bulk.density * 30
soil.data$poc.stock = soil.data$poc.percent * soil.data$bulk.density * 30
soil.data$tic.stock = soil.data$tic.percent * soil.data$bulk.density * 30

# write all soil data at quadrat level to csv
write.csv(soil.data, "Soil_Analysis/Clean_Data/Soil_Data_by_Quadrat_June2023.csv", row.names=F)

