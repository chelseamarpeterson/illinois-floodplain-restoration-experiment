path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(dplyr)

################################################################################
# load metadata

# treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# number of plots per treatment
trt.plot.strip.df = read.csv("Metadata/Treatments_Strips_Plots.csv")
plots = unique(trt.plot.strip.df$Plot)
n.p = length(plots)

# read in conversions, constants, and plot dimensions 
dim.df = read.csv("Metadata/Constants_Conversions_Dimensions.csv")
dim.list = list()
for (i in 1:nrow(dim.df)) { dim.list[[dim.df[i,"Name"]]] = dim.df[i,"Value"] }

################################################################################
# estimatic soil C stocks

# read in soil data
soil.data = read.csv("Soil_Analysis/Clean_Data/Soil_Data_by_Quadrat_June2023.csv", header=T)

# estimate SOC stocks (Mg/ha)
soil.data$tic.stock = soil.data$tic.percent * soil.data$bulk.density * dim.list[["auger.depth"]] # g/cm2 = Mg/ha
soil.data$maoc.stock = soil.data$maoc.percent * soil.data$bulk.density * dim.list[["auger.depth"]] # g/cm2 = Mg/ha
soil.data$poc.stock = soil.data$poc.percent * soil.data$bulk.density * dim.list[["auger.depth"]] # g/cm2 = Mg/ha
soil.data$soc.stock = soil.data$toc.percent * soil.data$bulk.density * dim.list[["auger.depth"]] # g/cm2 = Mg/ha
soil.data$tc.stock = soil.data$bulk.c.percent * soil.data$bulk.density * dim.list[["auger.depth"]] # g/cm2 = Mg/ha

# average results across plots
soil.aves = soil.data %>% 
            group_by(full.treatment.name, treatment, plot) %>%
            summarize(maoc.stock = mean(maoc.stock),
                      poc.stock = mean(poc.stock),
                      soc.stock = mean(soc.stock),
                      tic.stock = mean(tic.stock),
                      tc.stock = mean(tc.stock))

# update column names
colnames(soil.aves)[4:8] = c("MAOC","POC","SOC","TIC","TC")

################################################################################
# read in vegetation carbon stocks

# read in woody biomass/debris and understory C stock data
c.data = read.csv("Tree_Analysis/Clean_Data_By_Plot/All_Vegetation_C_Stocks_By_Plot.csv")

# join vegetation and soil c stocks
c.data = right_join(c.data, 
                    soil.aves, 
                    by=c("treatment","full.treatment.name","plot"))

# estimate aggregate C stocks in different pools
c.data = c.data %>%
         mutate(total.abg.wood.carbon = abg.live.stem.carbon + abg.live.tree.carbon,
                total.bg.wood.carbon = bg.live.stem.carbon + bg.live.tree.carbon,
                total.fwd.carbon = int.fwd.carbon + small.fwd.carbon,
                total.snag.carbon = snag.carbon.min + dead.stem.carbon.min,
                total.live.carbon = total.abg.wood.carbon + total.bg.wood.carbon + herbaceous.biomass.c.stock,
                total.dead.carbon = total.snag.carbon + total.cwd.carbon + total.fwd.carbon + herbaceous.litter.c.stock,
                total.veg.carbon = total.live.carbon + total.dead.carbon,
                total.ecosystem.carbon = TC + total.veg.carbon)

################################################################################
## estimating richness by in tree and understory layer

# read in tree species data
tree.C.df = read.csv("Tree_Analysis/Clean_Data_By_Species/WoodyBiomass_C_Stocks_By_Species.csv", header=T)

# combine genus and specific epithet
tree.C.df$spp = paste(tree.C.df$genus, tree.C.df$species, sep=" ")

# read in species cover data
cover.df = read.csv("Understory_Analysis/Clean_Data/Clean_Cover_Data_Sep2022.csv")
colnames(cover.df) = tolower(colnames(cover.df))

# make dataframe for total unique species
total.sp.df = data.frame(matrix(nrow=n.t*n.p, ncol=6))
colnames(total.sp.df) = c("treatment","plot","n.herb.only",
                          "n.tree.only","n.both.herb.tree","n.total")
k = 1
for (t in 1:n.t) {
  for (n in 1:n.p) {
    tree.id = which(tree.C.df$treatment == trt.letters[t] & tree.C.df$plot == plots[n])
    tree.df = tree.C.df[tree.id,]
    tree.sp = sort(unique(tree.df$spp))
    n.tree.sp = length(tree.sp)
    
    herb.id = which(cover.df$treatment == trt.letters[t] & cover.df$plot == plots[n])
    herb.df = cover.df[herb.id,]
    herb.sp = sort(unique(herb.df$spp))
    n.herb.sp = length(herb.sp)
    
    n.total.sp = length(union(tree.sp, herb.sp))
    n.shared.sp = length(intersect(tree.sp, herb.sp))
    total.sp.df[k, c("treatment","plot")] = c(trt.letters[t], as.integer(plots[n]))
    total.sp.df[k, "n.herb.only"] = n.herb.sp - n.shared.sp
    total.sp.df[k, "n.tree.only"] = n.tree.sp - n.shared.sp
    total.sp.df[k, "n.both.herb.tree"] = n.shared.sp
    total.sp.df[k, "n.total"] = n.total.sp
    k = k + 1
  }
}

# add richness values to biomass data frame
c.data$plot = as.numeric(c.data$plot)
total.sp.df$plot = as.numeric(total.sp.df$plot)
c.sp.data = right_join(c.data, total.sp.df, by=c("treatment","plot"))

# add column for treatment strip
trt.strip.plt.df = read.csv("Metadata/Treatments_Strips_Plots.csv")
colnames(trt.strip.plt.df) = tolower(colnames(trt.strip.plt.df))
c.sp.data = left_join(c.sp.data,
                      trt.strip.plt.df,
                      by=c("treatment","plot"))

# re-order columns
id.cols = c("treatment","full.treatment.name","strip","plot")
var.cols = colnames(c.sp.data)[4:41]
c.sp.data = c.sp.data[,c(id.cols,var.cols)]

# write data to file
write.csv(c.sp.data, "Tree_Analysis/Clean_Data_By_Plot/Clean_Veg_Soil_C_Stocks_Richness_by_Plot.csv", row.names=F)
