path_to_tree_folder = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/Floodplain-Experiment-Repo"
setwd(path_to_tree_folder)

library(allodb)
library(dplyr)
library(tidyr)
library(reshape2)

################################################################################
# load metadata

# treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# read in conversions, constants, and plot dimensions 
dim.df = read.csv("Metadata/Constants_Conversions_Dimensions.csv")
dim.list = list()
for (i in 1:nrow(dim.df)) { dim.list[[dim.df[i,"Name"]]] = dim.df[i,"Value"] }

# DBH data years
years = c(2022,2023)
n.y = length(years)

# list of genus families
family.df = read.csv("Tree_Analysis/Tree_Databases/Tree_Families.csv")
n.f = nrow(family.df)
genus.families = list()
for (i in 1:n.f) { genus.families[[family.df$Genus[i]]] = family.df$Family[i] }

################################################################################
# step 1: estimate above- and belowground biomass from diameter-at-breast height (DBH) data (cm) 

# read in all DBH data
dbh.all = read.csv("Tree_Analysis/Clean_Data_By_Species/DBH_Data_Clean_2022_2023.csv", header=T)

# make list for unique species
uni.spp = sort(unique(dbh.all$spp))
n.sp = length(uni.spp)

# read in allometric equation matrices
allo.df.jenkins = read.csv("Tree_Analysis/Tree_Databases/Jenkins2004.csv", header=T)
allo.df.chojnacky = read.csv("Tree_Analysis/Tree_Databases/Chojnacky2014.csv", header=T)

# update row names
rownames(allo.df.jenkins) = allo.df.jenkins$spp
rownames(allo.df.chojnacky) = allo.df.chojnacky$spp

# estimate biomass with allometric eqns in 1: Jenkins (2004) & 2: Chojnacky (2014)
new.mat = data.frame(matrix(nrow = length(dbh.all$plot), ncol = 6))
colnames(new.mat) = c("ab.jenkins","r.jenkins","bg.jenkins",
                      "ab.chojnacky","r.chojnacky","bg.chojnacky")
dbh.all = cbind(dbh.all, new.mat)
for (sp in uni.spp) {
  sp.ind = which(dbh.all$spp == sp)
  
  # Jenkins
  dbh.all[sp.ind,"ab.jenkins"] = exp(allo.df.jenkins[sp,"ab.b0"] + allo.df.jenkins[sp,"ab.b1"]*log(dbh.all[sp.ind,"dbh.cm"]))/dim.list[["kg.per.Mg"]] # bm = exp(b0 + b1*log(dbh [cm])) [kg -> Mg]
  dbh.all[sp.ind,"r.jenkins"] = exp(allo.df.jenkins[sp,"cr.b0"] + allo.df.jenkins[sp,"cr.b1"]/dbh.all[sp.ind,"dbh.cm"]) # ratio = exp(b0 + b1/(dbh [cm]))
  dbh.all[sp.ind,"bg.jenkins"] = dbh.all[sp.ind,"r.jenkins"]*dbh.all[sp.ind,"ab.jenkins"] 
  
  # Chojnacky
  dbh.all[sp.ind,"ab.chojnacky"] = exp(allo.df.chojnacky[sp,"ab.b0"] + allo.df.chojnacky[sp,"ab.b1"]*log(dbh.all[sp.ind,"dbh.cm"]))/dim.list[["kg.per.Mg"]] # bm = exp(b0 + b1*log(dbh [cm])) [kg -> Mg]
  dbh.all[sp.ind,"r.chojnacky"] = exp(allo.df.chojnacky[sp,"cr.b0"] + allo.df.chojnacky[sp,"cr.b1"]*log(dbh.all[sp.ind,"dbh.cm"])) # ratio = exp(b0 + b1*log(dbh [cm]))
  dbh.all[sp.ind,"bg.chojnacky"] = dbh.all[sp.ind,"r.chojnacky"]*dbh.all[sp.ind,"ab.chojnacky"]
}

# estimate biomass with allodb
dbh.all$ab.allodb = get_biomass(dbh = dbh.all$dbh.cm, 
                                genus = dbh.all$genus, 
                                species = dbh.all$species, 
                                coords = c(dim.list[["longitude"]], 
                                           dim.list[["latitude"]]))/dim.list[["kg.per.Mg"]] # [kg -> Mg]
dbh.all$bg.allodb = dbh.all$ab.allodb * dbh.all$r.chojnacky # using Chojnacky root:shoot ratio
mean(dbh.all$r.chojnacky)

# check for NAs
sum(is.na(dbh.all[,c("ab.jenkins","ab.chojnacky","r.jenkins","r.chojnacky","bg.jenkins","bg.chojnacky","ab.allodb","bg.allodb")]))

# sum above & belowground biomass by species within each plot
biomass.by.species = dbh.all %>% 
                     group_by(redo, year, treatment, full.treatment.name, strip, plot, spp, genus, species) %>% 
                     summarise(ab.jenkins = sum(ab.jenkins), 
                               ab.chojnacky = sum(ab.chojnacky), 
                               ab.allodb = sum(ab.allodb), 
                               bg.jenkins = sum(bg.jenkins), 
                               bg.chojnacky = sum(bg.chojnacky), 
                               bg.allodb = sum(bg.allodb))

# divide biomass by total plot area (Mg/ha)
bm.vars = c("ab.jenkins","ab.chojnacky","ab.allodb","bg.jenkins","bg.chojnacky","bg.allodb")
bm.vars.ha = c("ab.ha.jenkins","ab.ha.chojnacky","ab.ha.allodb","bg.ha.jenkins","bg.ha.chojnacky","bg.ha.allodb")
biomass.by.species[, bm.vars] = biomass.by.species[, bm.vars]/dim.list[["plot.area.ha"]]
colnames(biomass.by.species)[which(colnames(biomass.by.species) %in% bm.vars)] = bm.vars.ha

################################################################################
# step 2: estimate carbon stocks from biomass (Mg C/ha)

# clean wood database
wood.c.df = read.csv("Tree_Analysis/Tree_Databases/Doraisami_2021_Wood_C_Database.csv")
live.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "alive" & wood.c.df$growth.form == "tree"),]
live.stem.c.df = live.wood.c.df[which(live.wood.c.df$tissue == "stem"),]
live.stem.c.df = live.stem.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
live.stem.c.df$species = paste(live.stem.c.df$genus, live.stem.c.df$spp, sep= " ")

# get best-available wood C concentration for each species
biomass.by.species$taxon.c.content = rep(0, nrow(biomass.by.species))
biomass.by.species$binomial = paste(biomass.by.species$genus, biomass.by.species$species, sep=" ")
biomass.by.species$family = rep(NA, nrow(biomass.by.species))
for (i in 1:nrow(biomass.by.species)) {
  if (biomass.by.species$binomial[i] != "") {
    biomass.by.species$family[i] = genus.families[[biomass.by.species$genus[i]]]
    df.sp.id = which(live.stem.c.df$species == biomass.by.species$binomial[i])
    df.gen.id = which(live.stem.c.df$genus.resolved == biomass.by.species$genus[i]) 
    df.fam.id = which(live.stem.c.df$family.resolved == biomass.by.species$family[i]) 
    if (length(df.sp.id) > 0) {
      biomass.by.species[i,"taxon.c.content"] = mean(live.stem.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      biomass.by.species[i,"taxon.c.content"] = mean(live.stem.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      biomass.by.species[i,"taxon.c.content"] = mean(live.stem.c.df[df.fam.id,"tissue.c"])/100
    } else {
      biomass.by.species[i,"taxon.c.content"] = tree.C.content
    }
  } else {
    biomass.by.species[i,"taxon.c.content"] = tree.C.content
  }
}

# estimate C stock by species
id.vars = c("redo","year","treatment","full.treatment.name","strip","plot","spp","genus","species")
C.vars.ha = c("abC.ha.jenkins","abC.ha.chojnacky","abC.ha.allodb",
              "bgC.ha.jenkins","bgC.ha.chojnacky","bgC.ha.allodb")
C.stocks.by.species = biomass.by.species[, c(id.vars, bm.vars.ha)]
C.stocks.by.species[,bm.vars.ha] = biomass.by.species[,bm.vars.ha] * biomass.by.species$taxon.c.content
colnames(C.stocks.by.species)[which(colnames(C.stocks.by.species) %in% bm.vars.ha)] = C.vars.ha

# sum C stocks by plot (Mg/ha)
C.stocks.by.plot <- C.stocks.by.species %>% 
                    group_by(redo, year, treatment, full.treatment.name, strip, plot) %>% 
                    summarise(abC.ha.jenkins = sum(abC.ha.jenkins), 
                              abC.ha.chojnacky = sum(abC.ha.chojnacky), 
                              abC.ha.allodb = sum(abC.ha.allodb), 
                              bgC.ha.jenkins = sum(bgC.ha.jenkins), 
                              bgC.ha.chojnacky = sum(bgC.ha.chojnacky), 
                              bgC.ha.allodb = sum(bgC.ha.allodb))

# estimate total biomass C stocks per ha
C.stocks.by.plot = C.stocks.by.plot %>% 
                   group_by(redo, year, treatment, full.treatment.name, strip, plot) %>% 
                   summarise(abC.ha.jenkins = abC.ha.jenkins, 
                             abC.ha.chojnacky = abC.ha.chojnacky, 
                             abC.ha.allodb = abC.ha.allodb,
                             bgC.ha.jenkins = bgC.ha.jenkins, 
                             bgC.ha.chojnacky = bgC.ha.chojnacky, 
                             bgC.ha.allodb = bgC.ha.allodb,
                             bmC.ha.jenkins = abC.ha.jenkins + bgC.ha.jenkins,
                             bmC.ha.chojnacky = abC.ha.chojnacky + bgC.ha.chojnacky,
                             bmC.ha.allodb = abC.ha.allodb + bgC.ha.allodb)

################################################################################
# step 3: replace bad 2022 data with 2023 redo data

# redo plots from 2023
redo.ind = which(C.stocks.by.plot$redo == "Y")
redo.plot.df = C.stocks.by.plot[redo.ind,]
redo.plot.df$treatment_plot = paste(redo.plot.df$treatment, redo.plot.df$plot, sep="")
redo.trt.plots = unique(redo.plot.df$treatment_plot)

# original plots from 2022
C.stocks.by.plot$treatment_plot = paste(C.stocks.by.plot$treatment, C.stocks.by.plot$plot, sep="")
orig.ind = which((C.stocks.by.plot$redo == "N" & C.stocks.by.plot$year == 2022) & (C.stocks.by.plot$treatment_plot %in% redo.trt.plots))
C.stocks.by.plot[orig.ind,"treatment_plot"]

# remove C stocks from 2022 that were re-done in 2023 by plot
C.stocks.by.plot.redo = C.stocks.by.plot[-orig.ind,]
C.stocks.by.plot.redo[which(C.stocks.by.plot.redo$year == 2023),"year"] = 2022

# remove C stocks from 2022 that were re-done in 2023 by species
C.stocks.by.species$treatment_plot = paste(C.stocks.by.species$treatment, C.stocks.by.species$plot, sep="")
orig.ind.spp = which((C.stocks.by.species$redo == "N" & C.stocks.by.species$year == 2022) & (C.stocks.by.species$treatment_plot %in% redo.trt.plots))
C.stocks.by.species.redo = C.stocks.by.species[-orig.ind.spp,]
C.stocks.by.species.redo[which(C.stocks.by.species.redo$year == 2023),"year"] = 2022

# estimate total biomass C stocks per ha by species
C.stocks.by.species.redo = C.stocks.by.species.redo %>% 
                           group_by(treatment, full.treatment.name, strip, plot, spp, genus, species) %>%
                           summarise(abC.ha.jenkins = abC.ha.jenkins, 
                                     abC.ha.chojnacky = abC.ha.chojnacky, 
                                     abC.ha.allodb = abC.ha.allodb,
                                     bgC.ha.jenkins = bgC.ha.jenkins, 
                                     bgC.ha.chojnacky = bgC.ha.chojnacky, 
                                     bgC.ha.allodb = bgC.ha.allodb,
                                     bmC.ha.jenkins = abC.ha.jenkins + bgC.ha.jenkins,
                                     bmC.ha.chojnacky = abC.ha.chojnacky + bgC.ha.chojnacky,
                                     bmC.ha.allodb = abC.ha.allodb + bgC.ha.allodb)

# sort dataframes by treatment
plot.sort = sort(C.stocks.by.plot.redo$treatment_plot, index.return=T)
C.stocks.by.plot.redo = C.stocks.by.plot.redo[plot.sort$ix,-which(colnames(C.stocks.by.plot) %in% c("redo","year","treatment_plot"))]

# write C stocks results to file
write.csv(C.stocks.by.plot.redo, "Tree_Analysis/Clean_Data_By_Plot/Woody_Biomass_Carbon_Stocks_By_Plot.csv", row.names=F)
write.csv(C.stocks.by.species.redo, "Tree_Analysis/Clean_Data_By_Species/Woody_Biomass_Carbon_Stocks_By_Species.csv", row.names=F)
