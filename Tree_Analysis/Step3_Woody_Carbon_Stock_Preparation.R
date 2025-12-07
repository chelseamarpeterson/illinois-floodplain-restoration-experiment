path_to_repo = "C:/Users/Chels/OneDrive - University of Illinois - Urbana/Ch2_Floodplain_Experiment/Floodplain-Experiment-Repo"
setwd(path_to_repo)

library(tidyr)
library(dplyr)
library(reshape2)
library(allodb)

## script that integrates all biomass carbon stock data and writes the results to a CSV file

################################################################################
# name important variables and load wood C content database

# read in conversions, constants, and plot dimensions 
dim.df = read.csv("Metadata/Constants_Conversions_Dimensions.csv")
dim.list = list()
for (i in 1:nrow(dim.df)) { dim.list[[dim.df[i,"Name"]]] = dim.df[i,"Value"] }

# treatments
trt.df = read.csv("Metadata/Treatment_Letters_Names.csv")
trt.letters = trt.df[,"Treatment.letters"]
trt.names = trt.df[,"Treatment.names"]
n.t = nrow(trt.df)

# number of plots per treatment
trt.plot.strip.df = read.csv("Metadata/Treatments_Strips_Plots.csv")
colnames(trt.plot.strip.df) = tolower(colnames(trt.plot.strip.df))
n.p = length(unique(trt.plot.strip.df$plot))

# list of species families
family.df = read.csv("Tree_Analysis/Tree_Databases/Tree_Families.csv")
n.f = nrow(family.df)
genus.families = list()
for (i in 1:n.f) { genus.families[[family.df$Genus[i]]] = family.df$Family[i] }
families = family.df$Family

## read in woody C data 
wood.c.df = read.csv("Tree_Analysis/Tree_Databases/Doraisami_2021_Wood_C_Database.csv")

# isolate rows for dead material
dead.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "dead"),]
dead.wood.c.df = dead.wood.c.df[which(dead.wood.c.df$tissue == "stem"),]
dead.wood.c.df = dead.wood.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
dead.wood.c.df$species = paste(dead.wood.c.df$genus, dead.wood.c.df$spp, sep= " ")

# isolate rows for live material
live.wood.c.df = wood.c.df[which(wood.c.df$dead.alive == "alive" & wood.c.df$growth.form == "tree"),]
live.wood.c.df = live.wood.c.df[which(live.wood.c.df$position == "standing" & live.wood.c.df$tissue == "stem"),]
live.wood.c.df = live.wood.c.df %>% separate(binomial.resolved, c("genus","spp"), sep="_", remove=F)
live.wood.c.df$species = paste(live.wood.c.df$genus, live.wood.c.df$spp, sep= " ")

################################################################################
# clean coarse and fine woody debris data

## load dataframe
cwd.data = read.csv("Tree_Analysis/Raw_Data/CWD_June2023.csv")

# update column names
colnames(cwd.data) = c("treatment_plot","position","diameter.cm","species","fwd.count")

# split treatment_plot column
cwd.data = cwd.data %>% separate(treatment_plot, c("treatment","plot"), sep=1, remove=T)

# trim white space for species column
cwd.data$species = trimws(cwd.data$species)

# add column for full treatment name
cwd.data$full.treatment.name = 0
for (i in 1:n.t) { cwd.data$full.treatment.name[which(cwd.data$treatment == trt.letters[i])] = trt.names[i] }

## make data frame for fine woody debris (2.5-7.5 cm) counts (#)
fwd.counts = cwd.data[which(is.na(cwd.data$diameter.cm)), c("treatment","full.treatment.name","plot","fwd.count")]
row.names(fwd.counts) = seq(1, nrow(fwd.counts))
treatment.sort = sort(fwd.counts$treatment, index.return=T)
fwd.counts = fwd.counts[treatment.sort$ix,]

# estimate FWD volume (m^3/ha) and C storage (Mg/ha)
wd.vol.coef = (1/dim.list[["plot.length.ft"]]) * (dim.list[["C1"]]*(pi^2)/8)
fwd.counts$fwd.vol = wd.vol.coef * fwd.counts$fwd.count * ((dim.list[["median.int.fwd.diameter.cm"]] / dim.list[["cm.per.in"]])^2) #m3/ha
fwd.counts$int.fwd.carbon = fwd.counts$fwd.vol * dim.list[["cwd.wood.density"]] * dim.list[["fwd.C.content"]]  # m3/ha * g/cm3 * g/g = Mg/ha

## estimate CWD stocks by plot
cwd.diams = cwd.data[-which(is.na(cwd.data$diameter.cm)), -which(colnames(cwd.data) == "fwd.count")]
cwd.diams = cwd.diams %>% separate(species, c("genus","spp"), sep=" ", remove=F)

# add column for potential C content by taxon
cwd.diams$taxon.c.content = 0
cwd.diams$family = rep(NA, nrow(cwd.diams))
for (i in 1:nrow(cwd.diams)) {
  if (cwd.diams$species[i] != "") {
    cwd.diams$family[i] = genus.families[[cwd.diams$genus[i]]]
    df.sp.id = which(dead.wood.c.df$species == cwd.diams$species[i])
    df.gen.id = which(dead.wood.c.df$genus.resolved == cwd.diams$genus[i]) 
    df.fam.id = which(dead.wood.c.df$family.resolved == cwd.diams$family[i]) 
    if (length(df.sp.id) > 0) {
      cwd.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      cwd.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      cwd.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.fam.id,"tissue.c"])/100
    } else {
      cwd.diams[i,"taxon.c.content"] = dim.list[["cwd.C.content"]]
    }
  } else {
    cwd.diams[i,"taxon.c.content"] = dim.list[["cwd.C.content"]]
  }
}

# estimate CWD volume (m^3/ha) and carbon storage (Mg/ha)
cwd.diams$cwd.vol = wd.vol.coef * (cwd.diams$diameter/dim.list[["cm.per.in"]])^2 #m3/ha
cwd.diams$cwd.carbon = cwd.diams$cwd.vol * dim.list[["cwd.wood.density"]] * cwd.diams$taxon.c.content  # Mg/ha

# sum carbon stocks (Mg/ha) by plot
cwd.sum = cwd.diams %>% 
          group_by(treatment, full.treatment.name, plot) %>% 
          summarize(total.cwd.carbon = sum(cwd.carbon))

# add zeros to CWD c stocks
cwd.sum$plot = as.character(cwd.sum$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(cwd.sum$treatment == trt.letters[i] & cwd.sum$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = colnames(cwd.sum)
      zero.row[1, c("treatment","full.treatment.name","plot")] = c(trt.letters[i], trt.names[i], j)
      zero.row[1, "total.cwd.carbon"] = 0
      cwd.sum = rbind(cwd.sum, zero.row)
    }
  }
}
treatment.sort = sort(cwd.sum$treatment, index.return=T)
cwd.sum = cwd.sum[treatment.sort$ix,]

## estimate CWD stocks by species 

# sum carbon stocks (Mg/ha) by species
cwd.diams$species[which(cwd.diams$species == "")] = "Unknown"
cwd.sp.sum = cwd.diams %>% 
             group_by(treatment, full.treatment.name, plot, species) %>% 
             summarize(total.cwd.carbon = sum(cwd.carbon))

# reshape data frame
cwd.sp.wide = pivot_wider(cwd.sp.sum, 
                          id_cols = c(treatment, plot),
                          names_from = species, 
                          values_from = total.cwd.carbon)
cwd.sp.wide[is.na(cwd.sp.wide)] = 0

# add zeros to CWD stocks
cwd.spp = colnames(cwd.sp.wide)[3:6]
cwd.sp.wide$plot = as.character(cwd.sp.wide$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(cwd.sp.wide$treatment == trt.letters[i] & cwd.sp.wide$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=length(colnames(cwd.sp.wide)), nrow=1))
      colnames(zero.row) = colnames(cwd.sp.wide)
      zero.row[1, c("treatment","plot")] = c(trt.letters[i], j)
      zero.row[1, cwd.spp] = 0
      cwd.sp.wide = rbind(cwd.sp.wide, zero.row)
    }
  }
}
treatment.sort = sort(cwd.sp.wide$treatment, index.return=T)
cwd.sp.wide = cwd.sp.wide[treatment.sort$ix,]

# add strip number to CWD by species dataframe
cwd.sp.wide$plot = as.integer(cwd.sp.wide$plot)
cwd.sp.wide = left_join(trt.plot.strip.df, 
                        cwd.sp.wide, 
                        by=c("treatment","plot"))

# reshape dataframe
cwd.sp.melt = melt(cwd.sp.wide, 
                   id.vars=c("treatment","strip","plot"), 
                   variable.name="species",
                   value.name="total.cwd.carbon")

# write cleaned cwd data to csv
write.csv(cwd.sp.melt, "Tree_Analysis/Clean_Data_By_Species/CWD_Carbon_Stocks_By_Species.csv", row.names=F)

################################################################################
# clean snag data

# load dataframe
snag.dbh.data = read.table("Tree_Analysis/Raw_Data/DBH_Snags_June2023.csv", header=T, sep=",")

# update column names
colnames(snag.dbh.data) = c("redo","treatment_plot","live","species","dbh.cm","stem.count","notes")

# split plot name column
snag.dbh.data = snag.dbh.data %>% separate(treatment_plot, c("treatment","plot"), sep=1, remove=T)

# trim white space for species column
snag.dbh.data$species = trimws(snag.dbh.data$species)

# add full treatment name
snag.dbh.data$full.treatment.name = 0
for (i in 1:n.t) { snag.dbh.data$full.treatment.name[which(snag.dbh.data$treatment == trt.letters[i])] = trt.names[i] }

# isolate all dead trees
snag.data = snag.dbh.data[which(snag.dbh.data$live == "D"),]

# make dataframe for diameter measurements
snag.diams = snag.data[-which(snag.data$dbh == "<2.5"),]

# split species column
snag.diams = snag.diams %>% separate(species, c("genus","spp"), sep=" ", remove=F)

# calculate snag volume (m3/ha)
snag.diams$dbh.cm = as.numeric(snag.diams$dbh.cm) # cm
snag.diams$basal.area = (pi * (snag.diams$dbh.cm/2)^2) / dim.list[["plot.area.m2"]] #m2/ha = cm2/m2
snag.diams$vol.min = snag.diams$basal.area * (dim.list[["min.snag.height.ft"]] / dim.list[["ft.per.m"]]) #m3/ha

# calculate snag C storage
snag.diams$taxon.c.content = 0
snag.diams$family = rep(NA, nrow(snag.diams))
for (i in 1:nrow(snag.diams)) {
  if (snag.diams$species[i] != "Unknown") {
    snag.diams$family[i] = genus.families[[snag.diams$genus[i]]]
    df.sp.id = which(dead.wood.c.df$species == snag.diams$species[i])
    df.gen.id = which(dead.wood.c.df$genus.resolved == snag.diams$genus[i]) 
    df.fam.id = which(dead.wood.c.df$family.resolved == snag.diams$family[i]) 
    if (length(df.sp.id) > 0) {
      snag.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      snag.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      snag.diams[i,"taxon.c.content"] = mean(dead.wood.c.df[df.fam.id,"tissue.c"])/100
    } else {
      snag.diams[i,"taxon.c.content"] = dim.list[["snag.C.content"]]
    }
  } else {
    snag.diams[i,"taxon.c.content"] = dim.list[["snag.C.content"]]
  }
}
snag.diams$carbon.min = snag.diams$vol.min * dim.list[["snag.wood.density"]] * snag.diams$taxon.c.content  # Mg/ha

# calculate total basal area within each plot
snag.sum = snag.diams %>% 
           group_by(treatment, full.treatment.name, plot) %>% 
           summarize(snag.carbon.min = sum(carbon.min))

# add zeros to diameter sums
snag.sum$plot = as.character(snag.sum$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(snag.sum$treatment == trt.letters[i] & snag.sum$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = colnames(snag.sum)
      zero.row[1, c("treatment","full.treatment.name","plot")] = c(trt.letters[i], trt.names[i], j)
      zero.row[1, "snag.carbon.min"] = 0
      snag.sum = rbind(snag.sum, zero.row)
    }
  }
}
treatment.sort = sort(snag.sum$treatment, index.return=T)
snag.sum = snag.sum[treatment.sort$ix,]

# estimate volume and C stocks of dead stems
snag.counts = snag.data[which(snag.data$dbh == "<2.5"),]
snag.counts = snag.counts[,-which(colnames(snag.counts) %in% c("redo","live","dbh.cm","notes"))]
snag.counts$dead.stem.count = snag.counts$stem.count / dim.list[["plot.area.ha"]] # count/ha
snag.counts$dead.stem.volmin = snag.counts$dead.stem.count * pi * ((dim.list[["median.stem.diameter.cm"]]/dim.list[["cm.per.m"]])^2) * (dim.list[["min.snag.height.ft"]] / dim.list[["ft.per.m"]]) # m3/ha

# calculate dead stem C storage
snag.counts$taxon.c.content = 0 
snag.counts = snag.counts %>% separate(species, c("genus","spp"), sep = " ", remove=F)
snag.counts$family = rep(NA, nrow(snag.counts))
for (i in 1:nrow(snag.counts)) {
  if (snag.counts$species[i] != "Unknown") {
    snag.counts$family[i] = genus.families[[snag.counts$genus[i]]]
    df.sp.id = which(dead.wood.c.df$species == snag.counts$species[i])
    df.gen.id = which(dead.wood.c.df$genus.resolved == snag.counts$genus[i]) 
    df.fam.id = which(dead.wood.c.df$family.resolved == snag.counts$family[i]) 
    if (length(df.sp.id) > 0) {
      snag.counts[i,"taxon.c.content"] = mean(dead.wood.c.df[df.sp.id,"tissue.c"])/100
    } else if (length(df.gen.id > 0)) {
      snag.counts[i,"taxon.c.content"] = mean(dead.wood.c.df[df.gen.id,"tissue.c"])/100
    } else if (length(df.fam.id > 0)) {
      snag.counts[i,"taxon.c.content"] = mean(dead.wood.c.df[df.fam.id,"tissue.c"])/100
    } else {
      snag.counts[i,"taxon.c.content"] = dim.list[["snag.C.content"]]
    }
  } else {
    snag.counts[i,"taxon.c.content"] = dim.list[["snag.C.content"]]
  }
}
snag.counts$dead.stem.carbon.min = snag.counts$dead.stem.volmin * dim.list[["snag.wood.density"]] * snag.counts$taxon.c.content # Mg/ha

# sum snag counts over species
snag.count.sum = snag.counts %>% 
                 group_by(treatment, full.treatment.name, plot) %>% 
                 summarize(dead.stem.carbon.min = sum(dead.stem.carbon.min))

# add zeros to small snag count data
snag.count.sum$plot = as.character(snag.count.sum$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(snag.count.sum$treatment == trt.letters[i] & snag.count.sum$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = c("treatment","full.treatment.name","plot","dead.stem.carbon.min")
      zero.row[1, c("treatment","full.treatment.name","plot")] = c(trt.letters[i], trt.names[i], j)
      zero.row[1, "dead.stem.carbon.min"] = 0
      snag.count.sum = rbind(snag.count.sum, zero.row)
    }
  }
}
treatment.sort = sort(snag.count.sum$treatment, index.return=T)
snag.count.sum = snag.count.sum[treatment.sort$ix,]

# species-level estimates of snag C stocks
snag.sp.sum = snag.diams %>% 
              group_by(treatment, full.treatment.name, plot, species) %>% 
              summarize(snag.carbon.min = sum(carbon.min))

# reshape data frame
snag.sp.wide = pivot_wider(snag.sp.sum, 
                           id_cols = c(treatment, plot),
                           names_from = species, 
                           values_from = snag.carbon.min)
snag.sp.wide[is.na(snag.sp.wide)] = 0

# add zeros to diameter sums
snag.spp = colnames(snag.sp.wide)[3:13]
snag.sp.wide$plot = as.character(snag.sp.wide$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(snag.sp.wide$treatment == trt.letters[i] & snag.sp.wide$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=length(colnames(snag.sp.wide)), nrow=1))
      colnames(zero.row) = colnames(snag.sp.wide)
      zero.row[1, c("treatment","plot")] = c(trt.letters[i], j)
      zero.row[1, snag.spp] = 0
      snag.sp.wide = rbind(snag.sp.wide, zero.row)
    }
  }
}
treatment.sort = sort(snag.sp.wide$treatment, index.return=T)
snag.sp.wide = snag.sp.wide[treatment.sort$ix,]

# add strip number to CWD by species dataframe
snag.sp.wide$plot = as.integer(snag.sp.wide$plot)
snag.sp.wide = left_join(trt.plot.strip.df, 
                         snag.sp.wide, 
                         by=c("treatment","plot"))

# reshape dataframe
snag.melt = melt(snag.sp.wide, 
                 id.vars=c("treatment","strip","plot"), 
                 variable.name="species",
                 value.name="snag.carbon.min")
treatment.sort = sort(snag.melt$treatment, index.return=T)
snag.melt = snag.melt[treatment.sort$ix,]

# write cleaned snag data to csv
write.csv(snag.melt, "Tree_Analysis/Clean_Data_By_Species/Snag_Carbon_Stocks_By_Species.csv", row.names=F)

################################################################################
# clean live stem count data

# read in 2022 DBH data
dbh.data.2022 = read.csv('Tree_Analysis/Raw_Data/DBH_August2022.csv')

# update column names
colnames(dbh.data.2022) = c("treatment_plot","spp","dbh","stem.count")

# isolate stem counts
stem.data.2022 = dbh.data.2022[which(dbh.data.2022$dbh == "<3"), c("treatment_plot","spp","stem.count")]

# get full species names
allo.df = read.csv("Tree_Analysis/Tree_Databases/Chojnacky2014.csv")
stem.data.2022$species = 0
for (i in 1:nrow(stem.data.2022)) {
  spp.id = which(allo.df$spp == stem.data.2022$spp[i])
  stem.data.2022$species[i] = paste(allo.df[spp.id,"genus"], allo.df[spp.id,"species"], " ")
}

# split treatment_plot column
stem.data.2022 = stem.data.2022 %>% separate(treatment_plot, c("treatment","plot"), sep=1, remove=T)

# trim white space for species colum
stem.data.2022$species = trimws(stem.data.2022$species)

# read in 2023 DBH data
stem.data.2023 = snag.dbh.data[which(snag.dbh.data$dbh == "<2.5" & snag.dbh.data$live == "L"),]
stem.data.2023 = stem.data.2023[, c("treatment","plot","species","stem.count")]
stem.data.2023$species = trimws(stem.data.2023$species)

# make combined treatment-plot column
stem.data.2022$treatment_plot = paste(stem.data.2022$treatment, stem.data.2022$plot, sep="")
stem.data.2023$treatment_plot = paste(stem.data.2023$treatment, stem.data.2023$plot, sep="")

# combine 2022 and 2023 data
redo.trt.plots = unique(snag.dbh.data$treatment_plot[which(snag.dbh.data$redo == "Y")])
redo.ind = which(stem.data.2022$treatment_plot %in% redo.trt.plots) # no overlap, so safe to do simple combination
all.stem.data = rbind(stem.data.2022[,c("treatment","plot","species","stem.count")], 
                      stem.data.2023[,c("treatment","plot","species","stem.count")])
all.stem.data$species = trimws(all.stem.data$species)

# used allodb to estimate aboveground biomass for each species
colnames(all.stem.data)[which(colnames(all.stem.data) == "species")] = "spp"
all.stem.data = all.stem.data %>% separate(spp, c("genus","species"), sep=" ", remove=F)
all.stem.data$abg = get_biomass(dbh = rep(dim.list[["median.stem.diameter.cm"]], nrow(all.stem.data)),
                                genus = all.stem.data$genus, 
                                species = all.stem.data$spp, 
                                coords = c(dim.list[["longitude"]], dim.list[["latitude"]])) * all.stem.data$stem.count / dim.list[["kg.per.Mg"]] # [kg -> Mg]
all.stem.data$abg.ha = all.stem.data$abg / dim.list[["plot.area.ha"]] # Mg/ha

# get root ratios from Chojnacky and estimate belowground biomass
all.stem.data$r = exp(allo.df[1,"cr.b0"] + allo.df[1,"cr.b1"]*log(dim.list[["median.stem.diameter.cm"]])) # ratio = exp(b0 + b1*log(dbh [cm]))
all.stem.data$bg.ha = all.stem.data$abg.ha * all.stem.data$r # Mg/ha

# identify corresponding C density for each species
all.stem.data$c.content = 0
all.stem.data$level.c = 0
for (i in 1:nrow(all.stem.data)) {
  df.sp.id = which(live.wood.c.df$species == all.stem.data$species[i])
  df.gen.id = which(live.wood.c.df$genus.resolved == all.stem.data$genus[i])
  if (length(df.sp.id) == 0) {
    all.stem.data$c.content[i] = mean(live.wood.c.df[df.gen.id,"tissue.c"])/100
    all.stem.data$level.c[i] = "genus"
  } else {
    all.stem.data$c.content[i] = mean(live.wood.c.df[df.sp.id,"tissue.c"])/100
    all.stem.data$level.c[i] = "species"
  }
}

# estimate stem C content
all.stem.data$abg.live.stem.carbon = all.stem.data$abg.ha * all.stem.data$c.content # Mg C/ha
all.stem.data$bg.live.stem.carbon = all.stem.data$bg.ha * all.stem.data$c.content # Mg C/ha
all.stem.data$total.live.stem.carbon = all.stem.data$abg.live.stem.carbon + all.stem.data$bg.live.stem.carbon # Mg C/ha

# add full treatment name
all.stem.data$full.treatment.name = 0
for (i in 1:n.t) { all.stem.data$full.treatment.name[which(all.stem.data$treatment == trt.letters[i])] = trt.names[i] }

# sum data by plot
live.stem.sum = all.stem.data %>% 
                group_by(treatment, full.treatment.name, plot) %>% 
                summarize(abg.live.stem.carbon = sum(abg.live.stem.carbon),
                          bg.live.stem.carbon = sum(bg.live.stem.carbon),
                          total.live.stem.carbon = sum(total.live.stem.carbon))

# add zeros
live.stem.sum$plot = as.character(live.stem.sum$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(live.stem.sum$treatment == trt.letters[i] & live.stem.sum$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=6))
      colnames(zero.row) = colnames(live.stem.sum)
      zero.row[1, c("treatment","full.treatment.name","plot")] = c(trt.letters[i], trt.names[i], j)
      zero.row[1, c("abg.live.stem.carbon","bg.live.stem.carbon","total.live.stem.carbon")] = 0
      live.stem.sum = rbind(live.stem.sum, zero.row)
    }
  }
}
treatment.sort = sort(live.stem.sum$treatment, index.return=T)
live.stem.sum = live.stem.sum[treatment.sort$ix,]

# sum data by species
live.stem.sp.sum = all.stem.data %>% 
                   group_by(treatment, full.treatment.name, plot, spp) %>% 
                   summarize(total.live.stem.carbon = sum(total.live.stem.carbon))

# reshape data frame
stem.sp.wide = pivot_wider(live.stem.sp.sum, 
                           id_cols = c(treatment, plot),
                           names_from = spp, 
                           values_from = total.live.stem.carbon)
stem.sp.wide[is.na(stem.sp.wide)] = 0

# add zeros to diameter sums
live.stem.spp = colnames(stem.sp.wide)[3:8]
stem.sp.wide$plot = as.character(stem.sp.wide$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(stem.sp.wide$treatment == trt.letters[i] & stem.sp.wide$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(ncol=length(colnames(stem.sp.wide)), nrow=1))
      colnames(zero.row) = colnames(stem.sp.wide)
      zero.row[1, c("treatment","plot")] = c(trt.letters[i], j)
      zero.row[1, live.stem.spp] = 0
      stem.sp.wide = rbind(stem.sp.wide, zero.row)
    }
  }
}
treatment.sort = sort(stem.sp.wide$treatment, index.return=T)
stem.sp.wide = stem.sp.wide[treatment.sort$ix,]

# add strip number to CWD by species dataframe
stem.sp.wide$plot = as.integer(stem.sp.wide$plot)
stem.sp.wide = left_join(trt.plot.strip.df, 
                         stem.sp.wide, 
                         by=c("treatment","plot"))

# reshape dataframe
stem.sp.melt = melt(stem.sp.wide, 
                    id.vars=c("treatment","strip","plot"), 
                    variable.name="species",
                    value.name="total.live.stem.carbon")
treatment.sort = sort(stem.sp.melt$treatment, index.return=T)
stem.sp.melt = stem.sp.melt[treatment.sort$ix,]

# write cleaned live stem data to csv
write.csv(stem.sp.melt, "Tree_Analysis/Clean_Data_By_Species/Live_Stem_Carbon_Stocks_By_Species.csv", row.names=F)

################################################################################
## calculate hypothetical C stock if frax. pen. were alive

# isolate F. pennsylvanica diameter in snag dataframe
snag.frax.pen = snag.diams[which(snag.diams$species == "Fraxinus pennsylvanica"),]

# use allodb to calculate hypothetical aboveground biomass
snag.frax.pen$hypothetical.abg.bm.live = get_biomass(dbh = snag.frax.pen$dbh.cm, 
                                                     genus = snag.frax.pen$genus, 
                                                     species = snag.frax.pen$spp, 
                                                     coords = c(dim.list[["longitude"]], dim.list[["latitude"]])) / dim.list[["kg.per.Mg"]]
snag.frax.pen$hypothetical.abg.bm.ha = snag.frax.pen$hypothetical.abg.bm.live / dim.list[["plot.area.ha"]] # Mg/ha

# use Chojnacky allometric equations to calculate belowground to aboveground biomass ratio
snag.frax.pen$r = exp(allo.df[1,"cr.b0"] + allo.df[1,"cr.b1"] * log(dim.list[["median.int.fwd.diameter.cm"]])) # ratio = exp(b0 + b1*log(dbh [cm]))
snag.frax.pen$hypothetical.bg.bm.ha = snag.frax.pen$hypothetical.abg.bm.ha * snag.frax.pen$r # Mg/ha
snag.frax.pen$hypothetical.biomass.ha = snag.frax.pen$hypothetical.abg.bm.ha + snag.frax.pen$hypothetical.bg.bm.ha # Mg/ha

# estimate c content and multiply by biomass to get C stock
snag.frax.pen$taxon.c.content = 0
snag.frax.pen$level.c = 0
for (i in 1:nrow(snag.frax.pen)) {
  df.sp.id = which(live.wood.c.df$species == snag.frax.pen$species[i])
  df.gen.id = which(live.wood.c.df$genus.resolved == snag.frax.pen$genus[i])
  if (length(df.sp.id) == 0) {
    snag.frax.pen$taxon.c.content[i] = mean(live.wood.c.df[df.gen.id,"tissue.c"])/100
    snag.frax.pen$level.c[i] = "genus"
  } else {
    snag.frax.pen$taxon.c.content[i] = mean(live.wood.c.df[df.sp.id,"tissue.c"])/100
    snag.frax.pen$level.c[i] = "species"
  }
}
snag.frax.pen$hypothetical.c.stock = snag.frax.pen$hypothetical.biomass.ha * snag.frax.pen$taxon.c.content

# average by plot
snag.frax.sum = snag.frax.pen %>% 
                group_by(treatment, full.treatment.name, plot) %>% 
                summarize(snag.frax.live.carbon = sum(hypothetical.c.stock)) 

# fill in zeros
snag.frax.sum$plot = as.character(snag.frax.sum$plot)
for (i in 1:n.t) {
  for (j in 1:n.p) {
    treatment.plot.id = which(snag.frax.sum$treatment == trt.letters[i] & snag.frax.sum$plot == j)
    if (length(treatment.plot.id) == 0) {
      zero.row = data.frame(matrix(nrow=1, ncol=4))
      colnames(zero.row) = colnames(snag.frax.sum)
      zero.row[1, c("treatment","full.treatment.name","plot")] = c(trt.letters[i], trt.names[i], j)
      zero.row[1, "snag.frax.live.carbon"] = 0
      snag.frax.sum = rbind(snag.frax.sum, zero.row)
    }
  }
}
treatment.sort = sort(snag.frax.sum$treatment, index.return=T)
snag.frax.live = snag.frax.sum[treatment.sort$ix,]

# isolate F. pennsylvanica data in true snag C stock data, then average by plot
snag.frax.dead = snag.sp.wide[, c("treatment","plot","Fraxinus pennsylvanica")]
colnames(snag.frax.dead) = c("treatment","plot","snag.frax.dead.carbon")
snag.frax.dead$plot = as.character(snag.frax.dead$plot)
snag.frax.sum = left_join(snag.frax.live, 
                          snag.frax.dead,
                          by = c("treatment","plot"))

################################################################################
# read in calculated tree C stock data

# tree data
C.stock.data = read.csv("Tree_Analysis/Clean_Data_By_Plot/Woody_Biomass_C_Stocks_By_Plot.csv", header=T)
C.stock.data = C.stock.data[,c("treatment","full.treatment.name","plot","abC.ha3","bgC.ha3","bmC.ha3")]
C.stock.data$plot = as.character(C.stock.data$plot)

# update column names and sort data
colnames(C.stock.data)[4:6] = c("abg.live.tree.carbon","bg.live.tree.carbon","total.live.tree.carbon")
treatment.sort = sort(C.stock.data$treatment, index.return=T)
C.stock.data = C.stock.data[treatment.sort$ix,]

################################################################################
# read in C stocks for herbaceous litter, fine woody debris, and biomass

# read in understory data
understory.data = read.csv("Understory_Analysis/Clean_Data/CN_Stock_Summary_Sep2022.csv", header=T)

# convert g/m2 to Mg/ha [(g/m2) * (1 Mg/10^6 g) * (10^4 m2/ha)]
understory.data$c.mg.ha = understory.data$c.g.m2 / 100
understory.data$n.mg.ha = understory.data$n.g.m2 / 100

# reshape data frame
understory.wide = pivot_wider(understory.data, 
                              id_cols = c(treatment, plot),
                              names_from = type, 
                              values_from = c.mg.ha)

# update column names
colnames(understory.wide)[3:8] = c("herbaceous.biomass.c.stock","mixed.biomass.c.stock","p.arundinacea.c.stock",
                                   "h.japonicus.c.stock","herbaceous.litter.c.stock","small.fwd.carbon")
understory.wide$plot = as.character(understory.wide$plot)

## calculate C:N ratios of fine woody debris, herbaceous litter, and herbaceous biomass

# separate data by type
understory.fwd = understory.data[which(understory.data$type == "Fine.woody.debris"),]
understory.hl = understory.data[which(understory.data$type == "Herbaceous.litter"),]
understory.hb = understory.data[which(understory.data$type == "Herbaceous.biomass"),]

# FWD C:N ratio
understory.wide$fine.woody.debris.cn.ratio = understory.fwd$c.mg.ha / understory.fwd$n.mg.ha
understory.wide$fine.woody.debris.cn.ratio[which(is.na(understory.wide$fine.woody.debris.cn.ratio))] = 0

# HL C:N ratio
understory.wide$herbaceous.litter.cn.ratio = understory.hl$c.mg.ha / understory.hl$n.mg.ha
understory.wide$herbaceous.litter.cn.ratio[which(is.na(understory.wide$herbaceous.litter.cn.ratio))] = 0

# HN C:N ratioa
understory.wide$herbaceous.biomass.cn.ratio = understory.hb$c.mg.ha / understory.hb$n.mg.ha
understory.wide$herbaceous.biomass.cn.ratio[which(is.na(understory.wide$herbaceous.biomass.cn.ratio))] = 0

################################################################################
## combine all biomass C stock by plots dataframes 

# intermediate fine woody debris and coarse woody debris
all.bm.data = left_join(cwd.sum, fwd.counts[,c("treatment","full.treatment.name","plot","int.fwd.carbon")], 
                        by=c("treatment","full.treatment.name","plot"))

# large standing dead trees
all.bm.data = left_join(all.bm.data, snag.sum, 
                        by=c("treatment","full.treatment.name","plot"))

# small standing dead stems
all.bm.data = left_join(all.bm.data, snag.count.sum, 
                        by=c("treatment","full.treatment.name","plot"))

# large live trees
all.bm.data = left_join(all.bm.data, C.stock.data, 
                        by=c("treatment","full.treatment.name","plot"))

# small live stems
all.bm.data = left_join(all.bm.data, live.stem.sum[,c("treatment","full.treatment.name","plot",
                                                      "abg.live.stem.carbon","bg.live.stem.carbon","total.live.stem.carbon")], 
                        by=c("treatment","full.treatment.name","plot"))

# hypothetical C stocks of Frax pen
all.bm.data = left_join(all.bm.data, snag.frax.sum[,c("treatment","full.treatment.name","plot",
                                                      "snag.frax.live.carbon","snag.frax.dead.carbon")], 
                        by=c("treatment","full.treatment.name","plot"))

# understory C stocks and C:N ratios
all.bm.data = left_join(all.bm.data, understory.wide, 
                        by=c("treatment","plot"))

# write all C stock and C:N ratio to file
write.csv(all.bm.data, "Tree_Analysis/Clean_Data_By_Plot/All_Vegetation_C_Stocks_By_Plot.csv", row.names=F)