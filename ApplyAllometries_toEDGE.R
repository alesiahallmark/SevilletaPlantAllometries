### This script matches quad observations with known and inferred allometries to estimate biomass 
# Created and maintained by: Alesia Hallmark

# Load required libraries
library(ggplot2)
library(RColorBrewer)
library(lubridate)
library(plyr)
library(reshape2)
library(gridExtra)
library(cowplot)
library(knitr)
library(markdown)
library(xtable)
library(zoo)
library(tidyr)
library(data.table)
library(readxl)
library(XLConnect)

# Read in SEV297
prev.edge <- read.csv("~/Desktop/unsorted_rawdata/Sevilleta_quad_datasets/sev297_EDGE.csv", strip.white = T)
prev.edge$date <- as.Date(prev.edge$date, "%m/%d/%Y")

# Official sev297 dataset only goes to spring 2016. 
# Read in fall 2016 data from raw files
rawedgefiles <- list.files(path = "~/Desktop/EDGE/anpp/raw_data_sheets", pattern = "anpp.10", full.names = T)[list.files(path = "~/Desktop/EDGE/anpp/raw_data_sheets", pattern = "anpp.10", full.names = T) %in% list.files(path = "~/Desktop/EDGE/anpp/raw_data_sheets", pattern = "2016.xlsx", full.names = T)]

for (i in 1:length(rawedgefiles)) {
  # Read in files
  one.file <- as.data.frame(read_excel(rawedgefiles[i]))
  
  for (j in 1:nrow(one.file)) {
    if (is.na(one.file$site[j])) {
      one.file$site[j] <- one.file$site[j-1]}
    #if (is.na(one.file$block[j])) {
    #  one.file$block[j] <- one.file$block[j-1]}
    if (is.na(one.file$plot[j])) {
      one.file$plot[j] <- one.file$plot[j-1]}
    if (is.na(one.file$treatment[j])) {
      one.file$treatment[j] <- one.file$treatment[j-1]}
    if (is.na(one.file$quad[j])) {
      one.file$quad[j] <- one.file$quad[j-1]}
    if (is.na(one.file$spp[j])) {
      one.file$spp[j] <- one.file$spp[j-1]}}
  
  one.file$date <- as.Date(paste(
    as.numeric(strsplit(rawedgefiles[i], split = "[.]")[[1]][3]),
    as.numeric(strsplit(rawedgefiles[i], split = "[.]")[[1]][4]),
    as.numeric(strsplit(rawedgefiles[i], split = "[.]")[[1]][5])),
    format = "%m %d %Y")
  
  one.file$season <- 3
  
  if (i == 1) raw.edge <- one.file
  else raw.edge <- join(raw.edge, one.file, by=intersect(colnames(raw.edge), colnames(one.file)), type="full")
  
}

raw.edge$year <- year(raw.edge$date)
raw.edge$site <- as.factor(raw.edge$site)
raw.edge$block <- NULL
raw.edge$treatment <- as.factor(raw.edge$treatment)
raw.edge$spp <- as.factor(toupper(raw.edge$spp))
colnames(raw.edge)[colnames(raw.edge) %in% "spp"] <- "species"
colnames(raw.edge)[colnames(raw.edge) %in% "comments"] <- "comment"

raw.edge$site <- revalue(raw.edge$site, c(
  "b" = "EB", "g" = "EG"))
raw.edge$treatment <- revalue(raw.edge$treatment, c(
  "control" = "C", "full" = "D", "partial" = "E"))

raw.edge <- merge(raw.edge, unique(prev.edge[,c("treatment", "site", "block", "plot")]), all.x = T)


# Read in 2017 "QC'ed" file from Lauren (sent in email)
edge.2017 <- read.csv("~/Desktop/EDGE_2017_fromLauren.txt", header = T)
edge.2017$season[edge.2017$season == 1] <- 2
edge.2017$date[edge.2017$season == 2] <- "5 15 2017"
edge.2017$date[edge.2017$season == 3] <- "10 15 2017"
                                                 
edge.2017$date <- as.Date(edge.2017$date, "%m %d %Y")

# combine EDGE files
all.edge <- join(prev.edge, raw.edge, type = "full")
all.edge <- join(all.edge, edge.2017, type = "full")

all.edge$cover[all.edge$cover < 0] <- NA
all.edge$height[all.edge$height < 0] <- NA

all.edge$volume <- all.edge$cover * all.edge$height

all.edge$site <- revalue(all.edge$site, c(
  "EB" = "EDGE.blue",
  "EG" = "EDGE.black"))

all.edge$SiteCluster <- revalue(all.edge$site, c(
  "EB" = "BlueGrama",
  "EG" = "FivePoints"))

all.edge$treatment <- revalue(all.edge$treatment, c(
  "C" = "control",
  "D" = "delayed",
  "E" = "drought"))

# Fix some species names
colnames(all.edge)[colnames(all.edge) %in% "species"] <- "kartez"
all.edge$kartez <- revalue(all.edge$kartez, c(
  "LEDED" = "LEDE",
  "GACO5" = "OESU3",
  "SPWR" = "SPPO6",
  "SPOO6" = "SPPO6", 
  "MACA" = "MACA2", 
  "ERAD2" = "ERAB2"))

# Match site names with nearest met stations
all.edge$MetStation <- 
  revalue(all.edge$SiteCluster,
          c(DeepWell = 40, # Deep well
            CerroMontosa = 42, # Cerro Montosa PJ
            FivePoints = 49, # Five Points
            BlueGrama = 50)) # Blue grama

all.edge$season <- revalue(as.factor(all.edge$season), c(
  "1" = "winter", "2" = "spring", "3" = "fall"))

# Get rid of winter observations for now
all.edge <- all.edge[all.edge$season %in% c("spring", "fall"),]




# Read in allometry data
allos <- read.csv("~/Desktop/SevPlantBiomass_06Jul18_highlightsAPPLY.csv", header = T)
invar.allos <- allos[allos$Best.Invar.Model == "x" | allos$Best.Model == "x", c("kartez", "genus", "sp.epithet", "family", "LifeHistory", "PhotoPath", "FunctionalGroup", "season", "model.name", "Best.Model", "Best.Invar.Model", "main.eff.full", "mainclim.eff.full", "inter.eff.full")]


# Read in taxonomy data and merge
taxo <- read.csv("/Users/alesia/Documents/Project_SevPublicPlantRCode/SevilletaSpeciesList_AJH.csv", strip.white = T, na.strings = c("NA",""))
colnames(taxo) <- c("kartez", "family", "genus", "sp.epithet", "subspecies", "common.name", "native", "PhotoPath", "LifeHistory", "FunctionalGroup")

# Revalue a_p and g_f columns
taxo$LifeHistory <- revalue(taxo$LifeHistory, c(
  "a" = "annual",
  "p" = "perennial",
  "a/p" = "ann/peren"))
taxo$FunctionalGroup <- revalue(taxo$FunctionalGroup, c(
  "f" = "forb",
  "t" = "tree", 
  "g" = "grass",
  "s" = "shrub"))

# just call Phacelia an annual for now
taxo$LifeHistory[taxo$kartez %in% "PHACE"] <- "annual"

# read in climate data and merge
metdata <- read.csv("~/Desktop/SEVmetAJ_30Aug18.csv")
metdata$Sta <- as.factor(metdata$Sta)

# reshape into Year Season Site and reduce to columns of interest
springmetdata <- metdata[,c("Sta", "water_year", "precip.annual", "precip.spring", "GDD.spring", "spring6SPEI.thornth")]
colnames(springmetdata) <- c("MetStation", "Year", "annual.precip", "season.precip", "GDD", "SPEI6.thornth")
springmetdata$Season <- 2
monsoonmetdata <- metdata[,c("Sta", "water_year", "precip.annual", "precip.monsoon", "GDD.monsoon", "monsoon6SPEI.thornth")]
colnames(monsoonmetdata) <- c("MetStation", "Year", "annual.precip", "season.precip", "GDD", "SPEI6.thornth")
monsoonmetdata$Season <- 3

metdata <- rbind(springmetdata, monsoonmetdata)

metdata$MetStation <- as.factor(metdata$MetStation)




### Create temporary substitution allometries (sister species, other season, average of larger taxa)
## Replicate allometries for missing seasons
# fill dataset with all species-season combos
ex.sp.se <- expand.grid(
  kartez = as.factor(unique(invar.allos$kartez)), 
  season = as.factor(c("spring", "fall")))
ex.sp.se <- unique(merge(ex.sp.se, invar.allos[,c("kartez", "genus", "sp.epithet", "family", "LifeHistory", "PhotoPath", "FunctionalGroup")], all = T))

# Merge known allometries with full list of species-season combos
invar.allos <- merge(invar.allos, ex.sp.se, all = T)

# For each species, if one season doesn't have regression, fill with other season
allos.otherseason <- merge(invar.allos[is.na(invar.allos$model.name),c("kartez", "season")], invar.allos[!is.na(invar.allos$model.name), c("kartez", "genus", "sp.epithet", "family", "LifeHistory", "PhotoPath", "FunctionalGroup", "model.name", "Best.Model", "Best.Invar.Model", "main.eff.full", "mainclim.eff.full", "inter.eff.full")], all.x = T)
invar.allos <- merge(invar.allos[!is.na(invar.allos$model.name),], allos.otherseason, all = T)

# Find average slope for each plant functional group (PFG)
# Use only Volume model for substitutions
PFG.means <- ddply(allos[allos$model.name %in% c("Volume"),], c("LifeHistory", "PhotoPath", "FunctionalGroup", "season", "model.name"), function(x) data.frame(
  main.eff.full = mean(x$main.eff.full, na.rm = T)
))

# Apply other-season allometry from missing PFG-season combos
ex.PFG <- rbind(
  expand.grid(
    LifeHistory = c("annual"), 
    PhotoPath = c("C3", "C4"),
    FunctionalGroup = c("grass", "forb"),
    season = c("spring", "fall"),
    model.name = c("Volume")),
  expand.grid(
    LifeHistory = c("perennial"), 
    PhotoPath = c("C3", "C4"),
    FunctionalGroup = c("grass", "forb", "shrub"),
    season = c("spring", "fall"),
    model.name = c("Volume")),
  expand.grid(
    LifeHistory = c("perennial"), 
    PhotoPath = c("CAM"),
    FunctionalGroup = c("shrub"),
    season = c("spring", "fall"),
    model.name = c("Volume")))

PFG.means <- merge(PFG.means, ex.PFG, all = T)

# Fill in the missing PFG.means
# Give annual grasses the same (only) allometry
PFG.means$main.eff.full[PFG.means$LifeHistory %in% "annual" & PFG.means$FunctionalGroup %in% "grass" & is.na(PFG.means$main.eff.full)] <- 
  PFG.means$main.eff.full[PFG.means$LifeHistory %in% "annual" & PFG.means$FunctionalGroup %in% "grass" & !is.na(PFG.means$main.eff.full)]
# Give annual C4 forbs the same allometry
PFG.means$main.eff.full[PFG.means$LifeHistory %in% "annual" & PFG.means$PhotoPath %in% "C4" & PFG.means$FunctionalGroup %in% "forb" & is.na(PFG.means$main.eff.full)] <-
  PFG.means$main.eff.full[PFG.means$LifeHistory %in% "annual" & PFG.means$PhotoPath %in% "C4" & PFG.means$FunctionalGroup %in% "forb" & !is.na(PFG.means$main.eff.full)]
# Give perennial C4 shrubs the same allometry
PFG.means$main.eff.full[PFG.means$LifeHistory %in% "perennial" & PFG.means$PhotoPath %in% "C4" & PFG.means$FunctionalGroup %in% "shrub" & is.na(PFG.means$main.eff.full)] <-
  PFG.means$main.eff.full[PFG.means$LifeHistory %in% "perennial" & PFG.means$PhotoPath %in% "C4" & PFG.means$FunctionalGroup %in% "shrub" & !is.na(PFG.means$main.eff.full)]

# Find observations with no matching allometry
find.miss <- merge(all.edge, invar.allos, all.x = T)
find.miss <- find.miss[is.na(find.miss$main.eff.full),]
unique(find.miss$kartez)

substitution.allos <- expand.grid(
  kartez = unique(find.miss$kartez), 
  season = c("spring", "fall"))
substitution.allos <- unique(merge(substitution.allos, taxo[,c("kartez", "family", "genus", "sp.epithet", "PhotoPath", "LifeHistory", "FunctionalGroup")], all.x = T))

# Apply average slope for that functional group
substitution.allos <- merge(substitution.allos, PFG.means, by = c("LifeHistory", "PhotoPath", "FunctionalGroup", "season"), all.x = T)


# Apply a few special allometries - when we can do better than PFG mean
PFG.means.special <- ddply(allos[!allos$LifeHistory %in% "annual" & allos$genus %in% c("Aristida", "Astragalus", "Opuntia", "Phacelia", "Sporobolus", "Sphaeralcea") & allos$model.name %in% c("Volume"),], c("genus", "season", "model.name"), function(x) data.frame(
  main.eff.full = mean(x$main.eff.full, na.rm = T)
))

# Find groups that still don't have allometry and apply a substitute
substitution.allos[is.na(substitution.allos$main.eff.full),]

# Combine real and imagined allometries
all.invar.allos <- merge(invar.allos, substitution.allos, all = T)

# Merge quad data with taxonomy info
all.edge <- merge(all.edge, taxo[,c("kartez", "family", "genus", "sp.epithet", "PhotoPath", "LifeHistory", "FunctionalGroup")], all.x = T)

# Merge known allometries
quad.allo <- merge(all.edge, invar.allos, all.x = T)
# Merge substitutes: special genera allometries
quad.allo.subspec <- merge(quad.allo[is.na(quad.allo$main.eff.full), colnames(quad.allo)[!colnames(quad.allo) %in% c("model.name", "main.eff.full")]], PFG.means.special, all.x = T)
# Merge substitutes: PFG mean allometries
quad.allo.PFG <- merge(quad.allo.subspec[is.na(quad.allo.subspec$main.eff.full), colnames(quad.allo.subspec)[!colnames(quad.allo.subspec) %in% c("model.name", "main.eff.full")]], PFG.means, all.x = T)

summary(quad.allo.PFG[is.na(quad.allo.PFG$main.eff.full),])

quad.allo.all <- merge(quad.allo[!is.na(quad.allo$main.eff.full),], quad.allo.subspec[!is.na(quad.allo.subspec$main.eff.full),], all = T)
quad.allo.all <- merge(quad.allo.all, quad.allo.PFG, all = T)

# Calculate estimated biomass
# Because EDGE is a rainfall manipulation experiment, only invariant models should be applied
quad.allo.all <- quad.allo.all[quad.allo.all$model.name %in% c("Cover", "Volume"), c("year", "date", "season", "site", "kartez", "family", "genus", "sp.epithet", "LifeHistory", "PhotoPath", "FunctionalGroup", "treatment", "block", "plot", "quad", "obs", "cover", "height", "count", "volume", "SiteCluster", "MetStation", "model.name", "main.eff.full")]

quad.allo.all$biomass <- NA
quad.allo.all$biomass[quad.allo.all$model.name %in% "Volume"] <- 
  quad.allo.all$volume[quad.allo.all$model.name %in% "Volume"] * 
  quad.allo.all$main.eff.full[quad.allo.all$model.name %in% "Volume"]
quad.allo.all$biomass[quad.allo.all$model.name %in% "Cover"] <- 
  quad.allo.all$cover[quad.allo.all$model.name %in% "Cover"] * 
  quad.allo.all$main.eff.full[quad.allo.all$model.name %in% "Cover"]

quad.allo.all$biomass[quad.allo.all$kartez %in% "NONE"] <- NA

# Find illegal or outlier values

# Summarize biomass by quad and species
quad.sp.summ <- as.data.frame(data.table(quad.allo.all)[, list(
  cover = sum(cover * count, na.rm = T),
  volume = sum(volume * count, na.rm = T),
  biomass = sum(biomass * count, na.rm = T)),
  by = list(year, site, season, kartez, block, plot, treatment, quad, SiteCluster, MetStation, genus, sp.epithet, family, LifeHistory, PhotoPath, FunctionalGroup, total.quads)])

avg_dates <- as.data.frame(data.table(quad.allo.all)[, list(
  date = mean(as.Date(date), na.rm = T)),
  by = list(year, site, season)])

quad.sp.summ <- merge(quad.sp.summ, avg_dates, all.x = T)

# re-arrange column names
quad.sp.summ <- quad.sp.summ[,c("site", "year", "season", "date", "treatment", "block", "plot", "quad", "kartez", "genus", "sp.epithet", "family", "LifeHistory", "PhotoPath", "FunctionalGroup", "cover", "volume", "biomass")]

# Save for phenomass
write.csv(quad.sp.summ, "~/Desktop/unsorted_AJdata/EDGEbiomass_30Aug2018.csv", row.names = F)

# Count number of quads at each site
num.quads <- ddply(quad.sp.summ, c("site", "year", "treatment"), function(x) data.frame(
  total.quads = length(unique(paste(x$block, x$plot, x$quad)))
))
num.quads
quad.sp.summ <- merge(quad.sp.summ, num.quads, all.x = T)

# Add PFG column
quad.sp.summ$PFG <- paste(quad.sp.summ$LifeHistory, quad.sp.summ$PhotoPath, site.sp.summ$FunctionalGroup)
quad.sp.summ$PFG.PF <- paste(quad.sp.summ$PhotoPath, quad.sp.summ$FunctionalGroup)

# reshape data
quad.sp.wide <- spread(quad.sp.summ[,colnames(quad.sp.summ)[!colnames(quad.sp.summ) %in% c("date", "volume", "cover")],], key = season, value = biomass, fill = 0)
colnames(quad.sp.wide)[colnames(quad.sp.wide) %in% c("spring", "fall")] <- c("mass.spring", "mass.fall")
quad.sp.wide$mass.spring[quad.sp.wide$year == 2012] <- NA
quad.sp.wide$npp.fall <- quad.sp.wide$mass.fall - quad.sp.wide$mass.spring
quad.sp.wide$npp.fall[!is.na(quad.sp.wide$npp.fall) & quad.sp.wide$npp.fall < 0] <- 0

# Summarize biomass by site and species
biomean.site.sp <- as.data.frame(data.table(quad.sp.summ)[, list(
  volume.mean = sum(volume, na.rm = T) / total.quads,
  cover.mean = sum(cover, na.rm = T) / total.quads,
  biomass.mean = sum(biomass, na.rm = T) / total.quads),
  by = list(year, site, season, date, kartez, treatment, genus, sp.epithet, family, LifeHistory, PhotoPath, FunctionalGroup, total.quads)])

ggplot(quad.sp.summ, aes(x = date, y = biomass, group = kartez))


a <- c(1,2,3,4,5,6,7,8,9)
c(a, rep_len(0, 15-length(a)))

biomean.site.sp.wide <- as.data.frame(data.table(quad.sp.wide)[, list(
  mass.spring = sum(mass.spring, na.rm = T) / total.quads,
  mass.fall = sum(mass.fall, na.rm = T) / total.quads,
  npp.fall = sum(npp.fall, na.rm = T) / total.quads),
  by = list(year, site, kartez, treatment, genus, sp.epithet, family, LifeHistory, PhotoPath, FunctionalGroup, total.quads)])



site.PFG.total <- as.data.frame(data.table(site.sp.summ)[, list(
  volume = sum(volume, na.rm = T),
  cover = sum(cover, na.rm = T),
  biomass = sum(biomass, na.rm = T)),
  by = list(year, date, site, season, treatment, PFG.PF)])

site.PFG.wide.total <- as.data.frame(data.table(site.sp.wide)[, list(
  mass.spring = sum(mass.spring, na.rm = T),
  mass.fall = sum(mass.fall, na.rm = T),
  npp.fall = sum(npp.fall, na.rm = T)),
  by = list(year, site, treatment, PFG.PF)])

# Calculate mean and sd per species and per PFG



ggplot(site.sp.summ, aes(x = date, y = biomass, group = kartez)) + 
  geom_col(width = 100, aes(fill = kartez), show.legend = F) + 
  facet_grid(site~.)

ggplot(site.PFG.summ, aes(x = date, y = biomass, group = PFG.PF)) + 
  geom_col(width = 110, aes(fill = PFG.PF)) + 
  facet_grid(site~., scales = "free_y")


ggplot(site.PFG.summ, aes(x = date, y = biomass, group = treatment, colour = treatment)) + 
  geom_point() + geom_line() + 
  facet_grid(PFG.PF~site, scales = "free_y")

ggplot(site.PFG.wide[site.PFG.wide$PFG.PF %in% c("C4 grass", "C4 forb", "C3 forb", "C3 shrub", "CAM shrub"),], aes(x = year, y = npp.fall, group = treatment, colour = treatment)) + 
  geom_point() + geom_line() + 
  facet_grid(PFG.PF~site, scales = "free_y") +
  theme_classic(base_size = 18)







#####
#read in previous EDGE data and compare the values
old.edge.data <- read.csv("~/Desktop/EDGE/anpp/compiled_data/all_SEVEDGE_npp.csv")
colnames(old.edge.data)[colnames(old.edge.data) %in% c("swt", "fwt", "species", "path", "a_p", "g_f")] <- c("spring", "fall", "sp.epithet", "PhotoPath", "LifeHistory", "FunctionalGroup")

old.edge.data$site <- revalue(old.edge.data$site, c(
  "black" = "EDGE.black", "blue" = "EDGE.blue"))
old.edge.data$treatment <- revalue(old.edge.data$treatment, c(
  "Control" = "control", "DRT" = "delayed", "ESR" = "drought"))


# gather mass columns
old.edge.data <- gather(old.edge.data[,c("kartez", "site", "plot", "year", "block", "quad", "spring", "fall", "treatment", "family", "genus", "sp.epithet", "PhotoPath", "LifeHistory", "FunctionalGroup")], key = season, value = prev.biomass, spring, fall)


comp.edge <- merge(quad.sp.summ, old.edge.data, all = T)

comp.edge <- as.data.frame(data.table(comp.edge)[, list(
  biomass = sum(biomass, na.rm = T) / 40,
  prev.biomass = sum(prev.biomass, na.rm = T) / 40),
  by = list(year, season, date, site, treatment, kartez, LifeHistory, PhotoPath, FunctionalGroup)])


ggplot(comp.edge, aes(x = biomass, y = prev.biomass, group = kartez, colour = kartez)) + geom_point() + stat_smooth(method = "lm")


summary(comp.edge)
