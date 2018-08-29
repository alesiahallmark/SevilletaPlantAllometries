# This script reads in several datasets of Sevilleta quad measurements and combines them

# Load libraries
library(plyr)
library(knitr)
library(markdown)
library(tidyr)

# File location
all.quad.files <- list.files("~/Desktop/unsorted_rawdata/Sevilleta_quad_datasets", full.names = T)

for(i in 1:length(all.quad.files)) {
  # Read in files
  one.file <- read.csv(all.quad.files[i], strip.white = T)
  
  # add dataset column
  one.file$dataset <- as.factor(strsplit(strsplit(all.quad.files[i], "/")[[1]][length(strsplit(all.quad.files[i], "/")[[1]])], "_")[[1]][1])
  
  # make column names lowercase
  colnames(one.file) <- tolower(colnames(one.file))
  colnames(one.file)[colnames(one.file) == "species"] <- "kartez"
  if(any(one.file$dataset == "sev297")) { #change EDGE columns
    colnames(one.file)[colnames(one.file) == "plot"] <- "subplot"
    colnames(one.file)[colnames(one.file) == "block"] <- "plot"
  }
  
  # format date column
  one.file$month <- tstrsplit(as.character(one.file$date), split = "/")[[1]]
  one.file$day <- tstrsplit(as.character(one.file$date), split = "/")[[2]]
  one.file$date <- as.Date(paste(one.file$month, one.file$day, one.file$year), "%m %d %Y")
  
  # lengthen datasets that have "count" column
  one.file$count[is.na(one.file$count) | one.file$count < 1] <- 0
  one.file <- one.file[rep(row.names(one.file), one.file$count),]
  
  # only keep cool columns
  one.file <- one.file[,colnames(one.file) %in% c(
    "year", "season", "date", "dataset", 
    "site", "treatment", "web", "plot", "subplot", "transect", "quad",
    "kartez", "cover", "height")]
  
  # bind datasets together
  ifelse(i==1, 
         all.files <- one.file,
         all.files <- join(all.files, one.file, type = "full"))
}

all.files.two <- all.files

all.files.two$site <- as.factor(mapvalues(
  as.character(all.files.two$site), 
  from = "FALSE", to = "F", warn_missing = F))


# Make data match allquadNPP file
unique(all.files.two[all.files.two$site == "EB", c("site", "treatment", "web", "plot", "subplot", "quad", "transect")])


all.files.two$web <- as.character(all.files.two$web)
all.files.two$web[all.files.two$site %in% c("P", "EB", "EG")] <- as.character(all.files.two$plot[all.files.two$site %in% c("P", "EB", "EG")])
all.files.two$web <- as.factor(all.files.two$web)

all.files.two$plot <- as.character(all.files.two$plot)
all.files.two$plot[all.files.two$site %in% c("P")] <- as.character(all.files.two$transect[all.files.two$site %in% c("P")])
all.files.two$plot[all.files.two$site %in% c("EB", "EG")] <- as.character(all.files.two$subplot[all.files.two$site %in% c("EB", "EG")])
all.files.two$plot <- as.factor(all.files.two$plot)

all.files.two$subplot <- as.character(all.files.two$subplot)
all.files.two$subplot[all.files.two$site %in% c("B", "P", "G", "C", "F", "W", "M", "BURN09", "BURN03", "MG", "MS", "BG", "CG", "GG", "N", "EB", "EG")] <- all.files.two$quad[all.files.two$site %in% c("B", "P", "G", "C", "F", "W", "M", "BURN09", "BURN03", "MG", "MS", "BG", "CG", "GG", "N", "EB", "EG")]
all.files.two$subplot <- as.factor(all.files.two$subplot)


all.files.two <- all.files.two[,-which(colnames(all.files.two) %in% c("transect", "quad"))]

write.csv(all.files.two, "~alesia/Desktop/unsorted_AJdata/AJ_AllQuads.csv", row.names = F)
