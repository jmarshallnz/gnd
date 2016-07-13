# prep for patrick
library(dplyr)
source("code/read_metadata.R")

dat <- read.csv("no_error_abundance_with_ctrl.csv")
colSums(dat[,-1])
dat[,-1] <- sweep(dat[,-1], 2, colSums(dat[,-1]), FUN='/')

animal_cols <- names(dat)[-1]

animal <- sub(".*_([0-9ctrl]+)_.*", "\\1", animal_cols)
rep    <- sub("X([0-9]+)_.*", "\\1", animal_cols)

#' Read in metadata and update column names
names(dat)[-1] <- read_metadata(animal_cols=animal_cols)$label

write.csv(dat, "no_error_proportions_with_ctrl_for_patrick.csv", row.names=FALSE)

#' 
#' Patrick needs a nexus file...