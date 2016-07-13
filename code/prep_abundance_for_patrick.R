# prep for patrick
library(readxl)
library(dplyr)

dat <- read.csv("no_error_abundance_with_ctrl.csv")
colSums(dat[,-1])
dat[,-1] <- sweep(dat[,-1], 2, colSums(dat[,-1]), FUN='/')

animal_cols <- names(dat)[-1]

animal <- sub(".*_([0-9ctrl]+)_.*", "\\1", animal_cols)
rep    <- sub("X([0-9]+)_.*", "\\1", animal_cols)

#' Read in metadata
meta = read_excel("gnd_seqs_metadata.xlsx", sheet=2)

samp.meta = data.frame(sample=substring(animal_cols,2), stringsAsFactors = FALSE)

samp.meta = samp.meta %>% left_join(meta, by=c('sample'='Library Name')) %>% dplyr::select(sample, treatment=Treatment, source=`Description [Optional] `)
samp.meta$treatment = factor(samp.meta$treatment)
levels(samp.meta$treatment) = c("bf", "ct")
samp.meta$source = factor(samp.meta$source)
levels(samp.meta$source) = c("fec", "pob", "por", "pre")
samp.meta = samp.meta %>% mutate(animal = sub(".*_([0-9ctrl]+)_.*", "\\1", sample),
                                 label = paste0(animal, "_", source))

# relabel the ctrls
wch <- samp.meta$animal == "ctrl"
ctrl_labs = sub("([0-9ctrl]+)_.*", "\\1", samp.meta$sample)
samp.meta$label[wch] <- paste0("ctrl_", ctrl_labs[wch])

#' TODO: Update the column names...
names(dat)[-1] <- samp.meta$label

write.csv(dat, "no_error_proportions_with_ctrl_for_patrick.csv", row.names=FALSE)

#' 
#' Patrick needs a nexus file...