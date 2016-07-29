library(readxl)
library(dplyr)

#' Function for reading in the metadata on each sample
read_metadata <- function(file = "data/gnd_seqs_metadata.xlsx", animal_cols) {
  meta = read_excel(file, sheet=2)

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

  samp.meta
}
