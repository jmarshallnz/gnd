library(seqinr)
library(dplyr)
library(digest)

read_fasta_files <- function(md5s, fasta_paths, read_length = 284) {

  read_fasta_file <- function(file) {
    fa  <- read.fasta(file)
    seq <- lapply(fa, function(x) { toupper(paste(x,collapse='')) })
    md5 <- unlist(lapply(seq, function(x) { digest(x, serialize=FALSE)}))

    # get rid of duplicates
    names(fa) <- md5
    fa <- fa[!duplicated(md5)]

    fa
  }

  if (!is.list(fasta_paths)) {
    fasta_paths = as.list(fasta_paths)
  }
  fa_list <- lapply(fasta_paths, read_fasta_file)
  fa_all <- do.call(c, fa_list)

  # check we have sequences for all the MD5's requested
  if (sum(!md5s %in% names(fa_all)) != 0) {
    print(which(!md5s %in% names(fa_all)))
    stop("Error: We have MD5's that we don't have sequences for!")
  }

  # check the lengths
  wch <- which(lengths(fa_all) != read_length)
  if (length(wch) > 0) {
    warning(paste("There seems to be", length(wch), "sequences in the fasta file that are the wrong length (not", read_length, "). These will be removed."))
    fa_all <- fa_all[-wch]
  }

  # reshape into a data.frame
  fa <- data.frame(do.call(rbind, fa_all), stringsAsFactors = FALSE) %>%
    tibble::rownames_to_column("md5")

  # return object for further processing
  fa %>% filter(md5 %in% md5s)
}

read_fasta_paths <- function(path_meta, path_fasta, path_gnd_db, read_length=284) {

  # Read in the metadata (it has the abundances) and remove the non-functional groups
  fa_meta = read.table(path_meta, header=TRUE, sep="\t", stringsAsFactors = FALSE) %>%
    filter(functional != 'no')

  md5s = fa_meta$md5

  # TODO:  Update this so that a list of fasta files are sent in, and the ID used for mapping is always the MD5.
  # that way the label in the fasta file doesn't matter so much (or at all) - we can then take the identifier
  # from the abundance file if we want to do so later. We do need the metadata (abundance file) here because
  # in theory the fasta file may contain stuff that the abundance file doesn't. We could just pass in a set
  # of MD5's to match on though rather than the file - i.e. split loading of that out
  
  fasta_files = list(path_fasta, path_gnd_db)

  fa = read_fasta_files(md5s, fasta_files, read_length)

  final <- fa_meta %>% select(md5, serogroup) %>% inner_join(fa)

  final
}

read_fasta <- function(farm=1, qa=15, minTotal=10) {

  # Setup paths  
  base_dir <- "data/gnd2"
  farm_dir <- paste0("farm", farm)
  file_meta <- paste0("solexaQA",qa,"thou_minTotalGE", minTotal, ".txt")
  path_meta <- file.path(base_dir, farm_dir, file_meta)

  file_fa <- paste0("solexaQA",qa,"thou_nucleotideGE", minTotal, ".fa")
  path_fasta <- file.path(base_dir, farm_dir, file_fa)

  path_gnd_db <- file.path(base_dir, "gnd_DB_09012017.fas")

  read_fasta_paths(path_meta, path_fasta, path_gnd_db, read_length = 284)

}
