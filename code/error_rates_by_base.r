#' Compute the errors rates and compare them by base
source("code/read_fasta.R")

mapping = read.csv("error_mapping.csv")

fa15 = read_fasta()

#' compute difference maps and higlight them
map_source <- mapping %>% left_join(fa15 %>% select(-md5), by=c('sample' = 'serogroup'))
map_dest   <- mapping %>% left_join(fa15 %>% select(-md5), by=c('probable_source' = 'serogroup'))

# now check the difference between them and highlight it...
map_table <- matrix(0, 4, 4)
rownames(map_table) <- colnames(map_table) <- c("a","c","g","t")
map_table[1,-1] <- 1:3
map_table[2,-2] <- 4:6
map_table[3,-3] <- 7:9
map_table[4,-4] <- 10:12

diff <- matrix(NA, nrow(map_source), ncol(map_source)-2)
diff <- cbind(map_source[,1:2], diff)
for (i in 1:nrow(map_source)) {
  for (j in 3:ncol(map_source)) {
    diff[i,j] = map_table[map_dest[i,j], map_source[i,j]]
  }
}

map_diff = diff %>% arrange(probable_source)
map_matrix = as.matrix(map_diff[,-(1:2)])

colnames(map_matrix) = 1:284
rownames(map_matrix) = map_diff[,1]

# add in alternating background
for (i in 1:nrow(map_matrix)) {
  map_matrix[i, map_matrix[i,] == 0] <- 1-(as.numeric(map_diff[,2]) %% 2)[i] + 13
}

col=c("#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", 
      "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", 
      "grey70", "grey90")

actual_used <- unique(as.numeric(map_matrix))

pdf("figures/error_maps_by_base_pairs.pdf", width=10, height=20)
par(mar=c(4,6,2,6))
image(1:ncol(map_matrix), 1:nrow(map_matrix), t(map_matrix), col=col, xaxt="n", yaxt="n", xlab="", ylab="")
axis(2, 1:nrow(map_matrix), rownames(map_matrix), las=2, cex.axis=0.6)
axis(4, 1:nrow(map_matrix), map_diff[,2], las=2, cex.axis=0.6)

# do the legend
library(tidyr)
map_legend <- map_table %>% as.data.frame %>%
  mutate(source = rownames(map_table)) %>%
  gather("error", "colour", -source) %>%
  mutate(type = paste(source, "->", error)) %>%
  filter(colour %in% actual_used)
#legend(ncol(map_matrix)+40, nrow(map_matrix), legend=map_legend$type, fill=col[map_legend$col], xpd=TRUE, cex=0.6)
legend("topleft", legend=map_legend$type, fill=col[map_legend$col], xpd=TRUE, cex=0.6)
dev.off()

# purity computation
impurity <- function(x) {
  # TODO: Change to shannons entropy perhaps?
  p = prop.table(table(x))
  1 - sum(p^2)
}
# need to compute purity based on those that aren't errors
no_error_serogroups <- fa15 %>% filter(!(serogroup %in% map_diff$sample))
error_rates_by_base <- data.frame(errors = colSums(map_diff[,-c(1:2)] > 0), impurity = apply(no_error_serogroups[,1:284], 2, impurity), triple=c(rep(letters[c(2,3,1)],284/3),letters[1:2]))

# error rates by codon below here
error_rates_per_codon <- error_rates_by_base %>% group_by(triple) %>% summarise(errors = sum(errors)) %>% as.data.frame

pdf("figures/errors_by_codon.pdf", width=6, height=6)
barplot(error_rates_per_codon$errors, names=1:3, xlab = "Codon", ylab="Number of SNPs", ylim=c(0,120))
dev.off()

summary(glm(errors ~ impurity + triple, data=error_rates_by_base, family="quasipoisson"))
