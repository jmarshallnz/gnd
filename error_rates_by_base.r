#' Compute the errors rates and compare them by base

mapping = read.csv("error_mapping.csv")
#' compute difference maps and higlight them
map_source <- mapping %>% left_join(fa15 %>% select(-md5), by=c('sample' = 'serogroup'))
map_dest   <- mapping %>% left_join(fa15 %>% select(-md5), by=c('probable_source' = 'serogroup'))

# now check the difference between them and highlight it...
diff = map_source[,-(1:2)] != map_dest[,-(1:2)]
map_diff = cbind(map_source[,1:2], diff)

map_diff = map_diff %>% arrange(probable_source)
map_matrix = as.matrix(map_diff[,-(1:2)])
# add in an alternating thingee
map_matrix = map_matrix + 2*(as.numeric(map_diff[,2]) %% 2)
colnames(map_matrix) = 1:284
rownames(map_matrix) = map_diff[,1]
pdf("error_maps_by_base_pairs.pdf", width=10, height=20)
par(mar=c(4,6,2,6))
image(1:ncol(map_matrix), 1:nrow(map_matrix), t(map_matrix), col=c("white", "black", "grey80", "black"), xaxt="n", yaxt="n", xlab="", ylab="")
axis(2, 1:nrow(map_matrix), rownames(map_matrix), las=2, cex.axis=0.6)
axis(4, 1:nrow(map_matrix), map_diff[,2], las=2, cex.axis=0.6)
dev.off()

# purity computation
impurity <- function(x) {
  # TODO: Change to shannons entropy perhaps?
  p = prop.table(table(x))
  1 - sum(p^2)
}
# need to compute purity based on those that aren't errors
no_error_serogroups <- fa15 %>% filter(!(serogroup %in% map_diff$sample))
error_rates_by_base <- data.frame(errors = colSums(map_diff[,-c(1:2)]), impurity = apply(no_error_serogroups[,1:284], 2, impurity), triple=c(rep(letters[1:3],284/3),letters[1:2]))
write.csv(error_rates_by_base, "error_rates_by_base.csv", row.names=FALSE)

error_rates_by_base <- read.csv("error_rates_by_base.csv")
plot(jitter(errors) ~ impurity, col=triple, data=error_rates_by_base)

summary(glm(errors ~ impurity + triple, data=error_rates_by_base, family="quasipoisson"))
