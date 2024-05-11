library(ROCR)


args <- commandArgs(trailingOnly = TRUE)

df_surf_bp1 <- read.csv(args[1])
df_surf_bp2 <- read.csv(args[2])

dist_surf1 <- as.matrix(dist(df_surf_bp1[, c("x", "y", "z")]))
mean_list1 <- lapply(
  seq_len(nrow(df_surf_bp1)),
  \(i) mean(df_surf_bp1$color[dist_surf1[i, ] <= 6], na.rm = TRUE)
)
rm(dist_surf1)

dist_surf2 <- as.matrix(dist(df_surf_bp2[, c("x", "y", "z")]))
mean_list2 <- lapply(
  seq_len(nrow(df_surf_bp2)),
  \(i) mean(df_surf_bp2$color[dist_surf2[i, ] <= 6], na.rm = TRUE)
)
rm(dist_surf2)

get_auc_rocr <- function(means, bp) {
  pred <- prediction(means, bp)
  auc_rocr <- performance(pred, measure = "auc")
  return(1 - auc_rocr@y.values[[1]])
}

auc_rocr1 <- get_auc_rocr(unlist(mean_list1), df_surf_bp1$bp)
auc_rocr2 <- get_auc_rocr(unlist(mean_list2), df_surf_bp2$bp)

print(paste("The AUC of the first protein is:", round(auc_rocr1, 2)))
print(paste("The AUC of the second protein is:", round(auc_rocr2, 2)))
