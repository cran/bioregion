## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.width = 6, fig.height = 6)
# Packages --------------------------------------------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    library("bioregion")
    library("dplyr")
  })
})

options(tinytex.verbose = TRUE)

## -----------------------------------------------------------------------------
data("vegedf")
data("vegemat")

# Calculation of (dis)similarity matrices
vegedissim <- dissimilarity(vegemat, metric = c("Simpson"))
vegesim <- dissimilarity_to_similarity(vegedissim)

## -----------------------------------------------------------------------------
# Non hierarchical bioregionalization
vege_nhclu_kmeans <- nhclu_kmeans(vegedissim, n_clust = 3, index = "Simpson")
vege_nhclu_kmeans$cluster_info # 3

# Hierarchical bioregionalization
set.seed(1)
vege_hclu_hierarclust <- hclu_hierarclust(dissimilarity = vegedissim,
                                          method = "mcquitty", n_clust = 3,
                                          optimal_tree_method = "best")
vege_hclu_hierarclust$cluster_info # 3

# Network bioregionalization
set.seed(1)
vege_netclu_walktrap <- netclu_walktrap(vegesim,
                                        index = names(vegesim)[3])
vege_netclu_walktrap$cluster_info # 3

## -----------------------------------------------------------------------------
comp <- dplyr::left_join(vege_hclu_hierarclust$clusters,
                         vege_netclu_walktrap$clusters,
                         by = "ID")
colnames(comp) <- c("ID", "K_3_hclu", "K_3_netclu")
comp <- dplyr::left_join(comp,
                         vege_nhclu_kmeans$clusters,
                         by = "ID")
colnames(comp) <- c("ID", "K_3_hclu", "K_3_netclu", "K_3_nhclu")

head(comp)

## -----------------------------------------------------------------------------
hclu_vs_netclu <- compare_bioregionalizations(
  bioregionalizations = comp[, c("K_3_hclu", "K_3_netclu", "K_3_nhclu")],
  store_pairwise_membership = TRUE,
  cor_frequency = TRUE,
  store_confusion_matrix = TRUE)

str(hclu_vs_netclu)

## -----------------------------------------------------------------------------
nrow(hclu_vs_netclu$pairwise_membership) == nrow(comp)*(nrow(comp)-1)/2

## -----------------------------------------------------------------------------
comp[c(1, 9), ]

## -----------------------------------------------------------------------------
hclu_vs_netclu$pairwise_membership[8:10, ]

## -----------------------------------------------------------------------------
hclu_vs_netclu$freq_item_pw_membership[c(1, 8)]

## -----------------------------------------------------------------------------
table(hclu_vs_netclu$freq_item_pw_membership)

## -----------------------------------------------------------------------------
hclu_vs_netclu$confusion_matrix

## -----------------------------------------------------------------------------
hclu_vs_netclu$bioregionalization_freq_cor

