## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.width = 6, fig.height = 6)
# Packages --------------------------------------------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    library("bioregion")
    library("dplyr")
    library("ggplot2")
    library("sf")
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
vege_nhclu <- nhclu_kmeans(vegedissim, 
                           n_clust = 3, 
                           index = "Simpson",
                           seed = 1)
vege_nhclu$cluster_info 

# Hierarchical bioregionalization
set.seed(1)
vege_hclu <- hclu_hierarclust(dissimilarity = vegedissim,
                              index = "Simpson",
                              method = "average", 
                              n_clust = 3,
                              optimal_tree_method = "best",
                              verbose = FALSE)
vege_hclu$cluster_info

# Network bioregionalization
set.seed(1)
vege_netclu <- netclu_walktrap(vegesim,
                               index = "Simpson")
vege_netclu$cluster_info 

# Bipartite network bioregionalization
install_binaries(verbose = FALSE)
vege_netclubip <- netclu_infomap(vegedf,
                                 seed = 1, 
                                 bipartite = TRUE)
vege_netclubip$cluster_info


## -----------------------------------------------------------------------------
all_metrics <- site_species_metrics(
  bioregionalization = vege_netclubip,
  bioregion_metrics = c("Specificity", "NSpecificity", "Fidelity", 
                        "IndVal", "NIndVal", "Rho", "CoreTerms",
                        "Richness", "Rich_Endemics", "Prop_Endemics",
                        "MeanSim", "SdSim"), # You can also simply write "all"
  bioregionalization_metrics = c("P", "Silhouette"),
  data_type = "both",
  cluster_on = "both",
  comat = vegemat,
  similarity = vegesim,
  index = "Simpson",
  verbose = FALSE)

## -----------------------------------------------------------------------------
all_metrics

## -----------------------------------------------------------------------------
summary(all_metrics)

## -----------------------------------------------------------------------------
str(all_metrics)

## -----------------------------------------------------------------------------
nsb <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("Specificity", "NSpecificity",
                                                  "Fidelity", "IndVal", "NIndVal",
                                                  "Rho", 
                                                  "CoreTerms"),
                            bioregionalization_metrics = NULL,
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = NULL, # Name of similarity column
                            verbose = FALSE)

nsb


## -----------------------------------------------------------------------------
wsb <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("Specificity", "NSpecificity",
                                                  "Fidelity",
                                                  "IndVal", "NIndVal",
                                                  "Rho",
                                                  "CoreTerms"),
                            bioregionalization_metrics = NULL,
                            data_type = "abundance",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL, # Name of similarity column
                            index = NULL,
                            verbose = FALSE)

wsb


## -----------------------------------------------------------------------------
sim_metrics <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("Richness", "Rich_Endemics",
                                                  "Prop_Endemics"),
                            bioregionalization_metrics = NULL,
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = vegesim,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

sim_metrics

## -----------------------------------------------------------------------------
sim_metrics <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = c("MeanSim", "SdSim"),
                            bioregionalization_metrics = NULL,
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = vegesim,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

sim_metrics

## -----------------------------------------------------------------------------
gc <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = c("Specificity", "NSpecificity",
                                                  "Fidelity",
                                                  "IndVal", "NIndVal",
                                                  "Rho",
                                                  "CoreTerms"),
                            bioregionalization_metrics = "P",
                            data_type = "both",
                            cluster_on = "species",
                            comat = vegemat,
                            similarity = NULL,
                            index = NULL,
                            verbose = FALSE)

gc


## -----------------------------------------------------------------------------
sil_metrics <- site_species_metrics(bioregionalization = vege_nhclu,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "Silhouette",
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = vegesim,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

sil_metrics

## -----------------------------------------------------------------------------
p_occ_site <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "occurrence",
                            cluster_on = "species",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_occ_site

## -----------------------------------------------------------------------------
p_ab_site <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "abundance",
                            cluster_on = "species",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_ab_site

## -----------------------------------------------------------------------------
p_occ_sp <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "occurrence",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_occ_sp

## -----------------------------------------------------------------------------
p_ab_sp <- site_species_metrics(bioregionalization = vege_netclubip,
                            bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "abundance",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = "Simpson", # Name of similarity column
                            verbose = FALSE)

p_ab_sp

## -----------------------------------------------------------------------------
ps <- site_species_metrics(bioregionalization = vege_nhclu,
                           bioregion_metrics = NULL,
                            bioregionalization_metrics = "P",
                            data_type = "both",
                            cluster_on = "site",
                            comat = vegemat,
                            similarity = NULL,
                            index = NULL,
                            verbose = FALSE)

ps


## -----------------------------------------------------------------------------
bioregion_summary <- bioregion_metrics(bioregionalization = vege_nhclu,
                                       comat = vegemat)
bioregion_summary

## -----------------------------------------------------------------------------
# Spatial coherence
vegedissim <- dissimilarity(vegemat)
hclu <- nhclu_kmeans(dissimilarity = vegedissim, n_clust = 4)
vegemap <- map_bioregions(hclu, vegesf, write_clusters = TRUE, plot = FALSE)

bioregion_metrics(bioregionalization = hclu, comat = vegemat, map = vegemap,
col_bioregion = 2) 

## -----------------------------------------------------------------------------
ggplot(vegemap) +
  geom_sf(aes(fill = as.factor(K_4))) +
  scale_fill_viridis_d("Bioregion") +
  theme_bw() +
  theme(legend.position = "bottom")

