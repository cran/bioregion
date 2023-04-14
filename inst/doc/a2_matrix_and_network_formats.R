## ----setup, include=FALSE-----------------------------------------------------

knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE,
                      fig.width = 8, fig.height = 8)
# Packages --------------------------------------------------------------------
suppressPackageStartupMessages({
  suppressWarnings({
    library(bioregion)
  })
})

options(tinytex.verbose = TRUE)


## ----loading_vegedf-----------------------------------------------------------
data(vegedf)
head(vegedf)
dim(vegedf)
sum(!duplicated(vegedf[,1]))
sum(!duplicated(vegedf[,2]))

## ----loading_vegemat----------------------------------------------------------
data(vegemat)
vegemat[1:10,1:10]
dim(vegemat)

## ----mat_to_net_1-------------------------------------------------------------
net <- mat_to_net(vegemat, weight = TRUE, remove_zeroes = FALSE)

## ----print_net_1--------------------------------------------------------------
head(net)
dim(net)

## ----mat_to_net_2-------------------------------------------------------------
net <- mat_to_net(vegemat, weight = TRUE, remove_zeroes = TRUE)

## ----print_net_2--------------------------------------------------------------
head(net)
dim(net)

## ----net_to_mat_1-------------------------------------------------------------
mat <- net_to_mat(vegedf, weight = TRUE, squared = FALSE, symmetrical = FALSE, missing_value = 0)

## ----print_mat_1--------------------------------------------------------------
mat[1:5,1:5]
dim(mat)

## ----net_to_mat_2-------------------------------------------------------------
mat <- net_to_mat(vegedf, weight = TRUE, squared = TRUE, symmetrical = FALSE, missing_value = 0)

## ----print_mat_2--------------------------------------------------------------
mat[1:5,1:5]
dim(mat)

## ----net_to_mat_3-------------------------------------------------------------
temp <- data.frame(Site=c("35","36","36","38","39"), Species=c("36","35","37","37","39"), Abundance=c(1,2,3,4,0))
net <- rbind(temp,vegedf)
mat <- net_to_mat(net, weight = TRUE, squared = TRUE, symmetrical = FALSE, missing_value = -1)

## ----print_mat_3--------------------------------------------------------------
mat[1:5,1:5]

## ----net_to_mat_4-------------------------------------------------------------
mat <- net_to_mat(net, weight = TRUE, squared = TRUE, symmetrical = TRUE, missing_value = 0)

## ----print_mat_4--------------------------------------------------------------
mat[1:5,1:5]

