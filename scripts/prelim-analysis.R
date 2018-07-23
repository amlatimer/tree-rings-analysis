## This file loads and explores preliminary tree-ring data from Derek Young. 

## Load libraries and data
library(ggplot2)
library(sp)
library(lme4)

years <- read.csv("../data/years.csv")
trees <- read.csv("../data/trees.csv")
plots <- read.csv("../data/plots.csv")

## Initial look at the data

head(years)
head(trees)
length(unique(trees$tree.id)) # how many trees? (774)
min(years$year); max(years$year) # how many years (55)
length(unique(trees$plot.id)) # how many plots? (63)
table(trees$species) # how many trees of each species?

# Where are the plots?
treelocs <- SpatialPoints(coords=trees[,c("x", "y")], proj4string = CRS("+proj=albers"))
map()
plot(treelocs, pch=16, col=trees$species)

# Add cluster ID to trees data set 
trees <- merge(trees, plots[, c("plot.id", "cluster")], sort=FALSE, by = "plot.id", all=FALSE)

