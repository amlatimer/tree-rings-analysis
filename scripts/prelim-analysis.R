## This file loads and explores preliminary tree-ring data from Derek Young. 

## Load libraries and data
library(ggplot2)
library(sp)
library(lme4)
library(maps)

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

# fix problem with PSME id
levels(trees$species)
trees$species[trees$species==""] <- NA
trees$species[trees$species=="PSME "] <- "PSME"
trees$species <- droplevels(trees$species)
levels(trees$species)


# Add cluster ID to trees data set
trees <- merge(trees, plots[, c("plot.id", "cluster")], sort=FALSE, by = "plot.id", all=FALSE)
# Add cluster and plot ID and species to years data set
years <- merge(years, trees[, c("tree.id", "plot.id", "cluster", "species", "voronoi.area")], sort=FALSE, by = "tree.id", all=FALSE)

# how many trees have BAI data for each year?
with(years, table(bai.pres=!is.na(bai.comb), year)) # most years have values for most trees
with(years, table(bai.pres=!is.na(bai), year)) # about 2/3 of trees have "pith-based" BA for any given year. 
table(years$tree.id[!is.na(years$ba)])
table(years$cluster[!is.na(years$ba)])
table(years$cluster[!is.na(years$ba.comb)])

# Where are the plots?
treelocs <- SpatialPoints(coords=trees[,c("x", "y")], proj4string = CRS("+proj=albers"))
plot(treelocs, pch=16, col=trees$species)

# visulize some individual tree growth trends for trees with basal area measured from inside out, and from outside in. 
years_psme <- years[years$species=="PSME",]
unique(plots$plot.id)
p <- ggplot(years[years$plot.id=="ULC1" & years$species=="PSME",], aes(year, ba.comb)) + geom_line(aes(color = !is.na(ba))) + facet_wrap(~tree.id)
p
# No obvious difference between trees' basal area that's measured from core out (ba, bai) and basal area that's measured from outside in (ba.ext, bai.ext). Therefore, proceed with using the combined data for now (ba.comb, bai.comb). 


