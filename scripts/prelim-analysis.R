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

# Create relative growth rate (RGR) and log-RGR response variables
years$rgr.comb <- years$ba.comb / years$ba.prev.comb
hist(years$rgr.comb)
sum(is.na(years$rgr.comb)) / nrow(years)
sum(!is.na(years$rgr.comb))
sum(years$rgr.comb > 1.1, na.rm=T)
plot(rgr.comb~ba.comb, years)
plot(rgr.comb~ba.comb, years[years$rgr.comb<=1.1,]) # RGR is too heavily influenced by tree size at the smaller end to use
plot(bai.comb~ba.comb, years[years$rgr.comb<=1.1,]) 

  
# Standardize explanatory variables
vars_to_scale <- c("ba.comb", "ba.prev.comb", "ppt", "ppt.z", "tmean", "tmean.z")
for (i in 1:length(vars_to_scale)) years[,vars_to_scale[i]] <- scale(years[,vars_to_scale[i]])

# how many trees have BAI data for each year?
with(years, table(bai.pres=!is.na(bai.comb), year)) # most years have values for most trees
with(years, table(bai.pres=!is.na(bai), year)) # about 2/3 of trees have "pith-based" BA for any given year. 
table(years$tree.id[!is.na(years$ba)])
table(years$cluster[!is.na(years$ba)])
table(years$cluster[!is.na(years$ba.comb)])

# Where are the plots?
treelocs <- SpatialPoints(coords=trees[,c("x", "y")], proj4string = CRS("+proj=albers"))
plot(treelocs, pch=16, col=trees$species)

# Visulize some individual tree growth trends for trees with basal area measured from inside out, and from outside in. 
years_psme <- years[years$species=="PSME",]
unique(plots$plot.id)
p <- ggplot(years[years$plot.id=="ULC1" & years$species=="PSME",], aes(year, ba.comb)) + geom_line(aes(color = !is.na(ba))) + facet_wrap(~tree.id)
p
# No obvious difference between trees' basal area that's measured from core out (ba, bai) and basal area that's measured from outside in (ba.ext, bai.ext). Therefore, proceed with using the combined data for now (ba.comb, bai.comb). 


#### Q1: What are relationships between long-term precipitation and growth rate, interannual precipitation variation and growth rate, accounting for  temperature and tree size? 

# Specify which subset of data to use for this analysis 
d <- years[which(years$species=="PSME"),] # Douglas fir only



### Q2: For each tree, what is its sensitivity to variation in precipiation? 

# Specify which subset of data to use for this analysis 
d <- years[which(years$species=="PSME") ,]

# What is the history of precipitation variation in the region? 
ggplot(d, aes(year, ppt.z)) + geom_line(aes(col=plot.id)) + facet_wrap(~cluster) 
# Pretty consistent pattern across whole data set. Very consistent within clusters. Extreme dry years were 1976 and 1977. 2013-2014 and 1987-88 were next biggest dry spikes, with below-average precip extending through 1992. Weaker dry spell in 2007-2008. 
# Wettest years were 1982-83, followed by 2006 and 1969, 1995, 1998, 2006, 2011. 
ggplot(d, aes(year, ppt.z)) + geom_hist(aes(col=plot.id)) + facet_wrap(~cluster)


# Visual check: plot bai (growth) versus precipitation anomaly for each site for individual trees. 
# by plot
ggplot(d, aes(ppt.z)) + geom_histogram(bins=30) + facet_wrap(~cluster)
ggplot(d, aes(ppt)) + geom_histogram(bins=30) + facet_wrap(~cluster)

# by trees within a plot
# note for comparing response of individual tree to higher vs lower than average precip RWI makes most sense because it spreads out the values relative to the tree's growth rate variation. Sensitivity calculate on basis of RWI could then be compared to absolute mean growth rate. 
focal_plot <- "HARDIN B"
#ggplot(d[which(d$plot.id==focal_plot),], aes(ppt.z, bai.comb)) + geom_point() + facet_wrap(~tree.id)
ggplot(d[which(d$plot.id==focal_plot),], aes(ppt.z, rwi)) + geom_point() + facet_wrap(~tree.id)
ggplot(d[which(d$plot.id==focal_plot),], aes(ppt.z, bai.comb)) + geom_point() + facet_wrap(~tree.id)

# Does sensitivity change over time for plots as a whole? 
ggplot(d, aes(ppt.z, rwi)) + geom_point(aes(col = factor(year>1990))) + facet_wrap(~plot.id)
ggplot(d, aes(ppt.z, bai.comb)) + geom_point(aes(col = factor(year>1990))) + facet_wrap(~plot.id)
# not obviously, but there are some differences in mean and possibly sensitivity.




