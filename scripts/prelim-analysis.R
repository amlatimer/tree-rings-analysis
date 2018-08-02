#### Preliminary analysis of tree-ring data from Derek Young's PSME project. ####
## This script: 
# - can be run on filtered vs unfiltered data to check sensitivity of results to filtering.
# - loads, checks, and merges data 
# - creates metrics of individual tree-level growth and sensitivity of growth to precipitation
# - checks for overall trends in growth and sensitivity versus long-term site characteristics (average precip, temp, solar radiation)
# - checks for associations between sensitivity and growth rate, and between these and long-term site characteristics


## Load libraries and data ####
library(ggplot2)
library(sp)
library(lme4)
library(maps)
library(zoo)
library(MuMIn)

# Load filtered data, which removes "problematic" chronologies and segments that don't align well with the reference chronology without adding/deleting rings, or not (leaves those trees in). 
load_filtered <- FALSE

if (load_filtered) {
  years <- read.csv("../data/filtered/years.csv")
  trees <- read.csv("../data/filtered/trees.csv")
  plots <- read.csv("../data/filtered/plots.csv")
} else { 
  years <- read.csv("../data/unfiltered/years.csv")
  trees <- read.csv("../data/unfiltered/trees.csv")
  plots <- read.csv("../data/unfiltered/plots.csv")
}


## Clean up and merge the data sets ####

# Initial data check
head(years)
head(trees)
length(unique(trees$tree.id)) # how many trees? (776)
min(years$year); max(years$year) # how many years (55)
length(unique(trees$plot.id)) # how many plots? (63)
table(trees$species) # how many trees of each species?

# Fix problem with PSME id
levels(trees$species)
trees$species[trees$species==""] <- NA
trees$species[trees$species=="PSME "] <- "PSME"
trees$species <- droplevels(trees$species)
levels(trees$species)

# Add plot-level data to years data set 
years <- merge(years, plots[, c( "plot.id", "cluster", "ppt.norm", "tmean.norm", "rad.tot", "rad.03", "rad.06")])
head(years)

# Add measured dbh (size of tree when the plot was visited) to the years data set
years <- merge(years, trees[,c("tree.id", "dbh", "species")], by = "tree.id")

# Add cluster ID to trees data set
trees <- merge(trees, plots[, c("plot.id", "cluster")], sort=FALSE, by = "plot.id")

# Check how many trees have BAI data for each year
with(years, table(bai.pres=!is.na(bai.comb), year)) # most years have values for most trees
with(years, table(bai.pres=!is.na(bai), year)) # about 2/3 of trees have "pith-based" BA for any given year. 
table(years$tree.id[!is.na(years$ba)])
table(years$cluster[!is.na(years$ba)])
table(years$cluster[!is.na(years$ba.comb)])

# Display plots on a map
treelocs <- SpatialPoints(coords=trees[,c("x", "y")], proj4string = CRS("+proj=albers"))
plot(treelocs, pch=16, col=trees$species)

# Display growth rates of some trees vs time 
#ggplot(years[years$plot.id=="SSC5B",], aes(x=year, y=rwi)) + geom_line() + facet_wrap(~tree.id) # something wrong with this! 


#### Q1: What are relationships between long-term precipitation and growth rate, interannual precipitation variation and growth rate, accounting for  temperature and tree size? ####

# Specify which subset of data to use for this analysis 
d <- years[which(years$species=="PSME"),] # focus on Douglas fir only

# Standardize explanatory variables
vars_to_scale <- c("bai.comb", "ba.comb", "ba.prev.comb", "ppt.z", "tmean.z", "rad.tot", "ppt.norm", "tmean.norm", "dbh") 
for (i in 1:length(vars_to_scale)) d[,vars_to_scale[i]] <- scale(d[,vars_to_scale[i]])

# Fit model
plot(bai.comb~ppt.norm, d[d$year==1990,])
m1_q1 <- lm(bai.comb~ppt.norm + tmean.norm + ppt.z + tmean.z + rad.tot, data=d)
summary(m1_q1)


#### CHECKING POTENTIAL DATA PROBLEM ####
# Weirdly appears no variation among trees that can be explained by long-term climate. 
# Does this persist if we average up to the tree level? 
testdata_early <- aggregate(years$bai.comb[years$year >= 1990], by=list(years$tree.id[years$year >= 1990]), FUN=mean, na.rm=TRUE)
testdata_late <- aggregate(years$bai.comb[years$year < 1990], by=list(years$tree.id[years$year < 1990]), FUN=mean, na.rm=TRUE)
names(testdata_early) <- names(testdata_late) <- c("tree.id", "bai.comb.mean")
testdata_early <- merge(testdata_early, trees[,c("tree.id", "plot.id", "dbh", "species")], by = "tree.id")
testdata_late <- merge(testdata_late, trees[,c("tree.id", "plot.id", "dbh", "species")], by = "tree.id")
testdata_early <- merge(testdata_early, plots, by = "plot.id")
testdata_late <- merge(testdata_late, plots, by = "plot.id")
head(testdata)
cor(testdata)
summary(lm(bai.comb.mean~ppt.norm * rad.tot + dbh + species , testdata_early))
summary(lm(bai.comb.mean~ppt.norm * rad.tot + dbh + species , testdata_late))
summary(lm(bai.comb.mean~ppt.norm+tmean.norm+rad.tot , testdata))
plot(bai.comb.mean~ppt.norm , testdata)
plot(voronoi.area.mean ~ ppt.norm, testdata)
# No -- averaged up to tree level, we do see effect of overall precip on growth. 

# What if we go back to the whole data set and account for some of the variation using random effects? 
m <- lmer(bai.comb ~ ba.prev.comb + ppt.norm + rad.tot + ppt.z + (1|plot.id/tree.id), data = d)
summary(m)
r.squaredGLMM(m)
#  Compare to using RWI
m <- lm(rwi ~ tmean.z * ppt.z + ba.prev.comb*ppt.z, data = d)
summary(m)
r.squaredGLMM(m)

#### Q2: For each tree, what is its sensitivity to variation in precipiation? #### 

# Specify which subset of data to use for this analysis 
d <- years[which(years$species=="PSME") ,]

# What is the history of precipitation variation in the region? 
ggplot(d, aes(year, ppt.z)) + geom_line(aes(col=plot.id)) + facet_wrap(~cluster) 
# Pretty consistent pattern across whole data set. Very consistent within clusters. Extreme dry years were 1976 and 1977. 2013-2014 and 1987-88 were next biggest dry spikes, with below-average precip extending through 1992. Weaker dry spell in 2007-2008. 
# Wettest years were 1982-83, followed by 2006 and 1969, 1995, 1998, 2006, 2011. 
ggplot(d, aes(year, ppt.z)) + geom_histogram(aes(col=plot.id)) + facet_wrap(~cluster)
ggplot(d, aes(year, ppt)) + geom_histogram(aes(col=plot.id)) + facet_wrap(~cluster)
# There are distinct clumps of "high" vs "low" years, i.e. precipiation is somewhat bimodal. 
# When looking at absolute ppt values, most of the years in Plumas and Tahoe are wet, while most years in Yose and Sierra are dry. Clear N-S gradient in overal precipiation. 
# No obvious trend in variability over time 
cv <- function(x) {return(sd(x, na.rm=T)/mean(x, na.rm=T))}
plot(1964:2009,rollapply(years[which(d$tree.id=="1089"),"ppt"], 10, cv), col="blue", lwd=3, type="l", ylim=c(0.2,0.6), ylab="CV precipitation", xlab="year") # plumas
lines(1964:2009, rollapply(years[which(d$tree.id=="1411"),"ppt"], 10, cv),  col="yellow3",lwd=3) # sierra
#lines(1964:2009, rollapply(years[which(d$tree.id=="1172"),"ppt"], 10, cv), col="cyan3",lwd=3) # tahoe
#lines(1964:2009, rollapply(years[which(d$tree.id=="1427"),"ppt"], 10, cv), col="orange3",lwd=3) # yose
legend("topright", c("Sierra (S)", "Plumas (N)"), col=c("yellow3", "blue"), lwd=c(3, 3))
# Overall CV precipitation is higher in the sourh, and it stays high throughout. Whereas in the north (plumas), CV is low at the beginning and the end of the time series, and spikes during the 1970s drought and early 90s. 


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

## Calculate sensitivity metrics for each tree for whole time period ####
length(unique(years$tree.id)) == nrow(trees)
n_trees <- nrow(trees)
drought_years <- c(1976, 1977, 2013, 2014, 1987, 1988, 2013, 2014)
deluge_years <- c(1982, 1983, 2006, 1969, 1995, 1998, 2006, 2011)
sens_all <- sens_all_bai <- sens_all_early <- sens_all_late <- sens_wet <- sens_dry <- wet_dry_diff <- mean_bai <- sd_bai <- mean_ba <- rep(0, n_trees)

# loop over trees to calculate sensitivity and associated metrics
for (i in 1:n_trees) { 
  tempdata <- years[years$tree.id==trees$tree.id[i],]
  wetyears <- which(tempdata$ppt.z >= 0)
  sens_all[i] <- coef(lm(rwi~ppt.z, data=tempdata))[2]
  #sens_all_bai[i] <- coef(lm(bai.comb~ba.prev.comb + ppt.z, data=tempdata))[3]
  #sens_all_early[i] <- coef(lm(rwi~ppt.z, data=tempdata, subset=tempdata$year<1990))[2]
  #sens_all_late[i] <- coef(lm(rwi~ppt.z, data=tempdata, subset=tempdata$year>=1990))[2]
  sens_wet[i] <- coef(lm(rwi~ppt.z, data=tempdata[wetyears,]))[2]
  sens_dry[i] <- coef(lm(rwi~ppt.z, data=tempdata[-wetyears,]))[2]
  wet_dry_diff[i] <- mean(tempdata$rwi[tempdata$year %in% deluge_years], na.rm=T) - mean(tempdata$rwi[tempdata$year %in% drought_years], na.rm=T)
  mean_bai[i] <- mean(tempdata$bai.comb, na.rm=T)
  sd_bai[i] <- sd(tempdata$bai.comb, na.rm=T)
  mean_ba[i] <- mean(tempdata$ba.comb, na.rm=T)
}

sensdata <- data.frame(tree.id = trees$tree.id, sens_all, sens_wet, sens_dry, wet_dry_diff, mean_bai, sd_bai, mean_ba)
sensdata <- merge(sensdata, trees, by="tree.id")
sensdata <- merge(sensdata, plots, by="plot.id")

write.csv(sensdata, "../working-data/sensitivity-metrics.csv")

hist(sens_all)
plot(sens_dry~sens_all, ylim=c(-0.5, 0.5), xlim=c(-0.5, 0.5))
abline(c(0,1))
plot(sens_wet~sens_all)
abline(c(0,1))

sens_diff <- sens_wet-sens_dry
hist(sens_diff, nclass=20)
quantile(sens_diff, na.rm=T)
plot(sens_diff~sens_all, ylim=c(-0.5, 0.5), xlim=c(-0.5, 0.5)) # no strong association between difference in sensitivity between wet and dry, and overall sensitivity 

hist(wet_dry_diff)
plot(wet_dry_diff~sens_all)
cor(wet_dry_diff, sens_all, use="pairwise.complete") # Sensitivity from regression coef is strongly correlated with wet minus dry year difference in RWI

# plot absolute growth vs sensitivity
ggplot(sensdata, aes(y=sens_all, x=mean_bai)) + theme_light() + geom_point() + facet_wrap(~cluster) + stat_smooth(method="lm")

# plot sensitivity vs normal precipitation
sensdata$hirad <- as.numeric(sensdata$rad.tot > mean(sensdata$rad.tot))

ggplot(sensdata, aes(y=sens_all, x=ppt.norm, group=hirad, col=hirad)) + theme_light() + geom_point() + stat_smooth(method="lm") #+ facet_wrap(~cluster.x)

ggplot(sensdata[which(sensdata$cluster.x=="Plumas"),], aes(y=sens_all, x=ppt.norm, group=hirad, col=hirad)) + theme_light() + geom_point() + stat_smooth(method="lm") #+ facet_wrap(~cluster.x)



#########################################################
###### Plot PSME growth in refugium vs. non-refugium ####
#########################################################
# Adapted from Derek's script: dendro_weather_analysis_p2.plot.R

#are trees from refugium sites less sensitive?

# rename refuge plots

d.refinc <- years #trees.clim.rwi

# Set time interval to look at
d.refinc <- subset(d.refinc,(year > 1994) & (year < 2015))
#d.refinc <- subset(d.refinc,(year > 1960) & (year < 1990))

refuge.trees <- c("2589","2590","1222","2550","2551","2552")
refuge.plots <- c("SB01","RT02C","RS01")

d.refinc <- d.refinc[d.refinc$plot.id %in% refuge.plots,]

d.refinc$type <- ifelse(d.refinc$tree.id %in% refuge.trees,"refugium","standard")

## PSME only
d.refinc <- d.refinc[d.refinc$species == "PSME",]

# Fit model using RGR of basal area
m <- lm(bai.ba.comb ~ plot.id*ppt.z + plot.id*type + ppt.z*type ,data=d.refinc)
# Alternatively, fit model using RWI
#m <- lm(rwi ~ plot.id*ppt.z + plot.id*type + ppt.z*type ,data=d.refinc)



ppt.z.seq <- seq(from=-2,to=2,length.out=100)

# Make projections for multiple refugia areas separately
nd.r.1 <- data.frame(plot.id="SB01",ppt.z=ppt.z.seq,type="refugium")
nd.s.1 <- data.frame(plot.id="SB01",ppt.z=ppt.z.seq,type="standard")
nd.r.2 <- data.frame(plot.id="RT02C",ppt.z=ppt.z.seq,type="refugium")
nd.s.2 <- data.frame(plot.id="RT02C",ppt.z=ppt.z.seq,type="standard")
nd.r.3 <- data.frame(plot.id="RS01",ppt.z=ppt.z.seq,type="refugium")
nd.s.3 <- data.frame(plot.id="RS01",ppt.z=ppt.z.seq,type="standard")
newdata <- rbind(nd.r.1,nd.s.1,nd.r.2,nd.s.2,nd.r.3,nd.s.3)



pred <- as.data.frame(predict(m,newdata=newdata,interval="confidence"))

pred <- cbind(pred,newdata)


library(ggplot2)

pred$type <- factor(pred$type,levels=rev(levels(pred$type)))

pred$plot.id <- as.numeric(pred$plot.id)
pred$plot.id <- paste("P",pred$plot.id,sep="")


p <- ggplot(pred,aes(x=ppt.z,y=fit,colour=type,fill=type)) +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=lwr,ymax=upr),linetype=0,alpha=0.5) +
  facet_grid(. ~ plot.id,scales="free_x",space="free_x") +
  theme_bw(20) +
  labs(x="Precipitation anomaly",y="Relative growth rate") +
  scale_color_manual(values=c("#1A4D1B","#009999"),name="Type") +
  scale_fill_manual(values=c("#1A4D1B","#009999"),name="Type")
p

library(Cairo)
Cairo(file="test2.png",typ="png",width=960,height=600)
p
dev.off()


#### Pictures of tree growth for ESA slide show ####
par(mfrow=c(2,1))
plot(rwi~year, years[years$tree.id==1231,], pch=16, col="slateblue", xlab="Year", ylab="Ring width index", cex.axis=1.2, cex.lab=1.4)
plot(bai.comb/10000~year, years[years$tree.id==1231,], pch=16, col="darkgreen", xlab="Year", ylab="Growth (basal area increment)", cex.axis=1.2, cex.lab=1.4)
        