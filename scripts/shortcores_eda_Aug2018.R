# This script: 
# loads the data from short cores collected along elevational transect from Auburn Rec area to above Forest Hill in June and August 2015. 
# visualizes relationships between water content and starch/NSC content of the cores and characteristics of trees and environment
# Compares short-term sensitivity of physiological status (June vs August of a drought year) to longer-term sensitivity based on Derek's tree ring analysis. 


## Load data and set up ####

library(lme4); library(ggplot2); library(lattice)

setwd("/Volumes/GoogleDrive/My Drive/TreeDrought")
d = read.csv("JunAug2015_w_treerings.csv")
head(d)

## Visualize correlations in data ####

# Overall correlations
pairs(d[1:28,22:31])

# Water content vs tree size
plot(RWC_June~DBH, d, pch=16, col="cyan4", xlab="Tree size (DBH cm)", ylab="Water content (proportion)", cex.axis=1.2, cex.lab=1.5)
abline(coef(lm(RWC_June~DBH, d)), lwd=2, col="cyan4")
points(RWC_Aug~DBH, d, pch=16, col="orange3")
abline(coef(lm(RWC_Aug~DBH, d)),lwd=2, col="orange3")
legend("bottomright", c("June", "August"),pch=c(16, 16), col=c("cyan4", "orange3"), cex=1.5)

plot(dRWC~DBH, d, pch=16, col="black")

plot(RWC_June~plot.deficit, d, pch=16, col="cyan4")
abline(coef(lm(RWC_June~plot.deficit, d)), lwd=2, col="cyan4")
points(RWC_Aug~plot.deficit, d, pch=16, col="orange3")
abline(coef(lm(RWC_Aug~plot.deficit, d)),lwd=2, col="orange3")

plot(STARCH~plot.deficit, d, pch=16, col="cyan4")

plot(dRWC~plot.deficit, d, pch=16, col="orange3")
plot(STARCH~dRWC, d, pch=16, col="orange3") # weak hint of assn between water drop and starch

summary(lm(STARCH~dRWC, d))

plot(TNSC~plot.deficit, d, pch=16, col="cyan4")
abline(coef(lm(TNSC~plot.deficit, d)), lwd=2, col="cyan4")
summary(lm(TNSC~plot.deficit, d)) # trend, nonsignif

plot(TNSC~dRWC, d, pch=16, col="orange3")



plot(dRWC~plot.deficit, d, pch=16, col="black") # hint of greater water loss in higher water deficit
summary(lm(dRWC~plot.deficit, d)) # far from significant however

plot(RWC_June~plot.elev.m, d, pch=16, col="cyan4")
points(RWC_Aug~plot.elev.m, d, pch=16, col="orange3")

plot(dRWC~plot.elev.m, d, pch=16, col="black") # nothing

plot(RWC_June~plot.aet, d, pch=16, col="cyan4")
points(RWC_Aug~plot.aet, d, pch=16, col="orange3")

plot(dRWC~plot.aet, d, pch=16, col="black") # nothing


summary(lm(RWC_Aug~DBH*plot.elev.m, d))
summary(lm(dRWC~DBH*plot.elev.m, d))

summary(lm(RWC_Aug~DBH*STARCH, d))
summary(lm(dRWC~DBH*STARCH, d))
summary(lm(dKh~DBH*STARCH, d))

# Comparing water content against environmental factors
boxplot(RWC_Aug~PLOT, d, notch=T)
plot(RWC_Aug~plot.deficit, d)
plot(dKh.dt~plot.deficit, d, xlab="Annual water deficit (mm)", ylab="dKh/dt", cex.axis=1.2, cex.lab=1.5, pch=16, col="slateblue") # assn with water deficit
abline(coef(lm(dKh.dt~plot.deficit, d)), col="darkgray", lwd=3)
ggplot(d, aes(x=plot.deficit, y=dKh.dt)) + geom_point() + geom_smooth(method="lm") + theme_classic()
plot(dKh.dt~plot.elev.m, d) # no assn with elevation
boxplot(RWC_June~PLOT, d, notch=T)

# Comparing starch content against environmental factors
plot(STARCH~plot.deficit, d, xlab="Annual water deficit (mm)", ylab="Starch content", cex.axis=1.2, cex.lab=1.5, pch=16, col="slateblue") # assn with water deficit
abline(coef(lm(STARCH~plot.deficit, d)), col="darkgray", lwd=3)

# Total carbs against environment
plot(TNSC~plot.deficit, d, xlab="Annual water deficit (mm)", ylab="Total nonstructural carbs", cex.axis=1.2, cex.lab=1.5, pch=16, col="goldenrod") # assn with water deficit
summary(lm(TNSC~plot.deficit, d))
abline(coef(lm(TNSC~plot.deficit, d)), col="darkgray", lwd=3, lty=2)
# and vs tree size
plot(TNSC~DBH, d, xlab="Tree size (DBH cm)", ylab="Total nonstructural carbs", cex.axis=1.2, cex.lab=1.5, pch=16, col="goldenrod") # assn with water deficit
abline(coef(lm(TNSC~DBH, d)), col="darkgray", lwd=3)



cols = c("plot.elev.m","plot.deficit", "DBH", "BAF20_CALC","Density", "RWC_Aug", "RWC_June", "dKh", "dKh.dt","sd.rwi.20yr", "sd.rwi.all", "STARCH", "TNSC")
pairs(d[,cols])

# Lots of suggestive looking patterns here:

# Lower TNSC associated with greater reduction in Kh, and with lower water content. 
# Higher TNSC and water content assocaited with higher dbh and higher wood density. 
# Higher elev associated with higher starch
# Higher deficit associated with lower starch


# Comparing physiological status to wood density 

anova(lm(Density~PLOT, d)) # most variation in density within plot, about 1/3 among plots
plot(Density~plot.deficit, d, pch=16, col="brown")
plot(Density~plot.elev.m, d, pch=16, col="brown") # no association of density with defict or elevation

summary(lm(dKh.dt~Density, d))


plot(TNSC~Density, d)
plot(STARCH~Density, d)



plot(dKh~Density, d, pch=16, col="brown")
abline(coef(lm(dKh~Density, d)), lwd=3, col="gray")
# Higher density wood associated with less reduction in conductance

par(pty="s")
plot(dKh.dt~Density, d, pch=16, col="brown", xlab="Wood density", ylab="Reduction in conductance", cex.axis=1.2, cex.lab=1.5)
abline(coef(lm(dKh.dt~Density, d)), lwd=3, col="gray")



## comparing to long-term sensitivity
### NEED TO MERGE WITH SENSITIVITY METRICS FROM TREE RING ANALYSIS



p= qplot(Density, dKh.dt,  data=d, geom=c("point", "smooth"), method="lm") 
p = p+ theme_bw()
p
summary(lm(dKh.dt~Density, d))
summary(lmer(dKh.dt~Density+(1|PLOT), d))

p = ggplot(d, aes(x=dKh, y=Density)) + geom_point() + geom_smooth() + theme_bw()
p

p = ggplot(d, aes(y=dKh.dt, x=TNSC)) + geom_point() + geom_smooth() + theme_bw()
p
summary(lm(TNSC~RWC_Aug+dKh, d)) # TNSC (but not starch) associated with higher Aug water content


p = ggplot(d, aes(x=RWC_Aug, y=TNSC)) + geom_point() + geom_smooth() + theme_bw()
p
summary(lm(TNSC~RWC_Aug, d))
  # Total carbs correlated with august water content
summary(lm(TNSC~ dRWC, d)) # but not with change in water content

summary(lm(TNSC~dKh.dt, d)) # yes correlated with change in conductance June-Aug 

p = ggplot(d, aes(y=dKh, x=sd.rwi.20yr)) + geom_point() + geom_smooth() + theme_bw()
p
summary(lm(dKh~sd.rwi.20yr, d)) # no assn between sensitivity and dKh

p = ggplot(d, aes(y=dRWC, x=sd.rwi.20yr)) + geom_point() + geom_smooth() + theme_bw()
p
summary(lm(dRWC~sd.rwi.20yr+RWC_June, d)) # weak association between sensitivity and change in dRWC

p = ggplot(d, aes(y=STARCH, x=sd.rwi.20yr)) + geom_point() + geom_smooth() + theme_bw()
p
summary(lm(STARCH~sd.rwi.20yr, d)) # nothing

p = ggplot(d, aes(y=TNSC, x=sd.rwi.20yr)) + geom_point() + geom_smooth() + theme_bw()
p
summary(lm(TNSC~sd.rwi.20yr, d)) # nothing








s = shingle(d$DBH, intervals = matrix(c(27, 37, 57, 37, 57, 77), nrow=3, byrow=F))
s = shingle(d$DBH, intervals = matrix(c(27, 37, 37, 77), nrow=2, byrow=F))

xyplot(RWC_Aug~DBH|plot.elev.m, d) # Higher dbh correlated with higher RWC at all elevs
xyplot(dKh~DBH, d)
xyplot(dKh.dt~plot.deficit|s, d)
xyplot(dKh~plot.deficit, d)



xyplot(RWC_Aug~sd.rwi.20yr|plot.elev.m, d)
xyplot(dKh~sd.rwi.20yr|plot.elev.m, d)
xyplot(STARCH~sd.rwi.20yr|plot.elev.m, d)
xyplot(STARCH~sd.rwi.20yr|s, d)
# nothing much 

xyplot(STARCH~sd.rwi.20yr, d)
xyplot(RWC_Aug~sd.rwi.20yr, d)


xyplot(RWC_Aug~plot.aet|s, d)
xyplot(RWC_Aug~sd.rwi.20yr|s, d)

s = shingle(d$plot.deficit, intervals = matrix(c(397, 475, 475, 600), nrow=2, byrow=F))
xyplot(RWC_Aug~sd.rwi.20yr|s, d)

summary(lm(RWC_Aug~sd.rwi.20yr*plot.deficit, d)) # the model implied by the above plot
# No assn between deficit and sd.rwi 

xyplot(dKh~sd.rwi.20yr|s, d)
summary(lm(dKh~sd.rwi.20yr*plot.deficit, d)) # the model implied by the above plot

plot(dRWC~DBH, d)

dkh_resids = resid(lm(dRWC~DBH, d))
plot(dkh_resids~d$sd.rwi.20yr)
# no obvious interaction between sensitivity and tree size

names(d)
plot(d$dRWC, (d$RWC_Aug-d$RWC_June))


