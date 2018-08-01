# start with trees.clim.rwi
#            tree.norm
#            plot.norm

d <- trees.clim.rwi
d <- subset(d,(year > 1981) & (year < 2015))

refuge.trees <- c("2589","2590","1222","2550","2551","2552")
refuge.plots <- c("SB01","RT02C","RS01")
d.refuge <- d[d$tree.id %in% refuge.trees,]
d <- d[!(d$tree.id %in% refuge.trees),]

library(ggplot2)
library(viridis)

##############
##### Figure 1
##############

d.recent <- d[d$year %in% 1985:2014,]
# only include bai.ba values when ba was over 15 cm DBH
d.recent <- d.recent[d.recent$ba.comb > 17662,]
d.recent <- merge(d.recent,tree.norm[,c("tree.id","plot.id","cluster","rad.hr.03")],by="tree.id",all.x=TRUE)
d.recent$dbh <- as.numeric(d.recent$dbh)
d.recent$voronoi.prop <- d.recent$voronoi.area*100 / ((d.recent$dbh/2)^2)*3.14

d.recent <- d.recent[d.recent$species %in% c("PSME","PIPO"),]

m.all.sp <- lm(bai.ba.comb ~ ppt.z*tmean.normal*species + cluster*ppt.normal*species + ppt.z*ppt.normal*species+ppt.z*cluster*species+ppt.z*rad.hr.03.x*species + ppt.z*voronoi.prop*species + ppt.normal*voronoi.prop*species + tmean.z*species ,data=d.recent)

## predict growth in average and z=-1.5 for PSME and PIPO along ppt.normal gradient

clusters <- c("Plumas","Tahoe","Yose","Sierra")

diff.lowprecip.clusters <- data.frame()
diff.medprecip.clusters <- data.frame()

for(cluster in clusters) {
  
  # tahoe prediction sequence
  length.out <- 100
  min <- min(d.recent[d.recent$cluster==cluster,]$ppt.normal)
  max <- max(d.recent[d.recent$cluster==cluster,]$ppt.normal)
  ppt.seq <- seq(from=min,to=max,length.out=length.out)
  newdata.medprecip.psme <- data.frame(cluster=cluster,ppt.normal=ppt.seq,rad.hr.03.x=11,voronoi.prop=10,ppt.z=0,tmean.normal=12,tmean.z=0,species="PSME")
  newdata.lowprecip.psme <- data.frame(cluster=cluster,ppt.normal=ppt.seq,rad.hr.03.x=11,voronoi.prop=10,ppt.z=-1.5,tmean.normal=12,tmean.z=0,species="PSME")
  newdata.medprecip.pipo <- data.frame(cluster=cluster,ppt.normal=ppt.seq,rad.hr.03.x=11,voronoi.prop=10,ppt.z=0,tmean.normal=12,tmean.z=0,species="PIPO")
  newdata.lowprecip.pipo <- data.frame(cluster=cluster,ppt.normal=ppt.seq,rad.hr.03.x=11,voronoi.prop=10,ppt.z=-1.5,tmean.normal=12,tmean.z=0,species="PIPO")
  pred.lowprecip.psme <- as.data.frame(predict(m.all.sp,newdata=newdata.lowprecip.psme,se.fit=TRUE))
  pred.medprecip.psme <- as.data.frame(predict(m.all.sp,newdata=newdata.medprecip.psme,se.fit=TRUE))
  pred.lowprecip.pipo <- as.data.frame(predict(m.all.sp,newdata=newdata.lowprecip.pipo,se.fit=TRUE))
  pred.medprecip.pipo <- as.data.frame(predict(m.all.sp,newdata=newdata.medprecip.pipo,se.fit=TRUE))
  
  #bootstrap the difference between PSME and PIPO for each prediction point
  
 diff.lowprecip <- data.frame()
 diff.medprecip <- data.frame()
  
  
  for (i in 1:length.out) {
    
    #lowprecip
    mean.psme <- pred.lowprecip.psme[i,]$fit
    se.psme <- pred.lowprecip.psme[i,]$se.fit
    pred.psme <- rnorm(1000,mean.psme,se.psme)
    
    mean.pipo <- pred.lowprecip.pipo[i,]$fit
    se.pipo <- pred.lowprecip.pipo[i,]$se.fit
    pred.pipo <- rnorm(1000,mean.pipo,se.pipo)
    
    diff <- pred.psme - pred.pipo
    diff.sd <- sd(diff)
    diff.mean <- mean(diff)
    diff.lwr <- diff.mean - diff.sd*1.96
    diff.upr <- diff.mean + diff.sd*1.96
    
    diff.lowprecip <- rbind(diff.lowprecip,data.frame(diff.mean,diff.upr,diff.lwr,cluster))
    
    #medprecip
    mean.psme <- pred.medprecip.psme[i,]$fit
    se.psme <- pred.medprecip.psme[i,]$se.fit
    pred.psme <- rnorm(1000,mean.psme,se.psme)
    
    mean.pipo <- pred.medprecip.pipo[i,]$fit
    se.pipo <- pred.medprecip.pipo[i,]$se.fit
    pred.pipo <- rnorm(1000,mean.pipo,se.pipo)
    
    diff <- pred.psme - pred.pipo
    diff.sd <- sd(diff)
    diff.mean <- mean(diff)
    diff.lwr <- diff.mean - diff.sd*1.96
    diff.upr <- diff.mean + diff.sd*1.96
    
    diff.medprecip <- rbind(diff.medprecip,data.frame(diff.mean,diff.upr,diff.lwr,cluster))
  
  }
  
  
  diff.lowprecip$ppt.normal <- ppt.seq
  diff.medprecip$ppt.normal <- ppt.seq
  
  diff.lowprecip.clusters <- rbind(diff.lowprecip.clusters,diff.lowprecip)
  diff.medprecip.clusters <- rbind(diff.medprecip.clusters,diff.medprecip)
  
}


diff.lowprecip.clusters$ppt.z <- "low"
diff.medprecip.clusters$ppt.z <- "med"

diff.clusters <- rbind(diff.lowprecip.clusters,diff.medprecip.clusters)


ggplot(diff.clusters,aes(x=ppt.normal,y=diff.mean,colour=cluster,fill=cluster)) +
  theme_bw() +
  geom_line(size=1,aes(linetype=rev(ppt.z))) +
  geom_ribbon(data=diff.clusters[diff.clusters$ppt.z=="low",],aes(ymin=diff.lwr,ymax=diff.upr,),alpha=0,linetype=2) +
  geom_ribbon(data=diff.clusters[diff.clusters$ppt.z=="med",],aes(ymin=diff.lwr,ymax=diff.upr,),alpha=0.5,linetype=0)





#########################################
### Alternative approach: plot PSME-PIPO diff along axis of ppt.z for each cluster, given the mean ppt.norm of the cluster
#########################################


d.recent <- d[d$year %in% 1985:2014,]
# only include bai.ba values when ba was over 15 cm DBH
d.recent <- d.recent[d.recent$ba.comb > 17662,]
d.recent <- merge(d.recent,tree.norm[,c("tree.id","plot.id","cluster","rad.hr.03")],by="tree.id",all.x=TRUE)
d.recent$dbh <- as.numeric(d.recent$dbh)
d.recent$voronoi.prop <- d.recent$voronoi.area*100 / ((d.recent$dbh/2)^2)*3.14

d.recent <- d.recent[d.recent$species %in% c("PSME","PIPO"),]

m.all.sp <- lm(bai.ba.comb ~ cluster*ppt.normal*species + ppt.z*ppt.normal*species+ppt.z*cluster*species+ppt.z*rad.hr.03.x*species + ppt.z*voronoi.prop*species + ppt.normal*voronoi.prop*species + tmean.z*species ,data=d.recent)
m.all.sp <- lm(bai.ba.comb ~ cluster*species + ppt.z*species+ppt.z*cluster*species+ppt.z*rad.hr.03.x*species + ppt.z*voronoi.prop*species +voronoi.prop*species + tmean.z*species ,data=d.recent)


## predict growth in average and z=-1.5 for PSME and PIPO along ppt.normal gradient

clusters <- c("Plumas","Tahoe","Yose","Sierra")

diff.clusters <- data.frame()

for(cluster in clusters) {
  
  # tahoe prediction sequence
  length.out <- 100
  
  if(cluster=="Tahoe") ppt.normal <- 900
  if(cluster=="Plumas") ppt.normal <- 1600
  if(cluster=="Yose") ppt.normal <- 950
  if(cluster=="Sierra") ppt.normal <- 900
  
  ppt.z.seq <- seq(from=-2,to=2,length.out=length.out)
  newdata.psme <- data.frame(cluster=cluster,ppt.normal=ppt.normal,rad.hr.03.x=10,voronoi.prop=10,ppt.z=ppt.z.seq,tmean.normal=12,tmean.z=0,species="PSME")
  newdata.pipo <- data.frame(cluster=cluster,ppt.normal=ppt.normal,rad.hr.03.x=10,voronoi.prop=10,ppt.z=ppt.z.seq,tmean.normal=12,tmean.z=0,species="PIPO")
  pred.psme <- as.data.frame(predict(m.all.sp,newdata=newdata.psme,se.fit=TRUE))
  pred.pipo <- as.data.frame(predict(m.all.sp,newdata=newdata.pipo,se.fit=TRUE))
  
  #bootstrap the difference between PSME and PIPO for each prediction point
  
  diff <- data.frame()

  for (i in 1:length.out) {
    
    mean.psme <- pred.psme[i,]$fit
    se.psme <- pred.psme[i,]$se.fit
    sim.psme <- rnorm(1000,mean.psme,se.psme)
    
    mean.pipo <- pred.pipo[i,]$fit
    se.pipo <- pred.pipo[i,]$se.fit
    sim.pipo <- rnorm(1000,mean.pipo,se.pipo)
    
    diff.sim <- sim.psme - sim.pipo
    diff.sd <- sd(diff.sim)
    diff.mean <- mean(diff.sim)
    diff.lwr <- diff.mean - diff.sd
    diff.upr <- diff.mean + diff.sd
    
    diff <- rbind(diff,data.frame(diff.mean,diff.upr,diff.lwr,cluster))
    
  }
  
  
  diff$ppt.z <- ppt.z.seq
  
  diff.clusters <- rbind(diff.clusters,diff)
  
}


ggplot(diff.clusters,aes(x=ppt.z,y=diff.mean,color=cluster,fill=cluster)) +
  theme_bw() +
  geom_line(size=1) +
  geom_ribbon(aes(ymin=diff.lwr,ymax=diff.upr,),alpha=0.5,linetype=0) +
  scale_y_continuous(limits=c(-0.005,0.045))




#############################################################################
##### Plot PSME vs. PIPO growth vs. ppt.z (evaluate sensitivity vs. growth)
#############################################################################



d.recent <- d[d$year %in% 1985:2014,]
# only include bai.ba values when ba was over 15 cm DBH
d.recent <- d.recent[d.recent$ba.comb > 17662,]
d.recent <- merge(d.recent,tree.norm[,c("tree.id","plot.id","cluster","rad.hr.03")],by="tree.id",all.x=TRUE)
d.recent$dbh <- as.numeric(d.recent$dbh)
d.recent$voronoi.prop <- d.recent$voronoi.area*100 / ((d.recent$dbh/2)^2)*3.14

d.recent <- d.recent[d.recent$species %in% c("PSME","PIPO"),]


m.all.sp <- lm(bai.ba.comb ~ cluster*ppt.normal*species + ppt.z*ppt.normal*species+ppt.z*cluster*species+ppt.z*rad.hr.03.x*species + ppt.z*voronoi.prop*species + ppt.normal*voronoi.prop*species + tmean.z*species ,data=d.recent)
#m.all.sp <- lm(bai.ba.comb ~ cluster*species + ppt.z*species+ppt.z*cluster*species+ppt.z*rad.hr.03.x*species + ppt.z*voronoi.prop*species +voronoi.prop*species + tmean.z*species ,data=d.recent)





# dry end of each cluster, vary competition
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
  "ppt.normal" = c(rep(1600,4),rep(900,4),rep(950,4),rep(925,4)),
  "rad.hr.03.x" = 10,
  "voronoi.prop" = rep(c(10,40,70,100),4))

# dry end of each cluster, vary radiation
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "ppt.normal" = c(rep(1600,4),rep(900,4),rep(950,4),rep(925,4)),
                          "rad.hr.03.x" = rep(c(7,8,10,11),4),
                          "voronoi.prop" = 25)


# dry everywhere, vary radiation and competition (but not at same time)
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "ppt.normal" = c(rep(1600,4),rep(900,12)),
                          "rad.hr.03.x" = rep(c(7,11,9.5,9.5),4),
                          "voronoi.prop" = rep(c(25,25,10,100),4))



min(d.recent[d.recent$cluster=="Sierra",]$ppt.normal)


# precip gradient representative of each cluster
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "ppt.normal" = c(2000,1850,1700,1550,1700,1400,1100,900,1100,1050,1000,950,950,935,920,905),
                          "rad.hr.03.x" = 9.5,
                          "voronoi.prop" = 25)



# ppt of sierra everywhere
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "ppt.normal" = 905,
                          "rad.hr.03.x" = 9.5,
                          "voronoi.prop" = 25)




library(grid)
library(gridExtra)

plots <- list()

for(i in 1:nrow(param.combs)) {
  
  param.comb <- param.combs[i,]
  
  ppt.z.seq <- seq(from=-2,to=2,length.out=100)
  newdata.psme <- data.frame(param.comb,ppt.z=ppt.z.seq,tmean.normal=12,tmean.z=0,species="PSME")
  newdata.pipo <- data.frame(param.comb,ppt.z=ppt.z.seq,tmean.normal=12,tmean.z=0,species="PIPO")
  newdata <- rbind(newdata.psme,newdata.pipo)
  
  pred <- predict(m.all.sp,newdata=newdata,interval="confidence")
  
  pred.newdata <- data.frame(pred,newdata)
  
  library(gridExtra)
  
  plots[[i]] <- ggplot(pred.newdata,aes(x=ppt.z,y=fit,color=species,fill=species)) +
    theme_bw() +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=lwr,ymax=upr,),alpha=0.5,linetype=0) +
    scale_fill_manual(values=c("#0D6E23","#A180A9")) +
    scale_color_manual(values=c("#0D6E23","#A180A9")) +
    guides(fill=FALSE,color=FALSE) +
    scale_y_continuous(limits=c(-0.005,0.095)) +
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"lines")) +
    labs(x=NULL,y=NULL)
  
  if(i %in% c(2,3,4,6,7,8,10,11,12)) {
    plots[[i]] <- plots[[i]] + theme(axis.text.x=element_text(color="white"),axis.text.y=element_text(color="white"))
  }
  if(i %in% c(1,5,9)) {
    plots[[i]] <- plots[[i]] + theme(axis.text.x=element_text(color="white"))
  } 
  if(i %in% c(14,15,16)) {
    plots[[i]] <- plots[[i]] + theme(axis.text.y= element_text(color="white"))
  } 
  
  

}

library(Cairo)
Cairo(file="test.png",typ="png",width=900,height=800)

do.call(grid.arrange,plots)
dev.off()





#############################################################################
##### Plot PSME vs. PIPO growth vs. tmean.z (evaluate sensitivity vs. growth)
#############################################################################

m.all.sp <- lm(bai.ba.comb ~ cluster*tmean.normal*species + tmean.z*tmean.normal*species+tmean.z*cluster*species+tmean.z*rad.hr.03.x*species + tmean.z*voronoi.prop*species + tmean.normal*voronoi.prop*species + ppt.z*species ,data=d.recent)
#m.all.sp <- lm(bai.ba.comb ~ cluster*species + ppt.z*species+ppt.z*cluster*species+ppt.z*rad.hr.03.x*species + ppt.z*voronoi.prop*species +voronoi.prop*species + tmean.z*species ,data=d.recent)


# 
# 
# 
# # dry end of each cluster, vary competition
# param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
#                           "ppt.normal" = c(rep(1600,4),rep(900,4),rep(950,4),rep(925,4)),
#                           "rad.hr.03.x" = 10,
#                           "voronoi.prop" = rep(c(10,40,70,100),4))
# 
# # dry end of each cluster, vary radiation
# param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
#                           "ppt.normal" = c(rep(1600,4),rep(900,4),rep(950,4),rep(925,4)),
#                           "rad.hr.03.x" = rep(c(7,8,10,11),4),
#                           "voronoi.prop" = 25)


# hot everywhere, vary radiation and competition (but not at same time)
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "tmean.normal" = c(rep(15.5,8),rep(14.25,4),rep(11.5,4)),
                          "rad.hr.03.x" = rep(c(7,11,9.5,9.5),4),
                          "voronoi.prop" = rep(c(25,25,10,100),4))



# cold everywhere, vary radiation and competition (but not at same time)
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "tmean.normal" = c(rep(10.5,16)),
                          "rad.hr.03.x" = rep(c(7,11,9.5,9.5),4),
                          "voronoi.prop" = rep(c(25,25,10,100),4))




# comparable temperature everywhere, vary radiation and competition (but not at same time)
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "tmean.normal" = rep(12,16),
                          "rad.hr.03.x" = rep(c(7,11,9.5,9.5),4),
                          "voronoi.prop" = rep(c(25,25,10,100),4))




min(d.recent[d.recent$cluster=="Sierra",]$ppt.normal)


# temp gradient representative of each cluster
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "tmean.normal" = c(10.5,12,13.5,15.5,9.5,12,14,15.5,11,12.25,13.5,14.5,11,11.5,12,12.5),
                          "rad.hr.03.x" = 9.5,
                          "voronoi.prop" = 25)


# temp gradient representative of each cluster
param.combs <- data.frame("cluster"=c(rep("Plumas",4),rep("Tahoe",4),rep("Yose",4),rep("Sierra",4)),
                          "tmean.normal" = 15.5,
                          "rad.hr.03.x" = 9.5,
                          "voronoi.prop" = 25)



library(grid)

plots <- list()

for(i in 1:nrow(param.combs)) {
  
  param.comb <- param.combs[i,]
  
  tmean.z.seq <- seq(from=-2,to=2,length.out=100)
  newdata.psme <- data.frame(param.comb,tmean.z=tmean.z.seq,ppt.z=0,species="PSME")
  newdata.pipo <- data.frame(param.comb,tmean.z=tmean.z.seq,ppt.z=0,species="PIPO")
  newdata <- rbind(newdata.psme,newdata.pipo)
  
  pred <- predict(m.all.sp,newdata=newdata,interval="confidence")
  
  pred.newdata <- data.frame(pred,newdata)
  
  library(gridExtra)
  
  plots[[i]] <- ggplot(pred.newdata,aes(x=tmean.z,y=fit,color=species,fill=species)) +
    theme_bw() +
    geom_line(size=1) +
    geom_ribbon(aes(ymin=lwr,ymax=upr,),alpha=0.5,linetype=0) +
    scale_fill_manual(values=c("#0D6E23","#A180A9")) +
    scale_color_manual(values=c("#0D6E23","#A180A9")) +
    guides(fill=FALSE,color=FALSE) +
    scale_y_continuous(limits=c(-0.01,0.065)) +
    theme(plot.margin=unit(c(0.1,0.1,0.1,0.1),"lines")) +
    labs(x=NULL,y=NULL)
  
  if(i %in% c(2,3,4,6,7,8,10,11,12)) {
    plots[[i]] <- plots[[i]] + theme(axis.text.x=element_text(color="white"),axis.text.y=element_text(color="white"))
  }
  if(i %in% c(1,5,9)) {
    plots[[i]] <- plots[[i]] + theme(axis.text.x=element_text(color="white"))
  } 
  if(i %in% c(14,15,16)) {
    plots[[i]] <- plots[[i]] + theme(axis.text.y= element_text(color="white"))
  } 

  
  
  
}



library(Cairo)
Cairo(file="test2.png",typ="png",width=900,height=800)

do.call(grid.arrange,plots)
dev.off()





#########################################################
###### Plot PSME growth in refugium vs. non-refugium ####
#########################################################

#are trees from refugium sites less sensitive?

# rename refuge plots

d.refinc <- trees.clim.rwi
d.refinc <- subset(d.refinc,(year > 1981) & (year < 2015))

refuge.trees <- c("2589","2590","1222","2550","2551","2552")
refuge.plots <- c("SB01","RT02C","RS01")

d.refinc <- d.refinc[d.refinc$plot.id %in% refuge.plots,]

d.refinc$type <- ifelse(d.refinc$tree.id %in% refuge.trees,"refugium","standard")

## PSME only
d.refinc <- d.refinc[d.refinc$species == "PSME",]

m <- lm(bai.ba.comb ~ plot.id*ppt.z + plot.id*type + ppt.z*type ,data=d.refinc)


ppt.z.seq <- seq(from=-2,to=2,length.out=100)


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

library(Cairo)
Cairo(file="test2.png",typ="png",width=960,height=600)
p
dev.off()






##############
##### PSME and OTHER growth amount under average and dry climate
##############

fits <- data.frame()

for(tree in unique(d.recent$tree.id)) {
  
  d.tree <- subset(d.recent,tree.id==tree)
  if(sum(!(is.na(d.tree$rwi))) < 25) next() # skip if tree has less than 25 rings in the focal period
  
  m <- lm(bai.ba.comb ~ ppt.z, data=d.tree)
  
  newdata <- data.frame(ppt.z=0,ppt.z1=0,tmean.z=0,tmean.z1=0)
  pred <- predict(m,newdata,interval="confidence")
  med <- data.frame(tree.id=tree,pred,scenario="med")
  
  newdata <- data.frame(ppt.z=-1,ppt.z1=-1,tmean.z=0,tmean.z1=0)
  pred <- predict(m,newdata,interval="confidence")
  low <- data.frame(tree.id=tree,pred,scenario="low")

  fits <- rbind(fits,med,low)
}



# get species and plot columns
trees.fits <- merge(tree.norm,fits,by="tree.id")

# divide between PSME and OTHER
trees.fits$species <- ifelse(trees.fits$species=="PSME","PSME","OTHER")

# only include plots that have both groups
a <- table(trees.fits$species,trees.fits$plot.id)
a <- ifelse(a>0,1,0)
a <- colSums(a)
plots.keep <- names(a[which(a==2)])
trees.fits <- trees.fits[trees.fits$plot.id %in% plots.keep,]


# aggregate by species and plot and interval
plots.fits <- aggregate(trees.fits[,c("fit","lwr","upr")],by=list(trees.fits$plot.id,trees.fits$species,trees.fits$scenario),FUN=mean)
plots.fits$plot.id <- plots.fits$Group.1
plots.fits$species <- plots.fits$Group.2
plots.fits$scenario <- plots.fits$Group.3
plots.fits <- plots.fits[,-c(1,2,3)]

# merge with plot-level variables
plots.fits <- merge(plots.fits,plot.norm,by="plot.id")

# sort plot levels by cluster then increasing precip
plot.levels <- unique(plots.fits$plot.id[order(plots.fits$cluster,-plots.fits$ppt.normal)])
plots.fits$plot.id <- factor(plots.fits$plot.id,levels=plot.levels)

plots.fits$cluster <- factor(plots.fits$cluster,levels=c("Plumas","Tahoe","Yose","Sierra"))

levels(plots.fits$cluster) <- c("Plumas","Tahoe","Yosemite","Sierra")

levels(plots.fits$scenario) <- c("Average","Dry")

plots.fits$species.x <- factor(plots.fits$species.x,levels=c("PSME","OTHER"))

#plots.fits$plot.id <- paste("P",as.numeric(plots.fits$plot.id),sep="")

# plot these values
# plot these values
p <- ggplot(plots.fits,aes(x=plot.id,y=fit,colour=species.x,pch=scenario)) +
  geom_point(position=position_dodge(width=0.75),size=3) +
  geom_errorbar(aes(ymin=lwr,ymax=upr),width=0,size=0.5,position=position_dodge(width=0.75))+
  facet_grid(species.x ~ cluster,space="free_x",scales="free_x") +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
  scale_colour_manual(name="Species",values=c("#0D6E23","#725279"),labels=c("Doug-fir","Other")) +
  labs(y="Relative growth rate",x="Plot") +
  theme(text = element_text(size=16),strip.text.y=element_blank(),strip.background=element_blank())

library(Cairo)
Cairo(file="test2.png",typ="png",width=960,height=600)
p
dev.off()



p <- ggplot(plots.fits,aes(x=plot.id,y=1,fill=ppt.normal)) +
  geom_tile(stat="identity") +
  facet_grid(. ~ cluster,scales="free_x",space="free_x") +
  scale_fill_gradient(low="#993300",high="#99FFFF",name="Normal precipitation (mm)") +
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
  theme(legend.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(legend.position="none")

library(Cairo)
Cairo(file="test2.png",typ="png",width=960,height=600)
p
dev.off()









p <- ggplot(plots.fits[plots.fits$species.x=="PSME",],aes(x=plot.id,y=1,fill=ppt.normal)) +
  geom_tile(stat="identity") +
  facet_grid(. ~ cluster,scales="free_x",space="free_x") +
  scale_fill_gradient(low="#993300",high="#99FFFF",name="Normal precipitation (mm)") +
  theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5)) +
  theme(legend.background = element_rect(fill = "transparent",colour = NA),plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(legend.position="none")

library(Cairo)
Cairo(file="test2.png",typ="png",width=960,height=600,canvas="transparent")
p
dev.off()














##############
##### PSME and PIPO growth-climate relationship in plots with low and high PSME proportion
##############

d.recent <- d[d$year %in% 1985:2014,]
# only include bai.ba values when ba was over 15 cm DBH
d.recent <- d.recent[d.recent$ba.comb > 17662,]
d.recent <- merge(d.recent,tree.norm[,c("tree.id","plot.id","cluster","rad.hr.03")],by="tree.id",all.x=TRUE)
d.recent$dbh <- as.numeric(d.recent$dbh)
d.recent$voronoi.prop <- d.recent$voronoi.area*100 / ((d.recent$dbh/2)^2)*3.14

d.recent <- d.recent[d.recent$species %in% c("PSME","PIPO"),]
d.recent <- d.recent[d.recent$cluster=="Yose",]

psme.pipo.prop <- read.csv("Processed data/psme_pipo_prop.csv",header=TRUE)

psme.pipo.prop <- merge(psme.pipo.prop,plot.norm[,c("plot.id","cluster")],by="plot.id")

psme.pipo.prop <- psme.pipo.prop[psme.pipo.prop$cluster=="Yose",]


quant.33 <- quantile(psme.pipo.prop$psme.ratio,0.40)
quant.66 <- quantile(psme.pipo.prop$psme.ratio,0.60)

psme.pipo.prop$psme.pipo.prop <- "med"
psme.pipo.prop$psme.pipo.prop <- ifelse(psme.pipo.prop$psme.ratio > quant.66,"high",psme.pipo.prop$psme.pipo.prop)
psme.pipo.prop$psme.pipo.prop <- ifelse(psme.pipo.prop$psme.ratio < quant.33,"low",psme.pipo.prop$psme.pipo.prop)

psme.pipo.prop <- psme.pipo.prop[,c("plot.id","psme.ratio","psme.pipo.prop")]


d.recent <- merge(d.recent,psme.pipo.prop,by.x="plot.id.x",by.y="plot.id")




d.recent <- d.recent[d.recent$psme.pipo.prop != "med",]


m.psme <- lm(bai.ba.comb ~ ppt.z*psme.pipo.prop*ppt.normal,data=d.recent[d.recent$species=="PSME",])
m.pipo <- lm(bai.ba.comb ~ ppt.z*psme.pipo.prop*ppt.normal,data=d.recent[d.recent$species=="PIPO",])
m <- lm(bai.ba.comb ~ ppt.z*psme.pipo.prop*species,data=d.recent)


m.lowprop <- lm(bai.ba.comb ~ (ppt.z+species)^2,data=d.recent[d.recent$psme.pipo.prop=="low",])
m.highprop <- lm(bai.ba.comb ~ (ppt.z+species)^2,data=d.recent[d.recent$psme.pipo.prop=="high",])



ppt.z.seq <- seq(from=-2,to=2,length.out=100)

newdata.psme <- data.frame(ppt.z=ppt.z.seq,species="PSME",ppt.normal=1000)
newdata.pipo <- data.frame(ppt.z=ppt.z.seq,species="PIPO",ppt.normal=1000)

pred.psme <- as.data.frame(predict(m.lowprop,newdata.psme,interval="confidence"))
pred.pipo <- as.data.frame(predict(m.lowprop,newdata.pipo,interval="confidence"))

pred.psme$ppt.z <- ppt.z.seq
pred.pipo$ppt.z <- ppt.z.seq
pred.psme$species <- "PSME"
pred.pipo$species <- "PIPO"

pred <- rbind(pred.psme,pred.pipo)


ggplot(pred,aes(x=ppt.z,y=fit,color=species,fill=species)) +
  geom_line() +
  geom_ribbon(aes(ymin=lwr,ymax=upr),alpha=0.5,linetype=0) +
  scale_y_continuous(limits=c(0,0.05))


########################################
#### For plots, plot temp vs. precip ###
########################################

plot.norm

library(Cairo)

Cairo(file="test3.png",typ="png",width=600,height=600,canvas="transparent")



ggplot(plot.norm,aes(x=ppt.normal,y=tmean.normal,fill=cluster)) +
  geom_point(size=5,pch=21,color="black") +
  theme_bw(20) +
  scale_fill_viridis(discrete=TRUE,breaks=c("Plumas","Tahoe","Yose","Sierra"),labels=c("Plumas","Tahoe","Yosemite","Sierra"),name="Cluster") +
  labs(x="Annual precipitation (mm)",y="Annual temperature (deg C)")



dev.off()

