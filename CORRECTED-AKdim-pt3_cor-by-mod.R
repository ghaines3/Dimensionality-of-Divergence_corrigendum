# Correlations between Trait Suite LD1s

all.LD<-(merge(merge(LD.swim.df[c("FishID","Lake","swimming.LD1")],LD.def.df[c("FishID","defense.LD1")],by="FishID",all = TRUE, no.dups = T),
               LD.tro.df[c("FishID","trophic.LD1")],by.x ="FishID",by.y="FishID", all = TRUE, no.dups = T)%>%
           merge(LD.geo.df[c("CV1","FishID")],by.x="FishID",by.y = "FishID",all = TRUE))
rownames(all.LD)<-all.LD$FishID
## TABLE S8 USES lin.corr.pearson$r and lin.corr.spearman$r
lin.corr.pearson<-rcorr(as.matrix(all.LD[,c("CV1","defense.LD1","swimming.LD1","trophic.LD1")]), type = "pearson")
lin.corr.spearman<-rcorr(as.matrix(all.LD[,c("CV1","defense.LD1","swimming.LD1","trophic.LD1")]), type = "spearman")
# pairwise correlations plot (produces warnings b/c of different 
# sample sizes for different traits).
## FIGURE S1 PANEL A
LD.cor.plot<-ggpairs(all.LD, columns = c("CV1","defense.LD1","swimming.LD1","trophic.LD1"),
                     lower = list(continuous = wrap("points", alpha = 0.5,    size=0.25) )) +
  theme_linedraw()+theme(panel.grid = element_blank())

# 3d scatter plots (exploratory)
plot3d(x=all.LD$defense.LD1,y=all.LD$trophic.LD1,z=all.LD$swimming.LD1,
       xlab="Defense",ylab="Trophic",zlab="Swimming",box=F,sub="",size="4")
threeD.scatter(all.LD,"trophic.LD1",-4,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                         "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                         "sienna4","grey80","grey50"),c(3,4,5),c( 'Swimming LD1','Defense LD1','Trophic LD1'))
# correlations excluding G and Echo (because of extreme spine reduction)
ex.G.Echo<-which(all.LD$Lake!="G"&all.LD$Lake!="Echo")
lin.corr.pearson.ex.G.Echo<-rcorr(as.matrix(all.LD[ex.G.Echo,c("CV1","defense.LD1","swimming.LD1","trophic.LD1")]), type = "pearson")
lin.corr.spearman.ex.G.Echo<-rcorr(as.matrix(all.LD[ex.G.Echo,c("CV1","defense.LD1","swimming.LD1","trophic.LD1")]), type = "spearman")
# pairwise correlations plot excluding G and Echo
## FIGURE S1 PANEL B
LD.cor.plot.ex.G.Echo<-ggpairs(all.LD[ex.G.Echo,], columns = c("CV1","defense.LD1","swimming.LD1","trophic.LD1"), 
                               lower = list(continuous = wrap("points", alpha = 0.5,    size=0.25) )) +
  theme_linedraw()+theme(panel.grid = element_blank())

LD.cor.plot 
ggsave("CORRECTED_LD1_cor.jpg",width = 6,height = 6, dpi = 300)
LD.cor.plot.ex.G.Echo
ggsave("CORRECTED_LD1_cor-noGEcho.jpg",width = 6,height = 6, dpi = 300)




#
#
#
#
# lists trait names by module
def.trait.list<-c("LatPlNum_mean","adj.DS1L","adj.DS2L","adj.PS_mean","MidPl_mean")
swi.trait.list<-c("SL_mm","adj.BD","adj.CP")
tro.trait.list<-c("adj.BCL","adj.GW","adj.PW","RakerNum","adj.JL",
                  "adj.SnL","adj.ED","adj.HL")
two.d.coords<-two.d.array(AKmorph.df$coords)
coord.names<-colnames(two.d.coords)
#

#tests for multivariate normality
#defense (S180588 excluded in 1st line because this is the only 
# individual with the partial plate morph)
(mvn(na.omit(adj.data[which(adj.data$ID!="S180588"),def.trait.list]),mvnTest=c("royston"),multivariatePlot="qq"))
(mvn(na.omit(adj.data[,def.trait.list]),mvnTest=c("royston"),multivariatePlot="qq"))
#swimming
(mvn(na.omit(adj.data[,swi.trait.list]),mvnTest=c("royston"),multivariatePlot="qq"))
#trophic
(mvn(na.omit(adj.data[,tro.trait.list]),mvnTest=c("royston"),multivariatePlot="qq"))
#shape
(mvn(na.omit(two.d.coords),mvnTest=c("royston"),multivariatePlot="qq"))

# makes distance matrices to use as in input for mantel tests
adj.def.dist<-dist(adj.data[def.trait.list])
adj.swi.dist<-dist(adj.data[swi.trait.list])
adj.tro.dist<-dist(adj.data[tro.trait.list])
adj.geo.dist<-dist(merge(two.d.coords,adj.data["ID"],by.x="row.names",by.y="ID",all=T)[,-1])

###### mantel tests, starts with pearson's r, but also uses spearman's rho
# this chunk of code takes a while to run
set.seed(1)
adj.Mantel.def.swim <- mantel(adj.def.dist, adj.swi.dist, method = "pearson", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.def.tro <- mantel(adj.def.dist, adj.tro.dist, method = "pearson", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.def.geo <- mantel(adj.def.dist, adj.geo.dist, method = "pearson", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.def.swim_spear <- mantel(adj.def.dist, adj.swi.dist, method = "spearman", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.def.tro_spear <- mantel(adj.def.dist, adj.tro.dist, method = "spearman", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.def.geo_spear <- mantel(adj.def.dist, adj.geo.dist, method = "spearman", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.swim.tro <- mantel(adj.swi.dist, adj.tro.dist, method = "pearson", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.swim.geo <- mantel(adj.swi.dist, adj.geo.dist, method = "pearson", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.swim.tro_spear  <- mantel(adj.swi.dist, adj.tro.dist, method = "spearman", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.swim.geo_spear  <- mantel(adj.swi.dist, adj.geo.dist, method = "spearman", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.tro.geo <- mantel(adj.tro.dist, adj.geo.dist, method = "pearson", permutations = 499, na.rm = T)
set.seed(1)
adj.Mantel.tro.geo_spear  <- mantel(adj.tro.dist, adj.geo.dist, method = "spearman", permutations = 499, na.rm = T)

adj.Mantel.def.swim
adj.Mantel.def.tro
adj.Mantel.def.geo
adj.Mantel.swim.tro
adj.Mantel.swim.geo
adj.Mantel.tro.geo

adj.Mantel.def.swim_spear
adj.Mantel.def.tro_spear
adj.Mantel.def.geo_spear
adj.Mantel.swim.tro_spear
adj.Mantel.swim.geo_spear
adj.Mantel.tro.geo_spear

#stats table for mantel tests
# pearson
mantel.table.pearson<-data.frame(rbind(c("Mantel.def.swim",round(adj.Mantel.def.swim$statistic,3),adj.Mantel.def.swim$signif),
                                       c("Mantel.def.tro",round(adj.Mantel.def.tro$statistic,3),adj.Mantel.def.tro$signif),
                                       c("Mantel.def.geo",round(adj.Mantel.def.geo$statistic,3),adj.Mantel.def.geo$signif),
                                       c("Mantel.swim.tro",round(adj.Mantel.swim.tro$statistic,3),adj.Mantel.swim.tro$signif),
                                       c("Mantel.swim.geo",round(adj.Mantel.swim.geo$statistic,3),adj.Mantel.swim.geo$signif),
                                       c("Mantel.tro.geo",round(adj.Mantel.tro.geo$statistic,3),adj.Mantel.tro.geo$signif)))
colnames(mantel.table.pearson)<-c("trait suites","r","p")
# spearman
mantel.table.spearman<-data.frame(rbind(c("Mantel.def.swim",round(adj.Mantel.def.swim_spear$statistic,3),adj.Mantel.def.swim_spear$signif),
                                        c("Mantel.def.tro",round(adj.Mantel.def.tro_spear$statistic,3),adj.Mantel.def.tro_spear$signif),
                                        c("Mantel.def.geo",round(adj.Mantel.def.geo_spear$statistic,3),adj.Mantel.def.geo_spear$signif),
                                        c("Mantel.swim.tro",round(adj.Mantel.swim.tro_spear$statistic,3),adj.Mantel.swim.tro_spear$signif),
                                        c("Mantel.swim.geo",round(adj.Mantel.swim.geo_spear$statistic,3),adj.Mantel.swim.geo_spear$signif),
                                        c("Mantel.tro.geo",round(adj.Mantel.tro.geo_spear$statistic,3),adj.Mantel.tro.geo_spear$signif)))
colnames(mantel.table.spearman)<-c("trait suites","rho","p")
# Mantel results 
## TABLE S9
merge(mantel.table.pearson,mantel.table.spearman,by="trait suites")




#   2-BLOCK PLS CORRELATIONS between modules
#
# this small chunk merges scaled and centered traits into 1 df, 
# then splits by lake. Landmark coordinates not centered by units SD, 
# because they are already in the same units.
adj.dat.gpacoords<-merge(adj.data[,c("ID","Lake",def.trait.list,swi.trait.list,tro.trait.list)],
                         two.d.coords,by.x="ID",by.y="row.names",all=T)
split.adj.data<-split(adj.dat.gpacoords,adj.dat.gpacoords$Lake, drop = F)
split.adj.data.Geo<-split.adj.data[names(split.adj.data) %in% c("Walby", "Finger") == FALSE]

Def.Swim.pls<-lapply(split.adj.data,function(x){set.seed(1)
  two.b.pls(x[c(which(complete.cases(x[,c(def.trait.list,swi.trait.list)]))),def.trait.list],
            x[c(which(complete.cases(x[,c(def.trait.list,swi.trait.list)]))),swi.trait.list])})
Def.Tro.pls<-lapply(split.adj.data,function(x){set.seed(1)
  two.b.pls(x[c(which(complete.cases(x[,c(def.trait.list,tro.trait.list)]))),def.trait.list],
            x[c(which(complete.cases(x[,c(def.trait.list,tro.trait.list)]))),tro.trait.list])})
Def.Geo.pls<-lapply(split.adj.data.Geo,function(x){set.seed(1)
  two.b.pls(x[c(which(complete.cases(x[,c(def.trait.list,colnames(two.d.coords))]))),def.trait.list],
            x[c(which(complete.cases(x[,c(def.trait.list,colnames(two.d.coords))]))),c(colnames(two.d.coords))])})

Swim.Tro.pls<-lapply(split.adj.data,function(x){set.seed(1)
  two.b.pls(x[c(which(complete.cases(x[,c(swi.trait.list,tro.trait.list)]))),swi.trait.list],
            x[c(which(complete.cases(x[,c(swi.trait.list,tro.trait.list)]))),tro.trait.list])})
Swim.Geo.pls<-lapply(split.adj.data.Geo,function(x){set.seed(1)
  two.b.pls(x[c(which(complete.cases(x[,c(swi.trait.list,coord.names)]))),swi.trait.list],
            x[c(which(complete.cases(x[,c(swi.trait.list,coord.names)]))),coord.names])})

Geo.Tro.pls<-lapply(split.adj.data.Geo,function(x){set.seed(1)
  two.b.pls(x[c(which(complete.cases(x[,c(coord.names,tro.trait.list)]))),coord.names],
            x[c(which(complete.cases(x[,c(coord.names,tro.trait.list)]))),tro.trait.list])})

#stats tables for PLS
# w/fdr-adjusted p-values
def.swim.stats<-matrix(unlist(lapply(Def.Swim.pls,function(x){c(x[1:3])})),
                       ncol = 3,byrow=T,dimnames = list(names(Def.Swim.pls),c("r","p","Z")))
def.swim.stats[,"p"]<-p.adjust(def.swim.stats[,"p"],method = "fdr")
def.swim.stats<-round(def.swim.stats,3)
def.tro.stats<-matrix(unlist(lapply(Def.Tro.pls,function(x){c(x[1:3])})),
                      ncol = 3,byrow=T,dimnames = list(names(Def.Tro.pls),c("r","p","Z")))
def.tro.stats[,"p"]<-p.adjust(def.tro.stats[,"p"],method="fdr")
def.tro.stats<-round(def.tro.stats,3)
def.geo.stats<-matrix(unlist(lapply(Def.Geo.pls,function(x){c(x[1:3])})),
                      ncol = 3,byrow=T,dimnames = list(names(Def.Geo.pls),c("r","p","Z")))
def.geo.stats[,"p"]<-p.adjust(def.geo.stats[,"p"],method="fdr")
def.geo.stats<-round(def.geo.stats,3)
swim.tro.stats<-matrix(unlist(lapply(Swim.Tro.pls,function(x){c(x[1:3])})),
                       ncol = 3,byrow=T,dimnames = list(names(Swim.Tro.pls),c("r","p","Z")))
swim.tro.stats[,"p"]<-p.adjust(swim.tro.stats[,"p"],method="fdr")
swim.tro.stats<-round(swim.tro.stats,3)
swim.geo.stats<-matrix(unlist(lapply(Swim.Geo.pls,function(x){c(x[1:3])})),
                       ncol = 3,byrow=T,dimnames = list(names(Swim.Geo.pls),c("r","p","Z")))
swim.geo.stats[,"p"]<-p.adjust(swim.geo.stats[,"p"], method="fdr")
swim.geo.stats<-round(swim.geo.stats,3)
geo.tro.stats<-matrix(unlist(lapply(Geo.Tro.pls,function(x){c(x[1:3])})),
                      ncol = 3,byrow=T,dimnames = list(names(Geo.Tro.pls),c("r","p","Z")))
geo.tro.stats[,"p"]<-p.adjust(geo.tro.stats[,"p"], method="fdr")
geo.tro.stats<-round(geo.tro.stats,3)

# data for TABLE 3
write.csv(def.swim.stats,"CORRECTED_def-swim_PLS.csv")
write.csv(def.tro.stats,"CORRECTED_def-tro_PLS.csv")
write.csv(def.geo.stats,"CORRECTED_def-geo_PLS.csv")
write.csv(swim.tro.stats,"CORRECTED_swim-tro_PLS.csv")
write.csv(swim.geo.stats,"CORRECTED_swim-geo_PLS.csv")
write.csv(geo.tro.stats,"CORRECTED_geo-tro_PLS.csv")



#  Function to calculate angle between vectors. Used to compare PLS vectors from different pops
# vec.theta function calculates angle between vectors
# v1 and v2 are vectors
vec.theta<-function(v1, v2){
  acos((v1%*%v2)/(sqrt(v1%*%v1)*sqrt(v2%*%v2)))/(pi/180)
}
# pls.block.theta calculates pairwise angles between pls vectors of different populations
# pls.list is a list of pls objects (one for each population)
# block is "L" or "R" indicating whether angles are being calculated 
#  for the vectors of left or right pls blocks
pls.block.theta<-function(pls.list,block){
  mat=matrix(data=NA, nrow=length(pls.list),ncol=length(pls.list),
             dimnames=list(names(pls.list),names(pls.list)))
  if (block=="L"){block="left.pls.vectors"} 
  if (block=="R") {block="right.pls.vectors"}
  for (i in 1:length(pls.list)){
    for (j in 1:length(pls.list)){
      mat[i,j]<-if(i==j){NA} else (vec.theta((pls.list[[i]][[block]][,1]),
                                             (pls.list[[j]][[block]][,1])))}}
  print(as.matrix(mat))
}
# like eigenvectors, the signs of pls vectors are arbitrary. 
# But the sign of each block is NOT arbitrary with respect to the other block.
# The angle.check function ensures that the sum of the angles between populations A and B 
#  of the left and right block vectors is not greater than 180 degrees. If it is >= 180 degrees,
#  it subtracts both angles from 180, giving the angles as if the sign of both left and right vectors 
#  was reversed without changing their relationship to each other.
angle.check<-function(PW.theta){
  for (i in 1:dim(PW.theta[[1]])[1]){
    for (j in 1:dim(PW.theta[[1]])[2]){
      if (is.na(PW.theta[[1]][i,j])){PW.theta[[1]][i,j]<-NA} else
        if (PW.theta[[1]][i,j]+PW.theta[[2]][i,j]>=180) {
          PW.theta[[1]][i,j]<-180-PW.theta[[1]][i,j]
          PW.theta[[2]][i,j]<-180-PW.theta[[2]][i,j]
        }
    }
  }
  print(PW.theta)
}




# makes matrices of pairwise angles between pls axes of lakes 
# these are used below to make FIGURE S5
Def.Swim.PW.theta<-angle.check(list(L.def=pls.block.theta(Def.Swim.pls,"L"),R.swi=pls.block.theta(Def.Swim.pls,"R")))
Def.Tro.PW.theta<-angle.check(list(L.def=pls.block.theta(Def.Tro.pls,"L"),R.tro=pls.block.theta(Def.Tro.pls,"R")))
Def.Geo.PW.theta<-angle.check(list(L.def=pls.block.theta(Def.Geo.pls,"L"),R.geo=pls.block.theta(Def.Geo.pls,"R")))
Swim.Geo.PW.theta<-angle.check(list(L.swi=pls.block.theta(Swim.Geo.pls,"L"),R.geo=pls.block.theta(Swim.Geo.pls,"R")))
Swim.Tro.PW.theta<-angle.check(list(L.swi=pls.block.theta(Swim.Tro.pls,"L"),R.tro=pls.block.theta(Swim.Tro.pls,"R")))
Geo.Tro.PW.theta<-angle.check(list(L.geo=pls.block.theta(Geo.Tro.pls,"L"),R.tro=pls.block.theta(Geo.Tro.pls,"R")))
x.labs.full<-c("Co.","Ec.","En.","Fin.","G","Je.","Lo.",
               "Ru.","S.R.","Sp.","Te.","Wl.","Wt.","Wik")
x.labs.geo<-c("Co.","Ec.","En.","G","Je.","Lo.",
              "Ru.","S.R.","Sp.","Te.","Wt.","Wik")



#heatmap plots of angle matrices (FIGURE S5 B PANELS)
# pls.angle.plot generates a heatmap, with each cell represents
# an angle of a pairwise comparison for vectors of different lakes,
# and contains that angle value in degrees.
# theta.df is a list of two square matrices showing pairwise angles by block
# block.title is name of block represented by plot
# x.labs vector of population labels for x axis
pls.angle.plot<-function(theta.df,L.or.R,block.title,x.labs){
  dim=if (L.or.R=="L")1 else 2
  ggplot(data = melt(theta.df[[dim]]), aes(x=Var1, y=Var2, fill=value))+geom_tile()+theme_classic()+
    theme(axis.text.x=element_text(size=8))+scale_x_discrete(labels=x.labs)+
    scale_fill_continuous(low="turquoise1",high="firebrick3",limits = c(0, 156),na.value = "steelblue4", 
                          breaks=seq(min(0),max(156),45))+labs(x=NULL, y=NULL,fill="Angle",title=block.title)+
    geom_text(aes(label = round(value, 0)),size=2,color="grey25")
}
# PAIRWISE ANGLE HEATMAPS FOR FIGURE S5 B
# produces warnings because diagonals contain missing values, 
#  and deprecation of melt() in future versions
pls.DefSwi.map<-ggarrange(pls.angle.plot(Def.Swim.PW.theta,"L","Defense",x.labs.full),
                          pls.angle.plot(Def.Swim.PW.theta,"R","Swimming",x.labs.full),legend="none")
pls.DefTro.map<-ggarrange(pls.angle.plot(Def.Tro.PW.theta,"L","Defense",x.labs.full),
                          pls.angle.plot(Def.Tro.PW.theta,"R","Trophic",x.labs.full),legend="none")
pls.DefGeo.map<-ggarrange(pls.angle.plot(Def.Geo.PW.theta,"L","Defense",x.labs.geo),
                          pls.angle.plot(Def.Geo.PW.theta,"R","Shape",x.labs.geo),legend="none")
pls.SwiTro.map<-ggarrange(pls.angle.plot(Swim.Tro.PW.theta,"L","Swimming",x.labs.full),
                          pls.angle.plot(Swim.Tro.PW.theta,"R","Trophic",x.labs.full),legend="none")
pls.SwiGeo.map<-ggarrange(pls.angle.plot(Swim.Geo.PW.theta,"L","Swimming",x.labs.geo),
                          pls.angle.plot(Swim.Geo.PW.theta,"R","Shape",x.labs.geo),legend="none")
pls.GeoTro.map<-ggarrange(pls.angle.plot(Geo.Tro.PW.theta,"L","Shape",x.labs.geo),
                          pls.angle.plot(Geo.Tro.PW.theta,"R","Trophic",x.labs.geo),legend="none")
(pls.DefSwi.map/pls.DefTro.map/pls.DefGeo.map) # FIGURE S5 B1-3
ggsave("CORRECTED_PLSpairwise_theta1.jpg",width = 6,height = 8, dpi = 300)
(pls.SwiTro.map/pls.SwiGeo.map/pls.GeoTro.map) # FIGURE S5 B4-6
ggsave("CORRECTED_PLSpairwise_theta2.jpg",width = 6,height = 8, dpi = 300)

#angle dendrograms. degrees of angle used as distance metric
# FIGURE S5 A PANELS
def.swim.dendro<-((ggdendrogram(hclust(as.dist(Def.Swim.PW.theta$L.def)), rotate = T)+ggtitle("Def x Swi: Defense Block"))/
                    (ggdendrogram(hclust(as.dist(Def.Swim.PW.theta$R.swi)), rotate = T)+ggtitle("Def x Swi: Swimming Block")))
ggsave("CORRECTED_PLSdendro_def.swim.jpg",width = 3.5,height = 6, dpi = 300)
def.tro.dendro<-((ggdendrogram(hclust(as.dist(Def.Tro.PW.theta$L.def)), rotate = T)+ggtitle("Def x Tro: Defense Block"))/
                   (ggdendrogram(hclust(as.dist(Def.Tro.PW.theta$R.tro)), rotate = T)+ggtitle("Def x Tro: Trophic Block")))
ggsave("CORRECTED_PLSdendro_def.tro.jpg",width = 3.5,height = 6, dpi = 300)
def.geo.dendro<-((ggdendrogram(hclust(as.dist(Def.Geo.PW.theta$L.def)), rotate = T)+ggtitle("Def x Shape: Defense Block"))/
                   (ggdendrogram(hclust(as.dist(Def.Geo.PW.theta$R.geo)), rotate = T)+ggtitle("Def x Shape: Shape Block")))
ggsave("CORRECTED_PLSdendro_def.geo.jpg",width = 3.5,height = 6, dpi = 300)
swim.tro.dendro<-((ggdendrogram(hclust(as.dist(Swim.Tro.PW.theta$L.swi)), rotate = T)+ggtitle("Swi x Tro: Swimming Block"))/
                    (ggdendrogram(hclust(as.dist(Swim.Tro.PW.theta$R.tro)), rotate = T)+ggtitle("Swi x Tro: Trophic Block")))
ggsave("CORRECTED_PLSdendro_swim.tro.jpg",width = 3.5,height = 6, dpi = 300)
swim.geo.dendro<-((ggdendrogram(hclust(as.dist(Swim.Geo.PW.theta$L.swi)), rotate = T)+ggtitle("Swi x Shape: Swimming Block"))/
                    (ggdendrogram(hclust(as.dist(Swim.Geo.PW.theta$R.geo)), rotate = T)+ggtitle("Swi x Shape: Shape Block")))
ggsave("CORRECTED_PLSdendro_swim.geo.jpg",width = 3.5,height = 6, dpi = 300)
geo.tro.dendro<-((ggdendrogram(hclust(as.dist(Geo.Tro.PW.theta$L.geo)), rotate = T)+ggtitle("Shape x Tro: Shape Block"))/
                   (ggdendrogram(hclust(as.dist(Geo.Tro.PW.theta$R.tro)), rotate = T)+ggtitle("Shape x Tro: Trophic Block")))
ggsave("CORRECTED_PLSdendro_geo.tro.jpg",width = 3.5,height = 6, dpi = 300)


# heatmap plots of effect sizes of pls effect differences between lakes
# pls.Z.plot function creates heatmap
# melted is df of "melted" effect size matrix
# x.labs is vector of population labels for x axis
# plot.title is name of plot
pls.Z.plot<-function(melted,x.labs,plot.title){ggplot(data = melted, aes(x=Var1, y=Var2, fill=value))+geom_tile()+theme_classic()+
    theme(axis.text.x=element_text(size=8))+scale_x_discrete(labels=x.labs)+
    scale_fill_continuous(low="steelblue1",high="darkgreen",limits = c(0, 4.1),na.value = "palegreen3", 
                          breaks=seq(min(0),max(10),1))+labs(x=NULL, y=NULL,fill="Z-Score",title=plot.title)+
    geom_text(aes(label = round(value, 2)),size=1.8,color="grey10")
}
set.seed(1)
Def.Swim.comp.pls<-compare.pls(Def.Swim.pls)
Def.Swim.meltedZ<-reshape2::melt(Def.Swim.comp.pls$pairwise.z)
Def.Swim.meltedZ[which(Def.Swim.meltedZ$value==0),"value"]<-NA
def.x.swi_plsZ<-pls.Z.plot(Def.Swim.meltedZ,x.labs.full,"2B-PLS Def. x Swi.")

set.seed(1)
Def.Tro.comp.pls<-compare.pls(Def.Tro.pls)
Def.Tro.meltedZ<-reshape2::melt(Def.Tro.comp.pls$pairwise.z)
Def.Tro.meltedZ[which(Def.Tro.meltedZ$value==0),"value"]<-NA
def.x.tro_plsZ<-pls.Z.plot(Def.Tro.meltedZ,x.labs.full,"2B-PLS Def. x Tro.")

set.seed(1)
Def.Geo.comp.pls<-compare.pls(Def.Geo.pls)
Def.Geo.meltedZ<-reshape2::melt(Def.Geo.comp.pls$pairwise.z)
Def.Geo.meltedZ[which(Def.Geo.meltedZ$value==0),"value"]<-NA
def.x.geo_plsZ<-pls.Z.plot(Def.Geo.meltedZ,x.labs.geo,"2B-PLS Def. x Shape")

set.seed(1)
Swim.Geo.comp.pls<-compare.pls(Swim.Geo.pls)
Swim.Geo.meltedZ<-reshape2::melt(Swim.Geo.comp.pls$pairwise.z)
Swim.Geo.meltedZ[which(Swim.Geo.meltedZ$value==0),"value"]<-NA
swi.x.geo_plsZ<-pls.Z.plot(Swim.Geo.meltedZ,x.labs.geo,"2B-PLS Swi. x Shape")

set.seed(1)
Swim.Tro.comp.pls<-compare.pls(Swim.Tro.pls)
Swim.Tro.meltedZ<-reshape2::melt(Swim.Tro.comp.pls$pairwise.z)
Swim.Tro.meltedZ[which(Swim.Tro.meltedZ$value==0),"value"]<-NA
swi.x.tro_plsZ<-pls.Z.plot(Swim.Tro.meltedZ,x.labs.full,"2B-PLS Swi. x Tro.")

set.seed(1)
Geo.Tro.comp.pls<-compare.pls(Geo.Tro.pls)
Geo.Tro.meltedZ<-reshape2::melt(Geo.Tro.comp.pls$pairwise.z)
Geo.Tro.meltedZ[which(Geo.Tro.meltedZ$value==0),"value"]<-NA
geo.x.tro_plsZ<-pls.Z.plot(Geo.Tro.meltedZ,x.labs.geo,"2B-PLS Tro. x Shape")

# compound plot of effect size comparisons
# FIGURE S4
ggarrange(def.x.swi_plsZ,def.x.tro_plsZ,def.x.geo_plsZ,
          swi.x.geo_plsZ,swi.x.tro_plsZ,NULL,
          geo.x.tro_plsZ, NULL, NULL,ncol=3, nrow=3,common.legend = T)
ggsave("CORRECTED_PLSpairwiseZ.jpg",width = 10.5,height = 11, dpi = 300)

