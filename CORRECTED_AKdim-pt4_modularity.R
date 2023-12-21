# Preps dataframe for modularity tests within populations
# landmarks scaled by sd to avoid mixing covariance and correlation between modules in modularity tests
adj.geo<-data.frame("1.X"=((two.d.coords[,"1.X"]*1-mean(na.omit(two.d.coords[,"1.X"])))/sd(na.omit(two.d.coords[,"1.X"]))),
                    "1.Y"=((two.d.coords[,"1.Y"]*1-mean(na.omit(two.d.coords[,"1.Y"])))/sd(na.omit(two.d.coords[,"1.Y"]))),
                    "2.X"=((two.d.coords[,"2.X"]*1-mean(na.omit(two.d.coords[,"2.X"])))/sd(na.omit(two.d.coords[,"2.X"]))),
                    "2.Y"=((two.d.coords[,"2.Y"]*1-mean(na.omit(two.d.coords[,"2.Y"])))/sd(na.omit(two.d.coords[,"2.Y"]))),
                    "3.X"=((two.d.coords[,"3.X"]*1-mean(na.omit(two.d.coords[,"3.X"])))/sd(na.omit(two.d.coords[,"3.X"]))),
                    "3.Y"=((two.d.coords[,"3.Y"]*1-mean(na.omit(two.d.coords[,"3.Y"])))/sd(na.omit(two.d.coords[,"3.Y"]))),
                    "4.X"=((two.d.coords[,"4.X"]*1-mean(na.omit(two.d.coords[,"4.X"])))/sd(na.omit(two.d.coords[,"4.X"]))),
                    "4.Y"=((two.d.coords[,"4.Y"]*1-mean(na.omit(two.d.coords[,"4.Y"])))/sd(na.omit(two.d.coords[,"4.Y"]))),
                    "5.X"=((two.d.coords[,"5.X"]*1-mean(na.omit(two.d.coords[,"5.X"])))/sd(na.omit(two.d.coords[,"5.X"]))),
                    "5.Y"=((two.d.coords[,"5.Y"]*1-mean(na.omit(two.d.coords[,"5.Y"])))/sd(na.omit(two.d.coords[,"5.Y"]))),
                    "6.X"=((two.d.coords[,"6.X"]*1-mean(na.omit(two.d.coords[,"6.X"])))/sd(na.omit(two.d.coords[,"6.X"]))),
                    "6.Y"=((two.d.coords[,"6.Y"]*1-mean(na.omit(two.d.coords[,"6.Y"])))/sd(na.omit(two.d.coords[,"6.Y"]))),
                    "7.X"=((two.d.coords[,"7.X"]*1-mean(na.omit(two.d.coords[,"7.X"])))/sd(na.omit(two.d.coords[,"7.X"]))),
                    "7.Y"=((two.d.coords[,"7.Y"]*1-mean(na.omit(two.d.coords[,"7.Y"])))/sd(na.omit(two.d.coords[,"7.Y"]))),
                    "8.X"=((two.d.coords[,"8.X"]*1-mean(na.omit(two.d.coords[,"8.X"])))/sd(na.omit(two.d.coords[,"8.X"]))),
                    "8.Y"=((two.d.coords[,"8.Y"]*1-mean(na.omit(two.d.coords[,"8.Y"])))/sd(na.omit(two.d.coords[,"8.Y"]))),
                    "9.X"=((two.d.coords[,"9.X"]*1-mean(na.omit(two.d.coords[,"9.X"])))/sd(na.omit(two.d.coords[,"9.X"]))),
                    "9.Y"=((two.d.coords[,"9.Y"]*1-mean(na.omit(two.d.coords[,"9.Y"])))/sd(na.omit(two.d.coords[,"9.Y"]))),
                    "10.X"=((two.d.coords[,"10.X"]*1-mean(na.omit(two.d.coords[,"10.X"])))/sd(na.omit(two.d.coords[,"10.X"]))),
                    "10.Y"=((two.d.coords[,"10.Y"]*1-mean(na.omit(two.d.coords[,"10.Y"])))/sd(na.omit(two.d.coords[,"10.Y"]))),
                    "11.X"=((two.d.coords[,"11.X"]*1-mean(na.omit(two.d.coords[,"11.X"])))/sd(na.omit(two.d.coords[,"11.X"]))),
                    "11.Y"=((two.d.coords[,"11.Y"]*1-mean(na.omit(two.d.coords[,"11.Y"])))/sd(na.omit(two.d.coords[,"11.Y"]))),
                    "12.X"=((two.d.coords[,"12.X"]*1-mean(na.omit(two.d.coords[,"12.X"])))/sd(na.omit(two.d.coords[,"12.X"]))),
                    "12.Y"=((two.d.coords[,"12.Y"]*1-mean(na.omit(two.d.coords[,"12.Y"])))/sd(na.omit(two.d.coords[,"12.Y"]))),
                    "13.X"=((two.d.coords[,"13.X"]*1-mean(na.omit(two.d.coords[,"13.X"])))/sd(na.omit(two.d.coords[,"13.X"]))),
                    "13.Y"=((two.d.coords[,"13.Y"]*1-mean(na.omit(two.d.coords[,"13.Y"])))/sd(na.omit(two.d.coords[,"13.Y"]))),
                    "14.X"=((two.d.coords[,"14.X"]*1-mean(na.omit(two.d.coords[,"14.X"])))/sd(na.omit(two.d.coords[,"14.X"]))),
                    "14.Y"=((two.d.coords[,"14.Y"]*1-mean(na.omit(two.d.coords[,"14.Y"])))/sd(na.omit(two.d.coords[,"14.Y"]))),
                    "15.X"=((two.d.coords[,"15.X"]*1-mean(na.omit(two.d.coords[,"15.X"])))/sd(na.omit(two.d.coords[,"15.X"]))),
                    "15.Y"=((two.d.coords[,"15.Y"]*1-mean(na.omit(two.d.coords[,"15.Y"])))/sd(na.omit(two.d.coords[,"15.Y"]))),
                    "16.X"=((two.d.coords[,"16.X"]*1-mean(na.omit(two.d.coords[,"16.X"])))/sd(na.omit(two.d.coords[,"16.X"]))),
                    "16.Y"=((two.d.coords[,"16.Y"]*1-mean(na.omit(two.d.coords[,"16.Y"])))/sd(na.omit(two.d.coords[,"16.Y"]))),
                    "17.X"=((two.d.coords[,"17.X"]*1-mean(na.omit(two.d.coords[,"17.X"])))/sd(na.omit(two.d.coords[,"17.X"]))),
                    "17.Y"=((two.d.coords[,"17.Y"]*1-mean(na.omit(two.d.coords[,"17.Y"])))/sd(na.omit(two.d.coords[,"17.Y"]))),
                    "18.X"=((two.d.coords[,"18.X"]*1-mean(na.omit(two.d.coords[,"18.X"])))/sd(na.omit(two.d.coords[,"18.X"]))),
                    "18.Y"=((two.d.coords[,"18.Y"]*1-mean(na.omit(two.d.coords[,"18.Y"])))/sd(na.omit(two.d.coords[,"18.Y"]))),
                    "19.X"=((two.d.coords[,"19.X"]*1-mean(na.omit(two.d.coords[,"19.X"])))/sd(na.omit(two.d.coords[,"19.X"]))),
                    "19.Y"=((two.d.coords[,"19.Y"]*1-mean(na.omit(two.d.coords[,"19.Y"])))/sd(na.omit(two.d.coords[,"19.Y"]))),
                    "20.X"=((two.d.coords[,"20.X"]*1-mean(na.omit(two.d.coords[,"20.X"])))/sd(na.omit(two.d.coords[,"20.X"]))),
                    "20.Y"=((two.d.coords[,"20.Y"]*1-mean(na.omit(two.d.coords[,"20.Y"])))/sd(na.omit(two.d.coords[,"20.Y"]))),
                    "21.X"=((two.d.coords[,"21.X"]*1-mean(na.omit(two.d.coords[,"21.X"])))/sd(na.omit(two.d.coords[,"21.X"]))),
                    "21.Y"=((two.d.coords[,"21.Y"]*1-mean(na.omit(two.d.coords[,"21.Y"])))/sd(na.omit(two.d.coords[,"21.Y"]))))
colnames(adj.geo)<-coord.names
adj.data.geo<-merge(adj.data,adj.geo ,by.x ="ID",by.y="row.names",all=T,sort = F)
row.names(adj.data.geo)<-adj.data.geo$ID



#  Q3) Modularity tests
######
#splits trait data by lake
split.adj.data<-split(adj.data.geo,adj.data.geo$Lake, drop = F)
#splits trait data by lake including shape 
split.adj.data.Geo<-split.adj.data[names(split.adj.data) %in% c("Walby", "Finger") == FALSE]



# W gill rakers, W/O shape
split.adj.data.red<-lapply(split.adj.data, function(x){na.omit(x[,c(3:18)])})
# This bit separates each trait suite into a distinct module
## TABLE S11
mod.tot<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','b','b','b','c','c','c','c','c','c','c','c'),CI=T,iter=4999)
})
# The following lines test modularity by separating trait suites in a "one-against-all" fashion
## TABLE S12
mod.def<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','b','b','b','b','b','b','b','b','b','b','b'),CI=T,iter=4999)
})
mod.swi<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','b','b','b','a','a','a','a','a','a','a','a'),CI=T,iter=4999)
})
mod.tro<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','a','a','a','c','c','c','c','c','c','c','c'),CI=T,iter=4999)
})
mod.tot.results<-matrix(unlist(lapply(mod.tot,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.def.results<-matrix(unlist(lapply(mod.def,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.swi.results<-matrix(unlist(lapply(mod.swi,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.tro.results<-matrix(unlist(lapply(mod.tro,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
dimnames(mod.tot.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))### USED IN MAIN TEXT
mod.tot.results[,"p"]<-p.adjust(mod.tot.results[,"p"], method="fdr")
mod.tot.results<-round(mod.tot.results,3)
dimnames(mod.def.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))
mod.def.results[,"p"]<-p.adjust(mod.def.results[,"p"], method="fdr")
mod.def.results<-round(mod.def.results,3)
dimnames(mod.swi.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))
mod.swi.results[,"p"]<-p.adjust(mod.swi.results[,"p"], method="fdr")
mod.swi.results<-round(mod.swi.results,3)
dimnames(mod.tro.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))
mod.tro.results[,"p"]<-p.adjust(mod.tro.results[,"p"], method="fdr")
mod.tro.results<-round(mod.tro.results,3)

#Table S10
write.csv(mod.tot.results,"CORRECTED-def-swi-tro_Mod.csv")
#Table 4
write.csv(mod.def.results, "CORRECTED-def_Mod.csv")
write.csv(mod.swi.results, "CORRECTED-swi_Mod.csv")
write.csv(mod.tro.results, "CORRECTED-tro_Mod.csv")



# With Shape
split.adj.data.Geo.red<-lapply(split.adj.data.Geo,function(x){na.omit(x)[,-c(1:2)]})#removes ID & Lake cols
# This bit separates each trait suite into a distinct module
## TABLE S11
mod.Geo.tot<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','b','b','b','c','c','c','c','c','c','c','c',rep('d',42)),CI=T,iter=4999)
})
# The following lines test modularity by separating trait suites in a "one-against-all" fashion
## TABLE S12
mod.Geo.def<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','b','b','b','b','b','b','b','b','b','b','b',rep('b',42)),CI=T,iter=4999)
})
mod.Geo.swi<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','b','b','b','a','a','a','a','a','a','a','a',rep('a',42)),CI=T,iter=4999)
})
mod.Geo.tro<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','a','a','a','c','c','c','c','c','c','c','c',rep('a',42)),CI=T,iter=4999)
})
mod.Geo.geo<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x,c('a','a','a','a','a','a','a','a','a','a','a','a','a','a','a','a',rep('b',42)),CI=T,iter=4999)
})

# the following chunk puts results from above tests into tables
mod.Geo.tot.results<-matrix(unlist(lapply(mod.Geo.tot,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.def.results<-matrix(unlist(lapply(mod.Geo.def,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.swi.results<-matrix(unlist(lapply(mod.Geo.swi,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.tro.results<-matrix(unlist(lapply(mod.Geo.tro,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.geo.results<-matrix(unlist(lapply(mod.Geo.geo,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
dimnames(mod.Geo.tot.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
mod.Geo.tot.results[,"p"]<-p.adjust(mod.Geo.tot.results[,"p"], method="fdr")
mod.Geo.tot.results<-round(mod.Geo.tot.results,3)
dimnames(mod.Geo.def.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
mod.Geo.def.results[,"p"]<-p.adjust(mod.Geo.def.results[,"p"], method="fdr")
mod.Geo.def.results<-round(mod.Geo.def.results,3)
dimnames(mod.Geo.swi.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
mod.Geo.swi.results[,"p"]<-p.adjust(mod.Geo.swi.results[,"p"], method="fdr")
mod.Geo.swi.results<-round(mod.Geo.swi.results,3)
dimnames(mod.Geo.tro.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
mod.Geo.tro.results[,"p"]<-p.adjust(mod.Geo.tro.results[,"p"], method="fdr")
mod.Geo.tro.results<-round(mod.Geo.tro.results,3)
dimnames(mod.Geo.geo.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
mod.Geo.geo.results[,"p"]<-p.adjust(mod.Geo.geo.results[,"p"], method="fdr")
mod.Geo.geo.results<-round(mod.Geo.geo.results,3)

#Table S11
write.csv(mod.Geo.tot.results,"CORRECTED-def-swi-tro-geo_Mod.csv")
#Table S12
write.csv(mod.Geo.def.results,"CORRECTED-def-w.geo_Mod.csv")
write.csv(mod.Geo.swi.results,"CORRECTED-swi-w-geo_Mod.csv")
write.csv(mod.Geo.tro.results,"CORRECTED-tro-w-geo_Mod.csv")
write.csv(mod.Geo.geo.results,"CORRECTED-geo-w-geo_Mod.csv")





# W/O gill rakers or shape 
# Included for purpose of comparison, but not in any tables
mod.noGR.tot<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','b','b','b','c','c','c','c','c','c','c'),CI=T,iter=4999)
})
mod.noGR.def<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','b','b','b','b','b','b','b','b','b','b'),CI=T,iter=4999)
})
mod.noGR.swi<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','b','b','b','a','a','a','a','a','a','a'),CI=T,iter=4999)
})
mod.noGR.tro<-lapply(split.adj.data.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','a','a','a','c','c','c','c','c','c','c'),CI=T,iter=4999)
})
mod.noGR.tot.results<-matrix(unlist(lapply(mod.noGR.tot,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.noGR.def.results<-matrix(unlist(lapply(mod.noGR.def,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.noGR.swi.results<-matrix(unlist(lapply(mod.noGR.swi,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.noGR.tro.results<-matrix(unlist(lapply(mod.noGR.tro,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
dimnames(mod.noGR.tot.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))
dimnames(mod.noGR.def.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))
dimnames(mod.noGR.swi.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))
dimnames(mod.noGR.tro.results)<-list(levels(adj.data$Lake),c("CR","lowerCI","upperCI","p","Z"))



# W Shape, W/O gill rakers
# Included for purpose of comparison, but not in any tables
split.adj.data.Geo.red<-lapply(split.adj.data.Geo,function(x){na.omit(x)[,-c(1:2)]})
mod.Geo.noGR.tot<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','b','b','b','c','c','c','c','c','c','c',rep('d',42)),CI=T,iter=4999)
})
mod.Geo.noGR.def<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','b','b','b','b','b','b','b','b','b','b',rep('b',42)),CI=T,iter=4999)
})
mod.Geo.noGR.swi<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','b','b','b','a','a','a','a','a','a','a',rep('a',42)),CI=T,iter=4999)
})
mod.Geo.noGR.tro<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','a','a','a','c','c','c','c','c','c','c',rep('a',42)),CI=T,iter=4999)
})
mod.Geo.noGR.geo<-lapply(split.adj.data.Geo.red,function(x){set.seed(1)
  modularity.test(x[,-12],c('a','a','a','a','a','a','a','a','a','a','a','a','a','a','a',rep('b',42)),CI=T,iter=4999)
})
mod.Geo.noGR.tot.results<-matrix(unlist(lapply(mod.Geo.noGR.tot,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.noGR.def.results<-matrix(unlist(lapply(mod.Geo.noGR.def,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.noGR.swi.results<-matrix(unlist(lapply(mod.Geo.noGR.swi,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.noGR.tro.results<-matrix(unlist(lapply(mod.Geo.noGR.tro,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
mod.Geo.noGR.geo.results<-matrix(unlist(lapply(mod.Geo.noGR.geo,function(x){c(round(x$CR,3),round(x$CInterval[1],3),round(x$CInterval[2],3),x$P.value,round(x$Z,3))})),ncol = 5,byrow = T)
dimnames(mod.Geo.noGR.tot.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
dimnames(mod.Geo.noGR.def.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
dimnames(mod.Geo.noGR.swi.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
dimnames(mod.Geo.noGR.tro.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))
dimnames(mod.Geo.noGR.geo.results)<-list(levels(adj.data$Lake)[-c(4,12)],c("CR","lowerCI","upperCI","p","Z"))



#
#