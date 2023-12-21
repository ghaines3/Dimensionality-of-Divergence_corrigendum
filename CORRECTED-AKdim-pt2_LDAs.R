# Files required for body shape geometric morphometrics
#raw landmarks (semilandmark curve info supplied separately)
AKmorph <- readland.tps("Edited_8_semi.TPS", specID="imageID")
#slider file
slider <- as.matrix(read.csv("slider_8_semi.csv"))


# LD frequency plot function - plots freq. distribution of LD1 for each lake
# data= dataframe containing Lake vector and LD1 scores
# X = LD1 vector (e.g., data$LD1)
# Lake = Lake vector (e.g., data$Lake)
# Lakelab = title for lake legend ("Lake")
# title = plot title (as character string)
# LD.prop = % of trace explained by LD1 (as character string)

LD.freq<-function(data,X,Lake,Lakelab,title,LD.prop){data %>%
    mutate(Lake = fct_reorder(Lake, X, .fun='mean' )) %>%ggplot( aes(x = X, y = Lake, fill = ..x..)) + 
    scale_fill_viridis(option = "C") +
    xlab(paste("LD1",LD.prop))+ylab(Lakelab)+
    geom_density_ridges_gradient(scale = 1.85, rel_min_height = 0.01,jittered_points = TRUE, 
                                 point_shape = "|", point_size = 2, size = 0.25,
                                 position = position_points_jitter(height = 0)) +
    labs(title = title)+
    theme_minimal_hgrid()+theme(legend.position="none")
}

# function threeD.scatter makes a 3d scatter plot (for exploratory purposes, 
#  not included in paper)
# data = a dataframe including groupinf factor (e.g., "Lake") 
#   and at least 3 vectors of PC scores, LD scores, etc.
# anchor.axis = axis droplines should run along
# anchor.val = value along anchor.axis at which to start droplines
# col.vec = a vector of color names to apply to groups
# axis.index = vector of column indices of "data" that contain X, Y, and Z values
# axis.labs = X, Y, and Z axis labels
threeD.scatter<-function(data,anchor.axis,anchor.val,col.vec,axis.index,axis.labs){
  dat.2 <- replicate(2, data, simplify = F)
  dat.2[[2]][anchor.axis] <- anchor.val
  dat.3 <- dat.2 %>%bind_rows() %>%group2NA("FishID", "Lake")
  
  plot_ly(color = ~factor(Lake), showlegend = T, 
          colors = col.vec) %>%
    add_markers(data = data, x = ~data[,axis.index[1]], y = ~data[,axis.index[2]], z = ~data[,axis.index[3]],size = 100) %>%
    add_paths(data = dat.3,  x = ~dat.3[,axis.index[1]], y = ~dat.3[,axis.index[2]], z = ~dat.3[,axis.index[3]], opacity = 0.1,showlegend =F)%>%
    layout(scene=list(xaxis=list(title = axis.labs[[1]]),yaxis=list(title = axis.labs[[2]]), zaxis=list(title= axis.labs[[3]])))
  
}

# data = a dataframe including groupinf factor (e.g., "Lake") 
#   and at least 3 vectors of PC scores, LD scores, etc.
# axis.index = vector of column indices of "data" that contain X, Y, and Z values
# col.vec = a vector of color names to apply to groups
# vector.sum.df = dataframe summarizing magnitude and loadings of LDA vectors
#   row #2 is proportion of trace
# axis.labs = vector of X and Y axis names, pasted to proportion of trace for axis label,
#   (e.g., c("LD1 - ","LD2 - "))
# which.vec = indices of columns representing vectors on X and Y axes in vector.sum.df
#   (e.g., for X = LD1 and Y = LD2, c(1,2), for X = LD1 and Y = LD3, c(1,3))
twoD.scatter<-function(data,axis.index,col.vec,vector.sum.df,axis.labs,which.vec){
  ggplot(data = data, aes(x=data[,axis.index[1]],y=data[,axis.index[2]],color=Lake))+
    geom_point(alpha=.6,stroke=0,size=1.5)+theme_classic()+
    scale_color_manual(values=col.vec)+coord_fixed()+
    guides(color = guide_legend(override.aes = list(size = 4, alpha=1)))+
    labs(x=paste(axis.labs[1],round(vector.sum.df[2,which.vec[1]]/sum(vector.sum.df[2,])*100,1),"%"),
         y=paste(axis.labs[2],round(vector.sum.df[2,which.vec[2]]/sum(vector.sum.df[2,])*100,1),"%"))
}


### Maks a matrix of effective dimensionality calculations
dim.stats<-function(lda,euc.mat,lake.means){
  data.frame(dim.SVsq=estimate.ED.eig(prop.trace(lda$svd)),
             dim.matrix=unlist(estimate.ED(na.omit(euc.mat))),
             pop.mean.matrix=unlist(estimate.ED(subset(lake.means,select=-Lake))))
}

#LDAs#
#########
######    
#########

#### Defense ####
#Defense data frame#
Defense.df <- cbind.data.frame(ID=data$FishID,Lake= data$Lake, adj.PS_mean = data$adj.PS.mean, 
                               LatPlNum_mean = rowMeans(subset(data, select = c("LatPlNum_R", "LatPlNum_L")), na.rm = TRUE), 
                               MidPl_mean=rowMeans(subset(data, select = c("FirstPl_L", "LastPl_L","FirstPl_R","LastPl_R")), na.rm = TRUE), 
                               adj.DS1L= data$adj.DS1L,adj.DS2L = data$adj.DS2L)
rownames(Defense.df)<-Defense.df$ID


# LDAprep.A) This chunk splits Defense.df by Lake, then interpolates missing trait values
#   within lake, unlists the interpolated vector, and replaces the 
#   old vector with the interpolated one (Defense.df is already sorted by lake).
#   PS_mean and LatPlNum_mean are excluded from this because 
#   they are not continuous measurements.
Defense.df.split<-split(Defense.df,Defense.df$Lake, drop = F) 
NEW.adj.DS1L<-unlist(lapply(Defense.df.split,function(x){
  interp(data=x,y=x[,"adj.DS1L"],a=x[,"adj.PS_mean"],b=x[,"adj.DS2L"])$SecondaryNew.Vec}))
NEW.adj.DS2L<-unlist(lapply(Defense.df.split,function(x){
  interp(data=x,b=x[,"adj.DS1L"],a=x[,"adj.PS_mean"],y=x[,"adj.DS2L"])$SecondaryNew.Vec}))
NEW.adj.PS_mean<-unlist(lapply(Defense.df.split,function(x){
  if(var(na.omit(x$adj.PS_mean))==0)x$adj.PS_mean else(
    interp(data=x,b=x[,"adj.DS1L"],y=x[,"adj.PS_mean"],a=x[,"adj.DS2L"])$SecondaryNew.Vec)}))
Defense.df<-bind_rows(Defense.df.split) #this line recombines the list of dataframes for each 
#    lake into a single dataframe
Defense.df[,"adj.DS1L"]<-NEW.adj.DS1L
Defense.df[,"adj.DS2L"]<-NEW.adj.DS2L
Defense.df[,"adj.PS_mean"]<-NEW.adj.PS_mean
rownames(Defense.df)<-Defense.df$ID

# LDAprep.B) scaling to sd and mean-centering of traits. This is required to put the 
#   trait measurements into comparable units
Defense.df.euc<-data.frame(Defense.df["ID"],Defense.df["Lake"],
                           LatPlNum_mean=((Defense.df[,"LatPlNum_mean"]*1-mean(na.omit(Defense.df[,"LatPlNum_mean"])))/sd(na.omit(Defense.df[,"LatPlNum_mean"]))),
                           adj.DS1L=((Defense.df[,"adj.DS1L"]*1-mean(na.omit(Defense.df[,"adj.DS1L"])))/sd(na.omit(Defense.df[,"adj.DS1L"]))),
                           adj.DS2L=((Defense.df[,"adj.DS2L"]*1-mean(na.omit(Defense.df[,"adj.DS2L"])))/sd(na.omit(Defense.df[,"adj.DS2L"]))),
                           adj.PS_mean=((Defense.df[,"adj.PS_mean"]*1-mean(na.omit(Defense.df[,"adj.PS_mean"])))/sd(na.omit(Defense.df[,"adj.PS_mean"]))),
                           MidPl_mean=((Defense.df[,"MidPl_mean"]*1-mean(na.omit(Defense.df[,"MidPl_mean"])))/sd(na.omit(Defense.df[,"MidPl_mean"]))))
# makes numeric matrix for dimensionality calculation based on eigendecomp. (by individuals, then means)
Defense.df.euc.mat<-Defense.df.euc[,c("adj.PS_mean" ,"LatPlNum_mean","MidPl_mean","adj.DS1L", "adj.DS2L")]
def.means<-aggregate(data=na.omit(Defense.df.euc[,c("Lake","LatPlNum_mean","adj.DS1L","adj.DS2L", "adj.PS_mean", "MidPl_mean")]), .~Lake, FUN=mean)

########## LDA
Defense.lda<-lda(Lake ~ adj.PS_mean + LatPlNum_mean+
                   MidPl_mean+  adj.DS1L + adj.DS2L , data = Defense.df.euc, CV = F)
# yields matrices of lake assignments, probabilities, and LD scores for each individual
predictions.Defense <- predict(Defense.lda) 
# Dataframe of LD scoes fo individual fish
LD.def.df<-data.frame(FishID=na.omit(Defense.df.euc)$ID,Lake=na.omit(
  Defense.df.euc)$Lake, defense.LD1=predictions.Defense$x[,1], 
  defense.LD2=predictions.Defense$x[,2], defense.LD3=predictions.Defense$x[,3]) 
# LDA dimension summary 
## ROWS 3-7 IN TABLE S5
(lda.def.df<-round(lda.vector.df(Defense.lda),3))%>%write.csv("CORRECTED-lda-def_df.csv")
# Effective dimensions (Defense)
## PART OF TABLE 2
defense.dim.matrix<-round(dim.stats(Defense.lda,Defense.df.euc.mat,def.means),2)

### Defense LDA dim. Plots
# Defense LD1 frequency plot 
## FIGURE 5A
freq.def<-LD.freq(LD.def.df,LD.def.df$defense.LD1,LD.def.df$Lake,
                  "Lake",'LD1 of Defensive Traits by Lake',"89.2%")
# 3D scatter (EXPLORATORY PLOT, not in paper)
threeD.scatter(LD.def.df,"defense.LD3",-4,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                            "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                            "sienna4","grey80","grey50"),c(3,4,5),c('LD1-89.2%', 'LD2-8.1%','LD3-1.5%'))
# 2D scatter 
## FIGURE1 S1, A1-A2
def.plot.12<-twoD.scatter(LD.def.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "sienna4","grey80","grey50"),lda.def.df,c("LD1 - ","LD2 - "),c(1,2))
def.plot.13<-twoD.scatter(LD.def.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "sienna4","grey80","grey50"),lda.def.df,c("LD1 - ","LD3 - "),c(1,3))
def.plot.12/def.plot.13+ plot_layout(guides = 'collect')


## Same as above, but without G and Echo, 
# to demonstrate results without these 2 strongly divergent pops.
# makes numeric matrix for dimensionality calculation based on eigendecomp. (by individuals, then means)
Defense.df.euc.NoGEcho<-filter(Defense.df.euc,Lake!="G"&Lake!="Echo")
Defense.df.euc.mat.NoGEcho<-Defense.df.euc.NoGEcho[,c("adj.PS_mean" ,"LatPlNum_mean","MidPl_mean","adj.DS1L", "adj.DS2L")]
def.means.NoGEcho<-aggregate(data=na.omit(Defense.df.euc.NoGEcho[,c("Lake","LatPlNum_mean","adj.DS1L","adj.DS2L", "adj.PS_mean", "MidPl_mean")]), .~Lake, FUN=mean)

Defense.lda.NoGEcho<-lda(Lake ~ adj.PS_mean + LatPlNum_mean+
                           MidPl_mean+  adj.DS1L + adj.DS2L , data = Defense.df.euc.NoGEcho, CV = F)
# yields matrices of lake assignments, probabilities, and LD scores for each individual
predictions.Defense.NoGEcho <- predict(Defense.lda.NoGEcho) 
# Dataframe of LD scoes fo individual fish
LD.def.df.NoGEcho<-data.frame(FishID=na.omit(Defense.df.euc.NoGEcho)$ID,Lake=na.omit(
  Defense.df.euc.NoGEcho)$Lake, defense.LD1=predictions.Defense.NoGEcho$x[,1], 
  defense.LD2=predictions.Defense.NoGEcho$x[,2], defense.LD3=predictions.Defense.NoGEcho$x[,3]) 
# LDA dimension summary 
## ROWS 3-7 IN TABLE S5
(lda.def.df.NoGEcho<-round(lda.vector.df(Defense.lda.NoGEcho),3))
# Effective dimensions (Defense)
## not part of Table 2, but mentioned in discussion
defense.dim.matrix.NoGEcho<-round(dim.stats(Defense.lda.NoGEcho,Defense.df.euc.mat.NoGEcho,def.means.NoGEcho),2)

# 3D scatter (EXPLORATORY PLOT, not in paper)
threeD.scatter(LD.def.df.NoGEcho,"defense.LD3",-4,c("steelblue4","palegreen4","palegreen2",
                                                    "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                                    "sienna4","grey80","grey50"),c(3,4,5),c('LD1-78.7%', 'LD2-11.9%','LD3-5.9%'))



##################
#### Swimming ####

####LDA Swimming Trait####
#Swimming data frame#
Swimming.df <- cbind.data.frame(ID=data$FishID,Lake= data$Lake, 
                                SL_mm = log(data$SL_mm), 
                                adj.BD= data$adj.BD, 
                                adj.CP= data$adj.CP)
rownames(Swimming.df)<-Swimming.df$ID

# LDAprep.A) # SL data not interpolated because all 
#   individuals have SL data (and traits size corrected based on SL.
#   Body depth not interpolated because all 
#   individuals with SL data also have BD data
Swimming.df.split<-split(Swimming.df,Swimming.df$Lake, drop = F)
NEW.adj.CP<-unlist(lapply(Swimming.df.split,function(x){
  interp(data=x,y=x[,"adj.CP"],a=x[,"SL_mm"],b=x[,"adj.BD"])$SecondaryNew.Vec}))
Swimming.df<-bind_rows(Swimming.df.split)
Swimming.df[,"adj.CP"]<-NEW.adj.CP 
rownames(Swimming.df)<-Swimming.df$ID

# LDAprep.B) scaling to sd and mean-centering of traits
Swimming.df.euc<-data.frame(ID=Swimming.df["ID"],Lake=Swimming.df["Lake"],
                            SL_mm=((Swimming.df[,"SL_mm"]*1-mean(na.omit(Swimming.df[,"SL_mm"])))/sd(na.omit(Swimming.df[,"SL_mm"]))),
                            adj.CP=((Swimming.df[,"adj.CP"]*1-mean(na.omit(Swimming.df[,"adj.CP"])))/sd(na.omit(Swimming.df[,"adj.CP"]))),
                            adj.BD=((Swimming.df[,"adj.BD"]*1-mean(na.omit(Swimming.df[,"adj.BD"])))/sd(na.omit(Swimming.df[,"adj.BD"]))))
# makes numeric matrix for dimensionality based on eigendecomp. (by individuals, then means)
Swimming.df.euc.mat<-Swimming.df.euc[,c("SL_mm","adj.CP","adj.BD")]
swi.means<-aggregate(data=na.omit(Swimming.df.euc[,c("SL_mm","adj.CP","adj.BD")]), .~Lake, FUN=mean)

########## LDA Swimming
Swimming.lda <- lda(Lake ~ SL_mm + adj.BD + adj.CP, data = Swimming.df.euc)
predictions.Swimming <- predict(Swimming.lda)
LD.swim.df<-data.frame(FishID=na.omit(Swimming.df.euc)$ID,Lake=na.omit(Swimming.df.euc)$Lake, swimming.LD1=predictions.Swimming$x[,1],
                       swimming.LD2=predictions.Swimming$x[,2], swimming.LD3=predictions.Swimming$x[,3]) 
# LDA dimension summary 
## ROWS 3-5 IN TABLE S6
(lda.swi.df<-round(lda.vector.df(Swimming.lda),3))%>%write.csv("CORRECTED-lda-swi_df.csv")
# Effective dimensions (Swimming)
## PART OF TABLE 2
swimming.dim.matrix<-round(dim.stats(Swimming.lda,Swimming.df.euc.mat,swi.means),2)


### Swimming LDA dim. Plots
# Swimming LD1 frequency plot
## FIGURE 5B
freq.swi<-LD.freq(LD.swim.df,LD.swim.df$swimming.LD1,LD.swim.df$Lake,
                  "Lake","LD1 of Swimming Traits by Lake","67.1%")
# 3D scatter (EXPLORATORY PLOT, not in paper)
threeD.scatter(LD.swim.df,"swimming.LD3",-4,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "sienna4","grey80","grey50"),c(3,4,5),c('LD1-67.1%', 'LD2-26.4%','LD3-6.5%'))
# 2D scatter 
## FIGURE S1, B1-B2
swi.plot.12<-twoD.scatter(LD.swim.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "sienna4","grey80","grey50"),lda.swi.df,c("LD1 - ","LD2 - "),c(1,2))
swi.plot.13<-twoD.scatter(LD.swim.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                              "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                              "sienna4","grey80","grey50"),lda.swi.df,c("LD1 - ","LD3 - "),c(1,3))
swi.plot.12/swi.plot.13+ plot_layout(guides = 'collect')



##################
#### Trophic ####

####LDA Trophic Traits####
# Trophic data frame#
Trophic.df <- cbind.data.frame(ID=data$FishID,Lake= data$Lake, 
                               adj.BCL= data$adj.BCL, 
                               adj.GW= data$adj.GW, adj.PW= data$adj.PW, 
                               RakerNum= data$RakerNum, adj.JL= data$adj.JL,
                               adj.SnL= data$adj.SnL, adj.ED= data$adj.ED, 
                               adj.HL= data$adj.HL)
rownames(Trophic.df)<-Trophic.df$ID

# LDAprep.A)  RakerNum not interpolated because meristic (not continuous)
Trophic.df.split<-split(Trophic.df,Trophic.df$Lake, drop = F)
NEW.adj.GW<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.GW"],a=x[,"adj.BCL"],b=x[,"adj.JL"], c=x[,"adj.PW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.PW<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.PW"],a=x[,"adj.BCL"],b=x[,"adj.JL"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.JL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.JL"],a=x[,"adj.BCL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.BCL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.BCL"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.SnL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.SnL"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.BCL"],e=x[,"adj.ED"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.ED<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.ED"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.BCL"], f=x[,"adj.HL"])$SecondaryNew.Vec}))
NEW.adj.HL<-unlist(lapply(Trophic.df.split,function(x){
  interp(data=x,y=x[,"adj.HL"],a=x[,"adj.JL"],b=x[,"adj.PW"],c=x[,"adj.GW"],d=x[,"adj.SnL"],e=x[,"adj.ED"], f=x[,"adj.BCL"])$SecondaryNew.Vec}))
Trophic.df<-bind_rows(Trophic.df.split)
Trophic.df[,"adj.GW"]<-NEW.adj.GW
Trophic.df[,"adj.PW"]<-NEW.adj.PW
Trophic.df[,"adj.JL"]<-NEW.adj.JL
Trophic.df[,"adj.BCL"]<-NEW.adj.BCL
Trophic.df[,"adj.SnL"]<-NEW.adj.SnL
Trophic.df[,"adj.ED"]<-NEW.adj.ED
Trophic.df[,"adj.HL"]<-NEW.adj.HL

# LDAprep.B)
Trophic.df.euc<-data.frame(ID=Trophic.df["ID"],Lake=Trophic.df["Lake"],
                           adj.GW=((Trophic.df[,"adj.GW"]*1-mean(na.omit(Trophic.df[,"adj.GW"])))/sd(na.omit(Trophic.df[,"adj.GW"]))),
                           adj.PW=((Trophic.df[,"adj.PW"]*1-mean(na.omit(Trophic.df[,"adj.PW"])))/sd(na.omit(Trophic.df[,"adj.PW"]))),
                           RakerNum=((Trophic.df[,"RakerNum"]*1-mean(na.omit(Trophic.df[,"RakerNum"])))/sd(na.omit(Trophic.df[,"RakerNum"]))),
                           adj.BCL=((Trophic.df[,"adj.BCL"]*1-mean(na.omit(Trophic.df[,"adj.BCL"])))/sd(na.omit(Trophic.df[,"adj.BCL"]))),
                           adj.JL=((Trophic.df[,"adj.JL"]*1-mean(na.omit(Trophic.df[,"adj.JL"])))/sd(na.omit(Trophic.df[,"adj.JL"]))),
                           adj.SnL=((Trophic.df[,"adj.SnL"]*1-mean(na.omit(Trophic.df[,"adj.SnL"])))/sd(na.omit(Trophic.df[,"adj.SnL"]))),
                           adj.ED=((Trophic.df[,"adj.ED"]*1-mean(na.omit(Trophic.df[,"adj.ED"])))/sd(na.omit(Trophic.df[,"adj.ED"]))),
                           adj.HL=((Trophic.df[,"adj.HL"]*1-mean(na.omit(Trophic.df[,"adj.HL"])))/sd(na.omit(Trophic.df[,"adj.HL"]))))
# makes numeric matrix for dimensionality based on eigendecomp. (by individuals, then means)
#  "noGR" dataframes exclude gill rakers because inclusion of gill rakers 
#    substantially reduced sample size.
Trophic.df.euc.mat<-Trophic.df.euc[,c("adj.GW","adj.PW","RakerNum","adj.BCL",
                                      "adj.JL","adj.SnL","adj.ED","adj.HL")]
Trophic.df.euc.mat.noGR<-Trophic.df.euc[,c("adj.GW","adj.PW","adj.BCL",
                                           "adj.JL","adj.SnL","adj.ED","adj.HL")]
tro.means<-aggregate(data=na.omit(Trophic.df.euc[,c("Lake","adj.GW","adj.PW","RakerNum","adj.BCL",
                                                    "adj.JL","adj.SnL","adj.ED","adj.HL")]), .~Lake, FUN=mean)
tro.means.noGR<-aggregate(data=na.omit(Trophic.df.euc[,c("Lake","adj.GW","adj.PW","adj.BCL",
                                                         "adj.JL","adj.SnL","adj.ED","adj.HL")]), .~Lake, FUN=mean)

#LDA
Trophic.lda <- lda(Lake ~ adj.BCL + adj.GW + adj.PW + RakerNum + adj.JL +
                     adj.SnL + adj.ED + adj.HL, data = Trophic.df.euc)
Trophic.lda.noGR <- lda(Lake ~ adj.BCL + adj.GW + adj.PW  + adj.JL +
                          adj.SnL + adj.ED + adj.HL, data = Trophic.df.euc)
#         # note the following two lines are not replicated 
#           for noGR because these are not plotted 
predictions.Trophic <- predict(Trophic.lda)
LD.tro.df<-data.frame(FishID=na.omit(Trophic.df.euc)$ID,Lake=na.omit(Trophic.df.euc)$Lake, trophic.LD1=predictions.Trophic$x[,1],
                      trophic.LD2=predictions.Trophic$x[,2], trophic.LD3=predictions.Trophic$x[,3]) 

# LDA dim summary (Trophic)
## ROWS 3-10 IN TABLE S7
(lda.tro.df<-round(lda.vector.df(Trophic.lda),3))%>%write.csv("CORRECTED-lda-tro_df.csv")
(lda.tro.df.noGR<-lda.vector.df(Trophic.lda.noGR))%>%write.csv("CORRECTED-lda-tro_dfnoGR.csv")
# Effective Dimensions (Trophic)
## PART OF TABLE 2 
trophic.dim.matrix<-round(dim.stats(Trophic.lda,Trophic.df.euc.mat,tro.means),2)
trophic.dim.matrix.noGR<-round(dim.stats(Trophic.lda.noGR,Trophic.df.euc.mat.noGR,tro.means.noGR),2)


# Trophic LD1 frequency plot
## FIGURE 5C
freq.tro<-LD.freq(LD.tro.df,LD.tro.df$trophic.LD1,LD.tro.df$Lake,
                  "Lake","LD1 of Trophic Traits by Lake","44.9%")
# 3D scatter (EXPLORATORY PLOT, not in paper)
threeD.scatter(LD.tro.df,"trophic.LD3",-4,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                            "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                            "sienna4","grey80","grey50"),c(3,4,5),c('LD1-44.9%', 'LD2-22.9%','LD3-11.4%'))
# 2D scatter 
## FIGURE S1, C1-C2
tro.plot.12<-twoD.scatter(LD.tro.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "sienna4","grey80","grey50"),lda.tro.df,c("LD1 - ","LD2 - "),c(1,2))
tro.plot.13<-twoD.scatter(LD.tro.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "sienna4","grey80","grey50"),lda.tro.df,c("LD1 - ","LD3 - "),c(1,3))
tro.plot.12/tro.plot.13+ plot_layout(guides = 'collect')



#######  GEOMORPH  ############
######
#raw landmarks coordinates
AKmorph 
ID <- as.factor((unlist(strsplit(dimnames(AKmorph)[[3]], "_"))%>%matrix(ncol=3,byrow=T))[,3])
dimnames(AKmorph)<-list(NULL,NULL,as.vector(ID))
#slider file
slider 

#Using bending energy for semi-landmark superimposition and General Procrustes Analysis (GPA)
AKmorph.sup <- gpagen(AKmorph, curves= slider, ProcD=FALSE)
#Creating geomorph dataframe#
AKmorph.df <- geomorph.data.frame(coords = AKmorph.sup$coords, 
                                  Lake= data[match(ID,data$FishID),"Lake"], ID=ID)

# LDA & CVA (The same thing, but CVA results make visualizing shape change easier)
AKmorph.cva<-CVA(AKmorph.df$coords,AKmorph.df$Lake, weighting = T, plot = T, cv=T)
Geomorph.df<-data.frame(FishID=AKmorph.df$ID,Lake=AKmorph.df$Lake,AKmorph.cva$CVscores)
AKmorph.lda<-lda(AKmorph.df$Lake~two.d.array(AKmorph.df$coords))

LD.geo.df<-data.frame(FishID=Geomorph.df$FishID,Lake=Geomorph.df$Lake,CV1=Geomorph.df$CV.1,
                      CV2=Geomorph.df$CV.2, CV3=Geomorph.df$CV.3) 

# LDA dimension summary
(lda.geo.df<-lda.vector.df(AKmorph.lda))%>%write.csv("CORRECTED-lda-geo_df.csv")
geo.means<-aggregate(data=na.omit(data.frame(two.d.array(AKmorph.df$coords),Lake=AKmorph.df$Lake)), .~Lake, FUN=mean)

# Effective Dimensions (Shape)
## PART OF TABLE 2 
geo.dim.matrix<-round(dim.stats(
  AKmorph.lda,data.frame(two.d.array(AKmorph.df$coords)),geo.means),2)
#produces error that can be ignored because eigendecomp required for dimentionality est.
# yields negative values for the final two eigenvalues. These EVs are inconsequentially small
# (16-17 orders of magnitude smaller than the 1st eigenvalue)

#Visualize shape differences along CV1 and CV2 
## FIGURE S2
plotAllSpecimens(AKmorph.sup$coords)
plotOutliers(AKmorph.sup$coords)
CVposfive<-5*AKmorph.cva$CVvis[,1]+AKmorph.cva$Grandm
CVnegfive<-(-5)*AKmorph.cva$CVvis[,1]+AKmorph.cva$Grandm
deformGrid2d(CVposfive,CVnegfive,ngrid = 10,wireframe = c(1,2,7:21,1),
             cex1=.75,cex2=.75,col1 = "firebrick4",col2 = "steelblue2")
CV2posfive<-5*AKmorph.cva$CVvis[,2]+AKmorph.cva$Grandm
CV2negfive<-(-5)*AKmorph.cva$CVvis[,2]+AKmorph.cva$Grandm
deformGrid2d(CV2posfive,CV2negfive,ngrid = 10,wireframe = c(1,2,7:21,1),
             cex1=.5,cex2=.5,col1 = "firebrick4",col2 = "steelblue2")

# Geomorph CV1 frequency plot (shape excluded from Figure 5)
LD.freq(Geomorph.df,Geomorph.df$CV.1,Geomorph.df$Lake,
        "Lake","Shape CV1 by Lake", "35.2%")
# 3D scatter
threeD.scatter(LD.geo.df,"CV3",-4.5,c("steelblue4","firebrick3","palegreen4","steelblue2",
                                      "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                      "grey80","grey50"),c(3,4,5),c('LD1-35.2%', 'LD2-17.0%','LD3-11.9%'))
# 2D scatter 
## FIGURE S1 PANELS D1-D2
geo.plot.12<-twoD.scatter(LD.geo.df,c(3,4),c("steelblue4","firebrick3","palegreen4","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "grey80","grey50"),lda.geo.df,c("LD1 - ","LD2 - "),c(1,2))
geo.plot.13<-twoD.scatter(LD.geo.df,c(3,5),c("steelblue4","firebrick3","palegreen4","steelblue2",
                                             "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                             "grey80","grey50"),lda.geo.df,c("LD1 - ","LD3 - "),c(1,3))
geo.plot.12/geo.plot.13+ plot_layout(guides = 'collect')



#
#

# adj.data uses adjusted trait scaled to SD, adj.data.unscaled has unscaled adjusted data for later use in Env. analysis

##### CHECK TEXT FOR WHY "Pelvic_Total" excluded!!!!
adj.data<-merge(Defense.df.euc[,c("ID","Lake","LatPlNum_mean","adj.DS1L","adj.DS2L","adj.PS_mean","MidPl_mean")],
                merge(Swimming.df.euc[,c("ID","SL_mm","adj.BD","adj.CP")],
                      Trophic.df.euc[,c("ID","adj.BCL","adj.GW","adj.PW","RakerNum","adj.JL","adj.SnL","adj.ED","adj.HL")],
                      by= "ID",all= T), by = "ID", all = T)
rownames(adj.data)<-adj.data$ID

adj.data.mat<-subset(adj.data,select = -c(ID,Lake))
#same as line above but excludes GR
adj.data.mat.noGR<-subset(adj.data,select = -c(ID,Lake,RakerNum))

full.means<-aggregate(data=na.omit(subset(adj.data,select = -ID)), .~Lake, FUN=mean)
full.means.noGR<-aggregate(data=na.omit(subset(adj.data,select = -c(ID,RakerNum))), .~Lake, FUN=mean)


#LDA
total.lda <- lda(Lake ~ LatPlNum_mean+adj.DS1L+adj.DS2L+adj.PS_mean+MidPl_mean +SL_mm+adj.BD+adj.CP+
                   adj.BCL + adj.GW + adj.PW + RakerNum + adj.JL +adj.SnL + adj.ED + adj.HL, data = adj.data)
total.lda.noGR <- lda(Lake ~ LatPlNum_mean+adj.DS1L+adj.DS2L+adj.PS_mean+MidPl_mean +SL_mm+adj.BD+adj.CP+
                        adj.BCL + adj.GW + adj.PW +adj.JL +adj.SnL + adj.ED + adj.HL, data = adj.data)
predictions.total <- predict(total.lda)
LD.total.df<-data.frame(FishID=na.omit(adj.data)$ID,Lake=na.omit(adj.data)$Lake, total.LD1=predictions.total$x[,1],
                        total.LD2=predictions.total$x[,2], total.LD3=predictions.total$x[,3]) 
predictions.noGR.total <- predict(total.lda.noGR)
LD.total.noGR.df<-data.frame(ID=na.omit(subset(adj.data,select=-RakerNum))$ID,Lake=na.omit(subset(adj.data,select=-RakerNum))$Lake, 
                             total.LD1=predictions.noGR.total$x[,1],total.LD2=predictions.noGR.total$x[,2], 
                             total.LD3=predictions.noGR.total$x[,3])


# LDA dim summary (full/fullnoGR)
## TABLE S4
(lda.full.df<-round(lda.vector.df(total.lda),3))%>%write.csv("CORRECTED-lda-full_df.csv")
(lda.full.df.noGR<-round(lda.vector.df(total.lda.noGR),3))%>%write.csv("CORRECTED-lda-full_dfnoGR.csv")

# Effective Dimensions (full/fullnoGR)
full.dim.matrix<-round(dim.stats(total.lda,adj.data.mat,full.means),2)
full.dim.matrix.noGR<-dim.stats(total.lda.noGR,adj.data.mat.noGR,full.means.noGR)
# Table of effective dimensions for all trait suites
do.call(cbind,list(defense=round(defense.dim.matrix,2),trophic=round(trophic.dim.matrix,2),
                   swimming=round(swimming.dim.matrix,2),Morph=round(geo.dim.matrix,2),
                   total.noGeo=round(full.dim.matrix,2), total.noGeo.noGR=round(full.dim.matrix.noGR,2)))%>%
  write.csv("CORRECTED-dimMatrices.csv")



# Total LD1 frequency plot
## FIGURE 5D
freq.total<-LD.freq(LD.total.df,LD.total.df$total.LD1,LD.total.df$Lake,
                    "Lake","LD1 of D/Sw/T Traits by Lake","67.7%")
# Total LD 3D scatter plot
# 3D scatter
threeD.scatter(LD.total.df,"total.LD3",-5,c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                            "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                            "sienna4","grey80","grey50"),c(3,4,5),c('LD1-67.7%', 'LD2-11.6%','LD3-7.6%'))
# Total LD 2D scatter plot
## FIGURE S1, E1-E2 
total.plot.12<-twoD.scatter(LD.total.df,c(3,4),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                                 "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                                 "sienna4","grey80","grey50"),lda.full.df,c("LD1 - ","LD2 - "),c(1,2))
total.plot.13<-twoD.scatter(LD.total.df,c(3,5),c("steelblue4","firebrick3","palegreen4","palegreen2","steelblue2",
                                                 "firebrick1","turquoise4","turquoise1","tan3","tan1","darkgreen",
                                                 "sienna4","grey80","grey50"),lda.full.df,c("LD1 - ","LD3 - "),c(1,3))
total.plot.12/total.plot.13+ plot_layout(guides = 'collect')



# ALL COMPONENTS OF FIGURE S1
(def.plot.12+swi.plot.12+tro.plot.12+
    def.plot.13+swi.plot.13+tro.plot.13)+ 
  plot_layout(guides = 'collect',nrow=2,widths=c(1,1,1))
ggsave("CORRECTED-LDA_scatters1.jpg", width = 10, height = 6, dpi = 300)
geo.plot.12/geo.plot.13+ plot_layout(nrow=2,guides = 'collect')
ggsave("CORRECTED-LDA_scatters-geo.jpg", width = 4, height = 7, dpi = 300)
total.plot.12/total.plot.13+ plot_layout(guides = 'collect')
ggsave("CORRECTED-LDA_scatters-total.jpg", width = 6, height = 7, dpi = 300)


###### ALL COMPONENTS OF FIGURE 5
(freq.def+freq.swi+freq.tro)/freq.total
ggsave("CORRECTED-freq-compound.jpg", width = 16, height = 10, dpi = 300)

###############
###########
#   Dimensionality Rarefaction curves

# rarefy.dim functions be randomly (and cumulatively) drawing trait columns w/o replacement, 
# and repeatedly running LDAs and calculating effective dimensionality. 
# This is performed for permNum permutations

# dat = data frame of standardized (Z-scores) trait data for a trait suite
# trait.cols = the column indices with numerical trait data (not factors like ID, Lake, etc.)
# group = Lake column
# permNum = number of permutations desired
abind<-abind::abind
rarefy.dim<-function(dat,trait.cols,group,permNum){
  set.seed(1)
  perm.mat <-c()
  for (i in 1:permNum){
    perm <- sample(min(trait.cols):max(trait.cols), replace=F)
    perm.mat <- cbind(perm.mat, perm)
  }
  
  dim.list<-vector(mode="list")
  for (k in 2:length(min(trait.cols):max(trait.cols))) {
    ED<-matrix(nrow = permNum,ncol = 5) #creates matrix to be filled 
    # by ED calcs for each value of k
    for (i in 1:permNum) { #for each permutation for value of k, 
      # this loop creates a matrix of dimensionalities from lda vectors
      na.s=which(is.na(dat[min(trait.cols):max(trait.cols)]), 
                 arr.ind=TRUE)[,"row"]%>%unique()# identifies which lines in 
      # trait matrix have NAs
      dat.trait=dat[-na.s,perm.mat[1:k,i]] #creates vector of trait indices in 
      # trait matrix for trait number k and permutation i
      lda=if (length(na.s)==0) 
        lda(dat[,perm.mat[1:k,i]],group) else 
          lda(dat.trait,group[-na.s])
      
      EDrow<-c(estimate.ED.eig(prop.trace(lda$svd)),perm.num=i)
      ED[i,]<-EDrow
    }
    colnames(ED)<-names(EDrow)
    dim.list[[k]]<-data.frame(ED,dim.intrinsic=rep(k,permNum))
  }
  dim.list[[1]]<-cbind(matrix(rep(1,4*permNum),nrow = permNum,ncol=4),
                       perm.num=1:permNum,dim.intrinsic=rep(1,permNum))
  ED.rare<-abind(dim.list,along=1)
  ED.rare
}


def.dim.rare<-data.frame(rarefy.dim(Defense.df.euc,3:7,Defense.df.euc$Lake,50),TraitSuite="Defense")
swi.dim.rare<-data.frame(rarefy.dim(Swimming.df.euc,3:5,Swimming.df.euc$Lake,50),TraitSuite="Swimming")
tro.dim.rare<-data.frame(rarefy.dim(Trophic.df.euc,3:10,Trophic.df.euc$Lake,50),TraitSuite="Trophic")
# this procedure for shape produces a bunch of warnings, but that is only because there are two lakes
#  without shape data. These warnings can be ignored.
geo.dim.rare<-data.frame(rarefy.dim(as.data.frame(two.d.array(AKmorph.df$coords)),
                                    1:dim(two.d.array(AKmorph.df$coords))[2],AKmorph.df$Lake,50),TraitSuite="Shape")
full.dim.rare<-data.frame(rarefy.dim(adj.data,3:18,adj.data$Lake,50),TraitSuite="Full")
dim.rare.bySuite<-rbind(def.dim.rare,swi.dim.rare,tro.dim.rare,geo.dim.rare,full.dim.rare)


# creates line plot with loess regression line (thick line) for each trait suite and thin 
# lines for each permutation. rare.plot.zoom is same as rare.plot except zoomed in 
# and with loess reg line
# first plot produces warnings because many points are outside of plot dimensions
rare.plot.zoom<-ggplot(data=dim.rare.bySuite, aes(x=dim.intrinsic, y=n1, color=TraitSuite,
                                                  group = interaction(perm.num,TraitSuite)))+
  geom_smooth(aes(group=TraitSuite,fill=TraitSuite),linewidth=.9,method=loess,se=F,fullrange=T)+
  geom_line(alpha=.2,size=.25)+theme_classic()+xlim(1,5)+ylim(1,4.5)+xlab("traits")+ylab(expression("D"[n1]))+
  scale_color_manual(values=c("steelblue4","grey80","tan1","firebrick3","palegreen3"))+
  scale_fill_manual(values=c("steelblue4","grey80","tan1","firebrick3","palegreen3"))
rare.plot<-ggplot(data=dim.rare.bySuite, aes(x=dim.intrinsic, y=n1, color=TraitSuite,
                                             group = interaction(perm.num,TraitSuite)))+
  geom_line(alpha=.2,size=.15)+theme_classic()+xlab("traits")+
  ylab(expression("D"[n1]))+theme(legend.position="none")+
  scale_color_manual(values=c("steelblue4","grey80","tan1","firebrick3","palegreen3"))

## FIGURE 8
rare.plot+rare.plot.zoom  
ggsave("CORRECTED-Dim-rare.jpg", width = 6.5, height = 3.25, dpi = 300)