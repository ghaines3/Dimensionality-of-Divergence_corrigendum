####Library####
library(devtools)
install.packages("ggord")
install.packages("openair")
install.packages("factoextra")
install.packages("ggfortify")
install.packages("viridis")
install.packages("Hmisc")
install.packages("berryFunctions")
install.packages("plot3D")
remotes::install_github('fawda123/ggord')
install.packages("plotly")
install.packages("factoextra")
install.packages("remotes")
remotes::install_github("AkselA/R-ymse")
install.packages("ggalluvial")
install.packages("LOST")
install.packages("RVAideMemoire")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("pcaMethods")
install_github("vqv/ggbiplot")
install_github("gavinsimpson/ggvegan")
install.packages("ggdendro")
library(ggvegan)
library(LOST)
library(ymse)
library(factoextra)
library(ggord)
library(corrplot)
library(geomorph)
library(ggplot2)
library(dplyr)
library(stringr)
library(MASS)
library(shapes)
library(Momocs)
library(svd)
library(rgl)
library(ape)
library(gridExtra)
library(mvnormtest)
library(car)
library(readxl)
library(Morpho)
library(purrr)
library(forcats)
library(openair)
library(factoextra)
library(ggfortify)
library(ggridges)
library(viridis)
library(plyr)
library(Hmisc)
library(GGally)
library(data.table)
library(vegan)
library(berryFunctions)
library(plot3D)
library(plotly)
library(cowplot)
library(ggalluvial)
library(RVAideMemoire)
library(ggpubr)
library(gplots)
library(multigroup)
library(pracma)
library(patchwork)
library(ggrepel)
library(ggdendro)
library(MVN)

####Read data####
stickle_data <- read.csv("AKStickleData_1-1-21AmNat.csv") #non-defensive trait data
def_data <- select(read.csv("StickleDefensive_4-22-20AmNat.csv"),-Lake) #defensive trait data
data<-merge(stickle_data,def_data, by = "FishID")%>%
  mutate(FishID = as.character(FishID))%>%
  mutate(Lake = as.factor(Lake))#merges datasets by ID

#Vector object for Lake 
Lake <- data$Lake
Lake


# function to get coefficient for size adjustment
# trait - name of trait col to get coefficient(eg, data$DS1Length)
size.coef<-function(trait){
  ancova=lm(data=data, log(trait)~log(SL_mm)+Lake)
  coef=ancova$coefficients[2]
  coef
}

# spine length data contains zeroes that are excluded from size correction
DS1.coef<-lm(data=filter(data,DS1Length!=0),log(DS1Length)~log(SL_mm)+Lake)$coefficients[2]#DS1Length (DS1)
DS2.coef<-lm(data=filter(data,DS2Length!=0),log(DS2Length)~log(SL_mm)+Lake)$coefficients[2] #DS2Length (DS2)
PSL.coef<-lm(data=filter(data,PSLLength!=0),log(PSLLength)~log(SL_mm)+Lake)$coefficients[2] #PSLLength (PSL)
PSR.coef<-lm(data=filter(data,PSRLength!=0),log(PSRLength)~log(SL_mm)+Lake)$coefficients[2] #PSRLength (PSR)
BD.coef<-size.coef(data$BodyDepth) #BodyDepth (BD)
CP.coef<-size.coef(data$Caudal.peduncle) #Caudal.peduncle (CP)
BCL.coef<-size.coef(data$BuccalCavityLength) #BuccalCavityLength (BCL)
GW.coef<-size.coef(data$GapeWidth) #GapeWidth (GW)
PW.coef<-size.coef(data$PteroticWidth) #PteroticWidth (PW)
SnL.coef<-size.coef(data$Snout.Length) #Snout.Length (SnL)
ED.coef<-size.coef(data$Eye.Diameter) #Eye.Diameter (ED)
HL.coef<-size.coef(data$Head.Length) #Head.Length (HL)
JL.coef<-size.coef(data$Jaw.Length) #Jaw.Length (JL)


Coef<-matrix(c(DS1.coef, DS2.coef, PSL.coef, PSR.coef,BD.coef, CP.coef, BCL.coef, GW.coef, PW.coef, JL.coef, 
               SnL.coef, ED.coef, HL.coef),dimnames = list(
                 c("DorsalSpine1", "DorsalSpine2", "PelvicSpineLeft","PelvicSpineRight", 
                   "BodyDepth", "CaudalPeduncle", "BuccalCavityLength","GapeWidth", 
                   "PteroticWidth", "JawLength", "SnoutLength", "EyeDiameter", "HeadLength")))


data$adj.DS1L<-data$DS1Length*(mean(data$SL_mm)/data$SL_mm)^Coef["DorsalSpine1",]
data$adj.DS2L<-data$DS2Length*(mean(data$SL_mm)/data$SL_mm)^Coef["DorsalSpine2",]
data$adj.PSL<-data$PSLLength*(mean(data$SL_mm)/data$SL_mm)^Coef["PelvicSpineLeft",]
data$adj.PSR<-data$PSRLength*(mean(data$SL_mm)/data$SL_mm)^Coef["PelvicSpineRight",]
data$adj.PS.mean<-rowMeans(
  subset(data, select = c("adj.PSL", "adj.PSR")), 
  na.rm = TRUE) # This line averages L and R pelvic spine lengths
data$adj.BD<-data$BodyDepth*(mean(data$SL_mm)/data$SL_mm)^Coef["BodyDepth",]
data$adj.CP<-data$Caudal.peduncle*(mean(data$SL_mm)/data$SL_mm)^Coef["CaudalPeduncle",]
data$adj.BCL<-data$BuccalCavityLength*(mean(data$SL_mm)/data$SL_mm)^Coef["BuccalCavityLength",]
data$adj.GW<-data$GapeWidth*(mean(data$SL_mm)/data$SL_mm)^Coef["GapeWidth",]
data$adj.PW<-data$PteroticWidth*(mean(data$SL_mm)/data$SL_mm)^Coef["PteroticWidth",]
data$adj.JL<-data$Jaw.Length*(mean(data$SL_mm)/data$SL_mm)^Coef["JawLength",]
data$adj.SnL<-data$Snout.Length *(mean(data$SL_mm)/data$SL_mm)^Coef["SnoutLength",]
data$adj.ED<-data$Eye.Diameter*(mean(data$SL_mm)/data$SL_mm)^Coef["EyeDiameter",]
data$adj.HL<-data$Head.Length *(mean(data$SL_mm)/data$SL_mm)^Coef["HeadLength",]

data %>%
  group_by(Lake) %>%
  summarize(cor(SL_mm,adj.BD, na.rm=T))

# check that size corrction worked (SL_mm coefficient should be ~0)
lm(data=data, adj.DS1L~Lake+SL_mm)
lm(data=data, adj.DS2L~Lake+SL_mm)
lm(data=data, adj.PSL~Lake+SL_mm)
lm(data=data, adj.PSR~Lake+SL_mm)
lm(data=data, adj.BD~Lake+SL_mm)
lm(data=data, adj.CP~Lake+SL_mm)
lm(data=data, adj.BCL~Lake+SL_mm)
lm(data=data, adj.GW~Lake+SL_mm)
lm(data=data, adj.PW~Lake+SL_mm)
lm(data=data, adj.JL~Lake+SL_mm)
lm(data=data, adj.SnL~Lake+SL_mm)
lm(data=data, adj.ED~Lake+SL_mm)
lm(data=data, adj.HL~Lake+SL_mm)




### Interp() = function to estimate NAs
# this function uses regression against most highly correlated variable (as long as r>0.1) to estimate missing
# values. For cases in which a value for the most highly correlated variable is also missing, function defaults 
# to next most correlated variable
##
# Input: data = a dataframe; y = col to estimate variables for; a, b, and c = cols to test for best correlation
##
# Output: a list containing:
#   New.Vec = vector with missing values interpolated from most correlated variable
#   Any.na = remaining missing values; Arg.cor = Pearsons r for args a, b, & c; 
#   Best_Model = which col used for interpolation and linear regression terms
#   SecondaryNew.Vec = vector with missing values interpolated from most correlated & second most correlated variables
#   SecondarySum = if secondary col was required for interpolation, includes which col used
#   
interp<-function(data,y,a,b,c,d,e,f,g){
  if(var(na.omit(y))==0)stop("NAs not estimated, no variance in y")
  if(missing(b)){b<-a}
  if(missing(c)){c<-a}
  if(missing(d)){d<-a}
  if(missing(e)){e<-a}
  if(missing(f)){f<-a}
  if(missing(g)){g<-a}
  models<-rbind(a=(if(sd(na.omit(a))==0)NA else cor.test(y,a)$estimate),b=(if(sd(na.omit(b))==0)NA else cor.test(y,b)$estimate),
                c=(if(sd(na.omit(c))==0)NA else cor.test(y,c)$estimate),d=(if(sd(na.omit(d))==0)NA else cor.test(y,d)$estimate),
                e=(if(sd(na.omit(e))==0)NA else cor.test(y,e)$estimate),f=(if(sd(na.omit(f))==0)NA else cor.test(y,f)$estimate),
                g=(if(sd(na.omit(g))==0)NA else cor.test(y,g)$estimate))
  models.2<-models
  models.reg<-list(a=lm(y~a)[["coefficients"]][1:2],b=lm(y~b)[["coefficients"]][1:2],c=lm(y~c)[["coefficients"]][1:2],
                   d=lm(y~d)[["coefficients"]][1:2],e=lm(y~e)[["coefficients"]][1:2],f=lm(y~f)[["coefficients"]][1:2],g=lm(y~g)[["coefficients"]][1:2])
  bestfit<-if((max(abs(models),na.rm = T)<.1)==T)NA else models.reg[which.max(abs(models))][[1]]##
  models.2[which.max(abs(models))]<-NA
  models.2[duplicated(models)]<-NA
  nextbest<-if(length(which.max(abs(models.2)))==0)NA else if((max(abs(models.2),na.rm = T)<.1)==T)NA else models.reg[which.max(abs(models.2))][[1]]##
  vec<-if(length(which.max(abs(models)))==0)NA else if(which.max(abs(models))==1)a else if(which.max(abs(models))==2)b else if(which.max(abs(models))==3)c else 
    if(which.max(abs(models))==4)d else if(which.max(abs(models))==5)e else 
      if(which.max(abs(models))==6)f else if(which.max(abs(models))==7)g##
  y[which(is.na(y))]<-bestfit[1]+vec[which(is.na(y))]*bestfit[2]
  New.Vec<-y
  vec.2<-if(length(which.max(abs(models.2)))==0)NA else if(which.max(abs(models.2))==1)a else 
    if(which.max(abs(models.2))==2)b else if(which.max(abs(models.2))==3)c else if(which.max(abs(models.2))==4)d else 
      if(which.max(abs(models.2))==5)e else if(which.max(abs(models.2))==6)f else if(which.max(abs(models.2))==7)g
  y[which(is.na(y))]<-nextbest[1]+vec.2[which(is.na(y))]*nextbest[2]
  New.Vec.2<-y
  list(PrimaryNew.Vec= if((max(abs(models),na.rm = T)<.1)==T)c("no variables with r > 0.1 for interpolation") else New.Vec,
       Any.na=which(is.na(New.Vec)), Arg.cor=models,
       Best_Model=if((max(abs(models),na.rm = T)<.1)==T)c("no values estimated") else##
         if(which.max(abs(models))==1)c("a used for estimate",bestfit) else 
           if(which.max(abs(models))==2)c("b used for estimate",bestfit) else 
             if(which.max(abs(models))==3)c("c used for estimate",bestfit) else
               if(which.max(abs(models))==4)c("d used for estimate",bestfit) else 
                 if(which.max(abs(models))==5)c("e used for estimate",bestfit) else 
                   if(which.max(abs(models))==6)c("f used for estimate",bestfit) else
                     if(which.max(abs(models))==7)c("g used for estimate",bestfit),
       SecondaryNew.Vec=if(is.numeric(vec.2))New.Vec.2 else New.Vec,
       SecondarySum=if(length(which(is.na(New.Vec)))==0)"Only Primary cor. required" else 
         if(sum(na.omit(models.2))==0)c("no variables with r > 0.1 for interpolation") else
           if((max(abs(models.2),na.rm = T)<.1)==T)c("no variables with r > 0.1 for interpolation") else ## 
             list(New.Vec.2,Any.na=which(is.na(New.Vec.2)),
                  NextBest_Model=if(length(which.max(abs(models.2)))!=0)if(which.max(abs(models.2))==1)c("a used for estimate",nextbest) else 
                    if(which.max(abs(models.2))==2)c("b used for estimate",nextbest) else 
                      if(which.max(abs(models.2))==3)c("c used for estimate",nextbest) else
                        if(which.max(abs(models.2))==4)c("d used for estimate",nextbest) else 
                          if(which.max(abs(models.2))==5)c("e used for estimate",nextbest) else 
                            if(which.max(abs(models.2))==6)c("f used for estimate",nextbest) else
                              if(which.max(abs(models.2))==7)c("g used for estimate",nextbest)))
}

# calculates proportion of trace from lda
prop.trace<-function(svd){(svd^2/sum(svd^2))}

# makes dataframe of vector info from lda
lda.vector.df<-function(lda){
  round(data.frame(rbind("singular values"=lda$svd,"eigenvalues"=(lda$svd)^2,
                         lda$scaling)),4)
}




#### The following are to make sure that there is no allometric effect on traits in two ways:
#  first a trait~SL*Lake linear model, then spearmans's correlation to confirm results because traits are not continuous 
#linear regression
plate.allom.lm<-summary(lm(rowMeans(subset(data, select = c("LatPlNum_R", "LatPlNum_L")), na.rm = TRUE)~data$SL_mm*data$Lake))
p.adjust(plate.allom.lm$coefficients[-1,4], method = "fdr")# no significant effect of length
#spearman's rho
plate.spear<-by(data, data$Lake, FUN = function(x) cor.test(x$SL, rowMeans(x[,c("LatPlNum_R", "LatPlNum_L")]),method = "spearman",exact=F, data=x))
unlist(lapply(plate.spear,FUN= function(x) x$p.value))%>%p.adjust(method="fdr")# no significant effect of length

raker.allom.lm<-summary(lm(data$RakerNum~data$SL_mm*data$Lake))
p.adjust(raker.allom.lm$coefficients[-1,4], method = "fdr") # no significant effect of length
raker.spear<-by(data, data$Lake, FUN = function(x) cor.test(x$SL, x$RakerNum,method = "spearman",exact=F, data=x))
unlist(lapply(raker.spear,FUN= function(x) x$p.value))%>%p.adjust(method="fdr") # no significant effect of length
