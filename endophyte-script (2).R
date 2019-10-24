



# factors affecting endophytic fungal diversity and communities
###################################
###################################

## Pakages we need for the analysis

library(mvabund)
library(vegan)
library(boral)
library(MASS)
library(ggplot2)
library(gridExtra)
library(mgcv)
##################################
##################################

## Data input

EndoAbun= read.csv (file="matrix_endo.csv", 
                   header = T, row.names = 1)

EndoMetaData= read.csv (file="metadata_end.csv", header = T, row.names = 1)

# are the rownames matching?
rownames(EndoAbun) == rownames(EndoMetaData)

# setting the explaining predictors as factors

Tissue <- factor(EndoMetaData$SUBSTRATE, levels = c("Leaf", "Branch"))
Time <- factor(EndoMetaData$TIME, levels = c("1","2","3"))
Locality<- factor(EndoMetaData$LOCALITY, levels = c("Bisoton", "Hasan abad", "Khosro abad",
                                                   "Sarmil","Sorkhe dize"))
Temperature<-factor(EndoMetaData$TEMPRATURE, levels = c("25 Degree","4 Degree"))
IR= EndoMetaData$IR ### isolation rate
## Data distribution
hist(EndoMetaData$IR)
hist(log(EndoMetaData$IR))
boxplot(EndoMetaData$IR~EndoMetaData$SUBSTRATE)
boxplot(EndoMetaData$IR~EndoMetaData$TIME)
boxplot(EndoMetaData$IR~EndoMetaData$LOCALITY)
boxplot(EndoMetaData$IR~EndoMetaData$TEMPRATURE)

#########################
###### 1- Isolation rate
#########################

### GLM model for IR
### new model
Ir.m= glm(IR ~Temperature * Locality + Temperature *Time + Temperature * Tissue ,
           data = EndoMetaData)
par(mfrow= c(2,2))
plot(Ir.m)
AIC(Ir.m)
IR.m.anova = anova(Ir.m, test = "F")
IR.m.sammary = summary(Ir.m)

## don't know what kind of model or family is good for this kind of data non-integer 
#with lots of ZEROs
# plot the significant interactions
dev.off()

# ## Tissue * Time interaction
# timeplot= ggplot(data = EndoMetaData,aes(x=EndoMetaData$TIME,y=IR,
#                                          fill= factor(EndoMetaData$SUBSTRATE)))+
#   geom_bar(stat="identity",position="dodge",width=0.3)+
#   scale_fill_manual(values = c ("darkgray", "lightgray"),name="Tissue type",
#                     breaks=c("Leaf", "Branch"),
#                     labels=c("Leaf", "Branch"))+
#   xlab("Sampling time")+ylab("Isolation rate")+
#   scale_x_continuous(breaks= c(1,2,3),labels = c("May", "Jun","July"))+
#   theme_bw()+
#   theme(legend.position = "top", axis.text= element_text(size = rel(1.3)),
#         axis.title.y = element_text(size = rel(1.5)),
#         axis.title.x = element_text(size = rel(1.5)),
#         legend.text = element_text(size = rel(1.5)),
#         legend.title = element_text(size = rel(1.5)),
#         panel.grid.major = element_blank(),
#         axis.line = element_line(colour = "black", color = "black"))
# ## Tissue * Locality interaction
# 
ggplot(data = EndoMetaData,aes(x=Locality,y=IR, fill= factor(EndoMetaData$SUBSTRATE)))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  scale_fill_manual(values = c ("darkgray", "lightgray"),name="Tissue type",
                    breaks=c("Leaf", "Branch"),
                    labels=c("Leaf", "Branch"), guide = FALSE)+
  xlab("Sampling site")+ylab("Isolation rate")+
  theme_bw()+
  theme(legend.position = "top", axis.text= element_text(size = rel(1.3)),
        axis.title.y = element_text(size = rel(1.5)),
        axis.title.x = element_text(size = rel(1.5)),
        legend.text = element_text(size = rel(1.5)),
        legend.title = element_text(size = rel(1.5)),
        panel.grid.major = element_blank(),
              axis.line = element_line(colour = "black", color = "black"))

################################
### FIG 2
################################
# dev.off()
# png(file = "fig 2.png",  width = 2100, height = 2100, res= 300)
# grid.arrange(timeplot,locaplot,nrow=2, ncol=1 )
# dev.off()
########################################
# ### 2- Hill Diversities
########################################
# Remove zero observations for diversity and richness calculation
NotZero = EndoMetaData$SUCCESS > 0 #filter for zero-observation samples
EndoAbunZero = EndoAbun[NotZero,]
EndoMetaZero = EndoMetaData[NotZero,]
##aranging the factors with the new datasets
LocalityZ<- factor(EndoMetaZero$LOCALITY,levels = c("Bisoton", "Hasan abad", "Khosro abad",
                                                    "Sarmil","Sorkhe dize"))
TissueZ <- factor(EndoMetaZero$SUBSTRATE, levels = c("Leaf", "Branch"))
TimeZ <- factor(EndoMetaZero$TIME, levels = c("1","2","3"))
TemperatureZ<-factor(EndoMetaZero$TEMPRATURE, levels = c("25 Degree","4 Degree"))
#EndoRichness = specnumber(EndoAbunZero)
#hist(EndoRichness)

EndoHill = renyi(EndoAbunZero, scale=c(0,1,2), hill=T)
Endohill.1 = EndoHill$"0"#this is richness
Endohill.2 = EndoHill$"1"#antilogarithm of the Shannon representing the abundant species
Endohill.3 = EndoHill$"2"#inverse Simpson representing the very abundant species
hist(Endohill.1)
hist(Endohill.2)
hist(Endohill.3)

### MODELS
### I want to see if the time,loclaity and temperature have differensial effects on
# different tissues diversity

############### First hill
# EHill1.m= glm(Endohill.1~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
#                , data = EndoMetaZero, family = poisson(link = "log"))
# AIC(EHill1.m)
# par(mfrow= c(2,2))
# plot(EHill1.m)
# dev.off()
# 
# EHill1.m.anova= anova(EHill1.m, test = "Chisq")
# EHill1.m.summary= summary(EHill1.m)# why this doesn't show the significance?
########new model
EHill1.m.2= glm(Endohill.1~ TemperatureZ *LocalityZ+TemperatureZ *TimeZ+TemperatureZ * TissueZ
              , data = EndoMetaZero, family = poisson(link = "log"))
AIC(EHill1.m.2)
par(mfrow= c(2,2))
plot(EHill1.m.2)
dev.off()

EHill1.m.2anova= anova(EHill1.m.2, test = "Chisq")
EHill1.m.2summary= summary(EHill1.m.2)

#### discard the locality factor
EHill1.m.3= glm(Endohill.1~ TemperatureZ *TimeZ+TemperatureZ * TissueZ
                , data = EndoMetaZero, family = poisson(link = "log"))
AIC(EHill1.m.3)
par(mfrow= c(2,2))
plot(EHill1.m.3)
dev.off()

EHill1.m.3anova= anova(EHill1.m.3, test = "Chisq")
EHill1.m.3summary= summary(EHill1.m.3)

# See if glm.nb is a better fit
# EHill1.m.nb= glm.nb(Endohill.1~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
#               , data = EndoMetaZero, link = "log")
# AIC(EHill1.m.nb)
# par(mfrow= c(2,2))
# plot(EHill1.m.nb)
# dev.off()
## These models are not that different but due to the warnings that I get from glm.nb
## I choose the poisson glm model

############### second Hill
## use the same model orders here

# EHill2.m= glm(Endohill.2~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
#               ,data = EndoMetaZero, family =Gamma(link = "log"))
#AIC(EHill2.m)
# par(mfrow= c(2,2))
# plot(EHill2.m)
# dev.off()
# 
# EHill2.m.anova= anova(EHill2.m,test = "F")
# EHill2.m.summary= summary(EHill2.m)
##############new model
EHill2.m2= glm(Endohill.2~ TemperatureZ *LocalityZ+TemperatureZ *TimeZ+TemperatureZ * TissueZ 
              ,data = EndoMetaZero, family =Gamma(link = "log"))
AIC(EHill2.m2)
par(mfrow= c(2,2))
plot(EHill2.m2)
dev.off()

EHill2.m2.anova= anova(EHill2.m2,test = "F")
EHill2.m2.summary= summary(EHill2.m2)

#### discard the locality factor
EHill2.m3= glm(Endohill.2~ TemperatureZ *TimeZ+TemperatureZ * TissueZ 
               ,data = EndoMetaZero, family =Gamma(link = "log"))
AIC(EHill2.m3)
par(mfrow= c(2,2))
plot(EHill2.m3)
dev.off()

EHill2.m3.anova= anova(EHill2.m3,test = "F")
EHill2.m3.summary= summary(EHill2.m3)

## try the glm.nb
# Ehill2.m.nb= glm.nb(Endohill.2~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
#                     ,data = EndoMetaZero)
# warnings()

############# Third Hill
# EHill3.m= glm(Endohill.3~ TissueZ*LocalityZ+TissueZ*TimeZ+TissueZ*TemperatureZ 
#               ,data = EndoMetaZero, family =Gamma(link = "log"))
# AIC(EHill3.m)
# par(mfrow= c(2,2))
# plot(EHill3.m)
# dev.off()
# 
# EHill3.m.anova= anova(EHill3.m,test = "F")
# EHill3.m.summary= summary(EHill3.m)

########### new model
EHill3.m2= glm(Endohill.3~ TemperatureZ *TimeZ+TemperatureZ * TissueZ 
              ,data = EndoMetaZero, family =Gamma(link = "log"))
AIC(EHill3.m2)
par(mfrow= c(2,2))
plot(EHill3.m2)
dev.off()

EHill3.m2.anova= anova(EHill3.m2,test = "F")
EHill3.m2.summary= summary(EHill3.m2)


##############################
#### 3- Community composition
##############################

### Define CORE OTUs

## Summarize OTU observation
TotalCount = apply(EndoAbun,2,sum)

## The average observation of OTUs
MeanCount=apply(EndoAbun,2,function(vec) mean(vec[vec>0]))

## In how many samples is an OTU present?
TotalPresent = apply(EndoAbun,2,function(vec) sum(vec>0))

## The highest number of an OTU in a sample
MaximumCount=apply(EndoAbun,2,max)

## Plotting observation against abundance
plot(TotalPresent, MaximumCount, xlab="OTU Observation",
     ylab="OTU Maximum Abundance", pch=20)

plot(TotalPresent, log(MaximumCount), xlab="OTU Observation",
     ylab="log(OTU Maximum Abundance)", pch=20)

## Create a smoothed trendline
gam1 = gam(log(MaximumCount)~s(TotalPresent))

plot(gam1, residuals=T, shade=T, rug=F, cex=2.6,
     xlab="OTU Observation", ylab="log Mean Abundance") # , xaxp=c(0,150,15)

## keep core OTUs
OTUobserv = TotalPresent > 7
EndoCorAbun = EndoAbun[,OTUobserv]

### name of the Core OTUs
COREOTUS=colnames(EndoCorAbun)

#### Remove the Zero samples from the Core OTU abbundnace object and metadata
IsolSucc = apply(EndoCorAbun,1, sum)
NotZero2= IsolSucc>0
ECorAbunZero = EndoCorAbun[NotZero2,]
ECorMetaZero = EndoMetaData[NotZero2,]
row.names(ECorAbunZero)==row.names(ECorMetaZero)
corOTUs<-colnames(ECorAbunZero)
TissueC<-factor(ECorMetaZero$SUBSTRATE, levels = c("Leaf", "Branch"))
LocalityC<- factor(ECorMetaZero$LOCALITY,levels = c("Bisoton", "Hasan abad", "Khosro abad",
                                                    "Sarmil","Sorkhe dize"))
TimeC<- factor(ECorMetaZero$TIME, levels = c("1", "2","3"))
TemperatureC<- factor(ECorMetaZero$TEMPRATURE, levels = c("25 Degree","4 Degree"))

### Multispecies Model for Core OTUs
ECrMvabund= mvabund(ECorAbunZero)
plot(ECrMvabund)

EndoMV.m= manyglm (ECrMvabund~TissueC*LocalityC+TissueC*TimeC+TissueC*TemperatureC,
                   data = ECorMetaZero, family = "negative.binomial", show.residuals=T)
class(ECrMvabund)
plot.manyglm(EndoMV.m)
EndoMV.m.anova= anova.manyglm(EndoMV.m,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo")
EndoMV.m.sum= summary.manyglm(EndoMV.m, nBoot=100, test="LR",p.uni="adjusted", 
                              resamp="montecarlo")
## Which OTUs are significaantly affected
EnAnova <- as.data.frame(EndoMV.m.anova$uni.p)
otuTissueEf<-colnames(EnAnova)[EnAnova["TissueC",]<= 0.05]
otuLocEf<- colnames(EnAnova)[EnAnova["LocalityC",]<= 0.05]
otuTimEf<- colnames(EnAnova)[EnAnova["TimeC",]<= 0.05]
otuTempEf<- colnames(EnAnova)[EnAnova["TemperatureC",]<= 0.05]

### try to visualize these effects
## Tissue effects
OTUtissu<- c("Byssochlamys.spectabilis.","Gnomoniaceae.sp..66","Microsphaeriopsis.olivacea",
              "Penicillium.sp..A21","Preussia.sp..A31")
TissuABUN<- ECorAbunZero[OTUtissu]##Keeping only tissue affected OTUs
# get the mean valuse for each OTU in each tisse
Tissuemean <- aggregate(. ~ ECorMetaZero$SUBSTRATE, TissuABUN, mean)
#CReat a data frame of the mean valuse
TissuMeanfram<- as.data.frame(Tissuemean,optional=TRUE)
attr(TissuMeanfram, "row.names")<- c("Branch", "Leaf")

### creat a matrix of mean observation of OTUs affected by tissue for ploting
Tissudata<- data.matrix (TissuMeanfram[2:6],rownames.force = NA )

pdf(file = "Effect of Tissue on OTU observation.pdf", paper = "a4", width = 7, height = 4)
barplot(Tissudata, legend.text =TRUE, beside = TRUE,ylab= "mean observation per sample",
        names.arg= c("B. spectabilis", "Gnomoniaceae sp.", "M. olivacea","Penicillium sp.",
                     "Preussia sp."), axes= TRUE,ylim= c(0,1),  cex.names = 0.8,
              args.legend = list(x = "topright",bty= "n"), border = "Black" )
dev.off()

### Temprature effects                                                                                       
OTUtemp<- c ("Aspergillus.sp..A20","Aureobasidium.sp..A17","Byssochlamys.spectabilis."
             ,"Microsphaeriopsis.olivacea", "Preussia.sp..A31")
TempABUN<- ECorAbunZero[OTUtemp]##Keeping only Temp affected OTUs
# get the mean valuse for each OTU in each temp
Tempmean <- aggregate(. ~ ECorMetaZero$TEMPRATURE, TempABUN, mean)
#CReat a data frame of the mean valuse
TempMeanfram<- as.data.frame(Tempmean,optional=TRUE)
attr(TempMeanfram, "row.names")<- c("25 Degree", "4 Degree")

### creat a matrix of mean observation of OTUs affected by temp for ploting
Tempdata<- data.matrix (TempMeanfram[2:6],rownames.force = NA )
pdf(file = "Effect of Temprature on OTU observation.pdf", paper = "a4", width = 7, height = 4)
barplot(Tempdata,legend.text =TRUE, beside = TRUE,ylab= "mean observation per sample" ,
        names.arg= c ("A20","A17","B.spec","M.oliv",
                      "A31"), axes= TRUE,ylim= c(0,1.6), cex.names = 0.8,
        args.legend = list(x = "topleft",bty= "n"), border = "Black", 
        width = 0.5)
                      
dev.off()

## Locality effects
plot(ECorAbunZero$Byssochlamys.spectabilis.~ ECorMetaZero$LOCALITY)

# time effects 
plot(ECorAbunZero$Alternaria.sp..A25~ ECorMetaZero$TIME)
##########################################################
##### NEW GLMS MODEL FOR CORE OTUS WITH DIFFERENT INTERACTIONS
##########################################################
EndoMV.m2= manyglm (ECrMvabund~TemperatureC*LocalityC+TemperatureC*TimeC+TemperatureC*TissueC,
                   data = ECorMetaZero, family = "negative.binomial", show.residuals=T)

plot.manyglm(EndoMV.m2)
EndoMV.m2.anova= anova.manyglm(EndoMV.m2,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo")
EndoMV.m2.sum= summary.manyglm(EndoMV.m2, nBoot=100, test="LR",p.uni="adjusted", 
                              resamp="montecarlo")
## Which OTUs are significaantly affected
EnAnova2 <- as.data.frame(EndoMV.m2.anova$uni.p)
otuTissueEf2<-colnames(EnAnova2)[EnAnova2["TissueC",]<= 0.05]
otuLocEf2<- colnames(EnAnova2)[EnAnova2["LocalityC",]<= 0.05]
otuTimEf2<- colnames(EnAnova2)[EnAnova2["TimeC",]<= 0.05]
otuTempEf2<- colnames(EnAnova2)[EnAnova2["TemperatureC",]<= 0.05]
otutemplocEf<-colnames(EnAnova2)[EnAnova2["TemperatureC:LocalityC",]<= 0.05]
otutemtimEf<-colnames(EnAnova2)[EnAnova2["TemperatureC:TimeC",]<= 0.05]
otutemTISSEf<-colnames(EnAnova2)[EnAnova2["TemperatureC:TissueC",]<= 0.05]

####interaction plots
?boxplot()
boxplot(ECorAbunZero$Microsphaeriopsis.olivacea~TemperatureC)
boxplot(ECorAbunZero$Microsphaeriopsis.olivacea~TimeC)
# temp * Time interaction

# library(lattice)
bwplot(ECorAbunZero$Microsphaeriopsis.olivacea~LocalityC|TemperatureC,
       ylab = "M. olivacea", xlab = "Sampling site",
       scales=list(x=list(labels=c("Biseton","Mahi Dasht","Khosro abad","Sarmil","Sorkhe dize"))),
       par.settings = list(box.umbrella=list(col= "black"),
                            box.dot=list(col= FALSE),
                            box.rectangle = list(col= "black")
                           , plot.symbol = list(col= "black")))
# 
# barchart(ECorAbunZero$Cytospora.sp..AC35~TimeC|TemperatureC,
#        ylab = "Cytospora sp.", xlab = "Time of sampling",
#        scales=list(x=list(labels=c("May","Jun","July"))),
#        col= "gray", borders=FALSE)
# 
# bwplot(ECorAbunZero$Cytospora.sp..AC35~TimeC|TemperatureC,
#        ylab = "Cytospora sp.", xlab = "Time of sampling",
#        scales=list(x=list(labels=c("May","Jun","July"))),
#        par.settings = list(box.umbrella=list(col= "black"), 
#                            box.dot=list(col= FALSE), 
#                            box.rectangle = list(col= "black")
#                            , plot.symbol = list(col= "black")))
dev.off()
library(FSA)
sd(ECorAbunZero$Cytospora.sp..AC35)
### bar plot
Cytosporaplot= ggplot(data = ECorMetaZero,aes(x=ECorMetaZero$TIME,y=ECorAbunZero$Cytospora.sp..AC35,
                                         fill= factor(ECorMetaZero$TEMPRATURE)))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  scale_x_discrete(breaks=c("25 degree","4 degree"),labels=c("25 degree","4 degree"))+
  scale_fill_manual(values=c("gray27", "gray60"), name = "Incubation temperature", 
                    labels=c(expression(25*degree~C),expression(4*degree~C))) +
  xlab("Sampling time")+ylab("Cytospora sp.")+
  scale_x_continuous(breaks= c(1,2,3),labels = c("May", "June","July"))+
  theme_bw()+
  theme(legend.position = "top",axis.text= element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.title.x = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", color = "black"))

?geom_errorbar()
##new plot revision 2
ac35<-aggregate(ECorAbunZero$Cytospora.sp..AC35~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,mean)

ac35.se<-aggregate(ECorAbunZero$Cytospora.sp..AC35~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,se)
ac35.sd<-aggregate(ECorAbunZero$Cytospora.sp..AC35~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,sd)

ac35$time<-ac35$`ECorMetaZero$TIME`
ac35$`ECorMetaZero$TIME`<-NULL
ac35$temp<-ac35$`ECorMetaZero$TEMPRATURE`
ac35$`ECorMetaZero$TEMPRATURE`<-NULL
ac35$ac35mean<-ac35$`ECorAbunZero$Cytospora.sp..AC35`

ac35$`ECorAbunZero$Cytospora.sp..AC35`<-NULL
ac35$time<-as.factor(ac35$time)
ac35$se<-ac35.se$`ECorAbunZero$Cytospora.sp..AC35`
ac35$sd<-ac35.sd$`ECorAbunZero$Cytospora.sp..AC35`
ac35$se2<-ac35$se*0.25
#newplot fig 4
ac35.p<- ggplot(ac35, aes(x=time, y=ac35mean, fill=temp)) + 
  geom_bar(stat="identity",position="dodge",width=0.5,color="black") +
  geom_errorbar(aes(ymin=ac35mean, ymax=ac35mean+se2), width=0.2,
                position=position_dodge(0.5))+
  scale_x_discrete(breaks=c("25 degree","4 degree"),labels=c("25 degree","4 degree"))+
  scale_fill_manual(values=c("gray27", "gray60"), name = "Incubation temperature", 
                    labels=c(expression(25*degree~C),expression(4*degree~C))) +
  xlab("Sampling time")+ylab("Mean frequency of Cytospora sp. per sample")+
  scale_x_discrete(breaks= c(1,2,3),labels = c("May", "June","July"))+
  theme_bw()+
  theme(legend.position = "top",axis.text= element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.title.x = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", color = "black"))



olivaceaplot= ggplot(data = ECorMetaZero,aes(x=ECorMetaZero$TIME,y=ECorAbunZero$Microsphaeriopsis.olivacea,
                                              fill= factor(ECorMetaZero$TEMPRATURE)))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  scale_fill_manual(values=c("gray27", "gray60"), name = "Incubation temperature", 
                    labels=c("4 degree", "25 degree")) +guides(fill=FALSE)+
  xlab("Sampling time")+ylab("M. olivacea")+
  scale_x_continuous(breaks= c(1,2,3),labels = c("May", "June","July"))+
  theme_bw()+
  theme(axis.text= element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.title.x = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", color = "black"))
#temp*site interaction  for m.o
olivaceaplo2= ggplot(data = ECorMetaZero,aes(x=ECorMetaZero$LOCALITY,y=ECorAbunZero$Microsphaeriopsis.olivacea,
                                             fill= factor(ECorMetaZero$TEMPRATURE)))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  scale_x_discrete(breaks=c("4 degree","25 degree"),labels=c("4 degree","25 degree"))+
  scale_fill_manual(values=c("gray27", "gray60"), name = "Incubation temperature", 
                    labels=c(expression(4*degree~C),expression(25*degree~C)))+guides(fill=FALSE)+
  xlab("Sampling site")+ylab("M. olivacea")+
  scale_x_discrete(breaks= c("Bisoton","Hasan abad","Khosro abad","Sarmil","Sorkhe dize"),
                   labels = c("Biseton","Mahi Dasht","Khosro abad","Sarmil","Sorkhe dize"))+
  theme_bw()+
  theme(legend.position = "top",axis.text= element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.title.x = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", color = "black"))

############################################333
# fig 4 R2
##new plot data
ac35<-aggregate(ECorAbunZero$Cytospora.sp..AC35~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,mean)

ac35.se<-aggregate(ECorAbunZero$Cytospora.sp..AC35~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,se)
ac35.sd<-aggregate(ECorAbunZero$Cytospora.sp..AC35~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,sd)

ac35$time<-ac35$`ECorMetaZero$TIME`
ac35$`ECorMetaZero$TIME`<-NULL
ac35$temp<-ac35$`ECorMetaZero$TEMPRATURE`
ac35$`ECorMetaZero$TEMPRATURE`<-NULL
ac35$ac35mean<-ac35$`ECorAbunZero$Cytospora.sp..AC35`
ac35$`ECorAbunZero$Cytospora.sp..AC35`<-NULL
ac35$time<-as.factor(ac35$time)
ac35$se<-ac35.se$`ECorAbunZero$Cytospora.sp..AC35`
ac35$sd<-ac35.sd$`ECorAbunZero$Cytospora.sp..AC35`
ac35$se2<-ac35$se*0.25

#newplot fig 4
ac35.p<- ggplot(ac35, aes(x=time, y=ac35mean, fill=temp)) + 
  geom_bar(stat="identity",position="dodge",width=0.5,color="black") +
  geom_errorbar(aes(ymin=ac35mean, ymax=ac35mean+se2), width=0.2,
                position=position_dodge(0.5))+
  scale_x_discrete(breaks=c("25 degree","4 degree"),labels=c("25 degree","4 degree"))+
  scale_fill_manual(values=c("gray27", "gray60"), name = "Incubation temperature", 
                    labels=c(expression(25*degree~C),expression(4*degree~C))) +
  xlab("Sampling time")+ylab("Cytospora sp.(mean frequency per sample)")+
  scale_x_discrete(breaks= c(1,2,3),labels = c("May", "June","July"))+
  theme_bw()+
  theme(legend.position = "top",axis.text= element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.title.x = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", color = "black"))

## mo plot
mo.tt<-aggregate(ECorAbunZero$Microsphaeriopsis.olivacea~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,mean)
mo.tt.se<-aggregate(ECorAbunZero$Microsphaeriopsis.olivacea~ECorMetaZero$TIME+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,se)
mo.tt$se<-mo.tt.se$`ECorAbunZero$Microsphaeriopsis.olivacea`
mo.tt$time<-mo.tt$`ECorMetaZero$TIME`
mo.tt$`ECorMetaZero$TIME`<-NULL
mo.tt$temp<-mo.tt$`ECorMetaZero$TEMPRATURE`
mo.tt$`ECorMetaZero$TEMPRATURE`<-NULL
mo.tt$mean<-mo.tt$`ECorAbunZero$Microsphaeriopsis.olivacea`
mo.tt$`ECorAbunZero$Microsphaeriopsis.olivacea`<-NULL
mo.tt$se2<-mo.tt$se*0.5
#new time plot
MO.T.p<- ggplot(mo.tt, aes(x=time, y=mean, fill=temp)) + 
  geom_bar(stat="identity",position="dodge",width=0.5,color="black") +
  geom_errorbar(aes(ymin=mean, ymax=mean+se2), width=0.2,
                position=position_dodge(0.5))+
  scale_fill_manual(values=c("gray27", "gray60"), name = "Incubation temperature", 
                    labels=c("25 degree","4 degree")) +guides(fill=FALSE)+
  xlab("Sampling time")+ylab("M. olivacea(mean frequency per sample)")+
  scale_x_continuous(breaks= c(1,2,3),labels = c("May", "June","July"))+
  theme_bw()+
  theme(axis.text= element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.title.x = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", color = "black"))

#### site plot mo
mo.tt2<-aggregate(ECorAbunZero$Microsphaeriopsis.olivacea~ECorMetaZero$LOCALITY+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,mean)
mo.tt.se2<-aggregate(ECorAbunZero$Microsphaeriopsis.olivacea~ECorMetaZero$LOCALITY+ECorMetaZero$TEMPRATURE, data= ECorMetaZero,se)

mo.tt2$site<-mo.tt2$`ECorMetaZero$LOCALITY`
mo.tt2$`ECorMetaZero$LOCALITY`<-NULL
mo.tt2$temp<-mo.tt2$`ECorMetaZero$TEMPRATURE`
mo.tt2$`ECorMetaZero$TEMPRATURE`<-NULL
mo.tt2$mean<-mo.tt2$`ECorAbunZero$Microsphaeriopsis.olivacea`
mo.tt2$se<-mo.tt.se2$`ECorAbunZero$Microsphaeriopsis.olivacea`
mo.tt2$`ECorAbunZero$Microsphaeriopsis.olivacea`<-NULL
mo.tt2$se2<-mo.tt2$se*0.5

#new plot
MO.T.p2<- ggplot(mo.tt2, aes(x=site, y=mean, fill=temp)) + 
  geom_bar(stat="identity",position="dodge",width=0.5,color="black") +
  geom_errorbar(aes(ymin=mean, ymax=mean+se2), width=0.2,
                position=position_dodge(0.5))+
  scale_fill_manual(values=c("gray27", "gray60"), name = "Incubation temperature", 
                    labels=c("25 degree","4 degree")) +guides(fill=FALSE)+
  xlab("Sampling site")+ylab("M. olivacea(mean frequency per sample)")+
  scale_x_discrete(breaks= c("Bisoton","Hasan abad","Khosro abad","Sarmil","Sorkhe dize"),
                   labels = c("Biseton","Mahi Dasht","Khosro abad","Sarmil","Sorkhe dize"))+
  theme_bw()+
  theme(legend.position = "top",axis.text= element_text(size = rel(1)),
        axis.title.y = element_text(size = rel(1)),
        axis.title.x = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        legend.title = element_text(size = rel(1)),
        panel.grid.major = element_blank(),
        axis.line = element_line(colour = "black", color = "black"))


### export fig 4 r2
dev.off()
png(filename = "Fig4-r2.png", res= 300, width = 2100, height=3000)


grid.arrange(ac35.p, MO.T.p,MO.T.p2,nrow=3, ncol=1)

dev.off()



######################################################
## FIG 4 forest pathology
######################################################
library(gridExtra)
dev.off()
png(filename = "Fig3-final.png", res= 300, width = 2100, height=3000)


grid.arrange(Cytosporaplot, olivaceaplot,olivaceaplo2, nrow=3, ncol=1)
             
dev.off()


##########################################################
############ FIG 3-2
dev.off()
png(filename = "Fig3-2.png", res= 300, width = 2100, height=3000)
par(mfrow=c(3,1),xpd=TRUE, oma= c(0,0,2,0))
Cytosporaplot
olivaceaplot
olivaceaplo2
dev.off()
################################
#### 5- NMDS and similarities
################################
ENdoNMDS<-metaMDS(ECorAbunZero, distance = "bray", k= 2, trymax = 20)

#############################
## Plot NMDS for localities
#############################
dev.off()
plot(ENdoNMDS$points, xlab="NMDs 1", ylab="NMDS 2")
ordiplot(ENdoNMDS, type = "n", display = "sites",xlab="NMDS 1", ylab="NMDS 2",
         xlim= c(-2.5,2.5),ylim = c(-1.5,1.5))
points(ENdoNMDS$points, pch=20, col= "black", cex=0.5)

with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$LOCALITY,cex=.5, 
                       draw="polygon", col="blue",
                       alpha=100,kind="se",conf=0.95, 
                       show.groups=(c("Hasan abad"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$LOCALITY,cex=.5, 
                 draw="polygon", col=c("red"),
                   alpha=100,kind="se",conf=0.95, 
                   show.groups=(c("Khosro abad"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$LOCALITY,cex=.5, 
                draw="polygon", col=c("green"),
                  alpha=100,kind="se",conf=0.95, 
                show.groups=(c("Sarmil"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$LOCALITY,cex=.5, 
            draw="polygon", col=c("yellow"),
                      alpha=100,kind="se",conf=0.95, 
                   show.groups=(c("Sorkhe dize"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$LOCALITY,cex=1.5,
                              draw="polygon", col= "black",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("Bisoton"))))#red
mylegend = legend("topright", c("Bisoton", "Hasan abad", "Khosro abad",
                                "Sarmil","Sorkhe dize"), 
                  fill=c("black","blue","red","green","yellow"), border="white", bty="n")

## stress is not so good
## try 3D NMDS
NMDS.3<-metaMDS(ECorAbunZero, distance = "bray", k= 3, trymax = 20)
## stress is much better = 0.025
### Plot nmds1 &2
dev.off()
#pdf(file = "3D NMDS for localities.pdf", paper = "a4", width = 7, height = 4)
par(mfrow= c(1,3))
NMDS1.2=ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",xlab="NMDS 1",
                 ylab="NMDS 2"
         ,ylim = c(-2,2), xlim = c(-3,3))
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.8)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$LOCALITY, cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Hasan abad"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Khosro abad"))))
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sarmil"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("yellow"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sorkhe dize"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5,
                              draw="polygon", col= "black",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("Bisoton"))))#red
mylegend = legend("topright", c("Bisoton", "Hasan abad", "Khosro abad",
                                "Sarmil","Sorkhe dize"), 
                  fill=c("black","blue","red","green","yellow"), border="white", bty="n")
## plot nmds2&3
NMDS2.3=ordiplot(NMDS.3,choices=c(2,3), type = "n", display = "sites",xlab="NMDS 2", 
                 ylab="NMDS 3"
                 ,ylim = c(-1.5,1.5), xlim = c(-2,2))
points(NMDS.3$points[,2],NMDS.3$points[,3], pch=20, col= "black", cex=0.3)
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Hasan abad"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("black"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Khosro abad"))))
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sarmil"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("yellow"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sorkhe dize"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$LOCALITY,cex=1.5,
                              draw="polygon", col= "darkred",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("Bisoton"))))#red
### plot nmds 1&3
NMDS1.3=ordiplot(NMDS.3,choices=c(1,3), type = "n", display = "sites",xlab="NMDS 1", 
                 ylab="NMDS 3"
                 ,ylim = c(-1.5,1.5), xlim = c(-2,2))
points(NMDS.3$points[,1],NMDS.3$points[,3], pch=20, col= "black",cex=0.3)
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Hasan abad"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("black"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Khosro abad"))))
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sarmil"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("yellow"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sorkhe dize"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$LOCALITY,cex=1.5,
                              draw="polygon", col= "darkred",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("Bisoton"))))#red
mylegend = legend("topright", c("Bisoton", "Hasan abad", "Khosro abad",
                                "Sarmil","Sorkhe dize"), cex=0.5,
                  fill=c("red","blue","black","green","yellow"), 
                  border="white", bty="n")
dev.off()
################################
####### NMDS for time of sampling
################################
dev.off()
 plot(ENdoNMDS$points, xlab="dimension 1", ylab="dimension 2")
 ordiplot(ENdoNMDS, type = "n", display = "sites",xlab="NMDS 1", ylab="NMDS 2", 
         ylim = c(-2,2), xlim = c(-3,3))
points(ENdoNMDS$points, pch=20, col= as.numeric(ECorMetaZero$TIME))
#ordispider(ENdoNMDS,ECorMetaZero$TIME, col=c("grey"))

with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col="black",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("1"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS, ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("2"))))
with(ECorMetaZero,ordiellipse(ENdoNMDS,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("3"))))
legend("topright", c("May 2015","Jun 2015","July 2015"), 
                  fill=c("black","red","green"), 
                  border="white", bty="n")


##3D NMDS

## Plot nmds1 &2
dev.off()
#pdf(file = "3D NMDS time.pdf", paper = "a4", width = 7, height = 4)
#par(mfrow= c(1,3))
NMDS1.2=ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",xlab="NMDS 1", 
                 ylab="NMDS 2", xlim = c(-2,2)
                 )
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.75)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("1"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("2"))))
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("3"))))
legend("topright", c("May 2015","Jun 2015","July 2015"), 
       fill=c("blue","red","green"), 
       border="white", bty="n")
## plot nmds2&3
NMDS2.3=ordiplot(NMDS.3,choices=c(2,3), type = "n", display = "sites",xlab="Dimension 2", 
                 ylab="Dimension 3", ylim = c(-1,1), xlim = c(-1.5,1.5)
                 )
points(NMDS.3$points[,2],NMDS.3$points[,3], pch=20, col= "black", cex=0.3)
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("1"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("2"))))
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("3"))))
legend("topright", c("May 2015","Jun 2015","July 2015"), 
       fill=c("blue","red","green"), 
       border="white", bty="n")
## this shows the changes much better than the others

### plot nmds 1&3
NMDS1.3=ordiplot(NMDS.3,choices=c(1,3), type = "n", display = "sites",xlab="Dimension 1", 
                 ylab="Dimension 3"
                 ,ylim = c(-1.5,1.5), xlim = c(-2,2))
points(NMDS.3$points[,1],NMDS.3$points[,3], pch=20, col= "black",cex=0.3)
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("1"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("2"))))
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("green"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("3"))))
legend("topright", c("May 2015","Jun 2015","July 2015"), 
       fill=c("blue","red","green"), 
       border="white", bty="n")

dev.off()

############################
### 3D NMDS for tissue type

NMDS1.2=ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",xlab="NMDS 1", 
                 ylab="NMDS 2", xlim = c(-2,4)
)
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.75)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$SUBSTRATE, 
                              draw="polygon", col="green",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Leaf"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$SUBSTRATE, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Branch"))))

legend("topright", c("Leaf","Branch"), 
       fill=c("green","red"), 
       border="white", bty="n")
## plot nmds2&3
NMDS2.3=ordiplot(NMDS.3,choices=c(2,3), type = "n", display = "sites",xlab="NMDS 2", 
                 ylab="NMDS 3", ylim = c(-1,1), xlim = c(-1.5,1.5)
)
points(NMDS.3$points[,2],NMDS.3$points[,3], pch=20, col= "black", cex=0.75)
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$SUBSTRATE,cex=.5, 
                              draw="polygon", col="green",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Leaf"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$SUBSTRATE,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Branch"))))
legend("topright", c("Leaf","Branch"), 
       fill=c("green","red"), 
       border="white", bty="n")

### plot nmds 1&3
NMDS1.3=ordiplot(NMDS.3,choices=c(1,3), type = "n", display = "sites",xlab="NMDS 1", 
                 ylab="NMDS 3"
                 ,ylim = c(-1.5,1.5), xlim = c(-2,4))
points(NMDS.3$points[,1],NMDS.3$points[,3], pch=20, col= "black",cex=0.75)
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$SUBSTRATE,cex=.5, 
                              draw="polygon", col="green",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Leaf"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$SUBSTRATE,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Branch"))))
legend("topright", c("Leaf","Branch"), 
       fill=c("green","red"), 
       border="white", bty="n")

dev.off()

################
## NMDS for temprature

NMDS1.2=ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",xlab="NMDS 1", 
                 ylab="NMDS 2", xlim = c(-2,2)
)
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.75)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$TEMPRATURE, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("4 Degree"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$TEMPRATURE, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("25 Degree"))))

legend("topright", c("4 Degree","25 Degree"), 
       fill=c("blue","red"), 
       border="white", bty="n")
## plot nmds2&3
NMDS2.3=ordiplot(NMDS.3,choices=c(2,3), type = "n", display = "sites",xlab="NMDS 2", 
                 ylab="NMDS 3", xlim = c(-2,2)
)
points(NMDS.3$points[,2],NMDS.3$points[,3], pch=20, col= "black", cex=0.75)
with(ECorMetaZero,ordiellipse(NMDS2.3,ECorMetaZero$TEMPRATURE,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("4 Degree"))))
with(ECorMetaZero,ordiellipse(NMDS2.3, ECorMetaZero$TEMPRATURE,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("25 Degree"))))
legend("topright", c("4 Degree","25 Degree"), 
       fill=c("blue","red"), 
       border="white", bty="n")

### plot nmds 1&3
NMDS1.3=ordiplot(NMDS.3,choices=c(1,3), type = "n", display = "sites",xlab="NMDS 1", 
                 ylab="NMDS 3"
                 ,xlim = c(-2,2))
points(NMDS.3$points[,1],NMDS.3$points[,3], pch=20, col= "black",cex=0.75)
with(ECorMetaZero,ordiellipse(NMDS1.3,ECorMetaZero$TEMPRATURE,cex=.5, 
                              draw="polygon", col="blue",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("4 Degree"))))
with(ECorMetaZero,ordiellipse(NMDS1.3, ECorMetaZero$TEMPRATURE,cex=.5, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("25 Degree"))))
legend("topright", c(expression(4*degree~C),expression(25*degree~C)), 
       fill=c("blue","red"), 
       border="white", bty="n")

#########################################################
############# FIG 2
#########################################################
dev.off()

png(file = "fig 2-final.png",  width = 2100, height = 2100, res= 300)
par(mfrow= c(2,2) , mar=c(4,4,0.5,0.5), oma= c(2,0,2,0))
# tissue
ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites", 
                 ylab="NMDS 2", xlim = c(-2,4),ylim = c(-2,2)
)
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.75)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$SUBSTRATE, 
                              draw="polygon", col="forestgreen",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Leaf"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$SUBSTRATE, 
                              draw="polygon", col=c("brown"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Branch"))))

legend("topright", c("Leaf","Twig"), fill=c("forestgreen","brown"),
       border="white", bty="n",title = "A")
## location
ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",
         xlim = c(-2,4),ylim = c(-2,2))
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.75)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$LOCALITY, cex=.5, 
                              draw="polygon", col="cyan",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Hasan abad"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("forestgreen"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Khosro abad"))))
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("yellow"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sarmil"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5, 
                              draw="polygon", col=c("maroon3"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("Sorkhe dize"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$LOCALITY,cex=.5,
                              draw="polygon", col= "black",
                              alpha=100, kind="se", conf=0.95, 
                              show.groups=(c("Bisoton"))))#red
mylegend = legend("topright", c("Biseton", "Mahi Dasht", "Khosro abad",
                                "Sarmil","Sorkhe dize"), title = "B",
                  fill=c("black","cyan","forestgreen","yellow","maroon3"), border="white", bty="n")

### Time
ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",xlab="NMDS 1", 
                 ylab="NMDS 2", xlim = c(-2,4),ylim = c(-2,2)
)
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.75)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col="cyan",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("1"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("maroon3"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("2"))))
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$TIME,cex=.5, 
                              draw="polygon", col=c("forestgreen"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("3"))))
legend("topright", c("May 2015","June 2015","July 2015"), title = "C",
       fill=c("cyan","maroon3","forestgreen"), 
       border="white", bty="n")
### Temprature
ordiplot(NMDS.3,choices=c(1,2), type = "n", display = "sites",xlab="NMDS 1", 
         xlim = c(-2,4),ylim = c(-2,2)
)
points(NMDS.3$points[,1],NMDS.3$points[,2], pch=20, col= "black", cex= 0.75)
with(ECorMetaZero,ordiellipse(NMDS1.2,ECorMetaZero$TEMPRATURE, 
                              draw="polygon", col="cyan",
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("4 Degree"))))
with(ECorMetaZero,ordiellipse(NMDS1.2, ECorMetaZero$TEMPRATURE, 
                              draw="polygon", col=c("red"),
                              alpha=100,kind="se",conf=0.95, 
                              show.groups=(c("25 Degree"))))

legend("topright", c(expression(4*degree~C),expression(25*degree~C)), title = "D",
       fill=c("cyan","red"), 
       border="white", bty="n")
dev.off()
###############################################################
#### 6- community analysis at family level
###############################################################
## Reporting which factors are affecting the communities at family level
# Data input
### many OTUs were identified at famiy level so I chose this level to do this analysis
### the OTU abundances for each family was merged together
### out of 59 OTUs we were not able to identify 16 OTU and were grouped according to
#the assignment level

FmilyAbun= read.csv (file="matrix_family.csv", 
                    header = T, row.names = 1)
row.names(FmilyAbun)==row.names(EndoMetaData)
colnames(FmilyAbun)
FamilyMVABUND= mvabund(FmilyAbun)
plot(FamilyMVABUND)

FamilyMV.m= manyglm(FamilyMVABUND ~ Locality+ Time+ Temperature+ Tissue, data = EndoMetaData,
                    family = "negative.binomial", show.residuals=T)
FamilyMV.Anova= anova.manyglm(FamilyMV.m,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo" )
FamilyMV.summ= summary.manyglm(FamilyMV.m,nBoot=100, test="LR", p.uni="adjusted", 
                              resamp="montecarlo")
### which families are significantly affected?

FamilyAnova <- as.data.frame(FamilyMV.Anova$uni.p)
FmilyTissue<-colnames(FamilyAnova)[FamilyAnova["Tissue",]<= 0.05]
FmilyLoc<-colnames(FamilyAnova)[FamilyAnova["Locality",]<= 0.05]
FmilyTim<-colnames(FamilyAnova)[FamilyAnova["Time",]<= 0.05]
FmilyTemp<-colnames(FamilyAnova)[FamilyAnova["Temperature",]<= 0.05]


# get the sum for each family in each tisse
family.sum <- aggregate(. ~ EndoMetaData$SUBSTRATE,FmilyAbun , sum)
aggregate(. ~ EndoMetaData$TIME,FmilyAbun , sum)
###########################################################
####### Fig 2-paper 
###########################################################
dev.off()
#pdf(file = "PIE CHART L&B final.pdf", paper = "a4", width = 7, height = 3)
jpeg(file = "Fig 4persian.jpeg", width = 1500, height = 2500, res = 300)
par(mfrow=c(2,1))
# creat a Pie chart for leaf
L.slic<- c(1,1,2,25,1,2,1,9,6,7,2,3)#get the valus from family.sum and remove the zeros
L.lbls<- c("Coniochaetaceae","Cladosporiaceae","Dothioraceae","Gnomoniaceae",
           "Lasiosphaeriaceae","Montagnulaceae","Nectriaceae","Pleosporaceae",
           "Sporormiaceae","Trichocomaceae","Xylariaceae","Other families")
L.Percent<-round(L.slic/sum(L.slic)*100, digits=2)
L.lbls <- paste(L.lbls, L.Percent)
L.lbls<-paste(L.lbls,"%",sep="")
#png(file = "Fig 2-revised-leaf.png", width = 1000, height = 1000, res = 300)
pie(L.slic,labels =L.lbls, col = c("firebrick","indianred1","orchid1","magenta",
                                   "deeppink1","mediumblue","royalblue1","skyblue1","cyan",
                                  "yellow", "springgreen2","forestgreen"
                                   ) , main = "برگ", 
border = NA, radius = 0.7, cex=0.7,cex.main=0.7)

# creat a Pie chart for branch
B.slic<- c(6,45,157,6,37,48,85,15,60)#get the valus from family.sum and remove the zeros
B.lbls<- c("Coniochaetaceae","Dothioraceae","Montagnulaceae",
           "Cladosporiaceae","Pleosporaceae","Sporormiaceae","Trichocomaceae","Valsaceae",
           "Other families")
B.Percent<-round(B.slic/sum(B.slic)*100, digits=2)
B.lbls <- paste(B.lbls, B.Percent)
B.lbls<-paste(B.lbls,"%",sep="")
pie(B.slic,labels =B.lbls, col = c("firebrick","orchid1","mediumblue","indianred1",
                                   "skyblue1", "cyan","yellow","gold","forestgreen"),
    radius = 0.7, main = "شاخه",border = NA,cex=0.7,cex.main=0.7)
dev.off()
###################################################################
############## fig 5
###################################################################
                                                                                    
familtemp<- c ("Dothioraceae","Gnomoniaceae","Montagnulaceae","Pleosporaceae",
               "Sporormiaceae","Trichocomaceae")
familiTempABUN<- FmilyAbun[familtemp]##Keeping only Temp affected families
# get the mean valuse for each family in each temp
familiTempmean <- aggregate(. ~ EndoMetaData$TEMPRATURE, familiTempABUN, mean)
#CReat a data frame of the mean valuse
famiTMeanfram<- as.data.frame(familiTempmean,optional=TRUE)
attr(famiTMeanfram, "row.names")<- c("25 Degree", "4 Degree")

### creat a matrix of mean observation of families affected by temp for ploting
familiTdata<- data.matrix (famiTMeanfram[2:7],rownames.force = NA )
familiTdata2<-as.data.frame(familiTdata)
library(ggplot2)
dev.off()
png(file = "Fig 5-600.png", width = 3000, height = 2400, res = 600)
par(mar=c(6,4,1,0),oma=c(1,0,1,0))
bar.temp<-barplot(familiTdata,legend.text =c(expression(25*degree~C),expression(4*degree~C)),
        beside = TRUE,ylab= "Mean frequency per sample" ,xlab = "Fungal taxonomic families",
        names.arg= c ("Dothioraceae","Gnomoniaceae","Montagnulaceae","Pleosporaceae",
                      "Sporormiaceae","Trichocomaceae"), 
        axes= TRUE,ylim= c(0,0.7), 
        args.legend = list(x = "topright",bty= "n", cex=0.7), border = "Black", axis.lty=1,axisnames=T,
        width = 0.5)
error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
  arrows(x,y+upper, x, y-lower, angle=90, code=3,width= 0.5, length=length, ...)
}
error.bar(bar.temp,familiTdata,familiTdata2se)
dev.off()
?barplot()
??args.legend
# add error bars
#First I create a smell function that takes...in entry

# get the se valuse for each family in each temp
familiTempse <- aggregate(. ~ EndoMetaData$TEMPRATURE, familiTempABUN, se)
#CReat a data frame of the mean valuse
famiTsefram<- as.data.frame(familiTempse,optional=TRUE)
attr(famiTsefram, "row.names")<- c("25 Degree", "4 Degree")

### creat a matrix of mean observation of families affected by temp for ploting
familiTdatase<- data.matrix (famiTsefram[2:7],rownames.force = NA )
familiTdata2se<-familiTdatase*0.5

error.bar(bar.temp,familiTdata,familiTdata2se)

##################
## OTU isolation results in tissue and temperature
## OTU observation in each tissue type
TisOTUsum <- aggregate(. ~ EndoMetaData$SUBSTRATE, EndoAbun, sum)



###OTU observation in each isolation temperature
TempOTUsum<-aggregate(. ~ EndoMetaData$TEMPRATURE, EndoAbun, sum)
class(TempOTUsum)
library(plyr)
library(dplyr)
?tally
dev.off()


#################
######VEN diagram
library(VennDiagram)
EndoMetaData
EndoAbun
par(mfrow=c(1,2))
jpeg("fig6-per-organ.jpeg", width = 1000,height = 1000, res=300)
venn.organ <- draw.pairwise.venn(
  area1 = 22,
  area2 = 51,
  cross.area = 14,
  category = c("Leaf", "Twig"),
  fill = c("green", "gray"),cex =1.2,cat.cex = 1.2, 
  scaled=T,cat.pos = c(180,180),cat.dist = 0.05, lwd=c(0.5,0.5),col=c("green", "gray"))
grid.newpage()
?draw.pairwise.venn
venn.temp<-draw.pairwise.venn(
  area1 = 51,
  area2 = 15,
  cross.area = 7,
  category = c(expression(25*degree~C), expression(4*degree~C)),
  fill = c("red", "blue"), cex =1.2,cat.cex = 1.2, 
  scaled=T,cat.pos = c(0,0),cat.dist = 0.05,lwd=c(0.5,0.5),col=c("red", "blue"))
dev.off()



library(gridExtra)

# Fig 6
png(file = "Fig 6.png", width = 2000, height = 900, res = 300)
grid.arrange(gTree(children=venn.temp), gTree(children=venn.organ), ncol=2)

dev.off()

#another way for venn diagram
source("http://www.bioconductor.org/biocLite.R")    
biocLite("limma")
library(limma)

vennDiagram(vennCounts(EndoMetaData))
?vennDiagram

#subseting for site venn diagram
metadata2<-data.frame(out1=integer(446))

metadata2<-cbind(metadata2,site=EndoMetaData$LOCALITY)
metadata2<- cbind(metadata2,otu1=ifelse(EndoAbun$Alternaria.sp..A22 > 0, 1, 0))
metadata2<- cbind(metadata2,otu2=ifelse(EndoAbun$Alternaria.sp..Ac27 > 0, 1, 0))
metadata2<- cbind(metadata2,otu3=ifelse(EndoAbun$Alternaria.sp..A25 > 0, 1, 0))
metadata2<- cbind(metadata2,otu4=ifelse(EndoAbun$Alternaria.sp..A76 > 0, 1, 0))
metadata2<- cbind(metadata2,otu5=ifelse(EndoAbun$Alternaria.sp..A9 > 0, 1, 0))
metadata2<- cbind(metadata2,otu6=ifelse(EndoAbun$Arthrinium.marii > 0, 1, 0))
metadata2<- cbind(metadata2,otu7=ifelse(EndoAbun$Aspergillus.sp..A20 > 0, 1, 0))
metadata2<- cbind(metadata2,otu8=ifelse(EndoAbun$Aureobasidium.sp..A17 > 0, 1, 0))
metadata2<- cbind(metadata2,otu9=ifelse(EndoAbun$Aureobasidium.sp..A30 > 0, 1, 0))
metadata2<- cbind(metadata2,otu10=ifelse(EndoAbun$B25 > 0, 1, 0))
metadata2<- cbind(metadata2,otu11=ifelse(EndoAbun$Bc3 > 0, 1, 0))
metadata2<- cbind(metadata2,otu12=ifelse(EndoAbun$Biscogniauxia.mediterranea > 0, 1, 0))
metadata2<- cbind(metadata2,otu13=ifelse(EndoAbun$Byssochlamys.spectabilis. > 0, 1, 0))
metadata2<- cbind(metadata2,otu14=ifelse(EndoAbun$Chaetomiaceae.sp..A37 > 0, 1, 0))
metadata2<- cbind(metadata2,otu15=ifelse(EndoAbun$Cladosporium.herbarum.A8 > 0, 1, 0))
metadata2<- cbind(metadata2,otu16=ifelse(EndoAbun$Comoclathris.sedi > 0, 1, 0))
metadata2<- cbind(metadata2,otu17=ifelse(EndoAbun$Coniochaeta.sp..A85 > 0, 1, 0))
metadata2<- cbind(metadata2,otu18=ifelse(EndoAbun$Coniothyrium.sp..A41 > 0, 1, 0))
metadata2<- cbind(metadata2,otu19=ifelse(EndoAbun$Cytospora.sp..AC35 > 0, 1, 0))
metadata2<- cbind(metadata2,otu20=ifelse(EndoAbun$Cytospora.sp..C2 > 0, 1, 0))
metadata2<- cbind(metadata2,otu21=ifelse(EndoAbun$Diatrype.sp..C1 > 0, 1, 0))
metadata2<- cbind(metadata2,otu22=ifelse(EndoAbun$Diatrypella.sp..C6 > 0, 1, 0))
metadata2<- cbind(metadata2,otu23=ifelse(EndoAbun$Dikarya.sp..A38 > 0, 1, 0))
metadata2<- cbind(metadata2,otu24=ifelse(EndoAbun$Dothideomycetes.sp..A1 > 0, 1, 0))
metadata2<- cbind(metadata2,otu25=ifelse(EndoAbun$Dothideomycetes.sp..A79 > 0, 1, 0))
metadata2<- cbind(metadata2,otu26=ifelse(EndoAbun$Endoconidioma.populi.A39 > 0, 1, 0))
metadata2<- cbind(metadata2,otu27=ifelse(EndoAbun$Fusarium.sp..46 > 0, 1, 0))
metadata2<- cbind(metadata2,otu28=ifelse(EndoAbun$Gnomoniaceae.sp..66 > 0, 1, 0))
metadata2<- cbind(metadata2,otu29=ifelse(EndoAbun$Gnomoniaceae.sp..70 > 0, 1, 0))
metadata2<- cbind(metadata2,otu30=ifelse(EndoAbun$Humicola.sp..A52 > 0, 1, 0))
metadata2<- cbind(metadata2,otu31=ifelse(EndoAbun$Microsphaeriopsis.olivacea > 0, 1, 0))
metadata2<- cbind(metadata2,otu32=ifelse(EndoAbun$Penicillium.sp..A21 > 0, 1, 0))
metadata2<- cbind(metadata2,otu33=ifelse(EndoAbun$Penicillium.sp..A3 > 0, 1, 0))
metadata2<- cbind(metadata2,otu34=ifelse(EndoAbun$Pleosporaceae.sp..A5 > 0, 1, 0))
metadata2<- cbind(metadata2,ot35=ifelse(EndoAbun$Pleosporaceae.sp..Ac49 > 0, 1, 0))
metadata2<- cbind(metadata2,otu36=ifelse(EndoAbun$Pleosporaceae.sp..B27 > 0, 1, 0))
metadata2<- cbind(metadata2,otu37=ifelse(EndoAbun$Preussia.africana > 0, 1, 0))
metadata2<- cbind(metadata2,otu38=ifelse(EndoAbun$Preussia.australis > 0, 1, 0))
metadata2<- cbind(metadata2,otu39=ifelse(EndoAbun$Preussia.complex.sp..A36 > 0, 1, 0))
metadata2<- cbind(metadata2,otu40=ifelse(EndoAbun$Preussia.intermedia > 0, 1, 0))
metadata2<- cbind(metadata2,otu41=ifelse(EndoAbun$Preussia.sp..A31 > 0, 1, 0))
metadata2<- cbind(metadata2,otu42=ifelse(EndoAbun$Schizothecium.sp..B14 > 0, 1, 0))
metadata2<- cbind(metadata2,otu43=ifelse(EndoAbun$Sordariomycetes.sp..A13 > 0, 1, 0))
metadata2<- cbind(metadata2,otu44=ifelse(EndoAbun$Sordariomycetes.sp..A45 > 0, 1, 0))
metadata2<- cbind(metadata2,otu45=ifelse(EndoAbun$Sporormiaceae.sp..A29M > 0, 1, 0))
metadata2<- cbind(metadata2,otu46=ifelse(EndoAbun$Sporormiaceae.sp..L32 > 0, 1, 0))
metadata2<- cbind(metadata2,otu147=ifelse(EndoAbun$Unfidentified.Sordariomycetes.sp..A26 > 0, 1, 0))
metadata2<- cbind(metadata2,ot48=ifelse(EndoAbun$Unidentified.Ascomycota.sp..A86 > 0, 1, 0))
metadata2<- cbind(metadata2,otu49=ifelse(EndoAbun$Unidentified.Ascomycota.sp..Bc6 > 0, 1, 0))
metadata2<- cbind(metadata2,otu50=ifelse(EndoAbun$Unidentified.Dothideomycetes.sp..Bc4 > 0, 1, 0))
metadata2<- cbind(metadata2,otu51=ifelse(EndoAbun$Unidentified.Fungi.Ac44 > 0, 1, 0))
metadata2<- cbind(metadata2,otu52=ifelse(EndoAbun$Unidentified.Fungi.Ac52 > 0, 1, 0))
metadata2<- cbind(metadata2,otu53=ifelse(EndoAbun$Unidentified.Pezizomycotina.sp..B12 > 0, 1, 0))
metadata2<- cbind(metadata2,otu54=ifelse(EndoAbun$Unidentified.Pleosporales.sp..A65 > 0, 1, 0))
metadata2<- cbind(metadata2,otu55=ifelse(EndoAbun$Unidentified.Pleosporales.sp..A75 > 0, 1, 0))
metadata2<- cbind(metadata2,otu56=ifelse(EndoAbun$Unidentified.Xylariaceae.sp..B21 > 0, 1, 0))
metadata2<- cbind(metadata2,otu57=ifelse(EndoAbun$Ustilago.A14 > 0, 1, 0))
metadata2<- cbind(metadata2,otu58=ifelse(EndoAbun$Ustilago.A16 > 0, 1, 0))
metadata2<- cbind(metadata2,otu59=ifelse(EndoAbun$Valsaceae.sp..A32 > 0, 1, 0))

write.csv(metadata2, file = "metadata2.csv")
metadata1<-read.csv("metadata2.csv",header = T,row.names = 1)
site=metadata2$site
metadata1$site<-NULL
metadata3<-t((rowsum(metadata1)) group = EndoMetaData$LOCALITY, na.rm = T)

df1<-rowsum(metadata1,group = EndoMetaData$LOCALITY, na.rm = T)

metadataREV<-t(df1)
class(metadataREV)
metadataREV<-as.data.frame(metadataREV)

metadataREV<- cbind(metadataREV,Biseton=ifelse(metadataREV$Bisoton > 0, 1, 0))
metadataREV<- cbind(metadataREV,Mahidasht=ifelse(metadataREV$`Hasan abad`> 0, 1, 0))
metadataREV<- cbind(metadataREV,KhosroAbad=ifelse(metadataREV$`Khosro abad`> 0, 1, 0))
metadataREV<- cbind(metadataREV,SARMIL=ifelse(metadataREV$Sarmil> 0, 1, 0))
metadataREV<- cbind(metadataREV,SORKHEDIZE=ifelse(metadataREV$`Sorkhe dize`> 0, 1, 0))

metadataREV$Bisoton<-NULL
metadataREV$`Hasan abad`<-NULL
metadataREV$`Khosro abad`<-NULL
metadataREV$Sarmil<-NULL
metadataREV$`Sorkhe dize`<-NULL

venn.site<-vennDiagram(vennCounts(metadataREV), 
            circle.col = c("red", "darkblue", "green","black","darkviolet"),
            names = c("Biseton","Mahi Dasht", "Khosro Abad","Sarmil","Sorkhe Dize"))
dev.off()

# Fig 6- site
jpeg(file = "Fig 6site-per.jpeg", width = 1000, height = 1000, res = 300)
venn.site<-vennDiagram(vennCounts(metadataREV), cex=c(0.5,0.5,0.5),
                       circle.col = c("red", "darkblue", "green","black","darkviolet"),
                       names = c("Biseton","Mahi Dasht", "Khosro Abad","Sarmil","Sorkhe Dize"))
dev.off()
?vennDiagram
# subset for time venn diagram

metadata3<-data.frame(out1=integer(446))
metadata3$out1<-NULL

metadata3<-cbind(metadata3,time=EndoMetaData$TIME)
metadata3<- cbind(metadata3,otu1=ifelse(EndoAbun$Alternaria.sp..A22 > 0, 1, 0))
metadata3<- cbind(metadata3,otu2=ifelse(EndoAbun$Alternaria.sp..Ac27 > 0, 1, 0))
metadata3<- cbind(metadata3,otu3=ifelse(EndoAbun$Alternaria.sp..A25 > 0, 1, 0))
metadata3<- cbind(metadata3,otu4=ifelse(EndoAbun$Alternaria.sp..A76 > 0, 1, 0))
metadata3<- cbind(metadata3,otu5=ifelse(EndoAbun$Alternaria.sp..A9 > 0, 1, 0))
metadata3<- cbind(metadata3,otu6=ifelse(EndoAbun$Arthrinium.marii > 0, 1, 0))
metadata3<- cbind(metadata3,otu7=ifelse(EndoAbun$Aspergillus.sp..A20 > 0, 1, 0))
metadata3<- cbind(metadata3,otu8=ifelse(EndoAbun$Aureobasidium.sp..A17 > 0, 1, 0))
metadata3<- cbind(metadata3,otu9=ifelse(EndoAbun$Aureobasidium.sp..A30 > 0, 1, 0))
metadata3<- cbind(metadata3,otu10=ifelse(EndoAbun$B25 > 0, 1, 0))
metadata3<- cbind(metadata3,otu11=ifelse(EndoAbun$Bc3 > 0, 1, 0))
metadata3<- cbind(metadata3,otu12=ifelse(EndoAbun$Biscogniauxia.mediterranea > 0, 1, 0))
metadata3<- cbind(metadata3,otu13=ifelse(EndoAbun$Byssochlamys.spectabilis. > 0, 1, 0))
metadata3<- cbind(metadata3,otu14=ifelse(EndoAbun$Chaetomiaceae.sp..A37 > 0, 1, 0))
metadata3<- cbind(metadata3,otu15=ifelse(EndoAbun$Cladosporium.herbarum.A8 > 0, 1, 0))
metadata3<- cbind(metadata3,otu16=ifelse(EndoAbun$Comoclathris.sedi > 0, 1, 0))
metadata3<- cbind(metadata3,otu17=ifelse(EndoAbun$Coniochaeta.sp..A85 > 0, 1, 0))
metadata3<- cbind(metadata3,otu18=ifelse(EndoAbun$Coniothyrium.sp..A41 > 0, 1, 0))
metadata3<- cbind(metadata3,otu19=ifelse(EndoAbun$Cytospora.sp..AC35 > 0, 1, 0))
metadata3<- cbind(metadata3,otu20=ifelse(EndoAbun$Cytospora.sp..C2 > 0, 1, 0))
metadata3<- cbind(metadata3,otu21=ifelse(EndoAbun$Diatrype.sp..C1 > 0, 1, 0))
metadata3<- cbind(metadata3,otu22=ifelse(EndoAbun$Diatrypella.sp..C6 > 0, 1, 0))
metadata3<- cbind(metadata3,otu23=ifelse(EndoAbun$Dikarya.sp..A38 > 0, 1, 0))
metadata3<- cbind(metadata3,otu24=ifelse(EndoAbun$Dothideomycetes.sp..A1 > 0, 1, 0))
metadata3<- cbind(metadata3,otu25=ifelse(EndoAbun$Dothideomycetes.sp..A79 > 0, 1, 0))
metadata3<- cbind(metadata3,otu26=ifelse(EndoAbun$Endoconidioma.populi.A39 > 0, 1, 0))
metadata3<- cbind(metadata3,otu27=ifelse(EndoAbun$Fusarium.sp..46 > 0, 1, 0))
metadata3<- cbind(metadata3,otu28=ifelse(EndoAbun$Gnomoniaceae.sp..66 > 0, 1, 0))
metadata3<- cbind(metadata3,otu29=ifelse(EndoAbun$Gnomoniaceae.sp..70 > 0, 1, 0))
metadata3<- cbind(metadata3,otu30=ifelse(EndoAbun$Humicola.sp..A52 > 0, 1, 0))
metadata3<- cbind(metadata3,otu31=ifelse(EndoAbun$Microsphaeriopsis.olivacea > 0, 1, 0))
metadata3<- cbind(metadata3,otu32=ifelse(EndoAbun$Penicillium.sp..A21 > 0, 1, 0))
metadata3<- cbind(metadata3,otu33=ifelse(EndoAbun$Penicillium.sp..A3 > 0, 1, 0))
metadata3<- cbind(metadata3,otu34=ifelse(EndoAbun$Pleosporaceae.sp..A5 > 0, 1, 0))
metadata3<- cbind(metadata3,ot35=ifelse(EndoAbun$Pleosporaceae.sp..Ac49 > 0, 1, 0))
metadata3<- cbind(metadata3,otu36=ifelse(EndoAbun$Pleosporaceae.sp..B27 > 0, 1, 0))
metadata3<- cbind(metadata3,otu37=ifelse(EndoAbun$Preussia.africana > 0, 1, 0))
metadata3<- cbind(metadata3,otu38=ifelse(EndoAbun$Preussia.australis > 0, 1, 0))
metadata3<- cbind(metadata3,otu39=ifelse(EndoAbun$Preussia.complex.sp..A36 > 0, 1, 0))
metadata3<- cbind(metadata3,otu40=ifelse(EndoAbun$Preussia.intermedia > 0, 1, 0))
metadata3<- cbind(metadata3,otu41=ifelse(EndoAbun$Preussia.sp..A31 > 0, 1, 0))
metadata3<- cbind(metadata3,otu42=ifelse(EndoAbun$Schizothecium.sp..B14 > 0, 1, 0))
metadata3<- cbind(metadata3,otu43=ifelse(EndoAbun$Sordariomycetes.sp..A13 > 0, 1, 0))
metadata3<- cbind(metadata3,otu44=ifelse(EndoAbun$Sordariomycetes.sp..A45 > 0, 1, 0))
metadata3<- cbind(metadata3,otu45=ifelse(EndoAbun$Sporormiaceae.sp..A29M > 0, 1, 0))
metadata3<- cbind(metadata3,otu46=ifelse(EndoAbun$Sporormiaceae.sp..L32 > 0, 1, 0))
metadata3<- cbind(metadata3,otu147=ifelse(EndoAbun$Unfidentified.Sordariomycetes.sp..A26 > 0, 1, 0))
metadata3<- cbind(metadata3,ot48=ifelse(EndoAbun$Unidentified.Ascomycota.sp..A86 > 0, 1, 0))
metadata3<- cbind(metadata3,otu49=ifelse(EndoAbun$Unidentified.Ascomycota.sp..Bc6 > 0, 1, 0))
metadata3<- cbind(metadata3,otu50=ifelse(EndoAbun$Unidentified.Dothideomycetes.sp..Bc4 > 0, 1, 0))
metadata3<- cbind(metadata3,otu51=ifelse(EndoAbun$Unidentified.Fungi.Ac44 > 0, 1, 0))
metadata3<- cbind(metadata3,otu52=ifelse(EndoAbun$Unidentified.Fungi.Ac52 > 0, 1, 0))
metadata3<- cbind(metadata3,otu53=ifelse(EndoAbun$Unidentified.Pezizomycotina.sp..B12 > 0, 1, 0))
metadata3<- cbind(metadata3,otu54=ifelse(EndoAbun$Unidentified.Pleosporales.sp..A65 > 0, 1, 0))
metadata3<- cbind(metadata3,otu55=ifelse(EndoAbun$Unidentified.Pleosporales.sp..A75 > 0, 1, 0))
metadata3<- cbind(metadata3,otu56=ifelse(EndoAbun$Unidentified.Xylariaceae.sp..B21 > 0, 1, 0))
metadata3<- cbind(metadata3,otu57=ifelse(EndoAbun$Ustilago.A14 > 0, 1, 0))
metadata3<- cbind(metadata3,otu58=ifelse(EndoAbun$Ustilago.A16 > 0, 1, 0))
metadata3<- cbind(metadata3,otu59=ifelse(EndoAbun$Valsaceae.sp..A32 > 0, 1, 0))
metadata3$time<-NULL
write.csv(metadata3, file = "metadata3.csv")

metadata4<-read.csv("metadata3.csv",header = T)
df2<-rowsum(metadata4,group = EndoMetaData$TIME, na.rm = T)

metadata.time<-t(df2)
class(metadata.time)
metadata.time<-as.data.frame(metadata.time)

metadata.time<- cbind(metadata.time,May=ifelse(metadata.time$`1` > 0, 1, 0))
metadata.time<- cbind(metadata.time,June=ifelse(metadata.time$`2`> 0, 1, 0))
metadata.time<- cbind(metadata.time,July=ifelse(metadata.time$`3`> 0, 1, 0))

metadata.time$`1`<-NULL
metadata.time$`2`<-NULL
metadata.time$`3`<-NULL

jpeg(file = "Fig 6 time-per.jpeg", width = 1000, height = 1000, res = 300)
venn.time<-vennDiagram(vennCounts(metadata.time), cex=c(0.5,0.5,0.5),
                       circle.col = c("green","darkblue","red"))
dev.off()

