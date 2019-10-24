
##############################################
##############################################
# Greenhouse data analysis
##############################################
##############################################

library(vegan)
library(mvabund)
library(rjags)
library(MASS)
library(ggplot2)
library(ggtree)
library(tidyverse)
library(mvabund)
library(reshape)
library(scales)
library(agricolae)

##############################################
#https://rcompanion.org/handbook/G_09.html

#DATA INPUT (Cucumis sativus)

CUCUMBER.DATA<-read.csv("Cucumber.csv",header = TRUE)
CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)

CUCUMIS.DATA$Biomass<-CUCUMIS.DATA$Wdshoot+CUCUMIS.DATA$Wdroot

# Total biomass (g)
anova.BiomassCUCUMBER<-aov(Biomass~Treatment.code,data=CUCUMIS.DATA)
output.BiomassCUCUMBER <- duncan.test(anova.BiomassCUCUMBER,"Treatment.code", group=TRUE)
plot(output.BiomassCUCUMBER,horiz=TRUE,las=1)
print(output.BiomassCUCUMBER$groups)


CucumisBiomass<-aov(Biomass~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisBiomass.sum<-summary(CucumisBiomass)

ggplot(CUCUMIS.DATA, aes(x=Drought,y=Biomass, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge")+theme_bw()+ 
  ylab("Biomass(g)")+theme(legend.position="none",axis.text=element_text(size=12),
                           axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
                           strip.text = element_text(size = 12))


BiomassCUCUMBER<-read.csv(file = "BiomassCUCUMBER.csv")
ggplot(data = BiomassCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Biomass")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Shoot wet weight (g)

CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)

CucumisWshoot<-aov(Wshoot~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisWshoot.sum<-summary(CucumisWshoot)

ggplot(CUCUMIS.DATA, aes(x=Drought,y=Wshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))


WshootCUCUMBER<-read.csv(file = "WshootCUCUMBER.csv")
ggplot(data = WshootCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Wshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Shoot dry weight (g)

CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)

CucumisWdshoot<-aov(Wdshoot~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisWdshoot.sum<-summary(CucumisWdshoot)


ggplot(CUCUMIS.DATA, aes(x=Drought,y=Wdshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wdshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

WdshootCUCUMBER<-read.csv(file = "WdshootCUCUMBER.csv")
ggplot(data = WdshootCUCUMBER,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("Wdshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Root wet weight (g)

CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)
CucumisWroot<-aov(Wroot~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisWroot.sum<-summary(CucumisWroot)

ggplot(CUCUMIS.DATA, aes(x=Drought,y=Wroot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))
WrootCUCUMBER<-read.csv(file = "WrootCUCUMBER.csv")
ggplot(data = WrootCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Wroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Root dry weight (g)
CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)

CucumisWdroot<-aov(Wdroot~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisWdroot.sum<-summary(CucumisWdroot)

ggplot(CUCUMIS.DATA, aes(x=Drought,y=Wdroot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wdroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

WdrootCUCUMBER<-read.csv(file = "WdrootCUCUMBER.csv")
ggplot(data = WdrootCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Wdroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Total wet weight

CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)

CucumisTotalwetweight<-aov(Tww~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisS.sum<-summary(CucumisTotalwetweight)


ggplot(CUCUMIS.DATA, aes(x=Drought,y=Tww, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Total wet weight")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# chlorophyll concentration (%)

CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)

CucumisS<-aov(S~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisS.sum<-summary(CucumisS)

ggplot(CUCUMIS.DATA, aes(x=Drought,y=S, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("chlorophyll concentration")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

SCUCUMBER<-read.csv(file = "SCUCUMBER.csv")
ggplot(data =SCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("chlorophyll concentration")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Shoot Lenght (cm)

CUCUMIS.DATA<-read.csv("Cucumis.csv",header = TRUE)

CucumisLshoot<-aov(Lshoot~Fungi*Drought*Salinity,data = CUCUMIS.DATA)
CucumisS.sum<-summary(CucumisLshoot)


ggplot(CUCUMIS.DATA, aes(x=Drought,y=Lshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Lshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# POX

CUCUMBER.DATA<-read.csv("Cucumber.csv",header = TRUE)

CucumberPOX<-aov(POX~Fungi*Drought*Salinity,data = CUCUMBER.DATA)
CucumisPOX.sum<-summary(CucumberPOX)

POXCUCUMBER<-read.csv(file = "POXCUCUMBER.csv")
ggplot(data =POXCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("POX")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# CAT

CUCUMBER.DATA<-read.csv("Cucumber.csv",header = TRUE)

CucumberCAT<-aov(CAT~Fungi*Drought*Salinity,data = CUCUMBER.DATA)
CucumisCAT.sum<-summary(CucumberCAT)

CATCUCUMBER<-read.csv(file = "CATCUCUMBER.csv")
ggplot(data =CATCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("CAT")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


CATCUCUMBER<-read.csv(file = "CATCUCUMBER.csv")
ggplot(data =CATCUCUMBER,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("CAT")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# SOD

CUCUMBER.DATA<-read.csv("Cucumber.csv",header = TRUE)

CucumberSOD<-aov(SOD~Fungi*Drought*Salinity,data = CUCUMBER.DATA)
CucumisSOD.sum<-summary(CucumberSOD)

SODCUCUMBER<-read.csv(file = "SODCUCUMBER.csv")
ggplot(data =SODCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("SOD")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


SODCUCUMBER<-read.csv(file = "SODCUCUMBER.csv")
ggplot(data =SODCUCUMBER,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("SOD")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Protein

CUCUMBER.DATA<-read.csv("Cucumber.csv",header = TRUE)

CucumberProtein<-aov(Protein~Fungi*Drought*Salinity,data = CUCUMBER.DATA)
CucumisProtein.sum<-summary(CucumberProtein)


ggplot(CUCUMBER.DATA, aes(x=Drought,y=Protein, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Protein")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

ProteinCUCUMBER<-read.csv(file = "ProteinCUCUMBER.csv")
ggplot(data =ProteinCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Protein")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


ProteinCUCUMBER<-read.csv(file = "ProteinCUCUMBER.csv")
ggplot(data =ProteinCUCUMBER,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("Protein")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Proline

CUCUMBER.DATA<-read.csv("Cucumber.csv",header = TRUE)

CucumberProline<-aov(Proline~Fungi*Drought*Salinity,data = CUCUMBER.DATA)
CucumisS.sum<-summary(CucumberProline)

ProlineCUCUMBER<-read.csv(file = "ProlineCUCUMBER.csv")
ggplot(data =ProlineCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Proline")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


ProlineCUCUMBER<-read.csv(file = "ProlineCUCUMBER.csv")
ggplot(data =ProlineCUCUMBER,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Proline")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


#DATA INPUT (Solanum lycopersicum)
Solanum.DATA<-read.csv("Solanum.csv",header = TRUE)
Tomato.DATA<-read.csv("Tomato.csv",header = TRUE)

Solanum.DATA$Biomass<-Solanum.DATA$Wdshoot+Solanum.DATA$Wdroot

# Total biomass (g)
anova.BiomassCUCUMBER<-aov(Biomass~Treatment.code,data=CUCUMIS.DATA)
output.BiomassCUCUMBER <- duncan.test(anova.BiomassCUCUMBER,"Treatment.code", group=TRUE)
plot(output.BiomassCUCUMBER,horiz=TRUE,las=1)
print(output.BiomassCUCUMBER$groups)


SolanumBiomass<-aov(Biomass~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumBiomass.sum<-summary(SolanumBiomass)

ggplot(Solanum.DATA, aes(x=Drought,y=Biomass, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge")+theme_bw()+ 
  ylab("Biomass(g)")+theme(legend.position="none",axis.text=element_text(size=12),
                           axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
                           strip.text = element_text(size = 12))


BiomassSolanum<-read.csv(file = "BiomassSolanum.csv")
ggplot(data = BiomassSolanum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Biomass")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Shoot wet weight (g)
Solanum.DATA<-read.csv("Solanum.csv",header = TRUE)
SolanumWshoot<-aov(Wshoot~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumWshoot.sum<-summary(SolanumWshoot)

ggplot(Solanum.DATA, aes(x=Drought,y=Wshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))


WshootSolanum<-read.csv(file = "WshootSolanum.csv")
ggplot(data = WshootSolanum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Wshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Shoot dry weight (g)

Solanum.DATA<-read.csv("Solanum.csv",header = TRUE)

SolanumWdshoot<-aov(Wdshoot~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumWdshoot.sum<-summary(SolanumWdshoot)


ggplot(Solanum.DATA, aes(x=Drought,y=Wdshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wdshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# Root wet weight (g)
Solanum.DATA<-read.csv("Solanum.csv",header = TRUE)

SolanumWroot<-aov(Wroot~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumWroot.sum<-summary(SolanumWroot)

ggplot(Solanum.DATA, aes(x=Drought,y=Wroot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

WrootSolanum<-read.csv(file = "WrootSolanum.csv")
ggplot(data = WrootSolanum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Wroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

ggplot(data = WrootSolanum,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("Wroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Root dry weight (g)

Solanum.DATA<-read.csv("Solanum.csv",header = TRUE)

SolanumWdroot<-aov(Wdroot~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumWdroot.sum<-summary(SolanumWdroot)

ggplot(Solanum.DATA, aes(x=Drought,y=Wdroot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wdroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# Total wet weight
Solanum.DATA$TotalWW<-Solanum.DATA$Wshoot+Solanum.DATA$Wroot

SolanumTotalwetweight<-aov(TotalWW~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumTotalwetweight.sum<-summary(SolanumTotalwetweight)


ggplot(Solanum.DATA, aes(x=Drought,y=TotalWW, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Total wet weight")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# chlorophyll concentration (%)

Solanum.DATA<-read.csv("Solanum.csv",header = TRUE)

SolanumS<-aov(S~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumS.sum<-summary(SolanumS)

ggplot(Solanum.DATA, aes(x=Drought,y=S, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("chlorophyll concentration")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

SSolanum<-read.csv(file = "SSolanum.csv")
ggplot(data =SSolanum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("chlorophyll concentration")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Shoot Lenght (cm)

Solanum.DATA<-read.csv("Solanum.csv",header = TRUE)

SolanumLshoot<-aov(Lshoot~Fungi*Drought*Salinity,data = Solanum.DATA)
SolanumLshoot.sum<-summary(SolanumLshoot)


ggplot(Solanum.DATA, aes(x=Drought,y=Lshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Lshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

LshootSolanum<-read.csv(file = "LshootSolanum.csv")
ggplot(data =LshootSolanum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Lshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# POX

Tomato.DATA<-read.csv("Tomato.csv",header = TRUE)

TomatoPOX<-aov(POX~fungi*Drought*Salinity,data = Tomato.DATA)
TomatoPOX.sum<-summary(TomatoPOX)


POXTomato<-read.csv(file = "POXTomato.csv")
ggplot(data =POXTomato,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("POX")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

POXTomato<-read.csv(file = "POXTomato.csv")
ggplot(data =POXTomato,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("POX")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# CAT

Tomato.DATA<-read.csv("Tomato.csv",header = TRUE)

TomatoCAT<-aov(CAT~fungi*Drought*Salinity,data = Tomato.DATA)
TomatoCAT.sum<-summary(TomatoCAT)


CATTomato<-read.csv(file = "CATTomato.csv")
ggplot(data =CATTomato,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("CAT")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

CATTomato<-read.csv(file = "CATTomato.csv")
ggplot(data =CATTomato,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("CAT")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# SOD

Tomato.DATA<-read.csv("Tomato.csv",header = TRUE)

TomatoSOD<-aov(SOD~fungi*Drought*Salinity,data = Tomato.DATA)
TomatoSOD.sum<-summary(TomatoSOD)


SODTomato<-read.csv(file = "SODTomato.csv")
ggplot(data =SODTomato,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("POX")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

SODTomato<-read.csv(file = "SODTomato.csv")
ggplot(data =SODTomato,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("POX")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))
# Protein

Tomato.DATA<-read.csv("Tomato.csv",header = TRUE)

TomatoProtein<-aov(Protein~fungi*Drought*Salinity,data = Tomato.DATA)
TomatoProtein.sum<-summary(TomatoProtein)

ProteinTomato<-read.csv(file = "ProteinTomato.csv")
ggplot(data =ProteinTomato,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Protein")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


ProteinTomato<-read.csv(file = "ProteinTomato.csv")
ggplot(data =ProteinTomato,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("Protein")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Proline

Tomato.DATA<-read.csv("Tomato.csv",header = TRUE)
TomatoProline<-aov(Proline~fungi*Drought*Salinity,data = Tomato.DATA)
TomatoProline.sum<-summary(TomatoProline)


ProlineTomato<-read.csv(file = "ProlineTomato.csv")
ggplot(data =ProlineTomato,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Proline")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

#DATA INPUT (Hordeum vulgare)
Barly.DATA<-read.csv("Barly.csv",header = TRUE)
Hordeum.DATA<-read.csv("Hordeum.csv",header = TRUE)

Hordeum.DATA$Biomass<-Hordeum.DATA$Wdshoot+Hordeum.DATA$Wdroot

# Total biomass (g)
anova.BiomassCUCUMBER<-aov(Biomass~Treatment.code,data=CUCUMIS.DATA)
output.BiomassCUCUMBER <- duncan.test(anova.BiomassCUCUMBER,"Treatment.code", group=TRUE)
plot(output.BiomassCUCUMBER,horiz=TRUE,las=1)
print(output.BiomassCUCUMBER$groups)


HordeumBiomass<-aov(Biomass~Fungi*Drought*Salinity,data = Solanum.DATA)
HordeumBiomass.sum<-summary(HordeumBiomass)

ggplot(Hordeum.DATA, aes(x=Drought,y=Biomass, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge")+theme_bw()+ 
  ylab("Biomass(g)")+theme(legend.position="none",axis.text=element_text(size=12),
                           axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
                           strip.text = element_text(size = 12))


BiomassHordeum<-read.csv(file = "BiomassHordeum.csv")
ggplot(data = BiomassHordeum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Biomass")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Shoot wet weight (g)

Hordeum.DATA<-read.csv("Hordeum.csv",header = TRUE)

HordeumWshoot<-aov(Wshoot~Fungi*Drought*Salinity,data = Hordeum.DATA)
HordeumWshoot.sum<-summary(HordeumWshoot)

ggplot(Hordeum.DATA, aes(x=Drought,y=Wshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

WshootHordeum<-read.csv(file = "WshootHordeum.csv")
ggplot(data = WshootHordeum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Wshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Shoot dry weight (g)

Hordeum.DATA<-read.csv("Hordeum.csv",header = TRUE)

HordeumWdshoot<-aov(Wdshoot~Fungi*Drought*Salinity,data = Hordeum.DATA)
HordeumWdshoot.sum<-summary(HordeumWdshoot)


ggplot(Hordeum.DATA, aes(x=Drought,y=Wdshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wdshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# Root wet weight (g)
Hordeum.DATA<-read.csv("Hordeum.csv",header = TRUE)

HordeumWroot<-aov(Wroot~Fungi*Drought*Salinity,data = Hordeum.DATA)
HordeumWroot.sum<-summary(HordeumWroot)

ggplot(Hordeum.DATA, aes(x=Drought,y=Wroot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

WrootHordeum<-read.csv(file = "WrootHordeum.csv")
ggplot(data = WrootHordeum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Wroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Root dry weight (g)

Hordeum.DATA<-read.csv("Hordeum.csv",header = TRUE)

HordeumWdroot<-aov(Wdroot~Fungi*Drought*Salinity,data = Hordeum.DATA)
HordeumWdroot.sum<-summary(HordeumWdroot)

ggplot(Hordeum.DATA, aes(x=Drought,y=Wdroot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Wdroot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# Total wet weight

Hordeum.DATA<-read.csv("Hordeum.csv",header = TRUE)
Hordeum.DATA$TotalWW<-Hordeum.DATA$Wshoot+Hordeum.DATA$Wroot

HordeumTotalwetweight<-aov(TotalWW~Fungi*Drought*Salinity,data = Solanum.DATA)
HordeumTotalwetweight.sum<-summary(HordeumTotalwetweight)


ggplot(Hordeum.DATA, aes(x=Drought,y=TotalWW, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Total wet weight")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

# chlorophyll concentration (%)

HordeumS<-aov(S~Fungi*Drought*Salinity,data = Hordeum.DATA)
HordeumS.sum<-summary(HordeumS)

ggplot(Hordeum.DATA, aes(x=Drought,y=S, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("chlorophyll concentration")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

SHordeum<-read.csv(file = "SHordeum.csv")
ggplot(data =SHordeum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("chlorophyll concentration")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Shoot Lenght (cm)

Hordeum.DATA<-read.csv("Hordeum.csv",header = TRUE)

HordeumLshoot<-aov(Lshoot~Fungi*Drought*Salinity,data = Hordeum.DATA)
HordeumLshoot.sum<-summary(HordeumLshoot)


ggplot(Hordeum.DATA, aes(x=Drought,y=Lshoot, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Lshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

LshootHordeum<-read.csv(file = "LshootHordeum.csv")
ggplot(data =LshootHordeum,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Lshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# POX

Barly.DATA<-read.csv("Barly.csv",header = TRUE)

BarlyPOX<-aov(POX~Fungi*Drought*Salinity,data = Barly.DATA)
BarlyPOX.sum<-summary(BarlyPOX)

ggplot(Barly.DATA, aes(x=Drought,y=POX, fill=Fungi)) + 
  facet_wrap(~Salinity)+geom_boxplot(position = "dodge", width=0.5)+theme_bw()+ylab("Lshoot")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),legend.title = element_text(size=12),
        strip.text = element_text(size = 12))

POXBarly<-read.csv(file = "POXBarly.csv")
ggplot(data =POXBarly,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("POX")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))



# CAT

Barly.DATA<-read.csv("Barly.csv",header = TRUE)

BarlyCAT<-aov(CAT~Fungi*Drought*Salinity,data = Barly.DATA)
BarlyCAT.sum<-summary(BarlyCAT)


CATBarly<-read.csv(file = "CATBarly.csv")
ggplot(data =CATBarly,aes(x=Salinity,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Drought)+ theme_bw()+ ylab("CAT")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# SOD

Barly.DATA<-read.csv("Barly.csv",header = TRUE)

BarlySOD<-aov(SOD~Fungi*Drought*Salinity,data = Barly.DATA)
BarlySOD.sum<-summary(BarlySOD)


SODBarly<-read.csv(file = "SODBarly.csv")
ggplot(data =SODBarly,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("SOD")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))


# Protein

Barly.DATA<-read.csv("Barly.csv",header = TRUE)

BarlyProtein<-aov(Protein~Fungi*Drought*Salinity,data = Barly.DATA)
BarlyProtein.sum<-summary(BarlyProtein)


ProteinBarly<-read.csv(file = "ProteinBarly.csv")
ggplot(data =ProteinBarly,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Protein")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))

# Proline

Barly.DATA<-read.csv("Barly.csv",header = TRUE)
BarlyProline<-aov(Proline~Fungi*Drought*Salinity,data = Barly.DATA)
BarlyProline.sum<-summary(BarlyProline)


ProlineBarly<-read.csv(file = "ProlineBarly.csv")
ggplot(data =ProlineBarly,aes(x=Drought,y=Value, fill=Fungi ))+
  geom_bar(stat="identity",position="dodge",width=0.5)+
  facet_wrap(~Salinity)+ theme_bw()+ ylab("Proline")+
  theme(legend.position="top",axis.text=element_text(size=12),
        axis.title=element_text(size=14),legend.text=element_text(size=12),
        strip.text = element_text(size = 12))+ guides(fill=guide_legend(""))




