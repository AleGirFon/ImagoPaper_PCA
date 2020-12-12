### FINAL PLOTS IMAGO_PAPER
# This script is only valid once PCA_Analysis.R has been run
#install.packages("factoextra")
library(RColorBrewer)
library(ggplot2)
library(ggfortify)
library(dplyr)
library(hrbrthemes)
library(psych)
library(cowplot)
library(ggpubr)
library(factoextra)
library(gridExtra)
library(grid)

# Sup Figure 1.D --> CDB1 PCA8
# PCA 8 CDB1 --> add eigenvector plot


setwd("E:/Programming/R/Celia")
DB<-read.csv("E:/Programming/R/Celia/FinalCelia/segmented30042020_AGFormat.csv",header=TRUE,stringsAsFactors = FALSE)
head(DB)

summary(DB)
table(DB$Type)
table(DB$Stage)
table(DB$Genotype)
table(DB$Labels)

# Initial dim: 272040 entries

# Removing Empty entries from cell LABELS
DB<-DB[which(DB$Labels!=""),]
summary(DB)


# Removing L1/L2 Basal-sup

DB<-DB[-which(DB$Labels=="L2 basal sup"),]
DB<-droplevels(DB[-which(DB$Labels=="L1 basal sup"),])
table(DB$Labels)


#Pooling "pSMC" with "SMC"
DB$Labels[which(DB$Labels=="pSMC")]<-"SMC"
DB<-droplevels(DB)
table(DB$Labels)

#Adapting the color Palette to final samples (without L1/L2 Basal-sup)
ColorP<-read.csv("PaletteColor.csv", header=FALSE)
ColorP
ColorP<-droplevels(ColorP[-c(4,8,10),])
ColVec<-as.character(ColorP$V2)




# AfterFiltering dim: 77720 entries

#Describe data
dd<-ggplot(DB, aes(Value, color=Type)) 
dd+ geom_histogram(position ="dodge")
#moveme(names(DB), "Value last") #Moving the column
summary(DB)


#SUBSET DB by variables
VOL<-DB[which(as.character(DB$Type)=="Cell Volume"),]
dd<-ggplot(VOL, aes(Value, color=Labels, fill=Labels)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=ColVec) + 
  scale_color_manual(values=ColVec)+
  theme_classic() +ggtitle("1. Cell Volume Histogram All Stages")
pl
ggsave("CDB2/1_HIST_Volume_CDB2.png",pl,dpi=600)


AREA<-DB[which(as.character(DB$Type)=="Cell Area"),]
dd<-ggplot(AREA, aes(Value, color=Labels, fill=Labels)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=ColVec) + 
  scale_color_manual(values=ColVec)+
  theme_classic() +ggtitle("2. Cell Area Histogram All Stages")
pl
ggsave("CDB2/2_HIST_Area_CDB2.png",pl,dpi=600)



OBL<-DB[which(as.character(DB$Type)=="Cell Ellipticity (oblate)"),]
dd<-ggplot(OBL, aes(Value, color=Labels, fill=Labels)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=ColVec) + 
  scale_color_manual(values=ColVec)+
  theme_classic() +ggtitle("3. Cell Oblacity Histogram All Stages")
pl
ggsave("CDB2/3_HIST_Obl_CDB2.png",pl,dpi=600)


PRL<-DB[which(as.character(DB$Type)=="Cell Ellipticity (prolate)"),]
dd<-ggplot(PRL, aes(Value, color=Labels, fill=Labels)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=ColVec) + 
  scale_color_manual(values=ColVec)+
  theme_classic() +ggtitle("4. Cell Prolacity Histogram All Stages")
pl
ggsave("CDB2/4_HIST_Prl_CDB2.png",pl,dpi=600)


SPH<-DB[which(as.character(DB$Type)=="Cell Sphericity"),]
dd<-ggplot(SPH, aes(Value, color=Labels, fill=Labels)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=ColVec) + 
  scale_color_manual(values=ColVec)+
  theme_classic() +ggtitle("5. Cell Sphericity Histogram All Stages")
pl
ggsave("CDB2/5_HIST_Sphe_CDB2.png",pl,dpi=600)


## Now the value vectors:
Vol<-as.numeric(VOL$Value)
Area<-as.numeric(AREA$Value)
Oblate<-as.numeric(OBL$Value)
Prolate<-as.numeric(PRL$Value)
Sphericity<-as.numeric(SPH$Value)


#DB now will be reformated
DB<-data.frame(AREA[,1:3],AREA[,6:7],Oblate,Prolate,Vol,Area,Sphericity)

describeBy(DB, DB$Stage)
describeBy(DB, DB$Labels)


## Separation Col from WS&bot

table(DB$Genotype)


ColDB<-DB[which(DB$Genotype=="Col"),]
BotDB<-DB[which(DB$Genotype=="bot"),]
WsDB<-DB[which(DB$Genotype=="Ws4"),]

BowsDB<-rbind(BotDB,WsDB)
head(BowsDB)

# In case we don't remove any variable:
#1. Have a view to the data table:
head(ColDB)
head(BowsDB)
#2. See if variable correlate
pairs(ColDB[,6:10])
pairs(BowsDB[,6:10])

# We can see how Volume and Area are highly correlated. Oblate and Prolate are also correlated ("strange manner")
# Rest of variables are apparently independent to each other



AA<-DB
rm(DB)
# 
# CDB1:
#   
#   Sup Figure 1.D --> CDB1 PCA8
# PCA 8 CDB1 --> add eigenvector plot
# 



#2. PCA analysis
#2.1 Calculating Principal Components
AllVAr.PCA<-prcomp(ColDB[,6:10],scale=TRUE)
#2.2 Show Principal components stats summary
AllVAr.PCA
summary(AllVAr.PCA)

# #Weight table as a plot
# grid.arrange(
#   tableGrob(AllVAr.PCA$rotation)
# )


#2.3 Barplot of PCAs and percent of explained variance
#fviz_eig(AllVAr.PCA, barcol="black", barfill= brewer.pal(n = 5, name = "RdYlBu")) MULTI COLOR BARPLOT
fviz_eig(AllVAr.PCA, barcol="black", barfill= "grey")




#2.4 PCA plot and Vector contribution plot (variables as eigenvecs)

PCA.dot<-autoplot(AllVAr.PCA,data=ColDB,  colour="Labels", shape="Stage", fill="Labels", alpha=0.8)+
  scale_color_manual(values=as.character(ColorP$V2))+ 
  scale_fill_manual(values=as.character(ColorP$V2))+ 
  scale_shape_manual(values = c(21, 21, 21,22,23,24,25))+ theme_classic()

PCA.dot



PCA.vec<-fviz_pca_var(AllVAr.PCA,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec







### SEPARATED BY MAJOR STAGES (St0, St1, St2)  
NewColor<-c(as.character(ColVec),"#A0A0A0","#C0C0C0","#E0E0E0")
NC<-c()
for (i in 1: dim(ColDB)[1]){
  #NC[i]<-NA
  NC[i]<-as.character(ColDB$Stage[i])
  if ((ColDB$Stage[i]=="st0-I")|| (ColDB$Stage[i]=="st0-II")||(ColDB$Stage[i]=="st0-III")){
    NC[i]<-"St0"
  }
  if ((ColDB$Stage[i]=="st1-I")|| (ColDB$Stage[i]=="st1-II")){
    NC[i]<-"St1"
  }
  if ((ColDB$Stage[i]=="st2-I")|| (ColDB$Stage[i]=="st2-II")){
    NC[i]<-"St2"
  }
  
}

ColDB.stg<-data.frame(ColDB,NC) 
colnames(ColDB.stg)[11]<-"MajorStage"  
head(ColDB.stg)




### Stage 0   #######

ColDB.stg<-ColDB.stg[which(ColDB.stg$MajorStage=="St0"),]
Stg.pca<-prcomp(ColDB.stg[,6:10],scale=TRUE)


Stg.plot<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels",shape="Stage",frame=TRUE, main= "Stg0 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)
Stg.plot2<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels",shape="Stage",frame=TRUE,frame.colour="Stage", main= "Stg0 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)
Stg.plot3<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels", main= "Stg0 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)



pl<-Stg.plot+ scale_color_manual(values=as.character(NewColor))+
  #geom_point(aes(ColDB.stg$Labels %in% "CC"))+
  scale_fill_manual(values=as.character(NewColor))+
  scale_shape_manual(values = c(21, 22, 25,21,22,21,22))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg0_frameLabel.png",pl,dpi=600)




pl<-Stg.plot2+ scale_color_manual(values=as.character(NewColor))+
  #geom_point(aes(ColDB.stg$Labels %in% "CC"))+
  scale_fill_manual(values=as.character(NewColor))+
  scale_shape_manual(values = c(21, 22, 25,21,22,21,22))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg0_frameStage.png",pl,dpi=600)

pl<-Stg.plot3+ scale_color_manual(values=as.character(NewColor))+
  scale_fill_manual(values=as.character(NewColor))+
  facet_wrap(as.factor(ColDB.stg$Stage))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg0_sepPlotStage.png",pl,dpi=600)


PCA.vec<-fviz_pca_var(Stg.pca,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec

Stg.pca
summary(Stg.pca)



### Stage 1    #######

ColDB.stg<-data.frame(ColDB,NC) 
colnames(ColDB.stg)[11]<-"MajorStage"  
head(ColDB.stg)
ColDB.stg<-ColDB.stg[which(ColDB.stg$MajorStage=="St1"),]
Stg.pca<-prcomp(ColDB.stg[,6:10],scale=TRUE)

Stg.plot<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels",shape="Stage",frame=TRUE, main= "Stg1 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)
Stg.plot2<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels",shape="Stage",frame=TRUE,frame.colour="Stage", main= "Stg1 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)
Stg.plot3<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels", main= "Stg1 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)



pl<-Stg.plot+ scale_color_manual(values=as.character(NewColor))+
  #geom_point(aes(ColDB.stg$Labels %in% "CC"))+
  scale_fill_manual(values=as.character(NewColor))+
  scale_shape_manual(values = c(21, 22, 25,21,22,21,22))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg1_frameLabel.png",pl,dpi=600)




pl<-Stg.plot2+ scale_color_manual(values=as.character(NewColor))+
  #geom_point(aes(ColDB.stg$Labels %in% "CC"))+
  scale_fill_manual(values=as.character(NewColor))+
  scale_shape_manual(values = c(21, 22, 25,21,22,21,22))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg1_frameStage.png",pl,dpi=600)

pl<-Stg.plot3+ scale_color_manual(values=as.character(NewColor))+
  scale_fill_manual(values=as.character(NewColor))+
  facet_wrap(as.factor(ColDB.stg$Stage))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg1_sepPlotStage.png",pl,dpi=600)


PCA.vec<-fviz_pca_var(Stg.pca,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec

Stg.pca
summary(Stg.pca)




### Stage 2

ColDB.stg<-data.frame(ColDB,NC) 
colnames(ColDB.stg)[11]<-"MajorStage"  
head(ColDB.stg)
ColDB.stg<-ColDB.stg[which(ColDB.stg$MajorStage=="St2"),]
Stg.pca<-prcomp(ColDB.stg[,6:10],scale=TRUE)

Stg.plot<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels",shape="Stage",frame=TRUE, main= "Stg2 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)
Stg.plot2<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels",shape="Stage",frame=TRUE,frame.colour="Stage", main= "Stg2 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)
Stg.plot3<-autoplot(Stg.pca,data=ColDB.stg ,size=2.5, colour="Labels", fill="Labels", main= "Stg2 PCA -- Oblacity, Prolacity, Sphericity, Volume and Area", alpha=0.75)



pl<-Stg.plot+ scale_color_manual(values=as.character(NewColor))+
  #geom_point(aes(ColDB.stg$Labels %in% "CC"))+
  scale_fill_manual(values=as.character(NewColor))+
  scale_shape_manual(values = c(21, 22, 25,21,22,21,22))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg2_frameLabel.png",pl,dpi=600)




pl<-Stg.plot2+ scale_color_manual(values=as.character(NewColor))+
  #geom_point(aes(ColDB.stg$Labels %in% "CC"))+
  scale_fill_manual(values=as.character(NewColor))+
  scale_shape_manual(values = c(21, 22, 25,21,22,21,22))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg2_frameStage.png",pl,dpi=600)

pl<-Stg.plot3+ scale_color_manual(values=as.character(NewColor))+
  scale_fill_manual(values=as.character(NewColor))+
  facet_wrap(as.factor(ColDB.stg$Stage))+
  theme_classic()

pl

ggsave("FinalCelia/SupFig_1E_Stg2_sepPlotStage.png",pl,dpi=600)


PCA.vec<-fviz_pca_var(Stg.pca,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec




PCA.vec<-fviz_pca_var(Stg.pca,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec
Stg.pca
summary(Stg.pca)





# REMOVING L1-Sup and L2 basal --> ALL STAGES ONE PLOT
# Format 1 -- SVA
NewColor<-c(as.character(ColVec),"#A0A0A0","#C0C0C0","#E0E0E0")

NewColor<-c(as.character(ColorP$V2))
NewColor<-NewColor[-c(3,6)]


St<-ColDB[-which(as.character(ColDB$Labels)=="L1 basal"),]
St<-St[-which(as.character(St$Labels)=="L2 basal"),]
St.plot<-autoplot(prcomp(St[,8:10],scale=TRUE),data=St ,size=2, colour="Labels", fill="Labels",main= "St plots Volume - Area - Sphericity",alpha=0.7)

pl<-St.plot+
  scale_fill_manual(values=NewColor) +
  scale_color_manual(values=NewColor) +
  facet_wrap(as.factor(St$Stage)) + theme_classic()


pl



ggsave("FinalCelia/SupFig_2I_Reduced-Col0_sepPlotStage.png",pl,dpi=600)


PCA.vec<-fviz_pca_var(prcomp(St[,8:10],scale=TRUE),
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 


PCA.vec
NewColor<-c(as.character(ColVec),"#A0A0A0","#C0C0C0","#E0E0E0")








####################################################################################################################################################################


####################################################################################################################################################################
#########################                   WS vs bot                ##############################
####################################################################################################################################################################


####################################################################################################################################################################

#Stage 0-III

# OldPlots: ColVec<-c(as.character(ColorP$V2[1:7]),"#E0E0E0",as.character(ColorP$V2[8:10]))

table(BowsDB$Stage)
St0<-BowsDB[which(as.character(BowsDB$Stage)=="st0-III"),]

# # REMOVING OBLATE AND PROLATE
# St0<-subset(St0,select=-c(Oblate,Prolate))
St0.pca<-prcomp(scale(St0[,6:10]))
St0.pca
summary(St0.pca)
St0.plot<-autoplot(St0.pca,data=St0 ,size=2.8, colour="Labels", alpha=0.8, frame=TRUE,frame.colour="Stage",frame.type="norm")

pl<-St0.plot+ scale_color_manual(values=NewColor)+
  scale_fill_manual(values=rep("white",10))+ 
  facet_wrap( ~ Genotype)+ theme_classic() +
  ggtitle("9. PCA Obl-Prl-Sph-Vol-Area st0-III")
pl
ggsave("FinalCelia/SupFig_4B1_PCA_st0-III_BotVsWs.png",pl,dpi=600)


PCA.vec<-fviz_pca_var(St0.pca,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec
Bicolor<-c("lightsalmon","darkred")
St0.plot<-autoplot(St0.pca,data=St0 ,size=2.8, colour="Genotype", alpha=0.8, frame=TRUE,frame.colour="Genotype",frame.type="norm")
pl<-St0.plot+ scale_color_manual(values=Bicolor)+
  scale_fill_manual(values=Bicolor)+ 
  ggtitle("9. PCA Obl-Prl-Sph-Vol-Area st0-III") + theme_classic()
pl
ggsave("FinalCelia/SupFig_4B1_PCA_st0-III_BotVsWs_ALT.png",pl,dpi=600)


#Stage 1-I

table(BowsDB$Stage)
St1<-BowsDB[which(as.character(BowsDB$Stage)=="st1-I"),]
# 
# # REMOVING OBLATE AND PROLATE
# St0<-subset(St0,select=-c(Oblate,Prolate))
# St0<-subset(St0,select=-c(Oblate,Prolate))
St1.pca<-prcomp(scale(St1[,6:10]))
St1.pca
summary(St1.pca)
St1.plot<-autoplot(St1.pca,data=St1 ,size=2.8, colour="Labels",  alpha=0.8, frame=TRUE,frame.colour="Stage",frame.type="norm")

pl<-St1.plot+ scale_color_manual(values=NewColor)+
  scale_fill_manual(values=rep("white",10))+ 
  facet_wrap( ~ Genotype)+ theme_classic() +
  ggtitle("10. PCA Obl-Prl-Sph-Vol-Area st1-I")
pl
ggsave("FinalCelia/SupFig_4B2_PCA_st1-I_BotVsWs.png",pl,dpi=600)

PCA.vec<-fviz_pca_var(St1.pca,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec

St1.plot<-autoplot(St1.pca,data=St1 ,size=2.8, colour="Genotype", alpha=0.8, frame=TRUE,frame.colour="Genotype",frame.type="norm")
pl<-St1.plot+ scale_color_manual(values=Bicolor)+
  scale_fill_manual(values=Bicolor)+ 
  ggtitle("10. PCA Obl-Prl-Sph-Vol-Area st0-III") + theme_classic()
pl
ggsave("FinalCelia/SupFig_4B2_PCA_st1-I_BotVsWs_ALT.png",pl,dpi=600)




#Stage 1-II


St1II<-BowsDB[which(as.character(BowsDB$Stage)=="st1-II"),]

#St1II<-subset(St1II,select=-c(Oblate,Prolate))

St1II.pca<-prcomp(scale(St1II[,6:10]))
St1II.pca
summary(St1II.pca)
St1II.plot<-autoplot(St1II.pca,data=St1II ,size=2.8, colour="Labels",  alpha=0.8, frame=TRUE,frame.colour="Stage",frame.type="norm")

pl<-St1II.plot+ scale_color_manual(values=NewColor)+
  scale_fill_manual(values=rep("white",10))+ 
  facet_wrap( ~ Genotype)+ theme_classic() +
  ggtitle("11.  PCA Obl-Prl-Sph-Vol-Area st1-II")
pl

ggsave("FinalCelia/SupFig_4B3_PCA_st1-II_BotVsWs.png",pl,dpi=600)

PCA.vec<-fviz_pca_var(St1II.pca,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 

PCA.vec

St1II.plot<-autoplot(St1II.pca,data=St1II ,size=2.8, colour="Genotype", alpha=0.8, frame=TRUE,frame.colour="Genotype",frame.type="norm")
pl<-St1II.plot+ scale_color_manual(values=Bicolor)+
  scale_fill_manual(values=Bicolor)+ 
  ggtitle("10. PCA Obl-Prl-Sph-Vol-Area st0-III") + theme_classic()
pl
ggsave("FinalCelia/SupFig_4B3_PCA_st1-II_BotVsWs_ALT.png",pl,dpi=600)






##Mann-Whitney-Wilcoxon Test
#1. Volume SMC vs volume ALL in COL
#1.a) Making New Column on DB with SMC / no SMC var
NC<-c()

for (i in 1: dim(ColDB)[1]){
  
  NC[i]<-"Not_SMC"
  if ((ColDB$Labels[i]=="SMC")||(ColDB$Labels[i]=="pSMC")) {
    NC[i]<-"SMC"
  }
  
}


DB.SMC<-data.frame(ColDB,NC)

# 1.b) test all stages pooled
SMC.test<-wilcox.test(Vol ~ NC, data=DB.SMC) 
SMC.test
summary(SMC.test)
describeBy(DB.SMC$Vol, NC)

#1.c) Performing the test across stages
Stg<-names(table(ColDB$Stage))
for (stage in Stg){
  print(stage)
  DB.SMC2<-DB.SMC[which(as.character(DB.SMC$Stage)==stage),]
  
  SMC.test<-wilcox.test(Vol ~ NC, data=DB.SMC2) 
  print(SMC.test)
  print(describeBy(DB.SMC$Vol, NC))
  cat("\n\n")
  
}




### 2. L1.vol vs L2.vol by stage (importantly 0-1 and 0-II) for col


DBL<-data.frame()
NC<-c()

#2.a) Making New Column L1/L2 clasification
for (i in 1: dim(ColDB)[1]){
  
  NC[i]<-"Not_L"
  if ((ColDB$Labels[i]=="L1 apical")||(ColDB$Labels[i]=="L1 basal")||(ColDB$Labels[i]=="L1 basal sup")||(ColDB$Labels[i]=="L1 dome")) {
    NC[i]<-"L1"
  }
  if ((ColDB$Labels[i]=="L2 apical")||(ColDB$Labels[i]=="L2 basal")||(ColDB$Labels[i]=="L2 basal sup")) {
    NC[i]<-"L2"
  }
  
}

DBL<-data.frame(ColDB,NC)
DBL<-droplevels(DBL[which(as.character(DBL$NC)!="Not_L"),])


# 2.b) test all stages pooled
DBL.test<-wilcox.test(Vol ~ NC, data=DBL) 
DBL.test
summary(DBL.test)
describeBy(DBL$Vol, DBL$NC)





#2.c) Performing the test across stages
Stg<-names(table(DBL$Stage))

for (stage in Stg){
  print(stage)
  DBL2<-DBL[which(as.character(DBL$Stage)==stage),]
  
  if (length(names(table(DBL2$NC)))==2){
    DBL.test<-wilcox.test(Vol ~ NC, data=DBL2) 
    print(DBL.test)
    print(describeBy(DBL2$Vol, DBL2$NC))
    
  }else{
    print ("NOT ENOUGH DATA TO PERFORM THE TEST!")
  }
  cat("\n\n")
  
  
}






#3.   Cell Number (can you get it easily or do you have to compute it?) of  
#apical=(L1apical, L1dome, L2apical, SMC, CC) of bot against Ws for each stages 0-III, 1-I, 1-II  = 3 tests



#3.b) Labeling BowsDB with "apical" or "non_apical" tags
NC<-c()
for (i in 1: dim(BowsDB)[1]){
  
  NC[i]<-"Non_Apical"
  if ((BowsDB$Labels[i]=="L1 apical")||(BowsDB$Labels[i]=="L1 dome")||(BowsDB$Labels[i]=="L2 apical")||(BowsDB$Labels[i]=="SMC")||(BowsDB$Labels[i]=="CC")) {
    NC[i]<-"Apical"
  }
  
  
}
BowsDB<-rbind(BotDB,WsDB)
BowsDB<-data.frame(BowsDB,NC)
head(BowsDB)
colnames(BowsDB)<-c("Stack","Genotype","Stage","Labels" , "neighbourSMC", "Oblate", "Prolate" ,"Vol" ,"Area" ,"Sphericity" ,"POS")

#3.c)Getting Cell numbers by Genotype
WS4<- data.frame(BowsDB[which(BowsDB$Genotype=="Ws4"),] %>% group_by(Stack) %>% count(POS))
WS4<-data.frame(WS4,rep("Ws4",dim(WS4)[1]))
colnames(WS4)<-c("Stack","Position","CellNum","Genotype")

BOT<- data.frame(BowsDB[which(BowsDB$Genotype=="bot"),] %>% group_by(Stack) %>% count(POS))
BOT<- data.frame(BOT,rep("bot",dim(BOT)[1]))

colnames(BOT)<-c("Stack","Position","CellNum","Genotype")

head(WS4)
head(BOT)

CNDB<-rbind(WS4,BOT)

NewColor<-c("red3","peachpuff")

Comp<-list(c("Ws4","bot"),c("Ws4","bot"))
dd<-ggplot(CNDB, aes(Genotype,CellNum, fill=Genotype )) 
p1<-dd + geom_boxplot() +
  scale_fill_manual(values=NewColor) + 
  facet_wrap(~Position)+
  theme_classic() + geom_jitter(width = 0.2) + stat_compare_means()

p1

ggsave("FinalCelia/15_Boxplot_Wtest_CellNum_BotVsWs4_CDB2.png",p1,dpi=600)

dd<-ggplot(CNDB, aes(CellNum, fill=Genotype )) 
p2<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=NewColor) + 
  scale_color_manual(values=NewColor) +
  facet_wrap(~Position)+
  theme_classic()
p2


ggsave("FinalCelia/16_Hist_Wtest_CellNum_BotVsWs4_CDB2.png",p2,dpi=600)


#3.d) Do the test by stages:

Stg<-names(table(BowsDB$Stage))


for (stage in Stg){
  BowsDB.temp<-BowsDB[which(BowsDB$Stage==stage),]
  BowsDB.temp<-BowsDB.temp[which(BowsDB.temp$POS =="Apical"),]
  
  print(head(BowsDB.temp))
  WS4<- data.frame(BowsDB.temp[which(BowsDB.temp$Genotype=="Ws4"),] %>% group_by(Stack) %>% count(POS))
  WS4<-data.frame(WS4,rep("Ws4",dim(WS4)[1]))
  colnames(WS4)<-c("Stack","Position","CellNum","Genotype")
  
  BOT<- data.frame(BowsDB.temp[which(BowsDB.temp$Genotype=="bot"),] %>% group_by(Stack) %>% count(POS))
  BOT<- data.frame(BOT,rep("bot",dim(BOT)[1]))
  
  colnames(BOT)<-c("Stack","Position","CellNum","Genotype")
  print(stage)
  print(head(WS4))
  print(head(BOT))
  
  BowsDB.temp<-rbind(WS4,BOT)
  cat ("\n W. test for Cell Number in", stage, " --  Ws4 vs bot \n", sep=" ")
  
  BowsDB.test<-wilcox.test(CellNum ~ Genotype, data=BowsDB.temp, correct=FALSE) 
  print(BowsDB.test)
  print( describeBy(BowsDB.temp$CellNum , BowsDB.temp$Genotype))
  cat("\n\n")
  
}





#4. SMC Volume of bot against Ws for each stages 0-III, 1-I, 1-II= 3 tests

#4.a) Making New Column on DB with SMC / no SMC var
NC<-c()

DB.SMC<-BowsDB
for (i in 1: dim(DB.SMC)[1]){
  
  NC[i]<-"Not_SMC"
  if ((DB.SMC$Labels[i]=="SMC")||(DB.SMC$Labels[i]=="pSMC")) {
    NC[i]<-"SMC"
  }
  
}

TESTs<-c("st0-III", "st1-I","st1-II")

DB.SMC<-data.frame(DB.SMC,NC)
DB.SMC2<- DB.SMC[which(as.character(DB.SMC$NC)=="SMC"),]

SMC.test<-wilcox.test(Vol ~ Genotype, data=DB.SMC2) 
SMC.test
summary(SMC.test)
describeBy(DB.SMC2$Vol, DB.SMC2$Genotype)

#4.b) Performing the test

#very few samples:
table(DB.SMC2$Stage)

Stg<-names(table(DB.SMC2$Stage))
Stg



for (stage in Stg){
  print(stage)
  DB.temp<-DB.SMC2[which(as.character(DB.SMC2$Stage)==stage),]
  print(describeBy(DB.temp$Vol,DB.temp$Genotype))
  if (length(names(table(DB.temp$Genotype)))==2){
    SMC.test<-wilcox.test(Vol ~ Genotype, data=DB.temp) 
    print(SMC.test)
  }else{
    print ("NOT ENOUGH DATA TO PERFORM THE W. TEST! -- missing one group")
  }
  cat("\n\n")
}





#5. NON_SMC Volume of bot against Ws for each stages 0-III, 1-I, 1-II= 3 tests

#5.a) Making New Column on DB with SMC / no SMC var

NC<-c()

DB.SMC<-BowsDB
for (i in 1: dim(DB.SMC)[1]){
  
  NC[i]<-"Not_SMC"
  if ((DB.SMC$Labels[i]=="SMC")||(DB.SMC$Labels[i]=="pSMC")) {
    NC[i]<-"SMC"
  }
  
}


DB.SMC<- DB.SMC[which(as.character(DB.SMC$NC)=="Not_SMC"),]

