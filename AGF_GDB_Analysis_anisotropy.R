#install.packages("ggfortify")
library(ggplot2)
library(ggfortify)
library(factoextra)
setwd("E:/Programming/R/Celia/GDB_labeled/")
dir()

### READING THE DATA PER FOLDER
setwd ("st0-1/")
files<-dir()
st01<-do.call(rbind,lapply(files,read.csv))
st01<-cbind(rep("st0-1",nrow(st01)),st01)
colnames(st01)[1]<-"Stage"



setwd ("../st0-II/")
files<-dir()
st02<-do.call(rbind,lapply(files,read.csv))
st02<-cbind(rep("st0-II",nrow(st02)),st02)
colnames(st02)[1]<-"Stage"


setwd ("../st0-III/")
files<-dir()
st03<-do.call(rbind,lapply(files,read.csv))
st03<-cbind(rep("st0-III",nrow(st03)),st03)
colnames(st03)[1]<-"Stage"


setwd ("../st1-I/")
files<-dir()
st11<-do.call(rbind,lapply(files,read.csv))
st11<-cbind(rep("st1-I",nrow(st11)),st11)
colnames(st11)[1]<-"Stage"

setwd ("../st1-II/")
files<-dir()
st12<-do.call(rbind,lapply(files,read.csv))
st12<-cbind(rep("st1-II",nrow(st12)),st12)
colnames(st12)[1]<-"Stage"


setwd ("../st2-I/")
files<-dir()
st21<-do.call(rbind,lapply(files,read.csv))
st21<-cbind(rep("st2-I",nrow(st21)),st21)
colnames(st21)[1]<-"Stage"


setwd("E:/Programming/R/Celia/")
ColorP<-read.csv("PaletteColor.csv", header=FALSE)
ColorP
ColorP<-droplevels(ColorP[-c(4,6,7,8,10),])
ColVec<-as.character(ColorP$V2)

#MERGE ALL IN GDB + CREATING CELL_STAGE COLUMN
GDB<-rbind(st01,st02,st03,st11,st12,st21)
GDB$Cell_Stage<-paste(GDB$Cell.Type,"-",GDB$Stage,sep="")

#Columns 5 and 6 have too complex names
colnames(GDB)[5:6]<-c("Top_Surface_Area","Bottom_Surface_Area")


## Checking the data
summary(GDB)


### ADDING VALUES EXTRACTED FROM RATIOS
MID<-GDB$Max.anisotropy.length/GDB$Max.Mid.anisotropy
MIN<-GDB$Max.anisotropy.length/GDB$Max.Min.anisotropy

Sumall<-GDB$Max.anisotropy.length + MID + MIN

Data.backup<-GDB$Max.anisotropy.length
MAX<-GDB$Max.anisotropy.length / Sumall
MID<-MID/ Sumall
MIN <- MIN / Sumall


GDB<-data.frame(GDB$Cell.Index,GDB$Stage,GDB$Cell.Type,GDB$Cell.Label,GDB$Cell_Stage,MIN,MID,MAX,GDB$Volume,GDB$Top_Surface_Area,GDB$Bottom_Surface_Area,GDB$Top.Bottom.ratio)
colnames(GDB)<-c("Cell.Index","Stage","Cell.Type","Cell.Label","Cell_Stage","Min","Mid","Max","Volume","Top_Surface_Area","Bottom_Surface_Area","TB.Ratio")


GDB<-GDB[-583,]# This data point doesn't represent a real cell --> Thus, we deleted ith

### REMOVE SPACE IN THE NAMES --> Very important!
GDB$Cell.Type<-trimws(GDB$Cell.Type,which = c("left"))

### REMOVING GENERIC CELLS
table(GDB[which(as.character(GDB$Cell.Type) == "GenericCell"),3])
GDB.filt<-data.frame(GDB[which(as.character(GDB$Cell.Type) != "GenericCell"),])

GDB.filt<-droplevels(GDB.filt)
table(GDB.filt$Cell.Type)

#Removing Cells with Volume < 5
GDG.filt<-GDB.filt[which(GDB.filt$Volume >= 5),]






#### START OF THE ANALYSIS

setwd("E:/Programming/R/Celia/")
ColorP<-read.csv("PaletteColor.csv", header=FALSE)
ColorP
ColorP<-droplevels(ColorP[c(1,3,6,5,9),])
ColVec<-as.character(ColorP$V2)


GDB.pca<-prcomp(GDB.filt[,6:9],scale=TRUE)
GDB.pca
summary(GDB.pca)

GDB.plot<-autoplot(prcomp(GDB.filt[,6:9],scale=TRUE),data=GDB.filt, colour="Cell.Type", shape="Stage",fill="Cell.Type")
GDB.plot+
  scale_color_manual(values=ColVec) +
  scale_fill_manual(values=ColVec)  +
  scale_shape_manual(values = c(21, 21, 21,22,23,24,25))


#Plot Histograms

dd<-ggplot(GDB.filt, aes(Volume, color=Cell.Type, fill=Cell.Type)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=as.character(ColorP$V2)) +
  scale_color_manual(values=as.character(ColorP$V2)) +
  ggtitle("GDB Volume histogram - data distribution") 
pl


dd<-ggplot(GDB.filt, aes(Max, color=Cell.Type, fill=Cell.Type)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=as.character(ColorP$V2)) +
  scale_color_manual(values=as.character(ColorP$V2)) +
  ggtitle("GDB Max anisotropy histogram - data distribution") 
pl


dd<-ggplot(GDB.filt, aes(Mid, color=Cell.Type, fill=Cell.Type)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=as.character(ColorP$V2)) +
  scale_color_manual(values=as.character(ColorP$V2)) +
  ggtitle("GDB Mid anisotropy histogram - data distribution") 
pl


dd<-ggplot(GDB.filt, aes(Min, color=Cell.Type, fill=Cell.Type)) 
pl<-dd+ geom_histogram(position ="dodge") + 
  scale_fill_manual(values=as.character(ColorP$V2)) +
  scale_color_manual(values=as.character(ColorP$V2)) +
  ggtitle("GDB Min anisotropy histogram - data distribution") 
pl


## TWO DATA VARIABLES ("TOP" AND "BOTTOM" SURFACE AREAS) ARE MISSING IN SOME STAGES --> I WON'T USE THEM IN THIS ANALYSIS
#Not scaled data
GDB.plot<-autoplot(prcomp(GDB.filt[,6:9],scale=FALSE),data=GDB.filt, colour="Cell.Type", shape="Stage",fill="Cell.Type")
GDB.plot+ scale_color_manual(values=ColVec) +
  scale_fill_manual(values=ColVec)  +
  scale_shape_manual(values = c(21, 21, 21,22,23,24,25))

#Scaled Data
GDB.plot<-autoplot(prcomp(GDB.filt[,6:9],scale=TRUE),data=GDB.filt, colour="Cell.Type", shape="Stage",fill="Cell.Type")
GDB.plot+ scale_color_manual(values=ColVec) +
  scale_fill_manual(values=ColVec)  +
  scale_shape_manual(values = c(21, 21, 21,22,23,24,25))





#All variables
GDB.plot<-autoplot(prcomp(scale(GDB.filt[,6:12],scale=TRUE)),data=GDB.filt, colour="Cell.Type",shape="Stage", fill="Cell.Type",
                   loadings = TRUE, loadings.colour = 'black',
                   loadings.label = TRUE, loadings.label.size = 5)

GDB.plot+ scale_color_manual(values=ColorP$V2)+scale_shape_manual(values = c(21, 21, 21,22,23,24,25))



#### CREATING SMALL GDB == GDB WITHOUT REDUNDANT COLUMNS


GDB.small<-data.frame(GDB.filt[,1:3],GDB.filt$Cell_Stage,GDB.filt$Volume,GDB.filt$Max,GDB.filt$Mid,GDB.filt$Min)
colnames(GDB.small)<- colnames(GDB.filt[c(1:3,5,9,8,7,6)])


# PLOT With Volume
ColorVec<-as.character(ColorP$V2)
#ColorVec<-c(ColorVec[1],ColorVec[2],ColorVec[6],ColorVec[4],ColorVec[9])# SPECIAL VECTOR COLOR



GDB.small.plot<-autoplot(prcomp(scale(GDB.small[,5:8])),data=GDB.small, colour="Cell.Type",shape="Stage", fill="Cell.Type")
                         #loadings = TRUE, loadings.colour = 'black',
       ### EIGENVECTORS  #loadings.label = TRUE, loadings.label.size = 5, alpha=0.7)
  
  #format 1                       
GDB.small.plot+scale_fill_manual(values=ColorVec) +
  scale_color_manual(values=rep("black",5)) +
  scale_shape_manual(values = c(21, 21, 21,22,23,24,25))                    

 #format2 
GDB.small.plot+
  scale_fill_manual(values=ColorVec) +
  scale_color_manual(values=ColorVec) +
  scale_shape_manual(values = c(21, 21, 21,22,23,24,25))    

# PLOT Without Volume
GDB.small.plot<-autoplot(prcomp(scale(GDB.small[,6:8])),data=GDB.small, colour="Cell.Type",shape="Stage", fill="Cell.Type",alpha=0.7)
                         #loadings = TRUE, loadings.colour = 'black',
  ### EIGENVECTORS OFF      #loadings.label = TRUE, loadings.label.size = 5)


GDB.small.plot+
  scale_fill_manual(values=ColorVec) +
  scale_color_manual(values=ColorVec) +
  scale_shape_manual(values = c(21, 21, 21,22,23,24,25))  


GDB.small.plot<-autoplot(prcomp(scale(GDB.small[,6:8])),data=GDB.small, colour="Cell.Type", fill="Cell.Type",alpha=0.7)


GDB.small.plot+
  scale_fill_manual(values=ColorVec) +
  scale_color_manual(values=ColorVec) + 
  facet_wrap(~ GDB.small$Stage,ncol=3) + theme_classic()

# SCATTER PLOTS PAIR-WISE VARIABLES

gglot( GDB)

## REMOVING "L1" and "L3" CEll TYPES --> THEY DON't HAVE SURFACE MEASUREMENTS
GDB.small.ratio<-GDB.small[which(GDB.small$Cell.Type != " L1"),]
GDB.small.ratio<-GDB.small.ratio[which(GDB.small.ratio$Cell.Type != " L3"),]




# GDB -- Making new Cell stage plots based on SEPARATED CELL.TYPE --> 

CC<-GDB.filt[which(as.character(GDB.filt$Cell.Type)=="CC"),]
L1<-GDB.filt[which(as.character(GDB.filt$Cell.Type)=="L1"),]
L2<-GDB.filt[which(as.character(GDB.filt$Cell.Type)=="L2"),]
L3<-GDB.filt[which(as.character(GDB.filt$Cell.Type)=="L3"),]

SMC<-GDB.filt[which(as.character(GDB.filt$Cell.Type)=="pSMC"),]



CC.plot<-autoplot(prcomp(scale(CC[,c(8,11,12)])),data=CC, 
             colour="Stage", fill="Stage", frame=TRUE, frame.colour="Stage", main="CC")
CC.plot+
  theme_classic()



L1.plot<-autoplot(prcomp(scale(L1[,c(8,11,12)])),data=L1, 
                  colour="Stage", fill="Stage", frame=TRUE, frame.colour="Stage", main="L1")
L1.plot+
  theme_classic()


L2.plot<-autoplot(prcomp(scale(L2[,c(8,11,12)])),data=L2, 
                  colour="Stage", fill="Stage", frame=TRUE, frame.colour="Stage", main="L2")
L2.plot+
  theme_classic()


L3.plot<-autoplot(prcomp(scale(L3[,c(8,11,12)])),data=L3, 
                   colour="Stage", fill="Stage", frame=TRUE, frame.colour="Stage", main="L3")
L3.plot+
  theme_classic()

SMC.plot<-autoplot(prcomp(scale(SMC[,c(8,11,12)])),data=SMC, 
                  colour="Stage", fill="Stage", frame=TRUE, frame.colour="Stage", main="SMC")
SMC.plot+
  theme_classic()





# GDB -- Making new Cell stage plots based on SEPARATED STAGES --> MAX MID MIN
#v1
NewColor<-c(rep("steelblue1",4),rep("darkblue",2))
#v2
NewColor<-c("violet",rep("steelblue1",3),rep("darkblue",2))

St<-prcomp(GDB.filt[,6:8],scale=TRUE)

DB.stg<-cbind(GDB.filt,St$x)

Alphavalue<-which(DB.stg$Cell.Type=="pSMC")

Alphavalue<-ifelse(DB.stg$Cell.Type=="pSMC", 1, 0.6)

St
summary(St)

St.plot<-autoplot(St,data=GDB.filt ,size=2, colour="Cell.Type", fill="Cell.Type",main= "St plots Max - Mid - Min",alpha=0.6,frame= TRUE, frame.colour="Cell.Type")
St.plot

Stg.plot<-ggplot() +
  geom_point(data=DB.stg, aes(PC1,PC2, color= Cell.Type, fill=Cell.Type, alpha=0.7))+
  stat_ellipse(data=DB.stg, aes(x=PC1,y=PC2,fill=Cell.Type),geom="polygon",type="norm", level=0.95, alpha=0.3)+
  facet_wrap(GDB.filt$Stage) +
  scale_fill_manual(values=NewColor) +
  scale_color_manual(values=NewColor) +
  theme_classic()
Stg.plot
fviz_eig(St, barcol="black", barfill= "grey")

PCA.vec<-fviz_pca_var(St,
                      col.var = "contrib", 
                      gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),col.circle = "white",
                      repel = TRUE) 
PCA.vec


pl<-Stg.plot+facet_wrap(GDB.filt$Stage)+
  scale_fill_manual(values=NewColor) +
  scale_color_manual(values=NewColor) +
  facet_wrap(GDB.filt$Stage) +theme_classic()


pl
ggsave("GDB/3.png",pl,dpi=600)


# Adding violet to CC
NewColor<-c("violet",rep("steelblue1",3),rep("darkblue",2))

pl<-Stg.plot+
  scale_fill_manual(values=NewColor) +
  scale_color_manual(values=NewColor) +
  facet_wrap(GDB.filt$Stage) +theme_classic()


pl
ggsave("GDB/3.png",pl,dpi=600)


pl+geom_point(data=DB.stg, aes(PC1,PC2, color= Cell.Type, shape= Stage, fill=Cell.Type))+
stat_ellipse(data=DB.stg, aes(x=PC1,y=PC2,fill=Stage,color=Stage),geom="polygon", level=0.95,alpha=0.2)


### Original


NewColor<-c(rep("steelblue1",4),rep("darkblue",2))

St<-GDB.filt[,c(1:3,8,11,12)]
St<-droplevels(St)
St.plot<-autoplot(St,data=GDB.filt ,size=2, colour="Cell.Type", fill="Cell.Type",main= "St plots Max - Mid - Min",alpha=0.7)

pl<-St.plot+
  scale_fill_manual(values=NewColor) +
  scale_color_manual(values=NewColor) +
  facet_wrap(GDB.filt$Stage) 


pl
ggsave("GDB/3.png",pl,dpi=600)

