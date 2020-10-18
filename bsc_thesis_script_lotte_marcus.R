################################################################################
###Transcriptomics - alphasyn timeseries
################################################################################

###Set your work directory
setwd("~/Desktop/THESIS/IMAC/Verwerkte_data/") ###set this to the folder you want to use in R
workwd <- getwd()
filename <- "alphasynRIL"

###Load pre-made functions
#uses eQTL pipeline functions https://git.wur.nl/mark_sterken/eQTL_pipeline
#     transcriptomics functions https://git.wur.nl/mark_sterken/Transcriptomics.workbench
git_dir <- "~/Desktop/THESIS/IMAC/R_scripts/"   #the folder in which you place the NEMA_functions folder
source(paste(git_dir,"/NEMA_functions/Loader.R",sep=""))

################################################################################
###Dependencies
################################################################################

install <- FALSE ####first time you run it; set to true to install packages!
if(install){
  install.packages("tidyverse")
  install.packages("colorspace")
  install.packages("RColorBrewer")
  install.packages("BiocManager")
  BiocManager::install("limma")
  BiocManager::install("statmod")
  install.packages("gridExtra")
  install.packages("VennDiagram")
  install.packages("openxlsx")
  install.packages("rmarkdown")
}
###load
library("colorspace")
library("RColorBrewer")
library(limma)
library(gridExtra)
library("VennDiagram")
library(openxlsx)
library("rmarkdown")
library(tidyverse)

################################################################################
###Plotting theme, colours
################################################################################


###Set plotting theme
presentation <- theme(axis.text.x = element_text(size=10, face="bold", color="black"),
                      axis.text.y = element_text(size=10, face="bold", color="black"),
                      axis.title.x = element_text(size=12, face="bold", color="black"),
                      axis.title.y = element_text(size=12, face="bold", color="black"),
                      strip.text.x = element_text(size=12, face="bold", color="black"),
                      strip.text.y = element_text(size=12, face="bold", color="black"),
                      plot.title = element_text(size=14, face="bold"),
                      strip.background = element_rect(fill= "grey80",color="black"),
                      panel.background = element_rect(fill = "white",color="black"),
                      panel.grid.major = element_line(colour = "grey80"),
                      panel.grid.minor = element_blank(),
                      legend.position = "right")


blank_theme <- theme(plot.background = element_blank(),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank(),
                     axis.title.x = element_blank(),
                     axis.title.y = element_blank(),
                     axis.text.x = element_blank(),
                     axis.text.y = element_blank(),
                     axis.ticks = element_blank())



###Here you can set colours for plotting in theme using ggplot2
#display.brewer.all()
myColors <- c(brewer.pal(9,"Set1")[c(2,5)],brewer.pal(9,"Purples")[c(4,6,6)],"black","darkgrey","black",brewer.pal(12,"Paired")[c(1,7)])


names(myColors) <- c("CB4856","N2","aSlight","aS","trans","cis","notsig","RIL","SCH4856","NL5901")

colScale <- scale_colour_manual(name = "Treatment",values = myColors)
fillScale <- scale_fill_manual(name = "Treatment",values = myColors)

################################################################################
###load all the data
################################################################################


# ###Load normalized data
#     ###all quality is good
#     load(file = "./aSRIL_datafiles/Ruwe_data/obj_list.data.Rdata")
# 
# ###populations
#     load(file = "./aSRIL_datafiles/Ruwe_data/obj_popmap.Rdata")
#     load(file = "./aSRIL_datafiles/obj_popmrk.Rdata")
# 
# ###eQTL mapping
#     load(file="./aSRIL_datafiles/obj_aS.eQTL.table.Rdata")
#     load(file="./aSRIL_datafiles/obj_peak.aS.QTL.Rdata")
#     load(file="./aSRIL_datafiles/obj_aS.QTL.Rdata")
#     
# ###load qPCR data
#     load(file = "./aSRIL_datafiles/obj_qpcr_data.Rdata")
#     load(file = "./aSRIL_datafiles/obj_qpcr.QTL.Rdata")

### normalized data 
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_list.data.Rdata")

## Introgression Lines data
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_list.data.IL.Rdata")

### QTL & eQTL mapping
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_aS.eQTL.table.Rdata")
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_aS.QTL.Rdata")
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_peak.aS.QTL.Rdata")

##Populations##
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_popmap.Rdata")
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_popmrk.Rdata")

##qPCR data 
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_qpcr_data.Rdata")
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_qpcr.QTL.Rdata")

##excell sheet with strains containg the used introgression lines
strains_containing_IL <- read.xlsx(xlsxFile = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/strains_containing_IL.xlsx") 

################################################################################
###plot ILs
################################################################################

strain1 <- wur.pop.map[,53]; strain1[strain1==-1] <- 2
strain2 <- wur.pop.map[,69]

pdf("F1.pdf",width = 4,height=4)
plot.genotype.strain(wur.pop.marker,strain.map1 = strain1,strain.map2 = strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
dev.off()

pdf("P1.pdf",width = 8,height=4)   
par(mfrow=c(1,2))
plot.genotype.strain(wur.pop.marker,strain.map1 = strain1,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
plot.genotype.strain(wur.pop.marker,strain.map1 = strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
dev.off()

pdf("F2_after_selfing.pdf",width = 16,height=4)
par(mfrow=c(1,4))
plot.genotype.strain(wur.pop.marker,strain.map1 = rep(1,length(strain1)),both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
plot.genotype.strain(wur.pop.marker,strain.map1 = strain1*strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
plot.genotype.strain(wur.pop.marker,strain.map1 = strain1,strain.map2 = strain2,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
plot.genotype.strain(wur.pop.marker,strain.map1 = strain2,strain.map2 = strain1,both.chr = TRUE,col.lines=c("#1F78B4","grey","#FF7F00","darkgreen"))
dev.off()

################################################################################
###Polymorphic genes in a region
################################################################################

CB4856.DB[[1]]

aS.eQTL.table[1,]

dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
              gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
  dplyr::select(gene_sequence_name) %>%
  dplyr::filter(!duplicated(gene_sequence_name))


genes <- dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
                       gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
  dplyr::select(gene_sequence_name) %>%
  dplyr::filter(!duplicated(gene_sequence_name)) %>%
  unlist() %>%
  as.character()

filter(CB4856.DB[[2]],Sequence_name %in% genes)


################################################################################
###plot QTL profile
################################################################################

###ignore the warnings; you can select any gene
head(aS.eQTL.table)

data.plot <- prep.ggplot.QTL.profile(peak.aS.QTL,aS.QTL,"AGIWUR1000")
data.plot[[2]] <- mutate(data.plot[[2]], geno_strain=ifelse(genotype==-1,"CB4856","N2"))

sf2a <- ggplot(data.plot$QTL_profile,aes(x=qtl_bp,y=qtl_significance,alpha=0.2)) +
  geom_line(size=1.5,colour=brewer.pal(9,"Set1")[3]) + facet_grid(~qtl_chromosome,scales="free",space="free_x") + presentation + theme(legend.position = "none") +
  geom_abline(intercept=3.1,slope=0,linetype=2,size=1) + labs(x="QTL position (Mb)",y="significance (-log10(p))") +
  scale_x_continuous(breaks=c(5,10,15,20)*10^6,labels=c(5,10,15,20)) + ylim(0,5.5)

sf2b <- ggplot(data.plot[[2]],aes(x=geno_strain,y=trait_value)) +
  geom_jitter(height=0,width=0.25,aes(colour=geno_strain),alpha=0.2) + geom_boxplot(outlier.shape=NA,alpha=0.2,aes(fill=geno_strain)) +
  xlab("Genotype\nat marker") + ylab(data.plot[[2]]$trait[1]) + facet_grid(~Chromosome+Position) +
  presentation + colScale + fillScale + theme(legend.position = "none") +
  annotate("text",x=1.5,y=max(data.plot[[2]]$trait_value,na.rm=T),label=paste("italic(R)^{2}==",round(data.plot[[2]]$R_squared[1],digits=2),sep=""),parse=TRUE)

annotation.grobA <- title.grob <- textGrob(label = "A",x = unit(0, "lines"),y = unit(0, "lines"),hjust = 0, vjust = 0,gp = gpar(fontsize = 20,fontface="bold"))
annotation.grobB <- title.grob <- textGrob(label = "B",x = unit(0, "lines"),y = unit(0, "lines"),hjust = 0, vjust = 0,gp = gpar(fontsize = 20,fontface="bold"))

sf2a <- arrangeGrob(sf2a,top=annotation.grobA)
sf2b <- arrangeGrob(sf2b,top=annotation.grobB)

grid.arrange(sf2a,sf2b,widths=c(2,1))

pdf("Figure-eQTL_profile.pdf",width=8,height=4)
grid.arrange(sf2a,sf2b,widths=c(2,1))
dev.off()

################################################################################
###Select trans-bands ///// 04-05-20
################################################################################

transtable <- dplyr::filter(aS.eQTL.table,trans_band != "none") %>%
  dplyr::filter(!duplicated(trans_band)) %>%
  dplyr::select(trans_band) %>%
  tidyr::separate(trans_band,into=c("chromosome","rest"),sep=":",remove = FALSE) %>%
  tidyr::separate(rest,into=c("loc_left","loc_right"),sep="-") %>%
  dplyr::mutate(chromosome=gsub("chr","",chromosome), loc_right=gsub("Mb","",loc_right)) %>%
  dplyr::mutate(loc_left=1e6*as.numeric(loc_left),loc_right=1e6*as.numeric(loc_right))

output <- NULL
for(i in 1:nrow(transtable)){
  genes <- dplyr::filter(Agilent.Ce.V2,chromosome==transtable[i,2],!is.na(chromosome),
                         gene_bp_end > transtable[i,3], gene_bp_start < transtable[i,4]) %>%
    dplyr::select(gene_sequence_name) %>%
    dplyr::filter(!duplicated(gene_sequence_name)) %>%
    unlist() %>%
    as.character()
  
  tmp <- filter(CB4856.DB[[2]],Sequence_name %in% genes)
  
  output <- rbind(output,
                  cbind(transtable[i,],tmp)
  )
}

###check if is in list
#filter(output,trans_band %in% c("chrV:14.5-17.0Mb","chrV:5.5-6.0Mb"))
#filter(output,Sequence_name %in% filter(aS.eQTL.table, qtl_type == "cis")$gene_public_name)

###check if is not in a list
#filter(output,!trans_band %in% c("chrV:14.5-17.0Mb","chrV:5.5-6.0Mb"))

###check if larger or smaller then
#filter(output, loc_left < 1000000)

###check if not equal to
#filter(output, chromosome != "V")

################################################################################
###check gene expression at trans-band  ///// 04-05-20
################################################################################

transtable <- dplyr::filter(aS.eQTL.table,trans_band != "none") %>%
  dplyr::filter(!duplicated(trans_band)) %>%
  dplyr::select(trans_band) %>%
  tidyr::separate(trans_band,into=c("chromosome","rest"),sep=":",remove = FALSE) %>%
  tidyr::separate(rest,into=c("loc_left","loc_right"),sep="-") %>%
  dplyr::mutate(chromosome=gsub("chr","",chromosome), loc_right=gsub("Mb","",loc_right)) %>%
  dplyr::mutate(loc_left=1e6*as.numeric(loc_left),loc_right=1e6*as.numeric(loc_right))


for(i in 1:nrow(transtable)){
  
  data.plot <- filter(aS.eQTL.table,trans_band == transtable[i,1])
  
  ggplot(data.plot,aes(x=qtl_effect,y=qtl_significance,size=qtl_R2_sm)) +
    geom_point()
  print(
    ggplot(data.plot,aes(x=qtl_effect,y=PL_alpha_effect)) +
      geom_point() + geom_smooth(method="lm")
  )
}

CB4856.DB[[1]]

aS.eQTL.table[1,]

dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
              gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
  dplyr::select(gene_sequence_name) %>%
  dplyr::filter(!duplicated(gene_sequence_name))


genes <- dplyr::filter(Agilent.Ce.V2,chromosome=="V",!is.na(chromosome),
                       gene_bp_end > 15397385, gene_bp_start < 18222746) %>%
  dplyr::select(gene_sequence_name) %>%
  dplyr::filter(!duplicated(gene_sequence_name)) %>%
  unlist() %>%
  as.character()

filter(CB4856.DB[[2]],Sequence_name %in% genes)

################################################################################
###check gene expression at trans-band  ///// 04-05-20
################################################################################

filter(aS.eQTL.table,trait=="AGIWUR1000") %>%
  merge(list.data,by.x=1,by.y=2)

################################################################################
###check gene expression at trans-band  ///// 07-05-20
################################################################################

###IL data
load(file = "~/Desktop/THESIS/IMAC/Ruwe_data/aSRIL_datafiles/obj_list.data.IL.Rdata")
list.data.IL <- filter(list.data.IL,strain_type !="RIL")
save(list.data.IL,file="obj_list.data.IL.Rdata")

eQTL.table.DB[[1]] ###matrix
class(eQTL.table.DB[[1]])

head(eQTL.table.DB[[2]]) ###dataframe
class(eQTL.table.DB[[2]])

class(eQTL.table.DB)

###factor problem
a <- round(runif(10,0,10)); a
sum(a)
b <- factor(a)
sum(b)

as.numeric(b) ###mind, the indexes are returned, not the real numbers!

as.numeric(as.character(unlist(b)))

summary(eQTL.table.DB[[2]])


###Test if trans-bands are unique

head(aS.eQTL.table)

transtable <- dplyr::filter(aS.eQTL.table,trans_band != "none") %>%
  dplyr::filter(!duplicated(trans_band)) %>%
  dplyr::select(trans_band) %>%
  tidyr::separate(trans_band,into=c("chromosome","rest"),sep=":",remove = FALSE) %>%
  tidyr::separate(rest,into=c("loc_left","loc_right"),sep="-") %>%
  dplyr::mutate(chromosome=gsub("chr","",chromosome), loc_right=gsub("Mb","",loc_right)) %>%
  dplyr::mutate(loc_left=1e6*as.numeric(loc_left),loc_right=1e6*as.numeric(loc_right))

head(eQTL.table.DB[[9]])

###one trans-band
###correlation
i <- transtable[1,1]


data.test <- filter(aS.eQTL.table,trans_band==i)



data.plot <- filter(eQTL.table.DB[[9]], gene_WBID %in% data.test$gene_WBID)

ggplot(data.plot,aes(x=qtl_bp)) +
  geom_histogram(binwidth = 500000) + facet_grid(~qtl_chromosome)


data.test <- filter(aS.eQTL.table,trans_band!="none")

data.plot <- dplyr::select(data.test,gene_WBID,trans_band) %>%
  merge(eQTL.table.DB[[2]],by.x=1,by.y=14)

table(data.plot$trans_band)


ggplot(data.plot,aes(x=qtl_bp)) +
  geom_histogram(binwidth = 500000) + facet_grid(trans_band~qtl_chromosome)



########################################################################################################
###correlation
########################################################################################################

### -----> add extra column log2_ratio_N2*-1 = log2_ratio_N2_corrected
list.data.IL <- dplyr::mutate(list.data.IL, log2_ratio_N2_corrected=-1*log2_ratio_N2)


### -----> Before running the correlation!

for(i in transtable[,1])
  
  ##for(i in transtable[,1])   LOOP
  
########################################################################################################        
######## Correlation TRANSBAND 1 - chrI:12.5-13.0Mb (in transtable = 9)

for(i in transtable[,1]){   
  i <- transtable[9,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN210","WN211","WN212"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))

write.xlsx(data.plot, file = "IL_bonferroni_TB1.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.5)

write.xlsx(data.plot, file = "IL_bonferroni_TB1.xlsx")
###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}    

dev.off()


########################################################################################################        
######## Correlation TRANSBAND 2 - chrII:6.0-6.5Mb (in transtable = 4)

for(i in transtable[,1]){   
  i <- transtable[4,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN222","WN223"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB2.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.25)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}        

dev.off()
########################################################################################################        
######## Correlation TRANSBAND 3 - chrIV:14.5-15.0Mb (in transtable = 6)

for(i in transtable[,1]){   
  i <- transtable[6,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN260","WN261","WN262"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB3.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.34)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}       

dev.off()

########################################################################################################        
######## Correlation TRANSBAND 4 - chrV:11.0-12.5Mb (in transtable = 3)

for(i in transtable[,1]){   
  i <- transtable[3,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN268","WN269","WN270"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB4.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.35)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}        

dev.off()

########################################################################################################        
######## Correlation TRANSBAND 5 - chrV:14.5-17.0Mb (in transtable = 1) 

for(i in transtable[,1]){   
  i <- transtable[1,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN269","WN270"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB5.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.6)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}   

dev.off()

########################################################################################################        
######## Correlation TRANSBAND 6 - chrV:5.5-6.5Mb (in transtable = 8)

for(i in transtable[,1]){   
  i <- transtable[8,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN267","WN268","WN269"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB6.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.3)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

dev.off()

########################################################################################################        
######## Correlation TRANSBAND 7 - chrX:14.0-14.5Mb (in transtable = 7)

for(i in transtable[,1]){   
  i <- transtable[7,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN283","WN289","WN290"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB7.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.2)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}  

dev.off()

########################################################################################################        
######## Correlation TRANSBAND 8 - chrX:6.0-6.5Mb (in transtable = 2)

for(i in transtable[,1]){   
  i <- transtable[2,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN281","WN284","WN286"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB8.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.45)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}        

dev.off()

########################################################################################################        
######## Correlation TRANSBAND 9 - chrX:8.5-9.0Mb (in transtable = 5)

for(i in transtable[,1]){   
  i <- transtable[5,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN281","WN283","WN286"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB9.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.45)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}        


dev.off()

### -----> add extra column log2_ratio_N2*-1 = log2_ratio_N2_corrected
list.data.IL <- dplyr::mutate(list.data.IL, log2_ratio_N2_corrected=-1*log2_ratio_N2)

############################### CIS OR TRANS-EQTL; THAT'S THE QUESTION #################################
########################################################################################################        
######## CHECK IF TRANSBAND 2 ((chrII:6.0-6.5Mb (in transtable = 4)) IS A CIS OR TRANS_QTL --> FILTER OUT CHROMOSOME II 
for(i in transtable[,1]){   
  i <- transtable[4,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,qtl_chromosome != 'II', trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN222","WN223"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1,
                  dplyrr::filter(qtl_chromosome != 'II')) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  
  data.frame()

  ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))

write.xlsx(data.plot, file = "IL_bonferroni_TB2_cis_or_trans.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.25)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}        

dev.off()

########################################################################################################        
######## CHECK IF TRANSBAND 4 (chrV:11.0-12.5Mb (in transtable = 3)) IS A CIS OR TRANS_QTL --> FILTER OUT CHROMOSOME V 

for(i in transtable[,1]){   
  i <- transtable[3,1]
  
  data.test.eQTL <- filter(aS.eQTL.table,qtl_chromosome != 'V', trans_band==i)
  data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% c("WN261","WN262","WN261"))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}

data.test.eQTL <- filter(aS.eQTL.table,trans_band==i)
data.test.ILs <- filter(list.data.IL,SpotID %in% data.test.eQTL$trait)

### --> Correlation (I suppose bonferroni is better for this!)
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1,
                  dplyrr::filter(qtl_chromosome != 'V')) %>%
  group_by(strain,treatment) %>% #GROUP BY TREATMENT
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  
  data.frame()

ggplot(data.plot,aes(x=strain,y=correlation)) +
  geom_point() + facet_grid(~treatment)


### --> CORRELATION bonferroni
data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                  dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                  by.x = 1, by.y = 1) %>%
  group_by(strain,treatment) %>%
  summarise(correlation=cor(qtl_effect,log2_ratio_N2_corrected),
            pvalue=cor.test(qtl_effect,log2_ratio_N2_corrected)$p.value) %>%
  data.frame() %>%
  dplyr::mutate(bonf=p.adjust(pvalue,method="bonferroni"))


write.xlsx(data.plot, file = "IL_bonferroni_TB4_cis_trans.xlsx")

### --> filter correlation

filter(data.plot, correlation > 0.35)

###for the specific ILs

for(i in 1:nrow(strains_containing_IL)){   
  
  data.test.eQTL <- dplyr::filter(aS.eQTL.table,trans_band==strains_containing_IL[i,1])
  data.test.ILs <- dplyr::filter(list.data.IL,SpotID %in% data.test.eQTL$trait, tolower(strain) %in% tolower(strains_containing_IL[i,2]))
  
  data.plot <-merge(dplyr::select(data.test.eQTL,trait,qtl_effect,gene_WBID),
                    dplyr::select(data.test.ILs,SpotID,treatment,strain,log2_ratio_N2_corrected),
                    by.x = 1, by.y = 1)
  
  for(j in 1:length(unique(data.plot$strain))){
    print(ggplot(dplyr::filter(data.plot,strain==unique(data.plot$strain)[j]),
                 aes(x=qtl_effect,y=log2_ratio_N2_corrected,colour=treatment)) +
            geom_point() + geom_smooth(method = "lm") + labs(title = unique(data.plot$strain)[j]))
  }
}        

dev.off()




##########################################################################
########## RIL ANALYSIS IN DIFFERENT STUDIES     30-06-2020
##########################################################################


eQTL.table.DB[[1]]

plots <- as.list(NULL)
tables <- as.list(NULL)
more.than.expected <- NULL
tables2 <- as.list(NULL)
for(i in 2:12)
  {
  data.test <- filter(aS.eQTL.table,trans_band!="none")
  
  ##how many unique genes in aS study?
  unique.aS <- length(unique(data.test$gene_WBID))
  
  ##how many unique genes in compared study?
  unique.compared <- length(unique(eQTL.table.DB[[i]]$gene_WBID))
  
  data.plot <- dplyr::select(data.test,gene_WBID,trans_band) %>%
    merge(eQTL.table.DB[[i]],by.x=1,by.y=14)
  
  ##how many unique genes overlap?
  overlap <- length(unique(data.plot$gene_WBID))
  
  total_genes <- length(unique(Agilent.Ce.V2$gene_WBID))
  
  more.than.expected[i] <- phyper(q=overlap,m=unique.compared,n=total_genes-unique.compared,k=unique.aS,lower.tail=F)
  
  tables[[i]] <- table(data.plot$trans_band,data.plot$qtl_type)
  tables2[[i]] <- table(eQTL.table.DB[[i]]$qtl_type)
  
  plots[[i]] <- ggplot(data.plot,aes(x=qtl_bp,fill=qtl_type)) +
    geom_histogram(binwidth = 500000) + facet_grid(trans_band~qtl_chromosome) +
    labs(x=paste(eQTL.table.DB[[1]][i,c(1,3)],collapse="_"))

  write.xlsx(tables, file = "cistransstudies.xlsx")
  write.xlsx(tables2, file = "numberofeqtlperstudy.xlsx")
}  
  
  
###print plots to pdf
pdf("overlap_RIL_studies.pdf",width=8,height=8)
for(i in 2:12){
  print(plots[[i]])
}
dev.off()

####ENRICHEMENT ANALYSIS 

###NUMBER OF GENES PER DIFFERENT STUDY

#Study_1
length(unique(filter(eQTL.table.DB[[2]],qtl_type=="trans")$gene_public_name))

#Study_2
length(unique(filter(eQTL.table.DB[[3]],qtl_type=="trans")$gene_public_name))

#Study_3
length(unique(filter(eQTL.table.DB[[4]],qtl_type=="trans")$gene_public_name))

#Study_4
length(unique(filter(eQTL.table.DB[[5]],qtl_type=="trans")$gene_public_name))

#Study_5
length(unique(filter(eQTL.table.DB[[6]],qtl_type=="trans")$gene_public_name))

#Study_6
length(unique(filter(eQTL.table.DB[[7]],qtl_type=="trans")$gene_public_name))

#Study_7
length(unique(filter(eQTL.table.DB[[8]],qtl_type=="trans")$gene_public_name))

#Study_8
length(unique(filter(eQTL.table.DB[[9]],qtl_type=="trans")$gene_public_name))

#Study_9
length(unique(filter(eQTL.table.DB[[10]],qtl_type=="trans")$gene_public_name))

#Study_10
length(unique(filter(eQTL.table.DB[[11]],qtl_type=="trans")$gene_public_name))

#Study_11
length(unique(filter(eQTL.table.DB[[12]],qtl_type=="trans")$gene_public_name))

###Total number of genes ###

length(unique(Agilent.Ce.V2$gene_public_name))

phyper(3,326,77,18448,lower.tail = FALSE)
      