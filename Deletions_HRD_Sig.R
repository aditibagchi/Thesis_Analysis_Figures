##Main Objective is to generte a data frame with Insertions, Deletions, HRD signature and also large deletions. 
## Requires merging of multiple MAF files
##Large deletions are defined as any deletion >5 nuceolite in lenght. Calculated by sustracting varant end position from varaint start position.
##maftools::

##The recurrent(primary site) samples with matched primary and germline
##The recurrent samples called against gremline. 

RM_Germline <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/9_recurrent_germline_samples_GRCH38_crossmapped_to_GRCH37.maf.txt")

RM_Germline_Variants <- mafSummary(RM_Germline)
RM_Germline_Variants_Type <- RM_Germline_Variants$variant.type.summary

##View(RM_Germline_Variants_Type)

RM_Germline_Variants_Type$Name <- c("Tumor14-R", "Tumor3-R", "Tumor2-R","Tumor18-R","Tumor19-R","Tumor20-R","Tumor21-R","Tumor22-R","Tumor23-R")
RM_Germline_Variants_Type$Indel <- RM_Germline_Variants_Type$DEL + RM_Germline_Variants_Type$INS
RM_Germline_Variants_Type$IndelFR <- RM_Germline_Variants_Type$Indel/RM_Germline_Variants_Type$total
RM_Germline_Variants_Type$Sig3 <- c(0, 0.176, 0.13, 0.074, 0.264, 0.101, 0.098, 0.198, 0.066)
RM_Germline_Variants_Type$Sig8 <- c(0.0, 0.334, 0.265, 0.17, 0.0, 0.184, 0.304, 0.162, 0.151) 
RM_Germline_Variants_Type$HRD <- RM_Germline_Variants_Type$Sig3 + RM_Germline_Variants_Type$Sig8
##View(RM_Germline_Variants_Type)
write.csv(RM_Germline_Variants_Type, file = "RM_Germline_Variants_Type.csv")

#lenght of deletions
RM_Germline_DF <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/9_recurrent_germline_samples_GRCH38_crossmapped_to_GRCH37.maf.txt", 
                             header=TRUE, comment.char="#")
#View(RM_Germline_DF)
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
RM_Germline_DF_MB_XX <- subset(RM_Germline_DF, RM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-XX")
RM_Germline_DF_MB_XX_del <- subset(RM_Germline_DF_MB_XX, RM_Germline_DF_MB_15$Variant_Type == "DEL")
#View(RM_Germline_DF_MB_XX_del)
RM_Germline_DF_MB_XX_del$Del_Len <- RM_Germline_DF_MB_XX_del$End_Position - RM_Germline_DF_MB_XX_del$Start_Position
#View(RM_Germline_DF_MB_XX_del)
median(RM_Germline_DF_MB_XX_del$Del_Len)
mean(RM_Germline_DF_MB_XX_del$Del_Len)
RM_Germline_DF_MB_XX_del$LargeDel = ifelse(RM_Germline_DF_MB_XX_del$Del_Len>=5,TRUE,FALSE)
RM_Germline_DF_MB_XX_del_large = RM_Germline_DF_MB_XX_del[RM_Germline_DF_MB_XX_del$LargeDel==TRUE,]
#View(RM_Germline_DF_MB_XX_del_large)
write.csv(RM_Germline_DF_MB_XX_del_large, file = "RM_Germline_DF_MB_XX_del_large.csv")
write.csv(RM_Germline_DF_MB_XX_del, file = "RM_Germline_DF_MB_XX_del.csv")


##The primary samples with matched recurrent and germline  
##MB-REC-16 and MB-REC-06 was removed , as these two samples were reanalyzed at a later date 
PM_Germline <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/primary_germline_15_files.maf.txt")
PM_Germline_Variants <- mafSummary(PM_Germline)
PM_Germline_Variants_Type <- PM_Germline_Variants$variant.type.summary
View(PM_Germline_Variants_Type)
PM_Germline_Variants_Type <- PM_Germline_Variants_Type[-c(1,2,6),]
PM_Germline_Variants_Type$Name <- c("Tumor2","Tumor3","Tumor20","Tumor23","Tumor4", "Tumor9", "Tumor22", "Tumor17", "Tumor19","Tumor21","Tumor18", "Tumor8")
PM_Germline_Variants_Type$Indel <- PM_Germline_Variants_Type$DEL + PM_Germline_Variants_Type$INS
PM_Germline_Variants_Type$IndelFR <- PM_Germline_Variants_Type$Indel/PM_Germline_Variants_Type$total
PM_Germline_Variants_Type$Sig3 <- c(0.0, 0.0,0.101, 0.0, 0.164, 0.342, 0.0, 0.0, 0.116, 0.0, 0.201, 0.094)
PM_Germline_Variants_Type$Sig8 <- c(0.305, 0.475, 0.119, 0.328, 0.0, 0.0, 0.218, 0.15, 0.0, 0.293, 0.139, 0.0) 
PM_Germline_Variants_Type$HRD <- PM_Germline_Variants_Type$Sig3 + PM_Germline_Variants_Type$Sig8
View(PM_Germline_Variants_Type)
write.csv(PM_Germline_Variants_Type, file = "PM_Germline_Variants_Type.csv")


#lenght of deletions
PM_Germline_DF <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/primary_germline_15_files.maf.txt", 
                             header=TRUE, comment.char="#")
View(PM_Germline_DF)
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
PM_Germline_DF_MB_XX <- subset(PM_Germline_DF, PM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-XX")
PM_Germline_DF_MB_XX_del <- subset(PM_Germline_DF_MB_XX, PM_Germline_DF_MB_XX$Variant_Type == "DEL")
View(PM_Germline_DF_MB_XX_del)
PM_Germline_DF_MB_XX_del$Del_Len <- PM_Germline_DF_MB_XX_del$End_Position - PM_Germline_DF_MB_XX_del$Start_Position
PM_Germline_DF_MB_XX_del$LargeDel = ifelse(PM_Germline_DF_MB_XX_del$Del_Len>=5,TRUE,FALSE)
PM_Germline_DF_MB_XX_del_large = PM_Germline_DF_MB_XX_del[PM_Germline_DF_MB_XX_del$LargeDel==TRUE,]
View(PM_Germline_DF_MB_XX_del_large)
write.csv(PM_Germline_DF_MB_XX_del_large, file = "PM_Germline_DF_MB_XX_del_large.csv")
write.csv(PM_Germline_DF_MB_XX_del, file = "PM_Germline_DF_MB_XX_del.csv")

##MB_REC_16 primary 
PM_Germline_16 <- read.maf( maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/MB-REC-16_primary.crossmap.hg38_to_hg19.vep.maf.txt")
PM_Germline_Variants_16 <- mafSummary(PM_Germline_16)

PM_Germline_Variants_Type_16 <- PM_Germline_Variants_16$variant.type.summary
PM_Germline_Variants_Type_16$Name <- c("Tumor14")
PM_Germline_Variants_Type_16$Indel <- PM_Germline_Variants_Type_16$DEL + PM_Germline_Variants_Type_16$INS
PM_Germline_Variants_Type_16$IndelFR <- PM_Germline_Variants_Type_16$Indel/PM_Germline_Variants_Type_16$total
PM_Germline_Variants_Type_16$Sig3 <- c(0.0)
PM_Germline_Variants_Type_16$Sig8 <- c(0.0)
PM_Germline_Variants_Type_16$HRD <- PM_Germline_Variants_Type_16$Sig3 + PM_Germline_Variants_Type_16$Sig8
View(PM_Germline_Variants_Type_16)
write.csv(PM_Germline_Variants_Type_16, file = "PM_Germline_Variants_Type_16.csv")
PM_Germline_Variants_Type_16$Tumor_Sample_Barcode <- c("MB_REC_16")

#MB_REC_16-primary and MB_REC_06 both primary and Recurrent 
#lenght of deletions
PM_Germline_DF_16 <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/MB-REC-16_primary.crossmap.hg38_to_hg19.vep.maf.txt", header=TRUE, comment.char="#")
View(PM_Germline_DF_16)
PM_Germline_DF_16$Tumor_Sample_Barcode <- c("MB_REC_16")
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
PM_Germline_DF_16 <- subset(PM_Germline_DF_16, PM_Germline_DF_16$Tumor_Sample_Barcode == "MB-REC-16")
PM_Germline_DF_MB_16_del <- subset(PM_Germline_DF_16, PM_Germline_DF_16$Variant_Type == "DEL")
View(PM_Germline_DF_MB_16_del)
PM_Germline_DF_MB_16_del$Del_Len <- PM_Germline_DF_MB_16_del$End_Position - PM_Germline_DF_MB_16_del$Start_Position
PM_Germline_DF_MB_16_del$LargeDel = ifelse(PM_Germline_DF_MB_16_del$Del_Len>=5,TRUE,FALSE)
PM_Germline_DF_MB_16_del_large = PM_Germline_DF_MB_16_del[PM_Germline_DF_MB_16_del$LargeDel==TRUE,]
View(PM_Germline_DF_MB_16_del_large)
write.csv(PM_Germline_DF_MB_16_del_large, file = "PM_Germline_DF_MB_16_del_large.csv")
write.csv(PM_Germline_DF_MB_16_del, file = "PM_Germline_DF_MB_16_del.csv")


##MB_REC_06 Recurrent
RM_Germline_06 <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-06_recurrent_parental_germline_filtered.maf.txt")
RM_Germline_Variants_06 <- mafSummary(RM_Germline_06)

RM_Germline_Variants_Type_06 <- RM_Germline_Variants_06$variant.type.summary
RM_Germline_Variants_Type_06$Name <- c("Tumor13-R")
RM_Germline_Variants_Type_06$Indel <- RM_Germline_Variants_Type_06$DEL + RM_Germline_Variants_Type_06$INS
RM_Germline_Variants_Type_06$IndelFR <- RM_Germline_Variants_Type_06$Indel/RM_Germline_Variants_Type_06$total
RM_Germline_Variants_Type_06$Sig3 <- c(0.0)
RM_Germline_Variants_Type_06$Sig8 <- c(0.0)
RM_Germline_Variants_Type_06$HRD <- RM_Germline_Variants_Type_06$Sig3 + RM_Germline_Variants_Type_06$Sig8
View(RM_Germline_Variants_Type_06)
write.csv(RM_Germline_Variants_Type_06, file = "RM_Germline_Variants_Type_06.csv")
RM_Germline_Variants_Type_06$Tumor_Sample_Barcode <- c("MB_REC_06")


#lenght of deletions
RM_Germline_DF_MB_06 <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/MB-REC-06_recurrent_parental_germline_filtered.maf.txt", header=TRUE)
View(RM_Germline_DF_06)
RM_Germline_DF_MB_06$Tumor_Sample_Barcode <- c("MB-REC-06")
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.

RM_Germline_DF_MB_06_del <- subset(RM_Germline_DF_MB_06, RM_Germline_DF_MB_06$Variant_Type == "DEL")
View(RM_Germline_DF_MB_06_del)
RM_Germline_DF_MB_06_del$Del_Len <- RM_Germline_DF_MB_06_del$End_Position - RM_Germline_DF_MB_06_del$Start_Position
RM_Germline_DF_MB_06_del$LargeDel = ifelse(RM_Germline_DF_MB_06_del$Del_Len>=5,TRUE,FALSE)
RM_Germline_DF_MB_06_del_large = RM_Germline_DF_MB_06_del[RM_Germline_DF_MB_06_del$LargeDel==TRUE,]
View(RM_Germline_DF_MB_06_del_large)
write.csv(RM_Germline_DF_MB_06_del_large, file = "RM_Germline_DF_MB_06_del_large.csv")
write.csv(RM_Germline_DF_MB_06_del, file = "RM_Germline_DF_MB_06_del.csv")



##MB_REC_06 primary 
PM_Germline_06 <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-06_primary_parental_germline_filtered.maf.txt")
PM_Germline_Variants_06 <- mafSummary(PM_Germline_06)

PM_Germline_Variants_Type_06 <- PM_Germline_Variants_06$variant.type.summary
PM_Germline_Variants_Type_06$Name <- c("Tumor13")
PM_Germline_Variants_Type_06$Indel <- PM_Germline_Variants_Type_06$DEL + PM_Germline_Variants_Type_06$INS
PM_Germline_Variants_Type_06$IndelFR <- PM_Germline_Variants_Type_06$Indel/PM_Germline_Variants_Type_06$total
PM_Germline_Variants_Type_06$Sig3 <- c(0.0)
PM_Germline_Variants_Type_06$Sig8 <- c(0.0)
PM_Germline_Variants_Type_06$HRD <- PM_Germline_Variants_Type_06$Sig3 + PM_Germline_Variants_Type_06$Sig8
View(PM_Germline_Variants_Type_06)
write.csv(PM_Germline_Variants_Type_06, file = "PM_Germline_Variants_Type_06.csv")
PM_Germline_Variants_Type_06$Tumor_Sample_Barcode <- c("MB_REC_06")

#lenght of deletions
PM_Germline_DF_06 <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/MB-REC-06_primary_parental_germline_filtered.maf.txt", header=TRUE)
View(PM_Germline_DF_06)
PM_Germline_DF_06$Tumor_Sample_Barcode <- c("MB-REC-06")
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
PM_Germline_DF_MB_06 <- subset(PM_Germline_DF_06, PM_Germline_DF_06$Tumor_Sample_Barcode == "MB-REC-06")
PM_Germline_DF_MB_06_del <- subset(PM_Germline_DF_MB_06, PM_Germline_DF_MB_06$Variant_Type == "DEL")
View(PM_Germline_DF_MB_06_del)
PM_Germline_DF_MB_06_del$Del_Len <- PM_Germline_DF_MB_06_del$End_Position - PM_Germline_DF_MB_06_del$Start_Position
PM_Germline_DF_MB_06_del$LargeDel = ifelse(PM_Germline_DF_MB_06_del$Del_Len>=5,TRUE,FALSE)
PM_Germline_DF_MB_06_del_large = PM_Germline_DF_MB_06_del[PM_Germline_DF_MB_06_del$LargeDel==TRUE,]
View(PM_Germline_DF_MB_06_del_large)
write.csv(PM_Germline_DF_MB_06_del_large, file = "PM_Germline_DF_MB_06_del_large.csv")
write.csv(PM_Germline_DF_MB_06_del, file = "PM_Germline_DF_MB_06_del.csv")


##The recurrent(Metastatic site) samples with matched primary and germline 
RM_Germline_met <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/ICGC_germline_metastatic_reduced_maf_final.txt")

RM_Germline_met_Variants <- mafSummary(RM_Germline_met)
RM_Germline_met_Variants_Type <- RM_Germline_met_Variants$variant.type.summary

View(RM_Germline_met_Variants_Type)
RM_Germline_met_Variants_Type <- RM_Germline_met_Variants_Type[-c(1),]

RM_Germline_met_Variants_Type$Name <- c("Tumor4-R", "Tumor8-R","Tumor9-R","Tumor16","Tumor17-R")
RM_Germline_met_Variants_Type$Indel <- RM_Germline_met_Variants_Type$DEL + RM_Germline_met_Variants_Type$INS
RM_Germline_met_Variants_Type$IndelFR <- RM_Germline_met_Variants_Type$Indel/RM_Germline_met_Variants_Type$total
RM_Germline_met_Variants_Type$Sig3 <- c(0.076, 0.226, 0.224, 0.0, 0.0 )
RM_Germline_met_Variants_Type$Sig8 <- c(0.0284, 0.0, 0.355, 0.0, 0.15)
RM_Germline_met_Variants_Type$HRD <- RM_Germline_met_Variants_Type$Sig3 + RM_Germline_met_Variants_Type$Sig8
View(RM_Germline_met_Variants_Type)
write.csv(RM_Germline_met_Variants_Type, file = "RM_Germline_met_Variants_Type.csv")

#lenght of deletions
RM_Germline_met_DF <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_MAF_FILES/ICGC_germline_metastatic_reduced_maf_final.txt", 
                             header=TRUE, comment.char="#")
View(RM_Germline_met_DF)
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
RM_Germline_met_DF_MB_XX <- subset(RM_Germline_met_DF, RM_Germline_met_DF$Tumor_Sample_Barcode == "MB-REC-XX")
RM_Germline_met_DF_MB_XX_del <- subset(RM_Germline_met_DF_MB_XX , RM_Germline_met_DF_MB_XX$Variant_Type == "DEL")
View(RM_Germline_met_DF_MB_XX)
RM_Germline_met_DF_MB_XX_del$Del_Len <- RM_Germline_met_DF_MB_XX_del$End_Position - RM_Germline_met_DF_MB_XX_del$Start_Position
#View(RM_Germline_met_DF_DF_MB_XX_del)
RM_Germline_met_DF_MB_XX_del$LargeDel = ifelse(RM_Germline_met_DF_MB_XX_del$Del_Len>=5,TRUE,FALSE)
RM_Germline_met_DF_MB_XX_del_large = RM_Germline_met_DF_MB_XX_del[RM_Germline_met_DF_MB_XX_del$LargeDel==TRUE,]
View(RM_Germline_met_DF_MB_XX_del_large)
write.csv(RM_Germline_met_DF_MB_XX_del_large, file = "RM_Germline_met_DF_MB_XX_del_large.csv")
write.csv(RM_Germline_met_DF_MB_XX_del, file = "RM_Germline_met_DF_MB_XX_del.csv")


##The recurrent(primary&Metasstatic site) left over samples 
## From very first MAF file , a lot of samples were recalled and therefore those were removed, total 24 cases
RM <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/All_MAF_FILES/Primary_recurrent_reduced_maf_final.txt")

RM_Variants <- mafSummary(RM)
RM_Variants_Type <- RM_Variants$variant.type.summary
RM_Variants_Type <- RM_Variants_Type[-c(3,4,7,9,10,16,25,26,27),]
RM_Variants_Type <- RM_Variants_Type[-c(11),]

##View(RM_Variants_Type)

RM_Variants_Type$Name <- c("Tumor10", "Tumor1", "Tumor24","Tumor25","Tumor5",
                           "Tumor6","Tumor26","Tumor27","Tumor7","Tumor28", "Tumor30","Tumor31", "Tumor32", "Tumor33",
                           "Tumor34", "Tumor35", "Tumor36", "Tumor37", "Tumor38", "Tumor39", "Tumor40", "Tumor41","Tumor42", "Tumor43")
RM_Variants_Type$Indel <- RM_Variants_Type$DEL + RM_Variants_Type$INS
RM_Variants_Type$IndelFR <- RM_Variants_Type$Indel/RM_Variants_Type$total
RM_Variants_Type$Sig3 <- c(0.0, 0.0, 0.0, 0.0179, 0.21, 0.208, 0.327, 0.107, 0.14, 0.24, 0.352, 0.207, 0.157, 0.0,0.283,0.0, 0.224, 0.0, 0.084, 0.0, 0.0,0.0, 0.275, 0.0)
RM_Variants_Type$Sig8 <- c(0.0, 0.0, 0.0, 0.074, 0.182, 0.232, 0.254, 0.094, 0.175, 0.232, 0.0, 0.0, 0.147, 0.0, 0.0, 0.0, 0.16, 0.253, 0.154, 0.0, 0.0, 0.181, 0.0, 0.149) 
RM_Variants_Type$HRD <- RM_Variants_Type$Sig3 + RM_Variants_Type$Sig8
View(RM_Variants_Type)
write.csv(RM_Variants_Type, file = "RM_Variants_Type.csv")

#lenght of deletions
RM_Variants_DF <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/All_MAF_FILES/Primary_recurrent_reduced_maf_final.txt", 
                             header=TRUE, comment.char="#")
#View(RM_Variants_DF)
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
RM_Variants_DF_MB_XX <- subset(RM_Variants_DF, RM_Variants_DF$Tumor_Sample_Barcode == "MB-REC-XX")
RM_Variants_DF_MB_XX_del <- subset(RM_Variants_DF_MB_XX, RM_Variants_DF_MB_XX$Variant_Type == "DEL")
View(RM_Variants_DF_MB_XX_del)
RM_Variants_DF_MB_XX_del$Del_Len <- RM_Variants_DF_MB_XX_del$End_Position - RM_Variants_DF_MB_XX_del$Start_Position
#View(RM_Variants_DF_MB_01_del)

RM_Variants_DF_MB_XX_del$LargeDel = ifelse(RM_Variants_DF_MB_XX_del$Del_Len>=5,TRUE,FALSE)
RM_Variants_DF_MB_XX_del_large = RM_Variants_DF_MB_XX_del[RM_Variants_DF_MB_XX_del$LargeDel==TRUE,]
View(RM_Variants_DF_MB_XX_del_large)
write.csv(RM_Variants_DF_MB_XX_del_large, file = "RM_Variants_DF_MB_XX_del_large.csv")
write.csv(RM_Variants_DF_MB_XX_del, file = "RM_Variants_DF_MB_XX_del.csv")

##Combining all the data frames generated to create one DF, with large deletions, and HRD signatures. 
##Cleaning the data for final analysis
PM_Germline_Variants_Type_06 <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Thesis_Analysis_Figures/Deletions_ICGC_All_Samples/PM_Germline_Variants_Type_06.csv", header=TRUE)
##View(PM_Germline_Variants_Type_06)
PM_Germline_Variants_Type_06 <- PM_Germline_Variants_Type_06[,-1]
PM_Germline_Variants_Type_06$Type <- c("Primary")
PM_Germline_Variants_Type_06$largeDel <- c(3301)
M = PM_Germline_Variants_Type_06
PM_Germline_Variants_Type_16 <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Thesis_Analysis_Figures/Deletions_ICGC_All_Samples/PM_Germline_Variants_Type_16.csv")
##View(PM_Germline_Variants_Type_16)
PM_Germline_Variants_Type_16 <- PM_Germline_Variants_Type_16[,-1]
PM_Germline_Variants_Type_16$Type <- c("Primary")
PM_Germline_Variants_Type_16$largeDel <- c(3620)
N = PM_Germline_Variants_Type_16

PM_Germline_Variants_Type <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Thesis_Analysis_Figures/Deletions_ICGC_All_Samples/PM_Germline_Variants_Type.csv")
##View(PM_Germline_Variants_Type)
PM_Germline_Variants_Type  <- PM_Germline_Variants_Type[,-1]
PM_Germline_Variants_Type$Type  <- c("Primary")
PM_Germline_Variants_Type$largeDel <- c(84, 38,34, 55, 93, 32, 19, 52, 31, 43, 28,33)
O = PM_Germline_Variants_Type

RM_Germline_met_Variants_Type <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Thesis_Analysis_Figures/Deletions_ICGC_All_Samples/RM_Germline_met_Variants_Type.csv")
##View(RM_Germline_met_Variants_Type)
RM_Germline_met_Variants_Type <-RM_Germline_met_Variants_Type[,-1]
RM_Germline_met_Variants_Type$Type <- c("Recurrent")
RM_Germline_met_Variants_Type$LargeDel <- c(134, 130, 55, 34, 47)
 P = RM_Germline_met_Variants_Type


RM_Germline_Variants_Type_06 <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Thesis_Analysis_Figures/Deletions_ICGC_All_Samples/RM_Germline_Variants_Type_06.csv")
RM_Germline_Variants_Type_06 <- RM_Germline_Variants_Type_06[,-1]
RM_Germline_Variants_Type_06$Type <- c("Recurrent")
RM_Germline_Variants_Type_06$LargeDel <- c(555)

Q = RM_Germline_Variants_Type_06

RM_Germline_Variants_Type <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Thesis_Analysis_Figures/Deletions_ICGC_All_Samples/RM_Germline_Variants_Type.csv")

RM_Germline_Variants_Type$Sig3 <- c(0.208, 0.176, 0.130, 0.074, 0.264, 0.101, 0.098, 0.198, 0.066)
RM_Germline_Variants_Type$Sig8 <- c(0.371, 0.334, 0.265, 0.170, 0.000, 0.184, 0.304, 0.162, 0.151)
RM_Germline_Variants_Type$HRD <- RM_Germline_Variants_Type$Sig3 + RM_Germline_Variants_Type$Sig8
RM_Germline_Variants_Type$Type <- c("Recurrent")
RM_Germline_Variants_Type$LargeDel <- c(466, 237, 332, 252, 207, 193, 44, 139, 42)
RM_Germline_Variants_Type <- RM_Germline_Variants_Type[,-1]
 R = RM_Germline_Variants_Type

RM_Variants_Type <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Thesis_Analysis_Figures/Deletions_ICGC_All_Samples/RM_Variants_Type.csv")
RM_Variants_Type
RM_Variants_Type <- RM_Variants_Type[,-1]
RM_Variants_Type$Type <- c("Recurrent")
RM_Variants_Type$LargeDel <- c(171, 71, 99, 265, 225, 295, 138, 32, 142, 37, 112, 57, 175, 149, 28, 99, 171, 91, 38, 57, 50, 31, 53, 34)

 S = RM_Variants_Type
 
 

##joining all the 7 data frames to create one DF 
Final_HRD_DEl <- rbind.data.frame(PM_Germline_Variants_Type_06,PM_Germline_Variants_Type_16,PM_Germline_Variants_Type, RM_Germline_met_Variants_Type, RM_Germline_Variants_Type_06, RM_Germline_Variants_Type, RM_Variants_Type)
Y = rbind(M,N,O)
X = rbind(P,Q,R,S)
names(Y)
[1] "Tumor_Sample_Barcode" "DEL"                  "INS"                  "SNP"                 
[5] "total"                "Name"                 "Indel"                "IndelFR"             
[9] "Sig3"                 "Sig8"                 "HRD"                  "Type"                
[13] "largeDel"   

names(X)
[1] "Tumor_Sample_Barcode" "DEL"                  "INS"                  "SNP"                 
[5] "total"                "Name"                 "Indel"                "IndelFR"             
[9] "Sig3"                 "Sig8"                 "HRD"                  "Type"                
[13] "LargeDel" 

pylr::
Z = rename(Y, c("largeDel" = "LargeDel"))
Final_DF_HRD_DEl = rbind(X,Z)
write.csv(Final_DF_HRD_DEl, file = "Final_DF_HRD_DEl.csv")  
  ## Correlation relationship 


cor.test(Final_DF_HRD_DEl$LargeDel, Final_DF_HRD_DEl$HRD)
Final_DF_HRD_DEl$DelRate = Final_DF_HRD_DEl$LargeDel/Final_DF_HRD_DEl$DEL
cor.test(Final_DF_HRD_DEl$DelRate, Final_DF_HRD_DEl$HRD)
cor.test(Final_DF_HRD_DEl$DelRate, Final_DF_HRD_DEl$Sig3)
Final_DF_HRD_DEl$DelRate = Final_DF_HRD_DEl$DEL/Final_DF_HRD_DEl$INS
cor.test(Final_DF_HRD_DEl$DelRate, Final_DF_HRD_DEl$HRD)
cor.test(Final_DF_HRD_DEl$DelRate, Final_DF_HRD_DEl$Sig3)
cor.test(Final_DF_HRD_DEl$LargeDel, Final_DF_HRD_DEl$HRD)
