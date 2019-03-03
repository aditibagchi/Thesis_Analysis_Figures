##Main Objective is to generte adata frame with Insertions, Deletions, HRD signature and alos large deletions. 
## Requires merging of multiple MAF files

#maftools::
##The recurrent(primary site) samples with matched primary and germline 
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
View(RM_Germline_Variants_Type)
write.csv(RM_Germline_Variants_Type, file = "RM_Germline_Variants_Type.csv")

#lenght of deletions
RM_Germline_DF <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/9_recurrent_germline_samples_GRCH38_crossmapped_to_GRCH37.maf.txt", 
                             header=TRUE, comment.char="#")
#View(RM_Germline_DF)
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
RM_Germline_DF_MB_15 <- subset(RM_Germline_DF, RM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-15")
RM_Germline_DF_MB_15_del <- subset(RM_Germline_DF_MB_15, RM_Germline_DF_MB_15$Variant_Type == "DEL")
#View(RM_Germline_DF_MB_15_del)
RM_Germline_DF_MB_15_del$Del_Len <- RM_Germline_DF_MB_15_del$End_Position - RM_Germline_DF_MB_15_del$Start_Position
#View(RM_Germline_DF_MB_15_del)
median(RM_Germline_DF_MB_15_del$Del_Len)
mean(RM_Germline_DF_MB_15_del$Del_Len)
RM_Germline_DF_MB_15_del$LargeDel = ifelse(RM_Germline_DF_MB_15_del$Del_Len>=5,TRUE,FALSE)
RM_Germline_DF_MB_15_del_large = RM_Germline_DF_MB_15_del[RM_Germline_DF_MB_15_del$LargeDel==TRUE,]
#View(RM_Germline_DF_MB_15_del_large)
write.csv(RM_Germline_DF_MB_15_del_large, file = "RM_Germline_DF_MB_15_del_large.csv")
write.csv(RM_Germline_DF_MB_15_del, file = "RM_Germline_DF_MB_15_del.csv")






##The primary samples with matched recurrent and germline 
##MB-REC-16 and MB-REC-06 was removed , as these two samples were reanalyzed.
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
PM_Germline_DF_MB_14 <- subset(PM_Germline_DF, PM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-14")
PM_Germline_DF_MB_14_del <- subset(PM_Germline_DF_MB_14, PM_Germline_DF_MB_14$Variant_Type == "DEL")
View(PM_Germline_DF_MB_14_del)
PM_Germline_DF_MB_14_del$Del_Len <- PM_Germline_DF_MB_14_del$End_Position - PM_Germline_DF_MB_14_del$Start_Position
PM_Germline_DF_MB_14_del$LargeDel = ifelse(PM_Germline_DF_MB_14_del$Del_Len>=5,TRUE,FALSE)
PM_Germline_DF_MB_14_del_large = PM_Germline_DF_MB_14_del[PM_Germline_DF_MB_14_del$LargeDel==TRUE,]
View(PM_Germline_DF_MB_14_del_large)
write.csv(PM_Germline_DF_MB_14_del_large, file = "PM_Germline_DF_MB_14_del_large.csv")
write.csv(PM_Germline_DF_MB_14_del, file = "PM_Germline_DF_MB_14_del.csv")

##MB_REC_16 primary 
PM_Germline_16 <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-16_primary.crossmap.hg38_to_hg19.vep.maf.txt")
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

MB.REC.16_primary.crossmap.hg38_to_hg19.vep.maf <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-16_primary.crossmap.hg38_to_hg19.vep.maf.txt", header=FALSE, comment.char="#")
#lenght of deletions
PM_Germline_DF_16 <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-16_primary.crossmap.hg38_to_hg19.vep.maf.txt", header=TRUE, comment.char="#")
View(PM_Germline_DF_16)
PM_Germline_DF_16$Tumor_Sample_Barcode <- c("MB_REC_16")
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
PM_Germline_DF_MB_16 <- subset(PM_Germline_DF, PM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-16")
PM_Germline_DF_MB_16_del <- subset(PM_Germline_DF_MB_16, PM_Germline_DF_MB_16$Variant_Type == "DEL")
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
RM_Germline_DF_MB_06 <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-06_recurrent_parental_germline_filtered.maf.txt", header=TRUE)
View(RM_Germline_DF_06)
RM_Germline_DF_MB_06$Tumor_Sample_Barcode <- c("MB_REC_06")
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
PM_Germline_DF_06 <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-06_primary_parental_germline_filtered.maf.txt", header=TRUE)
View(PM_Germline_DF_06)
PM_Germline_DF_06$Tumor_Sample_Barcode <- c("MB_REC_06")
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
PM_Germline_DF_MB_06 <- subset(PM_Germline_DF, PM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-06")
PM_Germline_DF_MB_06_del <- subset(PM_Germline_DF_MB_06, PM_Germline_DF_MB_06$Variant_Type == "DEL")
View(PM_Germline_DF_MB_06_del)
PM_Germline_DF_MB_06_del$Del_Len <- PM_Germline_DF_MB_06_del$End_Position - PM_Germline_DF_MB_06_del$Start_Position
PM_Germline_DF_MB_06_del$LargeDel = ifelse(PM_Germline_DF_MB_06_del$Del_Len>=5,TRUE,FALSE)
PM_Germline_DF_MB_06_del_large = PM_Germline_DF_MB_06_del[PM_Germline_DF_MB_06_del$LargeDel==TRUE,]
View(PM_Germline_DF_MB_06_del_large)
write.csv(PM_Germline_DF_MB_06_del_large, file = "PM_Germline_DF_MB_06_del_large.csv")
write.csv(PM_Germline_DF_MB_06_del, file = "PM_Germline_DF_MB_06_del.csv")

##The recurrent(primary site) samples with matched primary and germline 
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
View(RM_Germline_Variants_Type)
write.csv(RM_Germline_Variants_Type, file = "RM_Germline_Variants_Type.csv")

#lenght of deletions
RM_Germline_DF <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/9_recurrent_germline_samples_GRCH38_crossmapped_to_GRCH37.maf.txt", 
                             header=TRUE, comment.char="#")
#View(RM_Germline_DF)
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
RM_Germline_DF_MB_15 <- subset(RM_Germline_DF, RM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-15")
RM_Germline_DF_MB_15_del <- subset(RM_Germline_DF_MB_15, RM_Germline_DF_MB_15$Variant_Type == "DEL")
#View(RM_Germline_DF_MB_15_del)
RM_Germline_DF_MB_15_del$Del_Len <- RM_Germline_DF_MB_15_del$End_Position - RM_Germline_DF_MB_15_del$Start_Position
#View(RM_Germline_DF_MB_15_del)
median(RM_Germline_DF_MB_15_del$Del_Len)
mean(RM_Germline_DF_MB_15_del$Del_Len)
RM_Germline_DF_MB_15_del$LargeDel = ifelse(RM_Germline_DF_MB_15_del$Del_Len>=5,TRUE,FALSE)
RM_Germline_DF_MB_15_del_large = RM_Germline_DF_MB_15_del[RM_Germline_DF_MB_15_del$LargeDel==TRUE,]
#View(RM_Germline_DF_MB_15_del_large)
write.csv(RM_Germline_DF_MB_15_del_large, file = "RM_Germline_DF_MB_15_del_large.csv")
write.csv(RM_Germline_DF_MB_15_del, file = "RM_Germline_DF_MB_15_del.csv")


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
RM_Germline_met_DF_MB_03 <- subset(RM_Germline_met_DF, RM_Germline_met_DF$Tumor_Sample_Barcode == "MB-REC-03")
RM_Germline_met_DF_MB_03_del <- subset(RM_Germline_met_DF_MB_03 , RM_Germline_met_DF_MB_03$Variant_Type == "DEL")
View(RM_Germline_met_DF_MB_03)
RM_Germline_met_DF_MB_03_del$Del_Len <- RM_Germline_met_DF_MB_03_del$End_Position - RM_Germline_met_DF_MB_03_del$Start_Position
#View(RM_Germline_met_DF_DF_MB_02_del)
RM_Germline_met_DF_MB_03_del$LargeDel = ifelse(RM_Germline_met_DF_MB_03_del$Del_Len>=5,TRUE,FALSE)
RM_Germline_met_DF_MB_03_del_large = RM_Germline_met_DF_MB_03_del[RM_Germline_met_DF_MB_03_del$LargeDel==TRUE,]
View(RM_Germline_met_DF_MB_03_del_large)
write.csv(RM_Germline_met_DF_MB_03_del_large, file = "RM_Germline_met_DF_MB_03_del_large.csv")
write.csv(RM_Germline_met_DF_MB_03_del, file = "RM_Germline_met_DF_MB_03_del.csv")


##The recurrent(primary&Metasstaticsite) left over samples 
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
RM_Variants_DF_MB_29 <- subset(RM_Variants_DF, RM_Variants_DF$Tumor_Sample_Barcode == "MB-REC-29")
RM_Variants_DF_MB_29_del <- subset(RM_Variants_DF_MB_29, RM_Variants_DF_MB_29$Variant_Type == "DEL")
View(RM_Variants_DF_MB_29_del)
RM_Variants_DF_MB_29_del$Del_Len <- RM_Variants_DF_MB_29_del$End_Position - RM_Variants_DF_MB_29_del$Start_Position
#View(RM_Variants_DF_MB_01_del)

RM_Variants_DF_MB_29_del$LargeDel = ifelse(RM_Variants_DF_MB_29_del$Del_Len>=5,TRUE,FALSE)
RM_Variants_DF_MB_29_del_large = RM_Variants_DF_MB_29_del[RM_Variants_DF_MB_29_del$LargeDel==TRUE,]
View(RM_Variants_DF_MB_29_del_large)
write.csv(RM_Variants_DF_MB_29_del_large, file = "RM_Variants_DF_MB_29_del_large.csv")
write.csv(RM_Variants_DF_MB_29_del, file = "RM_Variants_DF_MB_29_del.csv")




##The primary samples with matched recurrent and germline 
##MB-REC-16 and MB-REC-06 was removed , as these two samples were reanalyzed.
PM_Germline <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/primary_germline_15_files.maf.txt")
PM_Germline_Variants <- mafSummary(PM_Germline)
PM_Germline_Variants_Type <- PM_Germline_Variants$variant.type.summary
##View(PM_Germline_Variants_Type)
PM_Germline_Variants_Type <- PM_Germline_Variants_Type[-c(1,2,4),]
PM_Germline_Variants_Type$Name <- c("Tumor20","Tumor23","Tumor4","Tumor9","Tumor22","Tumor17", "Tumor19", "Tumor21", "Tumor18", "Tumor8")
PM_Germline_Variants_Type$Indel <- PM_Germline_Variants_Type$DEL + PM_Germline_Variants_Type$INS
PM_Germline_Variants_Type$IndelFR <- PM_Germline_Variants_Type$Indel/PM_Germline_Variants_Type$total
PM_Germline_Variants_Type$Sig3 <- c( 0.101, 0.0, 0.164, 0.342, 0.0, 0.0, 0.116, 0.0, 0.201, 0.094)
PM_Germline_Variants_Type$Sig8 <- c(0.119, 0.328, 0.0, 0.0, 0.218, 0.15, 0.0, 0.293, 0.139, 0.0) 
PM_Germline_Variants_Type$HRD <- PM_Germline_Variants_Type$Sig3 + PM_Germline_Variants_Type$Sig8
View(PM_Germline_Variants_Type)
write.csv(PM_Germline_Variants_Type, file = "PM_Germline_Variants_Type.csv")


#lenght of deletions
PM_Germline_DF <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/primary_germline_15_files.maf.txt", 
                             header=TRUE, comment.char="#")
View(PM_Germline_DF)
#To identify the lenght of deletions; substarct endposition from start postion; first subset one tumor and then substet all deletions within the tumor.
PM_Germline_DF_MB_10 <- subset(PM_Germline_DF, PM_Germline_DF$Tumor_Sample_Barcode == "MB-REC-10")
PM_Germline_DF_MB_10_del <- subset(PM_Germline_DF_MB_10, PM_Germline_DF_MB_10$Variant_Type == "DEL")
View(PM_Germline_DF_MB_10_del)
PM_Germline_DF_MB_10_del$Del_Len <- PM_Germline_DF_MB_10_del$End_Position - PM_Germline_DF_MB_10_del$Start_Position
PM_Germline_DF_MB_10_del$LargeDel = ifelse(PM_Germline_DF_MB_10_del$Del_Len>=5,TRUE,FALSE)
PM_Germline_DF_MB_10_del_large = PM_Germline_DF_MB_10_del[PM_Germline_DF_MB_10_del$LargeDel==TRUE,]
View(PM_Germline_DF_MB_10_del_large)
write.csv(PM_Germline_DF_MB_10_del_large, file = "PM_Germline_DF_MB_10_del_large.csv")
write.csv(PM_Germline_DF_MB_10_del, file = "PM_Germline_DF_MB_10_del.csv")





