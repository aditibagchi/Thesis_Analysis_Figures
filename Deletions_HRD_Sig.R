##Main Objective is to generte adata frame with Insertions, Deletions, HRD signature and alos large deletions. 
## Requires merging of multiple MAF files

#maftools::
##The recurrent samples with matched primary and germline 
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



