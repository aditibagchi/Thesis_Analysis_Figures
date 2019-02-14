Test <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/9_recurrent_germline_samples_GRCH38_crossmapped_to_GRCH37.maf.txt", header=TRUE, comment.char="#",stringsAsFactors = FALSE)
# View(Test)
TestX  <- read.delim("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/35_Exomes_pipeline_filtered.maf.txt", header=TRUE, comment.char="#",stringsAsFactors = FALSE)
# dim(Test)
# dim(TestX)
TestX = TestX[,which(colnames(TestX) %in% colnames(Test))]
# colnames(Test)==colnames(TestX)
# library(dplyr)

# m2 = dplyr::bind_rows(Test,TestX)
# str(Test$AF)
# str(TestX$AF)
m = rbind(Test,TestX)

# getwd()
merged.file='Test_Merged_MAF_File_To_Delete.maf'
write.table(m,file = merged.file,sep='\t',row.names=FALSE,quote=FALSE)
maf = maftools::read.maf(maf = merged.file)
str(maf)
class(maf)
