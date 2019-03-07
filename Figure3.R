ggplot2::
  ggthemes::
  scales::
  ggrepel::
  RColorBrewer::
  grid::
  gridBase::
  gridExtra::
  maftools::
  vcfR::
  ##Mutation signatures  for hypermutated cases 
 cols <- c("11" = "grey5", "18" = "gray46" , "3" = "deeppink4", "4" = "slategray1", "5" = "slateblue1", "6" = "navy", "8" = "lightpink3", "9" = "snow2", "U" = "grey70", "1" = "mediumseagreen", "12" = "cyan", "16" = "blue", "19" = "darkslateblue", "20"= "orange1", "10" = "gold", "15" = "darkseagreen4", "14" = "firebrick", "21" = "darkorchid1") 

##Plot2 with only recurrent cases ( MB-REC-44, MB-REc-26, MB-REC-16, MB-REC-06)
SignatureR <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Signature_Hyp.csv")
D <- ggplot(data= SignatureR, aes(x=SignatureR$Name, y=SignatureR$Fraction, fill= Signatures)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Contribution") + theme_minimal()+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor") + ggtitle("C")

## Lollipop plot for POLE (MB-REC-44, MB-REC-26)

mafX <- read.maf(maf = "/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/All_MAF_FILES/Primary_recurrent_reduced_maf_final.txt")
maf44 <- subsetMaf(mafX, tsb = "MB-REC-44")
maf44$Tumor_Sample_Barcode = c("Tumor10")
maf44 = maftools::read.maf(maf = maf44)
maf26 <- subsetMaf(mafX, tsb = "MB-REC-26")
maf26$Tumor_Sample_Barcode = c("Tumor1")
maf26 = maftools::read.maf(maf = maf26)


Y = lollipopPlot(maf26, gene = "POLE", labelPos = "all", showMutationRate = FALSE, labPosSize = 4, showDomainLabel = FALSE, pointSize = 3) 
X = lollipopPlot(maf44, gene = "POLE", labelPos = "all", showMutationRate = FALSE, labPosSize = 4, showDomainLabel = FALSE, pointSize = 3) 

M = rainfallPlot(maf26, pointSize = 0.25, fontSize = 10)
N = rainfallPlot(maf44, pointSize = 0.25, fontSize = 10)

##Density plots
#Mutationallelefrequency(MAF)
#Calculating MAF of POLE and also comparing it with other mutations.

#Method obtain the mutation allel frequency from VCF
# Specify your file:
MB_REC_44 ='/Volumes/projects_secondary/jewell/EGAD00001000946_4/primary_recurrent/raw_somatic_vcfs/MB-REC-44_tumor.sv.vcf'
# Load the VCF file:
vcf44 = read.vcfR(MB_REC_44)
# Load the function to get the chr, pos & depth:
get_tumor_ref_alt_depth=function(vcf44){
  # x=as.data.frame(x@gt)
  # list=strsplit(x=x$TUMOR,split = ':',fixed = TRUE)   
  m=as.data.frame(extract.gt(vcf44,element='AD'))
  Mx <-as.character(m$TUMOR)
  list=strsplit(Mx,',',fixed=TRUE)
  
  # m=matrix(unlist(list),ncol=10,byrow=TRUE) # does not work if lists not all same length
  d = as.data.frame(matrix(unlist(list), ncol=2,byrow = TRUE))
  d = as.data.frame(matrix(unlist(list), ncol=2,byrow = TRUE), stringsAsFactors = FALSE)
  colnames(d) = c('t_ref_count','t_alt_count')
  # colnames(m)=unlist(strsplit(x$FORMAT[1],':',TRUE))
  
  # Now get the chrom & pos
  chrpos=as.data.frame(vcf44@fix, stringsAsFactors = FALSE)
  d$CHROM=chrpos$CHROM
  d$POS=chrpos$POS
  return(d)
}


#  Create to dataframe with total depth and also "mutant allele frequency); 
# Convert d to a matrix first to convert the rows to numeric
MB44 <- data.matrix(d, rownames.force = NA)
#Then reconvert it to dataframe
MB44 <- as.data.frame(MB44, stringsAsFactors = FALSE)
#Then add columns by adding ref_depth and allel depth that will give the total depth for the varaint
MB44$Total_depth <- MB44$t_ref_count + MB44$t_alt_count
#Then add column by dividing allel_depth by total_depth to get Mutant allele frquency
MB44$MAF <- MB44$t_alt_count / MB44$Total_depth
#subsetting for total depth more than 10
MB44_Subset_Depth10X <- MB44[ which(MB44$Total_depth >10),]
#Subsetting for total MAF of more than 0.1
MB44_Subset_Depth10X_MAF <- MB44_Subset_Depth10X[ which(MB44_Subset_Depth10X$MAF > 0.49),]

# Data frame with al mutaions with  more MAF than 0.17
MB44SUB <- MB26_Subset_Depth10X[ which(MB44_Subset_Depth10X$MAF > 0.49),]
MB44SUB$NAME <- c("MAF>0.49")
#Data frame witk all mutations with less MAF than 0.17
MB44SUB1 <- MB44_Subset_Depth10X_MAF[ which(MB44_Subset_Depth10X_MAF$MAF < 0.49),]
MB26SUB1$NAME <- c("MAF<0.49")
#Now join the rows of MB26SUB and MB44SUB1 to create one Data frame
FinalMB44 <- rbind(MB44SUB, MB44SUB1)

##
MB_REC_26="/Volumes/projects_secondary/jewell/EGAD00001000946_4/primary_recurrent/raw_somatic_vcfs/MB-REC-26_tumor.sv.vcf"
MB_REC_26_Annotated= "/Volumes/projects_secondary/jewell/EGAD00001000946_4/primary_recurrent/annotated_vcfs/MB-REC-26_tumor_PASS.remap.sv.vep.vcf"
# Load the VCF file:
vcf26An = read.vcfR(MB_REC_26_Annotated)

# Load the function to get the chr, pos & depth:
get_tumor_ref_alt_depth=function(vcf26An){
  # x=as.data.frame(x@gt)
  # list=strsplit(x=x$TUMOR,split = ':',fixed = TRUE)   
  m=as.data.frame(extract.gt(vcf26An,element='AD'))
  Mx <-  as.character(m$TUMOR)
  list=strsplit(Mx,',',fixed=TRUE)
  
  # m=matrix(unlist(list),ncol=10,byrow=TRUE) # does not work if lists not all same length
  MB26 = as.data.frame(matrix(unlist(list), ncol=2,byrow = TRUE))
  MB26 = as.data.frame(matrix(unlist(list), ncol=2,byrow = TRUE), stringsAsFactors = FALSE)
  colnames(MB26) = c('t_ref_count','t_alt_count')
  # colnames(m)=unlist(strsplit(x$FORMAT[1],':',TRUE))
  
  # Now get the chrom & pos
  chrpos=as.data.frame(vcf26An@fix, stringsAsFactors = FALSE)
  MB26$CHROM=chrpos$CHROM
  MB26$POS=chrpos$POS
  
  
  return(MB26)
}

#  Create to dataframe with total depth and also "mutant allele frequency); 
# Convert d to a matrix first to convert the rows to numeric
MB26m <- data.matrix(MB26, rownames.force = NA)
#Then reconvert it to dataframe
MB26 <- as.data.frame(MB26m, stringsAsFactors = FALSE)
#Then add columns by adding ref_depth and allel depth that will give the total depth for the varaint
MB26$Total_depth <- MB26$t_ref_count + MB26$t_alt_count
#Then add column by dividing allel_depth by total_depth to get Mutant allele frquency
MB26$MAF <- MB26$t_alt_count / MB26$Total_depth
#subsetting for total depth more than 10
MB26_Subset_Depth10X <- MB26[ which(MB26$Total_depth >10),]
#Subsetting for total MAF of more than 0.17
MB26_Subset_Depth10X_MAF <- MB26_Subset_Depth10X[ which(MB26_Subset_Depth10X$MAF > 0.17),]
MB26_Subset_Depth10X_MAF$NAME <- c("MAF>0.17")
MB26_Subset_Depth10X_MAF_less <- MB26_Subset_Depth10X[ which(MB26_Subset_Depth10X$MAF < 0.17),]
MB26_Subset_Depth10X_MAF_less$NAME <- c("MAF<0.17")
FinalMB26An <- rbind(MB26_Subset_Depth10X_MAF_less, MB26_Subset_Depth10X_MAF)
mean(MB26_Subset_Depth10X_MAF$MAF)
[1] 0.2577228
median(MB26_Subset_Depth10X_MAF$MAF)
[1] 0.2307692





#Now to use ggplot2:: to create density plot ; The MAF for POLE mutation for position "132673703" is 0.
 B = ggplot(FinalMB44, aes(x=FinalMB44$MAF)) +
  geom_density() + theme_minimal() + geom_vline(aes(xintercept = median(FinalMB44$MAF)), linetype = "dashed") + 
  scale_fill_brewer ("blue")  + xlab("Mutation Allele Frequency")+
  ylab("Density") + theme(axis.title.x = element_text(face = "bold", color = "black", size = 12),axis.title.y = element_text(face = "bold", color = "black", size =12))

 
 C = ggplot(MB26, aes(x=MB26$MAF)) +
   geom_density() + theme_minimal() + geom_vline(aes(xintercept = median(MB26$MAF)), linetype = "dashed") + 
   scale_fill_brewer ("blue")  + xlab("Mutation Allele Frequency")+
   ylab("Density") + theme(axis.title.x = element_text(face = "bold", color = "black", size = 12),axis.title.y = element_text(face = "bold", color = "black", size =12))
  
 

##Arranging figures in one page
 grid.arrange(arrangeGrob(N,M, ncol = 2), arrangeGrob(X,Y,ncol = 2), arrangeGrob(B,C, ncol = 2),top = textGrob("Figure4", gp=gpar(fontsize=18, fon=3)))

 
