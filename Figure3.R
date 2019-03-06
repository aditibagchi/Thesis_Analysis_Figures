ggplot2::
  ggthemes::
  scales::
  ggrepel::
  RColorBrewer::
  grid::
  gridBase::
  gridExtra::
  maftools::
  
  ##Mutation signatures  for hypermutated cases 
 cols <- c("11" = "grey5", "18" = "gray46" , "3" = "deeppink4", "4" = "slategray1", "5" = "slateblue1", "6" = "navy", "8" = "lightpink3", "9" = "snow2", "U" = "grey70", "1" = "mediumseagreen", "12" = "cyan", "16" = "blue", "19" = "darkslateblue", "20"= "orange1", "10" = "gold", "15" = "darkseagreen4", "14" = "firebrick", "21" = "darkorchid1") 

##Plot2 with only recurrent cases ( MB-REC-44, MB-REc-26, MB-REC-16, MB-REC-06)
SignatureR <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Signature_Hyp.csv")
R <- ggplot(data= SignatureR, aes(x=SignatureR$Name, y=SignatureR$Fraction, fill= Signatures)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Contribution") + theme_minimal()+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor") + ggtitle("A")

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
Tools    
vcfR::
ggplot2::
  
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
  as.character(m$TUMOR)
  list=strsplit(m$TUMOR,',',fixed=TRUE)
  
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
# Get a df that has four columns:
# column 1 = tumor reference count
# column 2 = tumor alternate allele count
# column 3 = chromosome
# column 4 = position
df = get_tumor_ref_alt_depth(vcf44)

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
MB44_Subset_Depth10X_MAF <- MB44_Subset_Depth10X[ which(MB44_Subset_Depth10X$MAF > 0.2),]

# Data frame with al mutaions with  more MAF than 0.49
MB44SUB <- MB44_Subset_Depth10X[ which(MB44_Subset_Depth10X$MAF > 0.49),]
MB44SUB$NAME <- c("MAF>0.49")
#Data frame witk all mutations with less MAF than 0.49
MB44SUB1 <- MB44_Subset_Depth10X_MAF[ which(MB44_Subset_Depth10X_MAF$MAF < 0.49),]
MB44SUB1$NAME <- c("MAF<0.49")
#Now join the rows of MB44SUB and MB44SUB1 to create one Data frame
FinalMB44 <- rbind(MB44SUB, MB44SUB1)
#Now to use ggplot2:: to create density plot ; The MAF for POLE mutation for position "132673703" is 0.


##Arranging figures in one page
Figure2 <- grid.arrange(R, X, top = textGrob("Figure2", gp=gpar(fontsize=18, fon=3)))
Figure2
