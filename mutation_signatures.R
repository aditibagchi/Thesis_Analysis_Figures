ggplot2::
  ggthemes::
  scales::
  ggrepel::
  RColorBrewer::
  grid::
  gridBase::
  gridExtra::
  
##Mutation signatures 
  
  cols <- c("11" = "grey5", "18" = "gray46" , "3" = "deeppink4", "4" = "slategray1", "5" = "slateblue1", "6" = "navy", "8" = "lightpink3", "9" = "snow2", "U" = "grey70", "1" = "mediumseagreen", "12" = "cyan", "16" = "blue", "19" = "darkslateblue") 
  
  MB.REC.09 <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-09.csv")
M <- ggplot(data= MB.REC.09, aes(x=MB.REC.09$Type, y=MB.REC.09$Fraction, fill=MB.REC.09$Signature)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Subsitutions") + theme_minimal()+
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor")
MB.REC.07 <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-07.csv")
N <- ggplot(data= MB.REC.07, aes(x=MB.REC.07$Type, y=MB.REC.07$Fraction, fill=MB.REC.07$Signature)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Subsitutions") + theme_minimal()+
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor")
MB.REC.X <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/MB-REC-X.csv")
X <- ggplot(data= MB.REC.X, aes(x=MB.REC.X$Type, y=MB.REC.X$Fraction, fill= Signatures)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Contribution") + theme_minimal()+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor")
    