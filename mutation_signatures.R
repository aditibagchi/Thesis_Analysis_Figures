ggplot2::
  ggthemes::
  scales::
  ggrepel::
  RColorBrewer::
  grid::
  gridBase::
  gridExtra::
  
##Mutation signatures 
##Plot1 with cases matched primary and recurrent tumors

cols <- c("11" = "grey5", "18" = "gray46" , "3" = "deeppink4", "4" = "slategray1", "5" = "slateblue1", "6" = "navy", "8" = "lightpink3", "9" = "snow2", "U" = "grey70", "1" = "mediumseagreen", "12" = "cyan", "16" = "blue", "19" = "darkslateblue", "20"= "orange1") 
 
Signature <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Signatures_Matched.csv")
X <- ggplot(data= Signature, aes(x=Signature$Type, y=Signature$Fraction, fill= Signatures)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Contribution") + theme_minimal()+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor")
 S <- X+facet_wrap(~Name) + ggtitle("A")

 ##Plot2 with only recurrent cases 
 SignatureR <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Signatures_recurrent.csv")
 R <- ggplot(data= SignatureR, aes(x=SignatureR$Name, y=SignatureR$Fraction, fill= Signatures)) +
   geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Contribution") + theme_minimal()+
   scale_fill_manual(values = cols)+
   theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor") + ggtitle("B")
 ##Arranging figures in one page
 Figure2 <- grid.arrange(S, R, top = textGrob("Figure2", gp=gpar(fontsize=18, fon=3)))
 Figure2
 