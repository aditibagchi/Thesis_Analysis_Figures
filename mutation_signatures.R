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

  
  cols <- c("11" = "grey5", "18" = "gray46" , "3" = "deeppink4", "4" = "slategray1", "5" = "slateblue1", "6" = "navy", "8" = "lightpink3", "9" = "snow2", "U" = "grey70", "1" = "mediumseagreen", "12" = "cyan", "16" = "blue", "19" = "darkslateblue") 
 
Signature <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Signatures_Matched.csv")
X <- ggplot(data= Signature, aes(x=Signature$Type, y=Signature$Fraction, fill= Signatures)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Contribution") + theme_minimal()+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 90)) + xlab("Tumor")
 S <- X+facet_wrap(~Tumor)

 