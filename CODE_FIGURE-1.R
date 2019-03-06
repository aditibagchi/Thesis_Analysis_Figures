##Tumor Mutational Burden in Primary and Recurrent Medulloblastoma; Calculation of Outliers
Medulloblastoma <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/ALL_SAMPLES_TMD.csv", comment.char="#")
str(Medulloblastoma)
head(Medulloblastoma)
hist(Medulloblastoma$Total.Mutation.Burden)
##The histogram is a skewed one the majority of the data lies between 1-2; to geta better judegement of the outliers
## Use boxplot function to get a real visual of outliers. 
boxplot(Medulloblastoma$Total.Mutation.Burden, ylab = 'Mutation per Mb')
## 9 outliers noted 
## When is a value an outlier: *expert judgement, *Rule of thumb: Q1-1.5 *IQR 
                                                                  Q3 + 1.5*IQR  [ IQR ~ Interquartile range]
                                                                  
 ## Calculate IQR for the tumor mutational burdden
IQR(Medulloblastoma$Total.Mutation.Burden) 
##1.046815
quantile(Medulloblastoma$Total.Mutation.Burden)
##0%        25%        50%        75%       100% 
##0.2102102(Q1)  1.1003317(Q2)  1.4029851(Q3)  2.1471471(Q4) 39.4594595 

##Outliers  are >  2.14(Q3) + 1.5*1.04(IQR)
2.14 + (1.5*1.04)
3.7
##Any tumor with tumor mutational burden of more than 3.7 can be considered as out liers.

Medulloblastoma1 <- Medulloblastoma[ which(Medulloblastoma$Total.Mutation.Burden < 3.7),]
Medulloblastoma2 <- Medulloblastoma[ which(Medulloblastoma$Total.Mutation.Burden > 3.7),]
Medulloblastoma2

boxplot(Medulloblastoma1$Total.Mutation.Burden, ylab = "Mutaion per MB")
median(Medulloblastoma1$Total.Mutation.Burden)
[1] 1.377258
mean(Medulloblastoma1$Total.Mutation.Burden)
[1] 1.450055
## Now make TMD plot
## All tumors with Mutational burden more than 3.7 mutation per Mb are considered hypermutated. A total of 9 samples. 
## Two graphs are created for tumor mutational burden.

ggplot2::
ggthemes::
scales::
ggrepel::
RColorBrewer::
grid::
gridExtra:: 
lattice::  
  cols <- c("C>A" = "blue4", "C>G" = "grey66", "C>T" = "green4", "T>A" = "goldenrod2", "T>C" = "gray33", "T>G" = "red3")

X <- ggplot(Medulloblastoma, aes(Medulloblastoma$Name, Medulloblastoma$Total.Mutation.Burden))
Y <- X + geom_bar(stat = "identity", aes(fill = Type))+ 
  scale_y_continuous(name="Mutations per Mb", limits=c(0, 40)) + theme_minimal() + xlab("Tumor")+
  scale_fill_brewer(palette = "Set1") + theme(axis.text.x = element_text(color = "black", size = 6, angle = 90)) + ggtitle("A")

Subsitutions <- read.csv("/Volumes/G-DRIVE mobile/Data_Analysis/Thesis/Thesis/Subsitutions_Final.csv")
View(Subsitutions)
M <- ggplot(data= Subsitutions, aes(x=Subsitutions$Name, y=Subsitutions$Fraction, fill=Subsitution)) +
  geom_bar(stat="identity") +  scale_y_continuous(name="Fraction Subsitutions") + theme_minimal()+
  scale_fill_manual(values = cols)+
  theme(axis.text.x = element_text(color = "black", size = 6, angle = 90)) + xlab("Tumor") + ggtitle("C")

#leaving empty space
CNV_Data <- textGrob("CNV_Data_figure1D")
Oncoplot <- textGrob("Oncoplot_figure1B")

##laying out multiple plots in each page
Figure1 <- grid.arrange(Y,Oncoplot,M,CNV_Data,top = textGrob("Figure1", gp=gpar(fontsize=18, fon=3)), ncol=1)
Figure1
