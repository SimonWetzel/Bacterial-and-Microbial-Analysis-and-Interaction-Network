#Package requirements
install.packages("corrplot")
library("corrplot")
library(Hmisc)    # Computes correlation and p-value matrix

#Set directory

path <- "/home/simonw/Schreibtisch/Sequencing_Data/Fungiome/Crohn_Kohorte"


setwd("/home/simonw/Schreibtisch/Sequencing_Data/Fungiome/Crohn_Kohorte/")

#Microbiome
mydata = read.csv("sample_info_tab_excluded_taxa.csv",header=TRUE, sep="\t", )
#Mycobiome
mydata = read.csv("sample_info_crohn_AA83_excluded.csv.csv",header=TRUE, sep="\t", )

#clinical Data corr
mydata = read.csv("sample_info_tab_corr.csv",header=TRUE, sep="\t", )


path <- "/home/simonw/Schreibtisch"

setwd("/home/simonw/Schreibtisch")

mydata = read.csv("clinical_data_english_clinicaldata_corr_only_sig_ordered.csv",header=TRUE, sep=",", )


mydata = read.csv("clinical_data_english_corr_bacteria_fungi.csv",header=TRUE, sep=",", )
head(mydata)



#Correlation analysis

M = cor(mydata[sapply(mydata, is.numeric)],method = "spearman", use = "complete.obs")
head(M)

#Visualization

corrplot(M, method="circle",type = 'lower', tl.col="black",tl.cex=0.5)


#p-Values
M2 <- rcorr(as.matrix(mydata[sapply(mydata, is.numeric)]),type= "spearman")
M2$r
M2$P

# ++++++++++++++++++++++++++++
# flattenCorrMatrix
# ++++++++++++++++++++++++++++
# cormat : matrix of the correlation coefficients
# pmat : matrix of the correlation p-values
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}


library(Hmisc)
results= flattenCorrMatrix(M2$r, M2$P)
print(flattenCorrMatrix(M2$r, M2$P))
p <- (flattenCorrMatrix(M2$r, M2$P)["p"])
p <- as.numeric(unlist(p))
print(p)


write.csv((flattenCorrMatrix(M2$r, M2$P)), "Crohn.csv")

flattenCorrMatrix(M2$r, M2$P)
padj <- p.adjust(p, method = "BH", n = length(p))
padj <- array(padj)
results$padj<- padj
print(results)

write.csv(results, "Crohn.csv")


#Plotting Correlation results

corrplot(M2$r[1:6, 8:45], p.mat=M2$P[1:6, 8:45], method="color", addCoef.col="black")
corrplot(M2$r[15:41,1:14], p.mat=M2$P[15:41, 1:14], method="color", sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig',number.cex = 0.4,  tl.cex = 0.9, cl.cex =1, pch.cex =1.5)

corrplot(M2$r[8:22,1:7], p.mat=M2$P[8:22, 1:6], method="color",  addCoef.col="black")

corrplot(M2$r[8:39,40:50], p.mat=M2$P[8:39,40:50], method="color", sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig',number.cex = 0.1,  tl.cex = 0.5, cl.cex = .1, pch.cex =1)
corrplot(M2$r[45:55,15:44], p.mat=M2$P[45:55,15:44], method="color", sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig',number.cex = 0.1,  tl.cex = 1, cl.cex = 1, pch.cex =2)
corrplot(M2$r[15:55,15:55], p.mat=M2$P[15:55,15:55], method="color", sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig',number.cex = 0.1,  tl.cex = 0.8, cl.cex = 1, pch.cex =2,type = 'upper')


corrplot(M2$r[15:55,15:55], p.mat=M2$P[15:55,15:55], method="color", sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig',number.cex = 0.1,  tl.cex = 0.8, cl.cex = 1, pch.cex =2, type = 'upper', )
corrplot(M2$r[8:51,8:51], p.mat=M2$P[8:51, 8:51], method="color", sig.level = c(0.001, 0.01, 0.05), insig = 'label_sig',number.cex = 0.1,  tl.cex = 0.8, cl.cex = 0.8, pch.cex =0.7)

write.csv((flattenCorrMatrix(M2$r, M2$P)), "Crohn.csv")

#Only significant correlations
library(lares) 
corr_cross(mydata, # dataset
           max_pvalue = 0.05, # show only sig. correlations at selected level
           top = 50) # display top 10 correlations, any couples of variables  )
           
corr_var(mydata)
