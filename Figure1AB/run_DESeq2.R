library(ggplot2)
library(DESeq2)
library(readr)
library(dplyr)
library(lazyeval)

workpath <- "./"

df_count <- read.csv(file.path(workpath,"guide_count.csv"))


all_conditions <- c(
                    "C9B_high","C9B_presort","C9B_high","C9B_presort","C9B_high","C9B_presort",  
                    "dCas9_high","dCas9_presort","dCas9_high","dCas9_presort","dCas9_high","dCas9_presort", 
                    "KRAB_high","KRAB_presort","KRAB_high","KRAB_presort","KRAB_high","KRAB_presort", 
                    "VP64_high","VP64_presort","VP64_high","VP64_presort","VP64_high", "VP64_presort"
                    "VP64_presort","VP64_presort","VP64_presort","VP64_high","VP64_presort","VP64_high","VP64_presort"
                    )


outputname<-file.path(workpath,"guide_count_with_l2fc.csv")

df_count_out <- df_count
sgRNAcolumn <- colnames(df_count)[1]
samplenames <-colnames(df_count)[-c(1)]

samplecolumns_list=list(c(1,3,5),c(7,9,11),c(13,15,17),c(19,21,23))
controlcolumns_list=list(c(2,4,6),c(8,10,12),c(14,16,18),c(20,22,24))

for(i in 1:length(samplecolumns_list))
{
  samplecolumns=samplecolumns_list[[i]]
  controlcolumns=controlcolumns_list[[i]]
  conditions=all_conditions[append(samplecolumns,controlcolumns)]
  
  df_count_filtered <- df_count %>%
    select(sgRNAcolumn,samplenames[append(samplecolumns,controlcolumns)])
  
  sampleTable <- data.frame(condition = factor(conditions))
  row.names(sampleTable) <- colnames(df_count_filtered)[-c(1)]
  
  dds <- DESeqDataSetFromMatrix(as.data.frame(df_count_filtered), colData = sampleTable, ~condition, tidy = TRUE)
  dds <- DESeq(dds)
  a2=all_conditions[samplecolumns][1]
  a1=all_conditions[controlcolumns][1]
  res <- results(dds,cooksCutoff=FALSE,contrast=c("condition",a2,a1))
  
  
  originalLength=length(colnames(df_count_out))
  df_count_out=cbind(df_count_out, res$log2FoldChange)
  

  comparestring=paste(a2,"_vs_",a1, sep='')
  l2fccolumn=paste(comparestring,"__L2FC", sep='')
  
  colnames(df_count_out)<-c(colnames(df_count_out)[1:originalLength],l2fccolumn)
  
 
}

write.table(df_count_out, file = outputname, append = FALSE, quote = FALSE, sep = ",",
            eol = "\n", na = "NA", dec = ".", row.names = FALSE,
            col.names = TRUE)

