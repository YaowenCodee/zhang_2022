suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
rm(list=ls())
dataset <- read.csv("gain_H3K27me3_codinggene_651_volcana.csv",header = T)



cut_off_pvalue = 0.2
cut_off_logFC = 0.5

dataset$change = ifelse(dataset$pvalue < cut_off_pvalue & abs(dataset$log2FoldChange) >= cut_off_logFC, 
                        ifelse(dataset$log2FoldChange> cut_off_logFC ,'Up','Down'),
                        'Stable')
# ===================================
p <- ggplot(
  
  dataset, 
  aes(x = log2FoldChange, 
      y = -log10(pvalue), 
      colour=change)) +
  geom_point(alpha=1, size=3) +
  scale_color_manual(values=c("#546de5", "#A9A9A9","#DA070E"))+
  xlim(c(-1,1))+ ylim(c(-0.0002,3)) + 
  
  geom_vline(xintercept=c(-0.5,0.5),lwd=1,lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(cut_off_pvalue),lwd=1,lty=4,col="black",lwd=0.8) +
  
  
  labs(x="log2(fold change)",
       y="-log10 (p-value)")+
  theme_bw()+
  
  
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )

plot(p)

#dataset$label = ifelse(dataset$pvalue < cut_off_pvalue & dataset$log2FoldChange < -1, as.character(dataset$gene),"")
dataset$label = ifelse(
                       dataset$genename %in% "sog", 
                       as.character(dataset$gene),"")

p + geom_text_repel(data = dataset, aes(x = log2FoldChange, 
                                        y = -log10(pvalue), 
                                        label = label),
                    size = 2,box.padding = unit(0.5, "lines"),
                    point.padding = unit(0.8, "lines"), 
                    segment.color = "black", 
                    show.legend = FALSE)
ggsave(".pdf",width=7,height=4) 
