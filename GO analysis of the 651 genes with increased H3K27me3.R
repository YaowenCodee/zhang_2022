##### GO
# Create a list with genes from each sample
gene = lapply(peakAnnoRKR, function(i) as.data.frame(i)$geneId)
# Run GO enrichment analysis 
ego <- enrichGO(gene = entrez, 
                keytype = "ENTREZID", 
                OrgDb = org.Hs.eg.db, 
                ont = "BP", 
                pAdjustMethod = "BH", 
                qvalueCutoff = 0.05, 
                readable = TRUE)
# Dotplot visualization
dotplot(ego, showCategory=50)
# Multiple samples KEGG analysis
compKEGG <- compareCluster(geneCluster = gene, 
                           fun = "enrichKEGG",
                           organism = "human",
                           pvalueCutoff  = 0.05, 
                           pAdjustMethod = "BH")
dotplot(compKEGG, showCategory = 20, title = "KEGG Pathway Enrichment Analysis")