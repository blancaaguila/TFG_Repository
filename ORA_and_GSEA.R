
# From a set of differential expressed genes, perform ORA and GSEA analysis
###########################################################################################
library(tidyverse)
library(dplyr)
library(ggplot2)
library(clusterProfiler) 

setwd('/media/blanca/BlancaAguila/rna_final/')
a<- readxl::read_xlsx('/media/blanca/BlancaAguila/rna_final/expression log2_FC all cell lines.xlsx')
filt <- a
log2fc <- transform(filt, HH_primed_1=as.numeric(HH_primed_1),
                    HH_8i_1=as.numeric(HH_8i_1),
                    HH_gr8i_1=as.numeric(HH_gr8i_1),
                    H1_primed_2 = as.numeric(H1_primed_2),
                    H1_8i_2=as.numeric(H1_8i_2),
                    H1_gr8i_2=as.numeric(H1_gr8i_2),
                    D2_primed_3=as.numeric(D2_primed_3),
                    D2_8i_3=as.numeric(D2_8i_3),
                    D2_gr8i_3=as.numeric(D2_gr8i_3),
                    D2_primed_4=as.numeric(D2_primed_4),
                    D2_8i_4=as.numeric(D2_8i_4),
                    D2_gr8i_4=as.numeric(D2_gr8i_4))

genelist.naive <-as.vector(filt$H1_8i_2, mode = "any")
names(genelist.naive) <- filt$Entrez.Gene

# ORDERED USAGE
genelist.naive <- genelist.naive[order(-genelist.naive)]

gene.naive <- names(genelist.naive)[abs(genelist.naive) >= log2(2)]

library(clusterProfiler)
library(org.Hs.eg.db)
gene.df.up <- bitr(gene.naive, fromType = "ENTREZID",
                toType = c("SYMBOL", 'ENSEMBL'),
                OrgDb = org.Hs.eg.db)



###########################################################################################
#############  Heatmap; pluripotent vs differential gene markers 
###########################################################################################

n <- c('POU5F1', 'Nanog', 'Sox2', 'Klf2', 'Klf4', 'Klf5', 'ZFP42', #Naive
       'Esrrb', 'Dppa3','Tfcp2l1', 'Fgf4', 'Tbx3', 'Cdh1',
       'POU5F1', 'Nanog', 'Sox2', #both
       'Klf2', 'Klf4', 'Klf5', 'klf17', 'Dppa3', 'tbx3', 'Tfcp2l1',
       'POU5F1', 'Nanog', 'Sox2', #both
       'Klf4', 'Klf5', 'klf17', 'Dppa3', 'Dppa5', 'Tfcp2l1', 
       'POU5F1', 'Nanog', 'Sox2', #both
       'Klf4', 'klf17', 'Dppa3','dnmt3l', 'Tfcp2l1', 
       'CDH1')# naive table  
n <- toupper(n)
n <- unique(n)

p <- c('POU5F1', 'Sox2', 'Dnmt3b', 'Fgf5', 'Pou3f1', #Primed
       'Meis1', 'Otx2', 'Sox11', 'Gdf3', 
       'Otx2', 'deusP6', 
       'Otx2', 'zic2', 'dusP6',
       'CDH2',  # PAPER
       'Pax6', 'meis1', 'hoxb1', 'lhx5', 'otx1', 'neurog1', # LEGG ectoderm dev
       'hand1', 'dlx5', 'myf5', 'onecut1', #mesoderm dev
       'isl1','zfhx3', 'esx1' ) # endoderm and embryotic develpoment
p <- toupper(p)
p <- unique(p)

l <- c(n,p)
tab.l <- as.data.frame(l)

#write.table(tab.l,'markers.tsv',  quote=FALSE, sep='\t', col.names = NA)
l<-toupper(l)

keytypes(org.Hs.eg.db)
s <- bitr(l, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)


log2fc$Entrez.Gene <- as.character(log2fc$Entrez.Gene)
a <- left_join(s, log2fc, by = c('ENTREZID'='Entrez.Gene'))
a1 <- a[,c(6,7)]
rownames(a1) <- a$SYMBOL
a1 <- as.matrix(a1)

as.data.frame(a1[,1]<a1[,2])
as.data.frame( a1[,3]<a1[,4])
as.data.frame( a1[,5]<a1[,6])
pdf('heatmap_interested_rep_all_markers.pdf', height = 6, width = 5)
heatmap(a1, scale = "row", cexCol = 1, cexRow = 0.5)
dev.off()

###########################################################################################
# GO ONTOLOGIES
##########################################################################################
############ Over-expression Analysis 

ego.naive.all.up <- enrichGO(gene          = gene.naive,
                          universe      = names(genelist.naive),
                          OrgDb         = org.Hs.eg.db,
                          ont           = "all",
                          pAdjustMethod = "BH",
                          pvalueCutoff  = 0.01,
                          qvalueCutoff  = 0.05,
                          readable      = TRUE)


head(ego.naive.all.down)

############## PLOTTING
a <- dotplot(ego.naive.all.up, showCategory=10, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free") +theme(text=element_text(size=5), axis.text.y = element_text(size=9),
                                                                                    axis.text.x = element_text(size = 9))



############# GSE
ego3.all <- gseGO(geneList     = genelist.naive,
                  OrgDb        = org.Hs.eg.db,
                  ont          = "all",
                  nPerm        = 1000,
                  minGSSize    = 100,
                  maxGSSize    = 500,
                  pvalueCutoff = 0.05,
                  verbose      = FALSE)

dotplot(ego3.all, showCategory=10, split=".sign") + facet_grid(.~.sign)
### PLOTTTING
p1 <-dotplot(ego.naive.all.up, showCategory=20, split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free") +theme(text=element_text(size=5), axis.text.y = element_text(size=9),
                                                                                                               axis.text.x = element_text(size = 9)) + ggtitle("dotplotfor ORA") 
p2 <-dotplot(ego3.all, showCategory=20,split="ONTOLOGY") +facet_grid(ONTOLOGY~., scale="free") +theme(text=element_text(size=5), axis.text.y = element_text(size=9),
                                                                                                      axis.text.x = element_text(size = 9)) + ggtitle("dotplotfor GSEA")
library(cowplot)
#pdf("png_figures/GO_dotplots_all.png",w=2000,h=1500,res=200)
#plot_grid(p1, p2, ncol=2)
#dev.off()

########### saving tables 
#write.table(as.data.frame(ego.naive.all.up), file='Enrichment_Analyisis_GO_All.tsv', quote=FALSE, sep='\t', col.names = NA)
#write.table(as.data.frame(ego3.all), file='Gene_Set_Enrichment_Analyisis_GO_All.tsv', quote=FALSE, sep='\t', col.names = NA)
table(ego.naive.all.up$ONTOLOGY)
table(ego3.all$ONTOLOGY)
###############
keytypes(org.Hs.eg.db)
gene.df.up <- bitr(ego.naive.all.up$ID, fromType = "GO",
                   toType = c("ENTREZID","ENSEMBL", "ENSEMBLPROT",  "ENSEMBLTRANS"),
                   OrgDb = org.Hs.eg.db)

ego.naive.all.up[ego.naive.all.up$Description=="Ras GTPase binding",]
gene.df.up[gene.df.up$GO=='GO:0017016',]



###################################################################################################
########### GSEA KEGG
###################################################################################################

search_kegg_organism('hsa', by='kegg_code')
hum <- search_kegg_organism('Homo sapiens', by='scientific_name')
head(hum)

### ORA
kk.up <- enrichKEGG(gene         = gene.naive,
                 organism     = 'hsa',
                 pvalueCutoff = 0.01)


nrow(kk.up)
a<- emapplot(kk.up, layout = "kk",showCategory = nrow(kk.up)) 
a <- emapplot(kk.down, layout = "kk",showCategory = nrow(kk.down)) 
library(DOSE)

### GSE KEGG

kk2 <- gseKEGG(geneList     = genelist.naive,
               organism     = 'hsa',
               nPerm        = 1000,
               minGSSize    = 120,
               pvalueCutoff = 1*10^-2,
               verbose      = FALSE)
head(kk2)
nrow(kk2)
nrow(kk.up)

mapk <- kk2[kk2$Description=="MAPK signaling pathway",]
strsplit(mapk$core_enrichment, '/')
#write.table(as.data.frame(kk.up), file='Enrichment_Analyisis_KEGG_All.tsv', quote=FALSE, sep='\t', col.names = NA)
#write.table(as.data.frame(kk2), file='Gene_Set_Enrichment_Analyisis_KEGG_All.tsv', quote=FALSE, sep='\t', col.names = NA)


#pdf("trial_dot_kk.up.58.pdf", height = 50, width = 60)
#cnetplot(kk.up, showCategory = nrow(kk.up),foldChange=genelist.naive, circular = TRUE, colorEdge = TRUE)
#dev.off()
heatplot(kk2, foldChange=genelist.naive)
emapplot(kk2, showCategory = nrow(kk2))
### Plot
a <-dotplot(kk.up, showCategory=30) + ggtitle("dotplot for ORA ")
a1 <- dotplot(kk2, showCategory=30) +ggtitle("dotplot for GSEA")
#pdf("figures/heatmap_ora_keg.pdf", height =5, width = 10)
#final_productKEGG_Overrep
#dev.off()


# ORA PLOTS
edox <- setReadable(kk.up, 'org.Hs.eg.db', 'ENTREZID')

p4 <- cnetplot(edox,foldChange=genelist.naive, circular = TRUE, colorEdge = TRUE)
p3 <- cnetplot(edox, showCategory=30,foldChange=genelist.naive, circular = TRUE, colorEdge = TRUE)

p4 <- emapplot(kk.up,showCategory = nrow(kk.up), pie_scale=1.5,layout="kk") 
p5 <- emapplot(kk.up, pie_scale=1.5,layout="kk") 

# alternative to cnt plot
a<-upsetplot(kk.up) + ggtitle("KEGG ORA")


# GSEA PLOTS
a <- ridgeplot(kk2, showCategory = nrow(kk2))
a
library(enrichplot)
library(ggupset)
a <- upsetplot(kk2,nset=8)
head(kk2)

gseaplot2(kk2, geneSetID = 1:3, pvalue_table = T)
gseaplot2(kk2, geneSetID = 20, title = kk2$Description[20])

setwd('/media/blanca/BlancaAguila/rna_final/figures/gsea_enrichment_KEGG_59')

# GSEA for the 59 sets
for(i in 1:nrow(kk2)){
    g <- gseaplot2(kk2, geneSetID = i, title = kk2$Description[i])
    ggsave(plot = g, file=paste(kk2$Description[i],'_gsea_plot.pdf', sep = ''), width = 7, height = 5) ## change if downsampled
    
}
dev.off()



