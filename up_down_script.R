
# 
setwd('/media/blanca/BlancaAguila/rna_final/')
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
#### HEAT MAP WITH NAIVE GENE SIGNATURES
##########################################################################################

#pdf('Heatmap_rep.pdf', height = 6, width = 6)
heatmap(a1, scale = "row", cexCol = 1, cexRow = 0.8)
#dev.off()
pdf('mds_plot_interested_rep.pdf', height = 5, width = 10)
plotMDS(a1)
dev.off()
# 1 gene puede contener varios transcrito

########################
#############  link
########################
l0 <- c("KHDC1L",  "FAM151A" ,"HORMAD1", "ALPPL2",  "ZNF729",  "KHDC3L", #n
        "TRIM60",  "MEG8",  "OLAH",    "LYZ", "HYAL4",   "ALPP", #N
        "DUSP6",  "FAT3",    "THY1",    "STC1",    "KLHL4", #p
        "ZDHHC22", "NEFM",    "HMX2","PLA2G3",  "PTPRZ1",  "CYTL1", "SOX11") #p
# first
l1 <- c('POU5F1', 'Nanog', 'Sox2', 'Klf2', 'Klf4', 'Klf5', 'rex1', #Naive
       'Esrrb', 'Dppa3','Tfcp2l1', 'Fgf4', 'Tbx3', 'Cdh1', #Naive
       'POU5F1', 'Sox2', 'Dnmt3b', 'Fgf5', 'Pou3f1', #Primed
       'Meis1', 'Otx2', 'Sox11', 'Gdf3') # Primed
# third
l2 <- c('POU5F1', 'Nanog', 'Sox2', #both
       'Klf2', 'Klf4', 'Klf5', 'klf17', 'Dppa3', 'tbx3', 'Tfcp2l1', #N
       'Otx2', 'zic2', 'zic3', 'zic5') #P
# fourth
l3 <- c('POU5F1', 'Nanog', 'Sox2', #both
       'Klf4', 'Klf5', 'klf17', 'Dppa3', 'Dppa5', 'Tfcp2l1', #N
       'Otx2', 'zic2', 'zic3', 'deusP6') #P
# fifth
l4 <- c('POU5F1', 'Nanog', 'Sox2', #both
       'Klf4', 'klf17', 'Dppa3','dnmt3l', 'Tfcp2l1', #N
       'Otx2', 'zic2', 'dusP6') #P 

n <- c('POU5F1', 'Nanog', 'Sox2', 'Klf2', 'Klf4', 'Klf5', 'ZFP42', #Naive
       'Esrrb', 'Dppa3','Tfcp2l1', 'Fgf4', 'Tbx3', 'Cdh1',
       'POU5F1', 'Nanog', 'Sox2', #both
       'Klf2', 'Klf4', 'Klf5', 'klf17', 'Dppa3', 'tbx3', 'Tfcp2l1',
       'POU5F1', 'Nanog', 'Sox2', #both
       'Klf4', 'Klf5', 'klf17', 'Dppa3', 'Dppa5', 'Tfcp2l1', 
       'POU5F1', 'Nanog', 'Sox2', #both
       'Klf4', 'klf17', 'Dppa3','dnmt3l', 'Tfcp2l1', 
       'CDH1')# naive table  
length(n)
n <- toupper(n)
n <- unique(n)
n
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
length(p)
p
l <- c(n,p)
l
tab.l <- as.data.frame(l)
n
p
#write.table(tab.l,'markers.tsv',  quote=FALSE, sep='\t', col.names = NA)
l<-toupper(l)
unique(l)
length(l)
keytypes(org.Hs.eg.db)
s <- bitr(l, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
s
nrow(s)
library(tidyverse)
log2fc$Entrez.Gene <- as.character(log2fc$Entrez.Gene)
a <- left_join(s, log2fc, by = c('ENTREZID'='Entrez.Gene'))
nrow(a)
head(a)
a1 <- a[,c(6,7)]
rownames(a1) <- a$SYMBOL
a1 <- as.matrix(a1)
a1
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
pdf("png_figures/GO_dotplots_all.png",w=2000,h=1500,res=200)
plot_grid(p1, p2, ncol=2)
dev.off()

########### saving tables 
#write.table(as.data.frame(ego.naive.all.up), file='Enrichment_Analyisis_GO_All.tsv', quote=FALSE, sep='\t', col.names = NA)
#write.table(as.data.frame(ego3.all), file='Gene_Set_Enrichment_Analyisis_GO_All.tsv', quote=FALSE, sep='\t', col.names = NA)
ego.naive.all.up
ego3.all
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
a
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

head(kk2)
table(kk2$pvalue)
table(kk.up$pvalue)

pdf("trial_dot_kk.up.58.pdf", height = 50, width = 60)
cnetplot(kk.up, showCategory = nrow(kk.up),foldChange=genelist.naive, circular = TRUE, colorEdge = TRUE)
dev.off()
heatplot(kk2, foldChange=genelist.naive)
emapplot(kk2, showCategory = nrow(kk2))
### Plot
a<-dotplot(kk.up, showCategory=30) + ggtitle("dotplot for ORA ")
a1 <- dotplot(kk2, showCategory=30) +ggtitle("dotplot for GSEA")
pdf("figures/heatmap_ora_keg.pdf", height =5, width = 10)
final_productKEGG_Overrep
dev.off()


# ORA PLOTS
edox <- setReadable(kk.up, 'org.Hs.eg.db', 'ENTREZID')
head(edox)
nrow(kk.up)
p4 <- cnetplot(edox,foldChange=genelist.naive, circular = TRUE, colorEdge = TRUE)
p3 <- cnetplot(edox, showCategory=30,foldChange=genelist.naive, circular = TRUE, colorEdge = TRUE)

p4 <- emapplot(kk.up,showCategory = nrow(kk.up), pie_scale=1.5,layout="kk") 
p5 <- emapplot(kk.up, pie_scale=1.5,layout="kk") 

# alternative to cnt plot
a<-upsetplot(kk.up) + ggtitle("KEGG ORA")
a
head(kk.up)

# GSEA PLOTS
a <- ridgeplot(kk2, showCategory = nrow(kk2))
a
library(enrichplot)
library(ggupset)
a <- upsetplot(kk2,nset=8)
head(kk2)

head(kk2[,1:10])

gseaplot2(kk2, geneSetID = 1:3, pvalue_table = T)
gseaplot2(kk2, geneSetID = 20, title = kk2$Description[20])
kk2[20,]
nrow(kk2)
head(kk2)
head(kk2)
setwd('/media/blanca/BlancaAguila/rna_final/figures/gsea_enrichment_KEGG_59')

# GSEA for the 59 sets
for(i in 1:nrow(kk2)){
    g <- gseaplot2(kk2, geneSetID = i, title = kk2$Description[i])
    ggsave(plot = g, file=paste(kk2$Description[i],'_gsea_plot.pdf', sep = ''), width = 7, height = 5) ## change if downsampled
    
}
dev.off()

#final_productKEGG_Overrep <- heatplot(kk.up, foldChange = genelist.naive, showCategory 
 #                                     = 20000)+ ggplot2::coord_flip()



##### DAVIDS 


