# This file aims to annotate non.captured Hind III fragments. Later they will be used to classify the interactions, 
# in diverse categories. The categories can be visualized by a sankey Plot
################################################################################################
options(stringsAsFactors = T)
setwd('/home/blanca/Desktop/c') 
library('ape')
library(varhandle)
library(GenomicFeatures)
library(GenomicInteractions)
library(IRanges)
library(tidyr)
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(dplyr)

#################################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
#############################################################################################
#############################################################################################
### (A) ANOTATATE HINDIII fragments #########################################################
#############################################################################################
gff <- read.gff('Homo_sapiens.GRCh37.87_subset.gff3', GFF3=T)
str(gff)
gff$seqid <- unfactor(gff$seqid)
gff$source <- unfactor(gff$source)
gff$type <- unfactor(gff$type)

hind <- read.table('digest_and_probes_homo_sapiens_hg19_updated_16_02_2020.txt', header = T, stringsAsFactors = F)
hind[is.na(hind)] <- "NA"
hind <- hind[hind$annotation=='.',] # other-ends

# 836961-21876 = number of hind III other end
#############################################################################################
### TRANSCRIPT FILTERING
#############################################################################################

attributes <- str_split_fixed(string = gff$attributes, pattern = ';', n=11)
table(gsub(':.*','',attributes[,1]))
# ID=gene     ID=transcript Parent=transcript 
# 57905            196501           1476041 
table(gsub('=.*','',attributes[,2]))
#          Name  Parent 
#280277 1253669  196501 
gff.trans <- gff[gsub(':.*', '', attributes[,1])=='ID=transcript',]
nrow(gff.trans)

gff.gene <- gff[gsub(':.*', '', attributes[,1])=='ID=gene',]
nrow(gff.gene)

##### Make genomic ranges
r.hind.copy <- makeGRangesFromDataFrame(hind, keep.extra.columns = T, start.field = 'start', end.field = 'end', seqnames.field ='chr')
r.gff.trans <- makeGRangesFromDataFrame(gff.trans, keep.extra.columns = T)

#############################################################################################
### GENE BODY  (1993215)
#############################################################################################

#final.hind <- final.hind %>% mutate(.,FULL_GENE_BODY= ifelse((final.hind %in% gene.body.hind),'TRUE', 'FALSE'))

mixed.gene.body.pre <- mergeByOverlaps(r.hind.copy, r.gff.trans, type='any')
gene.body.hind.2 <- unique(mixed.gene.body.pre$r.hind.copy)
final.hind.2 <- r.hind.copy
final.hind.2 <- final.hind.2 %>% mutate(., GENE_BODY=ifelse((final.hind.2 %in% gene.body.hind.2), 'TRUE', 'FALSE'))
# TRUE: overlaping any kind, FALSE: non overlaping
#############################################################################################
#### TSS: strand + (start) strand - (end) 
#############################################################################################

# CALCULATE TSS FOR EACH OF THE STRANDS
tss.mas <-r.gff.trans %>%  #####aqui tenia tss.overlap
  filter(strand=='+') %>%
  mutate(end=start)

tss.menys <- r.gff.trans %>%  
  filter(strand=='-') %>%
  mutate(start=end)

tss <- c(tss.mas, tss.menys)

### TSS within 
tss.ovelap <- mergeByOverlaps(r.hind.copy, tss) # any, by default, however, sine tss only one position, it would be within tss 
tss.ovelap

tss.overlap.hind <- unique(tss.ovelap$r.hind.copy)
final.hind.2 <- final.hind.2 %>% mutate(., TSS=ifelse((final.hind.2 %in% tss.overlap.hind), 'TRUE',
                                                      ifelse((final.hind.2 %in% gene.body.hind.2), 'FALSE', 'NA')))

# TRUE: tss overlaping within, FALSE: tss not overlaping, but still a gene body, NA: intergenic

#############################################################################################
##### TSS CODING
#############################################################################################
#table(tss.ovelap$r.hind.copy %in% group.overlaps.1$r.hind.copy)

r.gff.gene <- makeGRangesFromDataFrame(gff.gene, keep.extra.columns = T)
r.gff.gene

tss.ovelap
####
#(1) 
t <- str_split_fixed(string = tss.ovelap$attributes, pattern = ';', n=11)
parent.gene <- gsub('Parent=gene:','', t[,2]) # select id 

#(2)
t1 <- str_split_fixed(string = gff.gene$attributes, pattern = ';', n=11)
head(t1)
id.gene.gff <- gsub('ID=gene:','',t1[,1]) #select id
biotype <- gsub('biotype=','',t1[,3]) # select biotype

# (3)
# id for tss.overlap (table)
id.for.tss.overlap <- as.data.frame(parent.gene)
colnames(id.for.tss.overlap) <- 'Gene_ID'
# table for id and biotype of gene (table)
id.for.gff.gene <- as.data.frame(id.gene.gff)
id.for.gff.gene[,2] <- biotype
colnames(id.for.gff.gene) <- c('Gene_ID', 'Biotype')

# (4) join
join.gene.and.tss <- left_join(id.for.tss.overlap, id.for.gff.gene, by=('Gene_ID'='Gene_ID'))
#table((s$Gene_ID==trial$Gene_ID))
#table(is.na(join.gene.and.tss$Biotype))
#table(trial$Biotype)

# (5) now i have the biotype, for all the tss that overlap with a hind III
#table(join.gene.and.tss$Gene_ID==parent.gene) # 94432 true means okay to do the below steps
tss.ovelap$Gene_Biotype <- join.gene.and.tss$Biotype
tss.ovelap

table(tss.ovelap$Gene_Biotype)
prot.coding <- tss.ovelap[tss.ovelap$Gene_Biotype=='protein_coding',]

prot.coding.hind <- unique(prot.coding$r.hind.copy)
tss.overlap.hind <- unique(tss.ovelap$r.hind.copy)
final.hind.2 <- final.hind.2 %>% mutate(., PROTEIN_CODING_TSS=ifelse((final.hind.2 %in% prot.coding.hind), 'TRUE',
                                                                     ifelse((final.hind.2 %in% tss.overlap.hind), 'FALSE', 'NA')))


final.hind.2
#table(final.hind$PROTEIN_CODING)

#############################################################################################
#### INTERGENIC  (358555)
#############################################################################################

non.overlaps <- filter_by_non_overlaps(r.hind.copy, r.gff.trans)

non.overlaps.hind <- unique(non.overlaps)
final.hind.2 <- final.hind.2%>% mutate(., FULL_INTERGENIC=ifelse((final.hind.2 %in% non.overlaps.hind), 'TRUE', 'FALSE'))
#table(final.hind$FULL_INTERGENIC)
final.hind.2

#############################################################################################
#### MIXED AND MIXED PROTEIN CODING 
#############################################################################################

# hay que restar de mixed los tss true

# (1) find which hind III fragments are mixed gene body
final.hind.2 <- final.hind.2 %>% mutate(., MIXED_GENE_BODY=ifelse((final.hind.2$GENE_BODY==TRUE & final.hind.2$TSS==FALSE), 'TRUE', 'FALSE'))

# (2) Look which ones code for a protein coding and which ones not
mix.gene <- final.hind.2[final.hind.2$MIXED_GENE_BODY==TRUE,]
# r.gff.gene

# (3) find overlaps
# s <- filter_by_overlaps(mix.gene, r.gff.trans) # just to see that if we filter, we have the same number, meaning that all of them overlap
mix.trans <- mergeByOverlaps(mix.gene, r.gff.trans)
mix.trans

# (4) create a column just with the parents id
s.mix <- str_split_fixed(string=mix.trans$attributes, pattern = ';', n=11)
parent.mix <- gsub('Parent=gene:','', s.mix[,2]) # select id 
mix.trans$parent_id <- parent.mix #adds a new column with just the parent id. 

# (5) take just de id genes,
# nrow(id.for.gff.gene) table with gene id and type

mix.trans.tab <- as.data.frame(parent.mix)
colnames(mix.trans.tab) <- 'Gene_ID'
nrow(mix.trans.tab)

join.gene.and.mix <- left_join(mix.trans.tab, id.for.gff.gene, by=('Gene_ID'='Gene_ID'))
#table((s$Gene_ID==trial$Gene_ID))
head(join.gene.and.mix)

# adds the biotype from the r-gff.gene table
mix.trans$gene_type <- join.gene.and.mix$Biotype
mix.trans

table(mix.trans$gene_type)
prot.cod.mix <- mix.trans[mix.trans$gene_type=='protein_coding',]

prot.cod.mix.hind <- unique(prot.cod.mix$mix.gene)

final.hind.2 <- final.hind.2 %>% mutate(., PROT_COD_MIX =ifelse((final.hind.2 %in% prot.cod.mix.hind), 'TRUE', 
                                                                ifelse((final.hind.2 %in% mix.gene), 'FALSE', 'NA')))

table(final.hind.2$PROT_COD_MIX)



####################################################################################################################

final.hind.2
a <- final.hind.2 %>% mutate(., CATEGORY=ifelse((final.hind.2$TSS==TRUE & final.hind.2$PROTEIN_CODING_TSS==TRUE), 'tss_prot_coding',
                                                ifelse((final.hind.2$TSS==TRUE & final.hind.2$PROTEIN_CODING_TSS==FALSE), 'tss_non_coding_p',
                                                       ifelse((final.hind.2$MIXED_GENE_BODY==TRUE & final.hind.2$PROT_COD_MIX==TRUE), 'gene_body_p_coding',
                                                              ifelse((final.hind.2$MIXED_GENE_BODY==TRUE & final.hind.2$PROT_COD_MIX==FALSE), 'gene_body_non_coding_p',
                                                                     ifelse((final.hind.2$FULL_INTERGENIC==TRUE), 'full_intergenic', 'no'))))))
a

table(a$CATEGORY)
category.tab <- as.data.frame(a)




write.table(category.tab, file = 'Annotated_Hind_III_Categories.txt', sep = '\t', row.names = F)


########################################################################################################################

b <- final.hind.2 %>% mutate(., CATEGORY=ifelse((final.hind.2$TSS==TRUE & final.hind.2$PROTEIN_CODING_TSS==TRUE), 'tss_coding',
                                                ifelse((final.hind.2$TSS==TRUE & final.hind.2$PROTEIN_CODING_TSS==FALSE), 'tss_noncoding',
                                                       ifelse((final.hind.2$MIXED_GENE_BODY==TRUE & final.hind.2$PROT_COD_MIX==TRUE), 'genebody_coding',
                                                              ifelse((final.hind.2$MIXED_GENE_BODY==TRUE & final.hind.2$PROT_COD_MIX==FALSE), 'genebody_noncoding',
                                                                     ifelse((final.hind.2$FULL_INTERGENIC==TRUE), 'intergenic', 'no'))))))

category.tab.2 <- as.data.frame(b)
head(category.tab.2)
table(category.tab.2$CATEGORY)
write.table(category.tab.2, file = 'Annotated_Hind_III_Categories.2.txt', sep = '\t', row.names = F)

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
#############################################################################################
#############################################################################################
######### (B) Integrate the annotated file and plot Sankey
#############################################################################################


#############################################################################################
##### Function joining interaction coordinates with the anotation of captured and uncaptured fragments
#############################################################################################
categorize.1 <- function(interactions.file, capture, uncapture){
  interactions.file <- left_join(interactions.file, capture,  by=c('chr_I'='chr', 'start_I'='start', 'end_I'='end'))
  colnames(interactions.file)[7] <- 'catured_left'
  interactions.file <- left_join(interactions.file, capture, by=c('chr_II'='chr', 'start_II'='start', 'end_II'='end'))
  colnames(interactions.file)[8] <- 'catured_right'
  interactions.file <- left_join(interactions.file, uncapture,  by=c('chr_II'='chr', 'start_II'='start', 'end_II'='end'))
  return(interactions.file)
}

#############################################################################################
###### CLASSIFY
#############################################################################################

classify.1 <- function(interactions.file){
  primed.inter.3 <- interactions.file
  # capture- capture
  primed.inter.3 <- primed.inter.3 %>% mutate(., deti=ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$catured_right=='non_promoter'),'a:i',
                                                             ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$catured_right=='non_coding' | primed.inter.3$catured_left=='non_coding' & primed.inter.3$catured_right=='non_promoter'), 'a:j',
                                                                    ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$catured_right=='coding'| primed.inter.3$catured_left=='coding' & primed.inter.3$catured_right=='non_promoter'), 'a:k',
                                                                           ifelse((primed.inter.3$catured_left=='non_coding' & primed.inter.3$catured_right=='non_coding'), 'b:j',
                                                                                  ifelse((primed.inter.3$catured_left=='coding' & primed.inter.3$catured_right=='coding'), 'c:k',
                                                                                         ifelse((primed.inter.3$catured_left=='non_coding' & primed.inter.3$catured_right == 'coding'  | primed.inter.3$catured_left=='coding' & primed.inter.3$catured_right=='non_coding'), 'b:k', 'NA')))))))
  # capture-noncapture
  primed.inter.3 <- primed.inter.3 %>% mutate(., deti2= ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$noncaptured_right=='tss_noncoding'), 'a:d',
                                                               ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$noncaptured_right=='tss_coding'), 'a:e',
                                                                      ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$noncaptured_right=='genebody_noncoding'), 'a:f',
                                                                             ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$noncaptured_right=='genebody_coding'), 'a:g',
                                                                                    ifelse((primed.inter.3$catured_left=='non_promoter' & primed.inter.3$noncaptured_right=='intergenic'), 'a:h',
                                                                                           ifelse((primed.inter.3$catured_left=='non_coding' & primed.inter.3$noncaptured_right=='tss_noncoding'), 'b:d',
                                                                                                  ifelse((primed.inter.3$catured_left=='non_coding' & primed.inter.3$noncaptured_right=='tss_coding'), 'b:e',
                                                                                                         ifelse((primed.inter.3$catured_left=='non_coding' & primed.inter.3$noncaptured_right=='genebody_noncoding'), 'b:f',
                                                                                                                ifelse((primed.inter.3$catured_left=='non_coding' & primed.inter.3$noncaptured_right=='genebody_coding'), 'b:g',
                                                                                                                       ifelse((primed.inter.3$catured_left=='non_coding' & primed.inter.3$noncaptured_right=='intergenic'), 'b:h',
                                                                                                                              ifelse((primed.inter.3$catured_left=='coding' & primed.inter.3$noncaptured_right=='tss_noncoding'), 'c:d',
                                                                                                                                     ifelse((primed.inter.3$catured_left=='coding' & primed.inter.3$noncaptured_right=='tss_coding'), 'c:e',
                                                                                                                                            ifelse((primed.inter.3$catured_left=='coding' & primed.inter.3$noncaptured_right=='genebody_noncoding'), 'c:f',
                                                                                                                                                   ifelse((primed.inter.3$catured_left=='coding' & primed.inter.3$noncaptured_right=='genebody_coding'), 'c:g',
                                                                                                                                                          ifelse((primed.inter.3$catured_left=='coding' & primed.inter.3$noncaptured_right=='intergenic'), 'c:h', 'NA'))))))))))))))))
  return(primed.inter.3)
}


#############################################################################################
# Prepear table for Sankey
#############################################################################################
table.create <- function(inter.file, total){
  primed.inter.3 <- inter.file
  plot.names <- c("Captured non-promoters",          "Captured non-coding promoters", 
                  "Captured coding promoters" ,        "Captured non-promoter", 
                  "Captured non-coding promoter"  ,    "Captured coding promoter"  ,  
                  "Non-captured non-coding promoter" , "Non-captured coding promoter",
                  "Non-captured non-coding gene body" ,"Non-captured coding gene body","Non-captured intergenic" )         
  plot.names <- as.data.frame(plot.names)
  original.names <- c( "non_promoter" ,      "non_coding"   ,      "coding"        ,     "non_promoter"  ,    
                       "non_coding"    ,     "coding"     ,        "tss_noncoding"  ,    "tss_coding"   ,     
                       "genebody_noncoding", "genebody_coding"   , "intergenic"       )
  plot.names$original.names <- original.names
  names.id <- c("a" , "b",  "c",  "i" , "j" , "k",  "d",  "e",  "f",  "g", "h")
  plot.names$names.id <- names.id
  
  
  cap <- as.data.frame(table(primed.inter.3$catured_left))
  cap <- cap %>% mutate(per=(Freq/total)*100)
  cap <- cap %>% mutate(., id = ifelse(cap$Var1=='non_promoter', 'a',
                                       ifelse(cap$Var1=='non_coding', 'b',
                                              ifelse(cap$Var1=='coding', 'c', 'na'))))
  
  cap.r <- as.data.frame(table(primed.inter.3$catured_right))
  nc.cap.r <- as.data.frame(table(primed.inter.3$noncaptured_right))
  r <- rbind(cap.r, nc.cap.r)
  r <- r %>% mutate(per=(Freq/total)*100)
  r <- r %>% mutate(., id= ifelse(r$Var1=='tss_noncoding', 'd',
                                  ifelse(r$Var1=='tss_coding', 'e',
                                         ifelse(r$Var1=='genebody_noncoding', 'f',
                                                ifelse(r$Var1=='genebody_coding', 'g',
                                                       ifelse(r$Var1=='non_promoter', 'i', 
                                                              ifelse(r$Var1=='non_coding', 'j',
                                                                     ifelse(r$Var1=='coding', 'k',
                                                                            ifelse(r$Var1=='intergenic','h','na')))))))))
  
  
  
  t <- rbind(cap, r)
  
  percentage.table <- left_join(plot.names, t, by=c('names.id'='id'))
  return(percentage.table)
}

#############################################################################################
# Change format function
#############################################################################################

# Chr_I, Start_I, End_I , Chr_II, Start_II, End_II
change.format <- function(data.file){
  #data.file <- read.table(new.file, header = T, stringsAsFactors = F )
  colnames(data.file) <-  c("chr_I","start_I","end_I","gene_I", #rename columns
                            "chr_II","start_II","end_II","gene_II","R_I","CS_I")
  
  data.file$R_II <- data.file$R_I
  data.file$CS_II <- data.file$CS_I
  data.file <- data.file[,c(1,2,3,4,9,10,5,6,7,8,11,12)]
  data.file <- data.file[,c(1,2,3,7,8,9)] # IF WE want to plot sankey
  return(data.file)
}
#############################################################################################
###### PLOT SANKEY
#############################################################################################

plot.sankey.inter <- function(interactions.file, percentage.table, libname){
  a <- as.data.frame(table(interactions.file$deti))
  b <- as.data.frame(table(interactions.file$deti2))
  d <- rbind(a,b)
  
  int <- str_split_fixed(string = d$Var1, pattern = ':', n=2)
  links <- cbind(int, d$Freq)
  colnames(links)<- c('sou', 'tar', 'value')
  rownames(links) <- NULL
  
  links <- transform(links, value=as.factor(value))
  specify_decimal <- function(x, k) trimws(format(round(x, k), nsmall=k))
  k <- 2
  links <- links %>% mutate(., source=ifelse(links$sou=='a', paste('Captured non-promoters (', specify_decimal(percentage.table$per[percentage.table$names.id=='a'], k), '%)'), 
                                             ifelse(links$sou=='b', paste('Captured non-coding promoters (',  specify_decimal(percentage.table$per[percentage.table$names.id=='b'], k), '%)'),
                                                    ifelse(links$sou=='c', paste('Captured coding promoters (', specify_decimal(percentage.table$per[percentage.table$names.id=='c'], k), '%)'), 'na'))))
  
  links <- links %>% mutate(., target=ifelse(links$tar=='d', paste('Non-captured non-coding promoter (', specify_decimal(percentage.table$per[percentage.table$names.id=='d'], k), '%)'),
                                             ifelse(links$tar=='e', paste('Non-captured coding promoter (',  specify_decimal(percentage.table$per[percentage.table$names.id=='e'], k), '%)'),
                                                    ifelse(links$tar=='f', paste('Non-captured non-coding gene body (',  specify_decimal(percentage.table$per[percentage.table$names.id=='f'], k), '%)'),
                                                           ifelse(links$tar=='g', paste('Non-captured coding gene body (',  specify_decimal(percentage.table$per[percentage.table$names.id=='g'], k), '%)'),
                                                                  ifelse(links$tar=='h', paste('Non-captured intergenic (',  specify_decimal(percentage.table$per[percentage.table$names.id=='h'], k), '%)'),
                                                                         ifelse(links$tar=='i', paste('Captured non-promoter (',  specify_decimal(percentage.table$per[percentage.table$names.id=='i'], k), '%)'),
                                                                                ifelse(links$tar=='j', paste('Captured non-coding promoter (',  specify_decimal(percentage.table$per[percentage.table$names.id=='j'], k), '%)'),
                                                                                       ifelse(links$tar=='k', paste('Captured coding promoter (', specify_decimal(percentage.table$per[percentage.table$names.id=='k'], k), '%)'), 'na')))))))))
  
  
  #From these flows we need to create a node data frame: it lists every entities involved in the flow
  nodes <- data.frame(
    name=c(as.character(links$source), as.character(links$target)) %>% 
      unique()
  )
  # With networkD3, connection must be provided using id, not using real name like in the links dataframe.. So we need to reformat it.
  links$IDsource <- match(links$source, nodes$name)-1 
  links$IDtarget <- match(links$target, nodes$name)-1
  
  my_color <- 'd3.scaleOrdinal() .domain([nodes]) .range(["#69b3a2"])'
  
  # Make the Network. I call my colour scale with the colourScale argument
  p <- sankeyNetwork(Links = links, Nodes = nodes, Source = "IDsource", Target = "IDtarget", 
                     Value = "value", NodeID = "name", fontSize = 10, colourScale = my_color)
  p %>% saveNetwork(file=paste(libname,'.html', sep = ''))
  #return(p)
}

##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
##########################################################################################################################################################################
#############################################################################################
############################################################################################# 
### (C) Iterate over all the interaction files 
#############################################################################################

#setwd('/media/blanca/JavierreLab/JavierreLab/CDK8i_project/correlations/')
library('ape')
library(varhandle)
library(GenomicFeatures)
library(GenomicInteractions)
library(IRanges)
library(tidyr)
library(GenomicRanges)
library(plyranges)
library(tidyverse)
library(dplyr)
library(networkD3)
library(webshot)
# install_phantomjs()
library(htmlwidgets)

#############################################################################################
## Read files
#############################################################################################
# Captured/non.captured Hind III annotated fragements
non.captured <-read.table('Annotated_Hind_III_Categories.2.txt', header = T)
captured <- read.table('captured_categories.txt', header = T)
head(non.captured)
head(captured)

# (1) Reorganize captured and uncaptured table such: 
# captured chr, start, end, category.1
# uncaptured chr, start, end, category.2
captured <- captured %>% separate(region, c('chr', 'start', 'end'), ':')
captured <- transform(captured, chr=as.factor(chr),
                      start=as.integer(start), 
                      end=as.integer(end))
captured <- as_tibble(captured)
captured <- captured[,c(1,2,3,5)]
colnames(captured) <- c("chr", "start", "end", "captured_category")

non.captured <- non.captured[,c(1,2,3,14)]
colnames(non.captured) <- c("chr", "start", "end", "noncaptured_right")


#setwd('/media/blanca/JavierreLab/BLANCA/lab_externo/sankey_plot/downsampled_corrected/') # so the plots are saved there
s <- dir('/home/blanca/Desktop/c/file_ibed/')
s <- str_split_fixed(s, pattern = '_', n=2)
s <- s[,1]
s

#### save html
path = '/home/blanca/Desktop/c/file_ibed/'
out.file<-""
file.names <- dir(path, pattern =".ibed")
datalist = list()
for(i in 1:length(file.names)){
  file <- read.table(paste(path,'/',file.names[i], sep = ''),header=TRUE, stringsAsFactors = F)
  file$i <- i 
  datalist[[i]] <- file
}

length(datalist)
as_tibble(datalist[[1]])

for(i in 1:length(s)){
  step.1 <- change.format(datalist[[i]]) # sankey_p53_30abril2020
  step.2 <- categorize.1(step.1, captured, non.captured) #26abirl2020
  step.3 <- classify.1(step.2) # 26abril2020
  tabi <- table.create(step.3, nrow(step.3)) #26abril2020
  plot.sankey.inter(step.3, tabi, s[i]) # sankey_p53_30abril2020
  message(paste('Plot', i, 'have worked'))
}






