### Classifies PCHi-C interactions by distance and by promoter-promoter vs promoter-other end
###################################################################################################################

library(stringr)
library(GenomicRanges)
# files alredy cleaned

# Get samples name in order
s <- dir('/media/blanca/JavierreLab/BLANCA/lab_externo/p53_chicago_ibed_files/chicago.ibed/hESCs/')
s <- str_split_fixed(s, pattern = '_', n=2)
s <- s[,1]
s

# duplicate reads and chicago score columns
change.format <- function(data.file){
  #data.file <- read.table(new.file, header = T, stringsAsFactors = F )
  colnames(data.file) <-  c("chr_I","start_I","end_I","gene_I", #rename columns
                            "chr_II","start_II","end_II","gene_II","R_I","CS_I")
  
  data.file$R_II <- data.file$R_I
  data.file$CS_II <- data.file$CS_I
  data.file <- data.file[,c(1,2,3,4,9,10,5,6,7,8,11,12)]
  return(data.file)
}
#n <- change.format('/media/blanca/JavierreLab/BLANCA/lab_externo/p53_chicago/chicago.ibed/ETO1KO_clean.ibed')
#n <- as_tibble(n)
#n

#write.table(n, file = paste('/media/blanca/JavierreLab/BLANCA/lab_externo/p53_chicago/chicago.ibed/', s[1], sep = ''), col.names=T )
# 
#new_data.anotated <- reanotation(n, "/media/blanca/JavierreLab/BLANCA/project/chicago_output/digest_and_probes_homo_sapiens_hg19.txt")

#new_data.anotated

###############################################################################################################
# By fragment 1 vs fragment 2 distance
###############################################################################################################
list.distances= list()
calculatedis <- function(file, libname){
  #file <- read.table(file, header = T)
  # Calculates mid distances
  file_left <- makeGRangesFromDataFrame(file, seqnames.field = "chr_I", start.field = "start_I", end.field = "end_I", keep.extra.columns = TRUE) #calculate distances from the left and the right side of the file
  file_right <- makeGRangesFromDataFrame(file, seqnames.field = "chr_II", start.field = "start_II", end.field = "end_II", keep.extra.columns = TRUE)
  file_inter <- GenomicInteractions(file_left, file_right) # rebind the range table
  file_dist <- calculateDistances(file_inter, method = 'midpoint')
  file_dist[is.na(file_dist)] <- 0 # to avoid errors with NA values
  # table(is.na(cjl1_dist))
  # the sum of this lenght (plus the trans) must equal the lenght of cjl1_dist
  # vector which takes the amount of interaction <1mb, between 1 and 5 mb ... 
  
  # constructing dataframe of distances
  dist_mb <- c(length(file_dist[file_dist<1000000 & file_dist > 0]),
               length(file_dist[1000000<file_dist & file_dist< 5000000]),
               length(file_dist[5000000<file_dist & file_dist< 10000000]),
               length(file_dist[10000000<file_dist & file_dist< 50000000]),
               length(file_dist[50000000<file_dist & file_dist< 100000000]),
               length(file_dist[file_dist>100000000]),
               length(file_dist[file_dist==0])) #trans interactions
  
  #to cheack if correct
  separate_vals <- sum(dist_mb) 
  if(separate_vals==length(file_dist)){
    message('Everyting is correct until this point')
  }
  
  # create dataframe
  dist_mb <- data.frame(dist_mb)
  col2 <-c('<1MB', '1-5MB', '5-10MB', '10-50MB', '50-100MB', '>100MB', 'Trans')
  dist_mb[,2] <- col2
  colnames(dist_mb) <-c('amount', 'size')
  #pdf(file = paste(libname, '.pdf', sep = '.'))
  return(dist_mb) # For the other plots
  p <- ggplot(dist_mb, aes(size, amount, fill=factor(size))) + 
    geom_text(aes(label=amount), vjust=-0.3, size=3.5)+
    ggtitle(libname) + #theme(text = element_text(size=9)) # when power point
    geom_bar(stat = 'identity') + scale_fill_brewer(breaks = c('<1MB', '1-5MB', '5-10MB', '10-50MB', '50-100MB', '>100MB', 'Trans')) +
    scale_x_discrete(limits = c('<1MB', '1-5MB', '5-10MB', '10-50MB', '50-100MB', '>100MB', 'Trans'))
 #ggsave(filename = paste(libname, '_distance.pdf', sep = ''), plot = p)
 #return(p) # when power point
  
}


#a<- calculatedis('/media/blanca/JavierreLab/BLANCA/lab_externo/p53_chicago/chicago.ibed/ETO1KO', s[1])


###################################################################################################################
#### P-P vs P-OE
###################################################################################################################

calcdistprominter <- function(file, libname){
  #file <- read.table(file, header = T)
  num_inter <- table(file=='.') # create a table with the amount of lines with at list a dot (TRUE)
  promoter_otherends <- num_inter[2] # Those lines with a dot are promoter other end interactions
  promoter_promoter <- nrow(file)-promoter_otherends
  if ((promoter_otherends+promoter_promoter) == nrow(file)){
    message('Amount of p-p vs p-oe is correct')
  }
  
  # perform a data frame so function can return it
  inter <- as.data.frame(promoter_otherends)
  inter[2,] <- promoter_promoter
  inter[,2] <- c('promoter_OE', 'promoter_promoter')
  rownames(inter) <- c('1', '2')
  colnames(inter) <- c('amount', 'type')
  return(inter) # for all plot
  p <- ggplot(inter, aes(type, amount, fill=type)) + 
    geom_text(aes(label=amount), vjust=-0.3, size=3.5)+
    ggtitle(libname) + #theme(text = element_text(size=9)) # when power point
    geom_bar(stat = 'identity') + scale_fill_brewer(palette="Accent")
  #ggsave(filename = paste(libname, '_distance1.pdf', sep = ''), plot = p)
  # return(p) # when power point
}

#a<- calcdistprominter('/media/blanca/JavierreLab/BLANCA/lab_externo/p53_chicago/chicago.ibed/ETO1KO', s[1])



library(cowplot)
library(ggplot2)
# CHANGE TO SAMPLED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
setwd('/media/blanca/JavierreLab/BLANCA/lab_externo/distances_chicago/sampled/')

#s <- dir('/media/blanca/JavierreLab/BLANCA/lab_externo/p53_chicago_ibed_files/downsampled_corrected.ibed/')
#s <- str_split_fixed(s, pattern = '_', n=2)
#s <- s[,1]
#s
s
path = '/media/blanca/JavierreLab/BLANCA/lab_externo/p53_chicago_ibed_files/chicago.ibed/hESCs/'
out.file<-""
file.names <- dir(path, pattern =".ibed")
datalist = list()
for(i in 1:length(file.names)){
  file <- read.table(paste(path,'/',file.names[i], sep = ''),header=TRUE, stringsAsFactors = F)
  file$i <- i 
  datalist[[i]] <- file
}
head(datalist[[1]])
length(datalist)
library(GenomicInteractions)
# (Esta echo sin reanotar) original 
for(i in 1:length(datalist)){
  step.1 <- change.format(datalist[[i]])
  calculatedis(step.1, s[i])
  calcdistprominter(step.1, s[i])
}

##################################
###############################################
# for poweer point
#  Before this make some changes in the function above 
dis1 = list()
dis2 =list()
for(i in 1:length(s)){
  step.1 <- change.format(datalist[[i]])
  dis1[[i]] <- calculatedis(step.1, s[i])
  dis2[[i]] <- calcdistprominter(step.1, s[i])
}

#############################
##################################################
### cdk8i primed PLOTS
dis1 = list()
dis2 =list()
for(i in 1:length(s)){
  step.1 <- change.format(datalist[[i]])
  #dis1[[i]] <- calculatedis(step.1, s[i])
  dis2[[i]] <- calcdistprominter(step.1, s[i])
}

# name dist
names(dis1) <- s
# TO PLOT
for(i in 1:length(dis1)){
  dis1[[i]][,3] <- s[i]
  colnames(dis1[[i]])[3] <- "rep"
}


tab.1 <- rbind(dis1[[1]], dis1[[2]], dis1[[3]], dis1[[4]])
tab.1
# saved in setwd 
pdf('All_hESC.pdf')
ggplot(tab.1, aes(size, amount)) + 
  #geom_text(aes(label=amount), vjust=1.6, position = position_dodge(0.9),size=3.5)+
  ggtitle("Distance") + #theme(text = element_text(size=9), legend.position = c(0.8, 0.2)) +# when power point
  geom_bar(aes(fill=rep),stat = 'identity', position = position_dodge()) +# scale_fill_brewer()+
  scale_x_discrete(limits = c('<1MB', '1-5MB', '5-10MB', '10-50MB', '50-100MB', '>100MB', 'Trans')) +
  theme_minimal() 

dev.off()
###########
names(dis2) <- s
dis2
for(i in 1:length(dis2)){
  dis2[[i]][,3] <- s[i]
  colnames(dis2[[i]])[3] <- "rep"
}
dis2
tab.2 <- rbind(dis2[[1]], dis2[[2]], dis2[[3]], dis2[[4]])
pdf('pp_pe_all.pdf')
ggplot(tab.2, aes(type, amount)) +
  geom_bar(aes(fill=rep), stat = 'identity', position = position_dodge())+
  theme_minimal() + ggtitle("PP/P-OE")
dev.off()

library(gridExtra)
# SAVE Sampled in an arrange
g <- arrangeGrob(dis1[[1]], dis1[[2]], dis1[[3]], dis1[[4]], nrow=2)
ggsave(file="ETO.png", g, width = 30, height = 20, units = 'cm')
g <- arrangeGrob(dis1[[5]], dis1[[6]], dis1[[7]], dis1[[8]], nrow=2)
ggsave(file="MOCK.png", g, width = 30, height = 20, units = 'cm')
g <- arrangeGrob(dis1[[9]], dis1[[10]], dis1[[11]], dis1[[12]], nrow=2)
ggsave(file="NUTLIN.png", g, width = 30, height = 20, units = 'cm')

g <- arrangeGrob(dis2[[1]], dis2[[2]], dis2[[3]], dis2[[4]], nrow=2)
ggsave(file="ETO_1.png", g, width = 30, height = 20, units = 'cm')
g <- arrangeGrob(dis2[[5]], dis2[[6]], dis2[[7]], dis2[[8]], nrow=2)
ggsave(file="MOCK_1.png", g, width = 30, height = 20, units = 'cm')
g <- arrangeGrob(dis2[[9]], dis2[[10]], dis2[[11]], dis2[[12]], nrow=2)
ggsave(file="NUTLIN_1.png", g, width = 30, height = 20, units = 'cm')

s
# SAVE Downsampled in an arrange
g <- arrangeGrob(dis1[[1]], dis1[[2]], nrow=2)
ggsave(file="MOCK.png", g, width = 30, height = 40, units = 'cm')
g <- arrangeGrob(dis1[[3]], dis1[[4]], nrow=2)
ggsave(file="NUT.png", g, width = 30, height = 40, units = 'cm')

g <- arrangeGrob(dis2[[1]], dis2[[2]], nrow=2)
ggsave(file="MOCK_1.png", g, width = 30, height = 40, units = 'cm')
g <- arrangeGrob(dis2[[3]], dis2[[4]], nrow=2)
ggsave(file="NUT_1.png", g, width = 30, height = 40, units = 'cm')
