# This file takes as an input a correlation matrix, from the original PCHI-C bam,
# and plots the correlation between the different samples. 
########################################################################################################################
library(viridis)
library(ggplot2)
library(gridExtra)
library(grid)
library(lattice)
#library(ggpubr)
#library(cowplot)
library(ggstatsplot)
library(rlist)
library(dplyr)
  
#######################
# call dispersion point table
#######################

setwd("/media/blanca/JavierreLab/BLANCA/lab_externo/correlations") #File locations
b<- read.delim('blanca.txt') # points to plot
b <- b[,3:14]
b
correlation_values <- read.table('cor_table.txt') # correlation values
c <- as.matrix(correlation_values)
c

##################

distances.calculated <- calculate_distances(b)
reference.table<- ref.tab(distances.calculated,ncol(b))
final <- plot1(distances.calculated, reference.table, c)
splotdistribution(distances.calculated)
########################################################################################
### CALCULATE DISTANCES : calculates diff on reads comparing each of the samples
#########################################################################################
calculate_distances <- function(filei){
  s <- ncol(filei)
  for(l in 1:s){ ##### CHANGE
    for(i in 1:s){  ##### CHANGE
      filei[,paste(colnames(filei)[l], colnames(filei)[i], sep="_")] <- abs(filei[,l]-filei[,i])
    }
  }
  return(filei)
}


##########################################################################################
### REORGANIZE MATRIX : into a data set (left colum comparisons, right correlation value of the comparison)
#########################################################################################
# not need it !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
reorg.matrix <- function(cor.mat){
  output <- vector(mode = "numeric", length = 0)
  output <- as.data.frame(output)
  for(j in 1:nrow(cor.mat)){
    for(i in 1:ncol(cor.mat)){
      output[paste(rownames(cor.mat)[j], colnames(cor.mat)[i], sep = '_'),] <- cor.mat[j,i]
    }
  }
  sub <- output[,1] #just picking variables
  sub_mat <- matrix(sub, nrow = 12,ncol = 12) ######### CHANGE (of original mat)
  if(table(sub_mat==cor.mat)){
    print('Matrix matches original values')
  }else{
    print('Something failed')
  }
  return(output)
}


############################################################################################################
##### REFERENCE TABLE
###########################################################################################################
# n col initials 
ref.tab <- function(dist_b,n){
  # creates a data frame indicating the position where the comparisons start until the end 
    # seq(ncol(initial.file)+1, ncols(distance.file), by=1)
  pos.com <- as.data.frame(seq(n+1, ncol(dist_b), by=1)) ####### CHANGE 
  colnames(pos.com) <- 'pos'
  pos.com
  # saves the comparisons names in a vector
  com.initial <- colnames(dist_b)[(n+1):ncol(dist_b)] ####  CHANGE
  # Split the comparion names of the initial dataset in a list 
  splited = list()
  for(i in 1:length(com.initial)){
    splited[i] <- strsplit(com.initial[i], '_')
  }
  # add them to the reference table
  for(i in 1:length(splited)){
    pos.com[i,2] <- splited[[i]][1]
    pos.com[i,3] <- splited[[i]][2]
  }
  colnames(pos.com)<- c('pos', 'comp1', 'comp2')
  # creates a data frame with the sample names and its index in the initial file
  ref <- data.frame(name = c("ETO1KO" , "ETO1WT"  ,"ETO2KO" , "ETO2WT" , "MOCK1KO", "MOCK1WT", 
                             "MOCK2KO", "MOCK2WT", "NUT1KO" , "NUT1WT","NUT2KO", "NUT2WT"),index = c(1,2,3,4,5,6,7,8,9,10,11,12), stringsAsFactors = F)
  
  ref <- as_tibble(ref)
  pos.com <- as_tibble(pos.com)
  # joins such that we have a data frame indicating
      # 1. position of the comparison in the initial datataset
      # 2. Sample name of the first comparison
      # 3. Sample name of the second compariso
      # 4. Position of the matrix for the first sample
      # 5. Position of the matrix for the second sample
  tf <- left_join(pos.com, ref, by=c('comp1'='name'))
  t3 <- left_join(pos.com, ref, by=c('comp2'='name'))
  t3
  tf[,5] <- t3[,4]
  colnames(tf) <- c('pos', 'comp1', 'comp2', 'pos.x', 'pos.y')
  return(tf)
  
}

#table(b1$ETO1KO==b1[,b2$pos.y[2]])
#####################################################################################################################
########## PLOT CORRELATION
#####################################################################################################################
#### IT has to start in 13

b2 ### error en el negth de la table (devemos restar 12)

#b2$pos.x[9]
#b2$pos.y[9]
# calling matrix
#c[b2$pos.x[9], b2$pos.y[9]]
plot1 <- function(file, re.tab, sub_mat){
  for(i in 1:ncol(file)){   ##### !!!
    G<- ggplot(file, aes(x=file[,re.tab$pos.x[i]],y=file[,re.tab$pos.y[i]], color=log(file[,re.tab$pos[i]]+1)))+ theme_minimal() +
      geom_point( position = position_jitterdodge(), size=0.5,alpha = 0.5) +
      scale_color_viridis(discrete = F,begin = 1, end = 0)+
      annotate(geom = 'text', -Inf, Inf,
               label=paste("atop(' R= '*", sub_mat[re.tab$pos.x[i],re.tab$pos.y[i]],",' R²= '*",sub_mat[re.tab$pos.x[i],re.tab$pos.y[i]]**2,")"), parse=T,
               hjust=0, vjust=1, size = 3)+
      labs(x=paste('Replicate', re.tab$comp1[i]), y=paste('Replicate', re.tab$comp2[i]))+
      geom_smooth(method = lm, se = F, colour='red', size=0.1) +
      theme(legend.title = element_blank())
    ggsave(plot = G, file=paste(re.tab$comp1[i],'_', re.tab$comp2[i],'.pdf', sep = '')) ## change if downsampled
  }
}


#################################################################################################################
######## PLOT DISTRIBUTION
###############################################################################################################

plotdistribution <- function(file){
  for(i in 13:ncol(file)){ ### chagne
    pdf(paste('dist_', colnames(file)[i], '.pdf', sep = ''))
    plot(table(file[,i]), main = paste('Distribution', colnames(file)[i], split=""), xlab = 'num of reads', ylab = 'counts')
    dev.off()
  }
}
