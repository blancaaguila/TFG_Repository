# Script which removes duplicates from a segmonk file (chicago output)
############################################################################################
clean_seqmonk <- function(file)
{
  ## This should be removed since when you install the package {tidyverse} should be
  ## installed
  
  if(!("dplyr" %in% (.packages())))
  {
    library(dplyr,quietly = T,warn.conflicts = F)
  }

 
  data <- read.table(file, header = F, stringsAsFactors = F) #reading seqmonk file
  new_data <- rbind(cbind(data[seq(1,nrow(data),2),],data[seq(2,nrow(data),2),]), #reorganising file from 1int every two lines to just one line
                    cbind(data[seq(2,nrow(data),2),],data[seq(1,nrow(data),2),])) #a copy of the data is created but reversed (remember laures pptx)
  
  new_data <- dplyr::as_tibble(new_data[order(as.numeric(rownames(new_data))),], #restructure into a tibble
                                 .name_repair = "minimal")
  
  colnames(new_data) <- c("chr_I","start_I","end_I","gene_I", "R_I","CS_I", #rename columns
                          "chr_II","start_II","end_II","gene_II","R_II","CS_II")
  
  new_data$gene_I <- stringr::str_replace_all(new_data$gene_I, "\\|",",") #changing all names to same format
  new_data$gene_II <- stringr::str_replace_all(new_data$gene_II, "\\|",",")
  
  new_data <- unique(new_data) #ensuring that there are no real duplicates that were perhaps missed during chicago pipeline 
  #i.e. dup1 = A->B dup2 = A->B
  
  new_data <- new_data %>% group_by(chr_I, start_I, end_I, chr_II, start_II, end_II) %>% filter(CS_I == max(CS_I)) #IMPORTANT LINE
  #THIS LINE OF CODE GROUPS DATA BY coordinates of fragments and filter for each group by chicago score, keeping the one with highest score
  
  if(identical(new_data$CS_I, new_data$CS_II)) #QUALITY CONTROL STAGE, the scores should be equal for the "on purpose repeated fragments"
  {
    message("Cleaning seqmonk correct")
    new_data <- dplyr::slice(dplyr::ungroup(new_data), seq(1,dplyr::n(),2)) #data has to be ungrouped to retain original format, and one of the two copies are removed
    
    message(paste0((nrow(data)/2) - nrow(new_data))," interactions removed because were duplicated") #tells you how many duplicates there were
    
    return(new_data)
  }
  else(
    stop("Cleaning seqmonk error")
  )
}
   
data
# way to run it 
cjl1 <- clean_seqmonk("/media/blanca/JavierreLab/BLANCA/project/chicago_output/CJL1_5_seqmonk.txt") 
cjl1

cjl2 <- clean_seqmonk("~/Desktop/project/chicago_output/CJL2_5_seqmonk.txt")    
nrow(cjl1)
# TO SAVE THE FILES
#write.table(cjl1, file = "CJL1_5_seqmonk_cleaned.txt", sep = "\t")
#write.table(cjl2, file = "CJL2_5_seqmonk_cleaned.txt", sep = "\t")


