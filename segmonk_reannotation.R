# This function reannotates "chicfile" (chicago output file), if it has the following format: 
# chr_I start_I   end_I gene_I            R_I  CS_I chr_II start_II end_II gene           R_II CS_II
# from a fraghindfile (which corresponds to all the fragments of hind3 of the genome)
###################################################################################################################

reanotation <- function(chicfile, fraghindfile){
  new_data <- read.table(fraghindfile,stringsAsFactors = F)
  hind3 <- dplyr::as_tibble(new_data[order(as.numeric(rownames(new_data))),], #restructure into a tibble
                               .name_repair = "minimal")
  colnames(hind3) <- c("chr","start","end","fragment", "gene") # rename columns colnames(new_data) <- c("chr","start_I","end_I","fragment_ID", "gene_I") # rename columns 
  print(str(hind3))
  hind3$gene[is.na(hind3$gene)] <- "non-annotated"
  chicfile_I <- chicfile[1:6] #must cheach the file format, and change this if necessary
  chicfile_II <- chicfile[7:12]
  left <- dplyr::left_join(chicfile_I, hind3, by=c('chr_I'='chr', 'start_I'='start', 'end_I'='end')) # we wannt all the rows of cjl1, we join with hind3 fragments
  right <- dplyr::left_join(chicfile_II, hind3, by=c('chr_II'='chr', 'start_II'='start', 'end_II'='end')) # we wannt all the rows of cjl1, we join with hind3 fragments
  left <- left[, c(1,2,3,8,5,6)] # Takes the columns of interest
  right <- right[, c(1,2,3,8,5,6)]
  colnames(left) <- c('chr_I', 'start_I', 'end_I', 'gene_I', 'R_I', 'CS_I') # renames left tibble
  chicreannotated <- cbind(left, right) # rebinds the file 
  chicreannotated <- dplyr::as_tibble(chicreannotated)
  return(chicreannotated)
  
}

#Use as file the output of running clean_segmonk.r
cjl1_annotated <- reanotation(cjl1, "~/Desktop/project/chicago_output/digest_and_probes_homo_sapiens_hg19.txt")
cjl2_annotated <- reanotation(cjl2, "~/Desktop/project/chicago_output/digest_and_probes_homo_sapiens_hg19.txt")
cjl2_annotated
# it works since nrow(cjl1) = nrow(cjl1_annotated)
# also we have cheac chr, start and end columns to see if they show the same numbers. 
table(is.na(cjl1_annotated$gene))

#write.table(cjl1_annotated, file = "CJL1_5_seqmonk_cleaned_annotated.txt", sep = "\t")
#write.table(cjl2_annotated, file = "CJL2_5_seqmonk_cleaned_annotated.txt", sep = "\t")

