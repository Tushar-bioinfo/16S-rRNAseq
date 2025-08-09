##### this is a seperate script which i used to get the specific answers as per the questions given in the asssignment 
ncol(seqtable_nochimeras)
dim(seqtable_nochimeras)  ### 20 230

top5 <- sort(colSums(seqtable_nochimeras), decreasing = TRUE)[1:5]
top5


taxa[names(top5), ]

top5 <- sort(colSums(seqtable_nochimeras), decreasing = TRUE)[1:5]
sum(top5) / sum(seqtable_nochimeras) * 100


##### Data cleaning  and changing some parameters

out <- filterAndTrim(FnFs,   
                     filterFs, 
                     FnRs,   
                     filterRs, 
                     truncLen=c(240, 160),
                     maxN=0,  
                     maxEE=c(1, 2),
                     truncQ=12,
                     minQ=4,      
                     rm.phix=TRUE, 
                     compress=TRUE, multithread=TRUE) 


## computing % seq loss 

out_df <- as.data.frame(out)
out_df$percent_loss <- (out_df$reads.in - out_df$reads.out) / out_df$reads.in * 100
out_df



#####


taxa <- assignTaxonomy(seqtable_nochimeras, "/Users/tusharsingh/Downloads/data_set/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE, minBoot = 50)

taxa <- assignTaxonomy(seqtable_nochimeras, "/Users/tusharsingh/Downloads/data_set/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE, minBoot = 80)

####

track_new <- track
track <- read.csv("/Users/tusharsingh/Documents/Comp_genomics/Assignment_10/track.csv")

all.equal(track, track_new)


#### comparing  taxa NAs

# Genus level NA count (original)
sum(is.na(taxa[, "Genus"]))  ## 111 with 80 min.boot and 93 with 50 min.boot

# Genus level NA count (low confidence)
sum(is.na(taxa_lowconf[, "Genus"]))

sum(is.na(taxa))


head(track_new)
