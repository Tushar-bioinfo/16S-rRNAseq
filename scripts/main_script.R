library(dada2)
library(vegan)
library("RColorBrewer")

path <- "/Users/tusharsingh/Downloads/data_set/MiSeq_SOP"


FnFs <- sort(list.files(path, pattern = "_R1_001.fastq",full.name = T))
FnRs <- sort(list.files(path, pattern = "_R2_001.fastq",full.name = T))

#sample_names <- sapply(FnFs, function (x) {strsplit(basename(x),"_")[[1]][1]})


sample_names <- sapply(strsplit(basename(FnFs), "_"), `[`, 1)
  

## checking quality of the reads 

plotQualityProfile(FnFs[1:2])
plotQualityProfile(FnRs[1:2])


########## Reading filtered files ######

filterFs <-  file.path(path, "filtered", paste0(sample_names,"_F_filt.fastq.gz"))
filterRs <-  file.path(path, "filtered", paste0(sample_names,"_R_filt.fastq.gz"))

filterRs


names(filterFs) <- sample_names
names(filterRs) <- sample_names
 


###### trimming the low quality reads   #####
out <- filterAndTrim(FnFs,   
                     filterFs, 
                     FnRs,   
                     filterRs, 
                     truncLen=c(240, 160),
                     maxN=0,  
                     maxEE=2,
                     truncQ=2,
                     rm.phix=TRUE, 
                     compress=TRUE, multithread=TRUE) 

out

#### learning error rate #######


errF <- learnErrors(filterFs, multithread=TRUE)
errF <- learnErrors(filterRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)




##### sample inference 

dadaFs <- dada(filterFs, err=errF, multithread=TRUE)
dadaRs <- dada(filterRs, err=errF, multithread=TRUE)


######## merged paired reads 


merged <-  mergePairs(dadaFs, filterFs, dadaRs, filterRs, verbose=TRUE)


#head(merged)

seqtable <- makeSequenceTable(merged)
dim(seqtable)
seqtable[1:3,1:3]


##### Removing chimeras 

seqtable_nochimeras <- removeBimeraDenovo(seqtable, method="consensus", multithread=TRUE, verbose=TRUE)


dim(seqtable_nochimeras)


getN <- function(x) sum(getUniques(x))

track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merged, getN), rowSums(seqtable_nochimeras))


colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchimera")
rownames(track) <- sample_names
head(track)



##### Assignin taxonomy by byesian 


taxa <- assignTaxonomy(seqtable_nochimeras, "/Users/tusharsingh/Downloads/data_set/silva_nr99_v138_wSpecies_train_set.fa.gz", multithread=TRUE)

head(taxa)



#### Evaluating accruracy 

unqs_mock <- seqtable_nochimeras["Mock",]
unqs_mock <- sort(unqs_mock[unqs_mock>0], decreasing=TRUE) # Drop ASVs absent in the Mock
cat("DADA2 inferred", length(unqs_mock), "sample sequences present in the Mock community.\n")


mock_ref <- getSequences(file.path(path, "HMP_MOCK.v35.fasta"))
match_ref <- sum(sapply(names(unqs_mock), function(x) any(grepl(x, mock_ref))))
cat("Of those,", sum(match_ref), "were exact matches to the expected reference sequences.\n")


##### Shannon diversity 


samples_out <- rownames(seqtable_nochimeras)
subject <- sapply(strsplit(samples_out, "D"), `[`, 1)
gender  <- substr(subject,1,1)
subject <- substr(subject,2,999)
day     <- as.integer(sapply(strsplit(samples_out, "D"), `[`, 2))
samdf   <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples_out
samdf$col <- rep('red', nrow(samdf))
samdf$col[samdf$When == 'Early'] <- 'navy'



##### SD ####
shannonDiversity <- function(df) 
{
  p <- t( apply(df, 1, function(x){x/sum(x)}) )
  H <- apply(p , 1, function(x){x <- x[x > 0]; -sum( log(x) * x )})
  return(H)
}

SimpsonDiversity <- function(df) 
{
  p <- t( apply(df, 1, function(x){x/sum(x)}) )
  H <- apply(p , 1, function(x){x <- x[x > 0];1-sum( (x ^ 2))})
  return(H)
}
shD <- shannonDiversity(seqtable_nochimeras)
siD <- SimpsonDiversity(seqtable_nochimeras)
par(mfrow = c(1,2))
plot(samdf$Day,shD, col = samdf$col, pch = 19 ,las = 1, 
     xlab = 'Days', ylab = 'Shannon Diversity')
plot(samdf$Day,siD, col = samdf$col, pch = 19 ,las = 1, 
     xlab = 'Days', ylab = 'Simpson Diversity')


## changing data to proportions 

ps.prop <- t(apply(seqtable_nochimeras,1,function(x){x/sum(x)}))

NMDS <- metaMDS(ps.prop[-20,],     # remove mock line
                distance = 'bray')



plot(NMDS$points[,1], NMDS$points[,2], col = samdf$col,
     pch = 19, las = 1, xlab = 'MDS1', ylab = 'MDS2')

dim(taxa)

dim(seqtable_nochimeras)



######## Taxonomic distribution plot 

idTop20   <- order(colSums(seqtable_nochimeras), decreasing = T)[1:20]
top20     <- ps.prop[-nrow(seqtable_nochimeras),idTop20] #remove Mock
top20Taxa <- data.frame( taxa[idTop20, ] )

colNeed <- data.frame(1:9, brewer.pal(9, "Set1"))
top20Taxa$category <- as.numeric(as.factor(top20Taxa$Family))
top20Taxa$category[is.na(top20Taxa$category)] <- 7

top20Taxa$color    <- colNeed[match(top20Taxa$category, colNeed[,1]),2]
samdfUse <- samdf[-nrow(samdf),]
top20prop <- top20[, order(top20Taxa$Family, decreasing = T)]
top20Taxa <- top20Taxa[order(top20Taxa$Family,decreasing = T), ]

par(mfrow = c(1,3))
barplot( t(top20prop[samdfUse$When == 'Early',]), col = top20Taxa$color,
         las = 2)
barplot( t(top20prop[samdfUse$When == 'Late',]), col = top20Taxa$color,
         las = 2)

##for legend
legendMatrix <- unique(cbind(top20Taxa$Family,
                             top20Taxa$color))
par(mar = c(0,0,0,0))
plot(NA,NA, xlab = '', ylab = '',
     xaxt = 'n', yaxt = 'n',xlim = c(0,1), ylim = c(0,1),
     bty = 'n')
legend('center', legend = legendMatrix[,1],
       col = legendMatrix[,2], pch = 15,
       bty = 'n')


#track
#write.csv(track,"/Users/tusharsingh/Documents/Comp_genomics/Assignment_10/track_new.csv")

#write.table("track_new.csv",sep="\t",row.names = T,quote =F,col.names = T)




